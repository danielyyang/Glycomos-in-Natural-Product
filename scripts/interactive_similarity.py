#!/usr/bin/env python3
"""
==========================================================================
  [EN] Glycan TED Interactive Explorer
       A local web-based tool for real-time glycan sequence similarity
       analysis with adjustable scoring weights. Runs a lightweight
       HTTP server on localhost — no external dependencies beyond stdlib.

  [CN] 糖链编辑距离 (TED) 交互式探索器
       基于本地 Web 的实时糖链序列相似性分析工具。
       通过滑块实时调整打分权重 (节点/边/拓扑), 即时查看 Top-N 结果。
       使用 Python 内置 http.server, 无需额外依赖。

  用法 (Usage):
       python scripts/interactive_similarity.py
       → 浏览器自动打开 http://localhost:8765

  [TEST DATA ONLY]
==========================================================================
"""
import os
import sys
import json
import time
import threading
import webbrowser
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import parse_qs, urlparse
from typing import List

# 确保项目根目录在 sys.path 中 (Ensure project root is in sys.path)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from lib.glycan_similarity import (
    parseGlycanSequence,
    computeSimilarity,
    computeTreeEditDistance,
    searchSimilar,
    ScoringWeights,
    GlycanNode,
)

# ── 配置 (Configuration) ─────────────────────────────────────────────
PORT = 8765
CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "reports", "GlycoNP_Saponin_DB.csv")

# ── 全局数据缓存 (Global data cache) ─────────────────────────────────
ALL_SEQUENCES: List[str] = []
SEQUENCE_NAMES: List[str] = []


def loadSequences():
    """加载所有皂苷糖链序列到内存 (Load all saponin sugar sequences into memory)."""
    global ALL_SEQUENCES, SEQUENCE_NAMES
    print("Loading saponin sequences...")
    df = pd.read_csv(CSV_PATH, low_memory=False, usecols=[
        "Consensus_Sugar_Sequence", "name", "identifier",
    ])
    df = df[df["Consensus_Sugar_Sequence"].notna()]
    ALL_SEQUENCES = df["Consensus_Sugar_Sequence"].tolist()
    names = []
    for _, row in df.iterrows():
        name = str(row.get("name", "")) if pd.notna(row.get("name")) else ""
        ident = str(row.get("identifier", "")) if pd.notna(row.get("identifier")) else ""
        names.append(name[:40] or ident[:20] or "Unknown")
    SEQUENCE_NAMES = names
    print(f"  Loaded {len(ALL_SEQUENCES):,} sequences.")


# ── HTML 交互界面 (Interactive HTML UI) ───────────────────────────────

HTML_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GlycoNP TED Explorer — Glycan Similarity</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: linear-gradient(135deg, #0a0e1a 0%, #111827 50%, #1a1a2e 100%);
    color: #e0e0e0; min-height: 100vh; padding: 24px;
  }
  h1 {
    text-align: center; font-size: 1.6em; color: #7ec8e3;
    margin-bottom: 6px; letter-spacing: 0.5px;
  }
  .subtitle { text-align: center; color: #6b7280; font-size: 0.85em; margin-bottom: 24px; }

  .container { max-width: 1100px; margin: 0 auto; }
  .panel {
    background: rgba(17, 24, 39, 0.85); border: 1px solid #2a3a4a;
    border-radius: 12px; padding: 20px; margin-bottom: 16px;
    backdrop-filter: blur(10px);
  }
  .panel-title {
    color: #7ec8e3; font-size: 1.05em; font-weight: 600;
    margin-bottom: 14px; padding-bottom: 8px; border-bottom: 1px solid #2a3a4a;
  }

  .input-row { display: flex; gap: 12px; align-items: flex-end; margin-bottom: 14px; }
  .input-group { flex: 1; }
  .input-group label { display: block; color: #9ca3af; font-size: 0.82em; margin-bottom: 4px; }
  .input-group input[type="text"] {
    width: 100%; padding: 10px 14px; background: #1a2332; border: 1px solid #374151;
    border-radius: 8px; color: #e0e0e0; font-size: 0.95em; font-family: 'Courier New', monospace;
    transition: border-color 0.2s;
  }
  .input-group input[type="text"]:focus { border-color: #7ec8e3; outline: none; }

  .slider-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 16px; }
  .slider-item { text-align: center; }
  .slider-item label {
    display: block; color: #9ca3af; font-size: 0.82em; margin-bottom: 6px;
  }
  .slider-item input[type="range"] {
    width: 100%; accent-color: #7ec8e3; cursor: pointer;
  }
  .slider-value {
    display: inline-block; background: #7ec8e3; color: #0a0e1a;
    padding: 2px 10px; border-radius: 12px; font-weight: 700;
    font-size: 0.9em; margin-top: 4px; min-width: 40px;
  }

  .btn {
    padding: 10px 28px; background: linear-gradient(135deg, #2563eb, #7c3aed);
    color: white; border: none; border-radius: 8px; cursor: pointer;
    font-size: 0.95em; font-weight: 600; transition: transform 0.15s, box-shadow 0.15s;
    letter-spacing: 0.3px;
  }
  .btn:hover { transform: translateY(-1px); box-shadow: 0 4px 12px rgba(37, 99, 235, 0.4); }
  .btn:active { transform: translateY(0); }
  .btn:disabled { opacity: 0.5; cursor: not-allowed; transform: none; }

  .status { color: #6b7280; font-size: 0.82em; margin-top: 8px; min-height: 18px; }
  .status.loading { color: #fbbf24; }
  .status.error { color: #f87171; }

  /* Results table */
  .results-table { width: 100%; border-collapse: collapse; font-size: 0.88em; margin-top: 10px; }
  .results-table th {
    background: #1a2332; color: #7ec8e3; padding: 8px 6px; text-align: left;
    border-bottom: 2px solid #374151; position: sticky; top: 0; font-size: 0.85em;
  }
  .results-table td {
    padding: 7px 6px; border-bottom: 1px solid #1e293b;
    vertical-align: top; word-break: break-word;
  }
  .results-table tr:hover td { background: rgba(126, 200, 227, 0.05); }

  .sim-bar {
    display: inline-block; height: 16px; border-radius: 3px; margin-right: 6px;
    vertical-align: middle; transition: width 0.3s;
  }
  .sim-score { font-weight: 700; font-family: 'Courier New', monospace; }
  .seq-cell { font-family: 'Courier New', monospace; font-size: 0.85em; color: #a0d0f0; }
  .name-cell { color: #9ca3af; font-size: 0.82em; }
  .rank-cell { color: #4ade80; font-weight: 700; text-align: center; }

  .examples { margin-top: 10px; }
  .example-tag {
    display: inline-block; padding: 3px 10px; margin: 3px 4px;
    background: #1a2332; border: 1px solid #374151; border-radius: 16px;
    color: #93c5fd; font-size: 0.78em; cursor: pointer;
    font-family: 'Courier New', monospace; transition: all 0.15s;
  }
  .example-tag:hover { background: #2563eb33; border-color: #2563eb; color: #60a5fa; }

  .info-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin-top: 12px; }
  .info-card {
    background: #1a2332; border: 1px solid #374151; border-radius: 8px;
    padding: 12px; text-align: center;
  }
  .info-card .num { font-size: 1.4em; font-weight: 700; color: #7ec8e3; }
  .info-card .label { font-size: 0.75em; color: #6b7280; margin-top: 2px; }

  @keyframes pulse { 0%,100% { opacity: 1; } 50% { opacity: 0.5; } }
  .loading-anim { animation: pulse 1.2s ease-in-out infinite; }
</style>
</head>
<body>
<div class="container">
  <h1>🧬 GlycoNP TED Explorer</h1>
  <p class="subtitle">Glycan Tree Edit Distance — Interactive Similarity Search with Real-Time Weight Tuning</p>

  <!-- Query Panel -->
  <div class="panel">
    <div class="panel-title">🔍 Query Sequence / 查询序列</div>
    <div class="input-row">
      <div class="input-group" style="flex:3">
        <label>Sugar Sequence (GlycoNP format) / 糖链序列</label>
        <input type="text" id="querySeq" value="D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA"
               placeholder="e.g. D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA">
      </div>
      <div class="input-group" style="flex:0 0 auto">
        <label>Top-N</label>
        <input type="text" id="topN" value="15" style="width:60px; text-align:center">
      </div>
      <div style="flex:0 0 auto">
        <button class="btn" id="searchBtn" onclick="doSearch()">⚡ Search</button>
      </div>
    </div>
    <div class="examples">
      <span style="color:#6b7280;font-size:0.78em">Examples / 示例: </span>
      <span class="example-tag" onclick="setQuery('D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA')">GlcA trisaccharide</span>
      <span class="example-tag" onclick="setQuery('L-Rha-(a1-2)-D-Glc')">Rha-Glc disaccharide</span>
      <span class="example-tag" onclick="setQuery('[D-Glc-(b1-3)]-D-Gal-(b1-4)-D-Xyl')">Branched triglycan</span>
      <span class="example-tag" onclick="setQuery('D-Glc')">Single Glc</span>
      <span class="example-tag" onclick="setQuery('L-Rha ; D-Glc-(a1-4)-D-Gal')">Multi-chain</span>
    </div>
    <div class="status" id="status"></div>
  </div>

  <!-- Weight Panel -->
  <div class="panel">
    <div class="panel-title">⚖️ Scoring Weights / 打分权重 (实时调整后点击 Search)</div>
    <div class="slider-grid">
      <div class="slider-item">
        <label>🧩 Node Weight / 节点 (单糖匹配)</label>
        <input type="range" id="wNode" min="0" max="30" value="10" step="1"
               oninput="updateSliderLabel('wNode','wNodeVal')">
        <div><span class="slider-value" id="wNodeVal">1.0</span></div>
        <div style="color:#6b7280;font-size:0.7em;margin-top:4px">
          Higher → sugar identity matters more<br>越高 → 糖种类差异影响越大
        </div>
      </div>
      <div class="slider-item">
        <label>🔗 Edge Weight / 边 (连接键匹配)</label>
        <input type="range" id="wEdge" min="0" max="30" value="10" step="1"
               oninput="updateSliderLabel('wEdge','wEdgeVal')">
        <div><span class="slider-value" id="wEdgeVal">1.0</span></div>
        <div style="color:#6b7280;font-size:0.7em;margin-top:4px">
          Higher → α/β and position differences penalized more<br>越高 → α/β 和连接位置差异惩罚越大
        </div>
      </div>
      <div class="slider-item">
        <label>🌲 Topology Weight / 拓扑 (分支结构)</label>
        <input type="range" id="wTopo" min="0" max="30" value="5" step="1"
               oninput="updateSliderLabel('wTopo','wTopoVal')">
        <div><span class="slider-value" id="wTopoVal">0.5</span></div>
        <div style="color:#6b7280;font-size:0.7em;margin-top:4px">
          Higher → branching pattern differences matter more<br>越高 → 分支拓扑差异影响越大
        </div>
      </div>
    </div>
  </div>

  <!-- Results Panel -->
  <div class="panel" id="resultsPanel" style="display:none">
    <div class="panel-title" id="resultsTitle">📊 Results / 结果</div>
    <div class="info-grid" id="infoGrid"></div>
    <table class="results-table">
      <thead>
        <tr>
          <th style="width:3%">#</th>
          <th style="width:18%">Similarity / 相似度</th>
          <th style="width:8%">Distance</th>
          <th style="width:38%">Sugar Sequence / 糖链序列</th>
          <th style="width:18%">Name / 名称</th>
        </tr>
      </thead>
      <tbody id="resultsBody"></tbody>
    </table>
  </div>
</div>

<script>
function updateSliderLabel(sliderId, labelId) {
  const v = document.getElementById(sliderId).value / 10;
  document.getElementById(labelId).textContent = v.toFixed(1);
}

function setQuery(seq) {
  document.getElementById('querySeq').value = seq;
}

async function doSearch() {
  const seq = document.getElementById('querySeq').value.trim();
  const topN = parseInt(document.getElementById('topN').value) || 15;
  const wNode = document.getElementById('wNode').value / 10;
  const wEdge = document.getElementById('wEdge').value / 10;
  const wTopo = document.getElementById('wTopo').value / 10;

  if (!seq) { setStatus('Please enter a sequence', 'error'); return; }

  const btn = document.getElementById('searchBtn');
  btn.disabled = true; btn.textContent = '⏳ Searching...';
  setStatus('Computing TED across database...', 'loading');

  try {
    const params = new URLSearchParams({
      q: seq, n: topN, wn: wNode, we: wEdge, wt: wTopo
    });
    const resp = await fetch('/api/search?' + params.toString());
    const data = await resp.json();

    if (data.error) { setStatus(data.error, 'error'); return; }

    renderResults(data, seq);
    setStatus(`Found ${data.results.length} matches in ${data.elapsed_ms}ms (searched ${data.total_candidates} sequences)`, '');
  } catch (e) {
    setStatus('Network error: ' + e.message, 'error');
  } finally {
    btn.disabled = false; btn.textContent = '⚡ Search';
  }
}

function setStatus(msg, cls) {
  const el = document.getElementById('status');
  el.textContent = msg; el.className = 'status ' + (cls || '');
}

function renderResults(data, querySeq) {
  const panel = document.getElementById('resultsPanel');
  panel.style.display = 'block';

  // Info cards
  const grid = document.getElementById('infoGrid');
  const qInfo = data.query_info || {};
  grid.innerHTML = `
    <div class="info-card"><div class="num">${qInfo.node_count || '?'}</div><div class="label">Query Nodes / 查询节点数</div></div>
    <div class="info-card"><div class="num">${qInfo.depth || '?'}</div><div class="label">Query Depth / 查询深度</div></div>
    <div class="info-card"><div class="num">${data.total_candidates || '?'}</div><div class="label">Database Size / 数据库大小</div></div>
  `;

  // Table body
  const tbody = document.getElementById('resultsBody');
  tbody.innerHTML = '';
  data.results.forEach((r, i) => {
    const pct = Math.round(r.similarity * 100);
    const hue = Math.round(r.similarity * 120); // 0=red, 120=green
    const barColor = `hsl(${hue}, 70%, 50%)`;
    const row = document.createElement('tr');
    row.innerHTML = `
      <td class="rank-cell">${i + 1}</td>
      <td>
        <span class="sim-bar" style="width:${pct}px; background:${barColor}"></span>
        <span class="sim-score" style="color:${barColor}">${(r.similarity * 100).toFixed(1)}%</span>
      </td>
      <td style="color:#9ca3af;font-family:monospace">${r.distance.toFixed(2)}</td>
      <td class="seq-cell">${escapeHtml(r.sequence)}</td>
      <td class="name-cell">${escapeHtml(r.name || '')}</td>
    `;
    tbody.appendChild(row);
  });

  document.getElementById('resultsTitle').textContent =
    `📊 Top ${data.results.length} Results / 最相似 ${data.results.length} 条 (Weights: Node=${data.weights.node}, Edge=${data.weights.edge}, Topo=${data.weights.topology})`;
}

function escapeHtml(text) {
  const d = document.createElement('div');
  d.textContent = text;
  return d.innerHTML;
}

// Auto-search on Enter
document.getElementById('querySeq').addEventListener('keydown', e => {
  if (e.key === 'Enter') doSearch();
});
</script>
</body>
</html>"""


# ── HTTP 请求处理器 (HTTP Request Handler) ────────────────────────────

class TedHandler(BaseHTTPRequestHandler):
    """处理 Web UI 和 API 请求 (Handle Web UI and API requests)."""

    def do_GET(self):
        """处理 GET 请求: 页面 or API (Handle GET: page or API)."""
        parsed = urlparse(self.path)

        if parsed.path == "/" or parsed.path == "/index.html":
            self._servePage()
        elif parsed.path == "/api/search":
            self._handleSearch(parsed.query)
        else:
            self.send_error(404)

    def _servePage(self):
        """返回交互式 HTML 页面 (Serve the interactive HTML page)."""
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(HTML_PAGE.encode("utf-8"))

    def _handleSearch(self, queryString: str):
        """执行 TED 检索并返回 JSON (Run TED search and return JSON).

        参数 (Parameters):
            q: 查询序列 (query sequence)
            n: Top-N 结果数
            wn: 节点权重 (node weight)
            we: 边权重 (edge weight)
            wt: 拓扑权重 (topology weight)
        """
        params = parse_qs(queryString)
        querySeq = params.get("q", [""])[0]
        topN = int(params.get("n", ["15"])[0])
        weightNode = float(params.get("wn", ["1.0"])[0])
        weightEdge = float(params.get("we", ["1.0"])[0])
        weightTopo = float(params.get("wt", ["0.5"])[0])

        # 验证输入 (Validate input)
        if not querySeq:
            self._sendJson({"error": "No query sequence provided"})
            return

        queryTree = parseGlycanSequence(querySeq)
        if queryTree is None:
            self._sendJson({"error": f"Failed to parse query: {querySeq}"})
            return

        weights = ScoringWeights(
            nodeWeight=weightNode,
            edgeWeight=weightEdge,
            topologyWeight=weightTopo,
        )

        # 执行检索 (Run search)
        t0 = time.time()
        results = searchSimilar(querySeq, ALL_SEQUENCES, topN=topN, weights=weights)
        elapsed = time.time() - t0

        # 附加名称 (Attach names)
        for r in results:
            idx = r.get("index", -1)
            r["name"] = SEQUENCE_NAMES[idx] if 0 <= idx < len(SEQUENCE_NAMES) else ""

        response = {
            "query": querySeq,
            "query_info": {
                "node_count": queryTree.nodeCount(),
                "depth": queryTree.depth(),
                "leaf_count": queryTree.leafCount(),
            },
            "weights": {
                "node": weightNode,
                "edge": weightEdge,
                "topology": weightTopo,
            },
            "total_candidates": len(ALL_SEQUENCES),
            "elapsed_ms": round(elapsed * 1000),
            "results": results,
        }
        self._sendJson(response)

    def _sendJson(self, data: dict):
        """发送 JSON 响应 (Send JSON response)."""
        self.send_response(200)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(json.dumps(data, ensure_ascii=False).encode("utf-8"))

    def log_message(self, format, *args):
        """静默日志 (Suppress access logs for cleaner output)."""
        pass


# ── 主入口 (Main entry point) ─────────────────────────────────────────

def main():
    print("=" * 60)
    print("  GlycoNP TED Explorer — Interactive Similarity Search")
    print("=" * 60)

    loadSequences()

    # 注意: 搜索 26K 序列可能需要几秒钟
    # Note: searching 26K sequences may take a few seconds per query
    # 对于大规模使用, 建议采样子集 (For large-scale use, subsample)
    print(f"\n  Server starting on http://localhost:{PORT}")
    print(f"  Press Ctrl+C to stop.\n")

    # 延迟打开浏览器 (Delayed browser open)
    def openBrowser():
        time.sleep(1.0)
        webbrowser.open(f"http://localhost:{PORT}")
    threading.Thread(target=openBrowser, daemon=True).start()

    server = HTTPServer(("localhost", PORT), TedHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n  Server stopped.")
        server.server_close()


if __name__ == "__main__":
    main()
