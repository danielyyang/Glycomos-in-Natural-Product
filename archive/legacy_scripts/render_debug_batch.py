import pandas as pd
import os
import sys
from PIL import Image, ImageDraw, ImageFont
import textwrap

sys.stdout.reconfigure(encoding='utf-8')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

try:
    from lib.visualizer import StructureVisualizer
except ImportError:
    print("❌ 无法导入 lib.visualizer")
    sys.exit(1)

def create_data_card(image_path, row_data):
    """
    在原有的 RDKit 图像底部拼接一个信息面板，写入核心数据。
    """
    try:
        # 打开原图 (RDKit 生成的图)
        mol_img = Image.open(image_path)
        img_w, img_h = mol_img.size
        
        # 提取要展示的信息，处理空值
        identifier = str(row_data.get('identifier', 'N/A'))
        name = str(row_data.get('name', 'N/A'))
        if name == 'nan' or name == '': name = 'Unknown'
        
        seq = str(row_data.get('Sugar_Sequence', 'N/A'))
        mods = str(row_data.get('Sugar_Functional_Group', 'None'))
        organisms = str(row_data.get('organisms', 'Unknown'))
        dois = str(row_data.get('dois', 'No Reference'))
        
        # 格式化文本 (处理太长导致出界的问题)
        text_lines = [
            f"ID: {identifier} | Name: {name}",
            f"Sequence: {textwrap.shorten(seq, width=80, placeholder='...')}",
        ]
        if mods != 'None' and mods != 'nan':
            text_lines.append(f"Modifications: {textwrap.shorten(mods, width=80, placeholder='...')}")
            
        text_lines.append(f"Organism: {textwrap.shorten(organisms, width=80, placeholder='...')}")
        text_lines.append(f"DOI: {textwrap.shorten(dois, width=80, placeholder='...')}")
        
        # 计算文本面板的高度 (每行预留 30 像素)
        panel_h = len(text_lines) * 30 + 40 
        
        # 创建一张新的空白大图 (包含原图和底部文本面板)
        new_img = Image.new('RGB', (img_w, img_h + panel_h), 'white')
        new_img.paste(mol_img, (0, 0))
        
        draw = ImageDraw.Draw(new_img)
        
        # 尝试加载默认字体，如果失败则使用系统默认内置字体
        try:
            # 你可以替换为你系统中的真实字体路径，如 arial.ttf
            font = ImageFont.truetype("arial.ttf", 16) 
        except IOError:
            font = ImageFont.load_default()
            
        # 将文本逐行画在面板上
        y_text = img_h + 10
        for line in text_lines:
            draw.text((20, y_text), line, font=font, fill=(0, 0, 0))
            y_text += 30
            
        # 覆盖保存为带名片的新图
        new_img.save(image_path)
        
    except Exception as e:
        print(f"⚠️ 生成数据名片失败 {row_data.get('identifier')}: {e}")

def render_debug_images(csv_path, output_dir="reports/debug_images", limit=100):
    """主调度引擎"""
    print(f"🎨 启动 Phase 7 数据名片渲染引擎...")
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(csv_path).head(limit)
    visualizer = StructureVisualizer()
    
    success_count = 0
    for idx, row in df.iterrows():
        identifier = str(row['identifier'])
        smiles = str(row['canonical_smiles'])
        if pd.isna(smiles) or smiles == 'nan': continue
            
        # 先按常规画出分子图
        img_path = os.path.join(output_dir, f"{identifier}_card.png")
        try:
            visualizer.analyze_glycolipid(smiles, output_path=img_path)
            
            # 🌟 核心：画完图后，调用拼接函数盖上文本面板！
            if os.path.exists(img_path):
                create_data_card(img_path, row)
                
            success_count += 1
            if success_count % 20 == 0: print(f"   => 已生成 {success_count} 张数据名片...")
        except Exception as e:
            pass
            
    print(f"\n✅ 数据名片批次完成！成功生成 {success_count} 张。请查阅 {output_dir}")

if __name__ == "__main__":
    render_debug_images("reports/Coconut_Sugar_Phase6.csv", limit=50)
