import sqlite3
conn = sqlite3.connect(r'data\chembl_36\chembl_36_sqlite\chembl_36.db')
cur = conn.cursor()

# List tables
cur.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
tables = [r[0] for r in cur.fetchall()]
print(f"Total tables: {len(tables)}")
print("Tables:", tables[:30])

# Check compound_structures table
for t in ['compound_structures', 'molecule_dictionary', 'compound_properties']:
    if t in tables:
        cur.execute(f"PRAGMA table_info({t})")
        cols = [r[1] for r in cur.fetchall()]
        print(f"\n{t} columns: {cols}")
        cur.execute(f"SELECT COUNT(*) FROM {t}")
        print(f"  Rows: {cur.fetchone()[0]:,}")

conn.close()
