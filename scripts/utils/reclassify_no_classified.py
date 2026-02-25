import pandas as pd
import os
import shutil

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    csv_file = os.path.join(base_dir, "output", "Sugar_Sort_Refined.CSV")
    excel_file = os.path.join(base_dir, "data", "processed", "Coconut_Sugars_Unique.xlsx")
    img_dir1 = os.path.join(base_dir, "images", "No_Classificed")
    img_dir2 = os.path.join(base_dir, "images", "No Classified")
    
    img_dir = img_dir1 if os.path.exists(img_dir1) else (img_dir2 if os.path.exists(img_dir2) else None)
    
    print("Reading CSV...")
    df = pd.read_csv(csv_file)
    
    # Identify the real sheet name from unique sheet sources
    sheets = df['Sheet_Source'].unique()
    target_sheet = next((s for s in sheets if 'No' in str(s) and ('Classifi' in str(s) or 'classifi' in str(s))), None)
    
    if not target_sheet:
        print("Could not find No Classified sheet among:", sheets)
        target_sheet = 'No_Classificed' if 'No_Classificed' in sheets else 'No Classified'
        
    print(f"Target Sheet for reclassifying: {target_sheet}")
    
    # Build scaffold_to_class mapping based on frequency of other classes
    mapping_df = df[df['Sheet_Source'] != target_sheet]
    # Filter known specific sheets if any
    valid_map_df = mapping_df[mapping_df['Sheet_Source'].notna()]
    mapping = {}
    if not valid_map_df.empty:
        mapping = valid_map_df.groupby('Aglycan_Scaffold_ID')['Sheet_Source'].agg(lambda x: x.value_counts().index[0]).to_dict()
    
    # Reclassify No Classified subset
    no_class_df = df[df['Sheet_Source'] == target_sheet].copy()
    
    def get_class(row):
        scaf_id = row['Aglycan_Scaffold_ID']
        hint = row.get('Scaffold_Class_Hint')
        if pd.notna(scaf_id) and scaf_id in mapping and scaf_id != 'Linear_Or_No_Scaffold' and scaf_id != 'Polycyclic':
            return mapping[scaf_id]
        elif pd.notna(hint) and str(hint).strip() != '':
            return str(hint).strip()
        elif scaf_id == 'Polycyclic':
            return 'Polycyclic_Organic_Framework'
        else:
            return 'Unknown'
            
    no_class_df['Classification'] = no_class_df.apply(get_class, axis=1)
    
    print("Reclassification summary:")
    print(no_class_df['Classification'].value_counts())
    
    # Save to Excel
    print("Saving to New_Classfied sheet in Excel...")
    try:
        from openpyxl import load_workbook
        with pd.ExcelWriter(excel_file, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
            no_class_df.to_excel(writer, sheet_name='New_Classfied', index=False)
        print("Appended New_Classfied sheet successfully.")
    except Exception as e:
        print(f"Could not append to existing Excel: {e}")
        backup_out = os.path.join(base_dir, "output", "New_Classfied_Data.xlsx")
        no_class_df.to_excel(backup_out, sheet_name='New_Classfied', index=False)
        print(f"Saved to alternative location: {backup_out}")
        
    # Process Images
    if img_dir and os.path.exists(img_dir):
        print(f"Organizing images in {img_dir}...")
        # map identifier to Classification
        id_to_class = dict(zip(no_class_df['identifier'].astype(str), no_class_df['Classification']))
        
        moved_count = 0
        for fname in os.listdir(img_dir):
            fpath = os.path.join(img_dir, fname)
            if not os.path.isfile(fpath):
                continue
            
            file_class = 'Unknown'
            # Look for identifier prefix
            for ident, cls in id_to_class.items():
                safe_id = "".join([c for c in ident if c.isalnum() or c in ('-','_')]).strip()
                if fname.startswith(safe_id):
                    file_class = cls
                    break
            
            # Clean up class directory name
            safe_class_name = "".join([c for c in file_class if c.isalnum() or c in (' ', '-', '_')]).strip()
            if not safe_class_name:
                safe_class_name = "Unknown"
                
            class_dir = os.path.join(img_dir, safe_class_name)
            os.makedirs(class_dir, exist_ok=True)
            
            shutil.move(fpath, os.path.join(class_dir, fname))
            moved_count += 1
            
        print(f"Moved {moved_count} images based on new classification.")
    else:
        print("No image directory found for", target_sheet)
        print("If you need images generated, run split_glycosides.py or generate_test_images.py first.")
        
    print("✅ 重分类过程完成 (Reclassification complete).")

if __name__ == '__main__':
    main()
