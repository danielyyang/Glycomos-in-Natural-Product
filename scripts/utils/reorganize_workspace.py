import os
import shutil

def reorganize():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    # Ensure directories exist
    dirs_to_make = ["log", "reports", "output", "data", "lib", "scripts", "docs", "images", "tests"]
    for d in dirs_to_make:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)
        
    # Standard mapping for loose files
    file_mapping = {
        "debug_chirality_log.txt": "log",
        "debug_chirality_v2_log.txt": "log",
        "debug_rs_log.txt": "log",
        "isomers_log.txt": "log",
        "process_log.txt": "log",
        "rs_library_log.txt": "log",
        "rs_sigs_clean.txt": "log",
        "test_manual_output.txt": "log",
        "test_mods_fail.txt": "log",
        "test_output.txt": "log",
        "test_refine_log.txt": "log",
        "verify_log.txt": "log",
        
        "verification_results.txt": "reports",
        "verification_results_v2.txt": "reports",
        "sheet_info.txt": "reports",
    }
    
    moved_count = 0
    for filename, target_folder in file_mapping.items():
        src = os.path.join(base_dir, filename)
        dest = os.path.join(base_dir, target_folder, filename)
        if os.path.exists(src):
            shutil.move(src, dest)
            print(f"Moved {filename} -> {target_folder}/")
            moved_count += 1
            
    print(f"\n====================\nOrganization Complete. Moved {moved_count} files.")

if __name__ == "__main__":
    reorganize()
