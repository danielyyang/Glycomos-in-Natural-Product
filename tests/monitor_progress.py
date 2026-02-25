
import os
import time
import pandas as pd
import tkinter as tk
from tkinter import ttk, messagebox
import threading

IMAGE_BASE_DIR = r"d:\Lipid Database\images"
CSV_FILE = r"d:\Lipid Database\data\processed\Coconut_Sugar_Sort.csv"

class ProgressMonitorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Glycoside Image Generation Progress")
        self.root.geometry("400x250")
        
        self.total_images = 0
        self.current_count = 0
        self.is_running = True
        
        # Load Total Count
        self.load_total_count()
        
        # UI Elements
        self.lbl_status = ttk.Label(root, text="Scanning Images...", font=("Arial", 12))
        self.lbl_status.pack(pady=10)
        
        self.progress_bar = ttk.Progressbar(root, orient="horizontal", length=300, mode="determinate")
        self.progress_bar.pack(pady=10)
        
        self.lbl_count = ttk.Label(root, text=f"0 / {self.total_images}", font=("Arial", 10))
        self.lbl_count.pack(pady=5)
        
        self.lbl_percent = ttk.Label(root, text="0.0%", font=("Arial", 12, "bold"))
        self.lbl_percent.pack(pady=5)
        
        self.btn_refresh = ttk.Button(root, text="Force Refresh", command=self.update_progress)
        self.btn_refresh.pack(pady=10)
        
        # Start monitoring thread
        self.monitor_thread = threading.Thread(target=self.monitor_loop, daemon=True)
        self.monitor_thread.start()
        
    def load_total_count(self):
        try:
            if os.path.exists(CSV_FILE):
                # Count rows where Is_Glycoside is True
                # Reading only needed column for speed
                df = pd.read_csv(CSV_FILE, usecols=['Is_Glycoside'])
                self.total_images = len(df[df['Is_Glycoside'] == True])
            else:
                self.total_images = 93255 # Fallback based on known count
        except Exception as e:
            print(f"Error loading CSV: {e}")
            self.total_images = 93255

    def count_images(self):
        count = 0
        for root, dirs, files in os.walk(IMAGE_BASE_DIR):
            count += len([f for f in files if f.endswith('.png')])
        return count

    def update_progress(self):
        self.current_count = self.count_images()
        
        # Update UI in main thread
        self.root.after(0, self._update_ui)
        
    def _update_ui(self):
        self.lbl_count.config(text=f"{self.current_count} / {self.total_images}")
        
        if self.total_images > 0:
            percent = (self.current_count / self.total_images) * 100
            self.progress_bar['value'] = percent
            self.lbl_percent.config(text=f"{percent:.1f}%")
            
            if self.current_count >= self.total_images:
                self.lbl_status.config(text="Generation Complete!", foreground="green")
            else:
                self.lbl_status.config(text="Generating...", foreground="blue")

    def monitor_loop(self):
        while self.is_running:
            self.update_progress()
            time.sleep(2) # Check every 2 seconds

if __name__ == "__main__":
    if not os.path.exists(IMAGE_BASE_DIR):
        os.makedirs(IMAGE_BASE_DIR, exist_ok=True)
        
    root = tk.Tk()
    app = ProgressMonitorApp(root)
    root.mainloop()
