#!/usr/bin/env python3
import subprocess, sys, os, re
from pathlib import Path

# Configuration
R_PATH = None  # Set to r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" or leave None
SKIP = ["_04_second_layer_calculation.py"]

def find_r():
    if R_PATH and os.path.exists(R_PATH): return R_PATH
    if sys.platform == 'win32':
        for base in [r"C:\Program Files\R", r"C:\Program Files (x86)\R"]:
            if os.path.exists(base):
                for v in Path(base).glob("R-*/bin/Rscript.exe"):
                    return str(v)
    try:
        subprocess.run(['Rscript', '--version'], capture_output=True, timeout=5)
        return 'Rscript'
    except: pass
    return None

def find_scripts(dir="."):
    scripts = []
    for f in Path(dir).rglob("*"):
        if f.is_file() and f.name not in SKIP:
            m = re.search(r'_(\d+)_', f.name)
            if m and f.suffix in ['.py', '.R', '.r']:
                scripts.append((int(m.group(1)), f.suffix == '.py', f.name, str(f)))
    return sorted(scripts, key=lambda x: (x[0], x[2]))

def run(cmd, path):
    try:
        subprocess.run(cmd, check=True, cwd=os.path.dirname(path) or '.')
        return True
    except: return False

def main():
    r = find_r()
    scripts = find_scripts(os.path.dirname(os.path.abspath(__file__)))
    
    if not scripts:
        print("No scripts found")
        return 1
    
    print(f"Running {len(scripts)} script(s)...\n")
    failed = []
    
    for num, is_py, name, path in scripts:
        print(f"[{num:02d}] {name}")
        if is_py:
            if not run([sys.executable, path], path): failed.append(name)
        elif r:
            if not run([r, path], path): failed.append(name)
        else:
            print("  -> Skipped (R not found)")
        if not run([sys.executable, path], path): failed.append(name)
    
    print(f"\n{len(scripts) - len(failed)}/{len(scripts)} successful")
    return 1 if failed else 0

if __name__ == "__main__":
    sys.exit(main())