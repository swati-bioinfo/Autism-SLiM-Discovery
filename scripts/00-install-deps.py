"""Check which packages are installed, install missing ones."""
import subprocess
import sys

required = ["requests", "biopython", "pandas", "numpy"]
missing = []

for pkg in required:
    result = subprocess.run(
        [sys.executable, "-m", "pip", "show", pkg],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        missing.append(pkg)

if missing:
    print(f"Installing missing: {', '.join(missing)}")
    subprocess.run([sys.executable, "-m", "pip", "install", *missing])
else:
    print("All required packages already installed.")
