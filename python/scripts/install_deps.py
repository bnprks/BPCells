
import tempfile
from pathlib import Path
import sys
import subprocess

# Run the install script for each dependency. Do builds in a temporary directory

file_dir = Path(__file__).parent.resolve()

temp_dir = tempfile.TemporaryDirectory()

if len(sys.argv) > 1:
    sys.argv[1] = Path(sys.argv[1]).resolve()

def run_script(script_name):
    args = ["bash", file_dir / script_name]
    if len(sys.argv) > 1:
        args.append(sys.argv[1])
    subprocess.run(args, check=True, stdout=sys.stdout, stderr=sys.stderr, cwd=temp_dir.name)

print("Running builds in ", temp_dir.name)
run_script("install_eigen.sh")
run_script("install_hdf5.sh")
run_script("install_highway.sh")