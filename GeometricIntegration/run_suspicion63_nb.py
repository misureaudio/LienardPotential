"""
Headless execution + verification of suspicion63_implicit_midpoint_blowup.ipynb
Run: D:/Source/hermes-dir/.venv/Scripts/python.exe run_suspicion63_nb.py
"""
import json
import os
import sys
import tempfile

import nbformat
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
NB_PATH = os.path.join(HERE, "suspicion63_implicit_midpoint_blowup.ipynb")
VENV_PY = r"D:\Source\hermes-dir\.venv\Scripts\python.exe"

tmp = tempfile.mkdtemp(prefix="jupyter_spec_")
kernels_dir = os.path.join(tmp, "kernels", "suspicion63-venv")
os.makedirs(kernels_dir, exist_ok=True)
with open(os.path.join(kernels_dir, "kernel.json"), "w") as fh:
    json.dump({
        "argv": [VENV_PY, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "Python 3 (project venv)",
        "language": "python",
    }, fh)

os.environ["JUPYTER_DATA_DIR"] = tmp
os.environ["MPLBACKEND"] = "Agg"

nb = nbformat.read(NB_PATH, as_version=4)
client = NotebookClient(nb, timeout=1200, kernel_name="suspicion63-venv")
client.execute()

nbformat.write(nb, NB_PATH)
print("executed and re-saved %s" % NB_PATH)

# ---- verification ------------------------------------------------------
nb = nbformat.read(NB_PATH, as_version=4)
nbformat.validate(nb)

n_code = 0
n_err = 0
n_png = 0
n_shows = 0
for c in nb.cells:
    if c.cell_type != "code":
        continue
    n_code += 1
    n_shows += c.source.count("plt.show()") + c.source.count("fig.show()")
    for o in c.get("outputs", []):
        if o.get("output_type") == "error":
            n_err += 1
            print("ERROR in cell %s:" % c.get("id"))
            print(o.get("ename"), o.get("evalue"))
        if o.get("output_type") == "display_data" and o.get("data", {}).get("image/png"):
            n_png += 1
        if o.get("output_type") == "execute_result" and o.get("data", {}).get("image/png"):
            n_png += 1

print("cells: %d total, %d code" % (len(nb.cells), n_code))
print("error outputs: %d" % n_err)
print("embedded PNGs: %d (plt.show() calls: %d)" % (n_png, n_shows))

ok = n_err == 0 and n_png == n_shows
print("VERIFICATION: %s" % ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
