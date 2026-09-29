"""
Headless execution + verification of suspicion63_implicit_midpoint_blowup_v2w.ipynb

Run: D:/Source/hermes-dir/.venv/Scripts/python.exe run_suspicion63_v2w.py
Paths are derived from this script's own directory; the kernel uses the venv python
that is running this script (sys.executable), so it works on any drive.
"""
import json, os, sys, tempfile
import nbformat
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
NB_PATH = os.path.join(HERE, "suspicion63_implicit_midpoint_blowup_v2w.ipynb")
VENV_PY = sys.executable

tmp = tempfile.mkdtemp(prefix="jupyter_v2w_")
kernels_dir = os.path.join(tmp, "kernels", "suspicion63-v2w-venv")
os.makedirs(kernels_dir, exist_ok=True)
with open(os.path.join(kernels_dir, "kernel.json"), "w") as fh:
    json.dump({"argv": [VENV_PY, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
               "display_name": "Python 3 (project venv)", "language": "python"}, fh)

os.environ["JUPYTER_DATA_DIR"] = tmp
os.environ["MPLBACKEND"] = "Agg"

nb = nbformat.read(NB_PATH, as_version=4)
client = NotebookClient(nb, timeout=1200, kernel_name="suspicion63-v2w-venv")
client.execute()
nbformat.write(nb, NB_PATH)
print("executed and re-saved %s" % NB_PATH)

nb = nbformat.read(NB_PATH, as_version=4)
nbformat.validate(nb)
n_code = n_err = n_png = n_shows = 0
for c in nb.cells:
    if c.cell_type != "code":
        continue
    n_code += 1
    n_shows += c.source.count("plt.show()") + c.source.count("fig.show()")
    for o in c.get("outputs", []):
        if o.get("output_type") == "error":
            n_err += 1
            print("ERROR in cell %s: %s %s" % (c.get("id"), o.get("ename"), o.get("evalue")))
        if o.get("output_type") == "display_data" and o.get("data", {}).get("image/png"):
            n_png += 1
        if o.get("output_type") == "execute_result" and o.get("data", {}).get("image/png"):
            n_png += 1
print("cells: %d total, %d code" % (len(nb.cells), n_code))
print("error outputs: %d" % n_err)
print("embedded PNGs: %d (plt.show() calls: %d)" % (n_png, n_shows))
ok = n_err == 0 and n_png == n_shows
print("VERIFICATION: %s" % ("PASS" if ok else "FAIL"))

# Print key numeric outputs for the record
print("\n===== KEY OUTPUTS =====")
for c in nb.cells:
    if c.cell_type != "code":
        continue
    for o in c.get("outputs", []):
        if o.get("output_type") in ("stream", "execute_result"):
            txt = o.get("text") or o.get("data", {}).get("text/plain")
            if txt:
                t = "".join(txt) if isinstance(txt, list) else txt
                if any(k in t for k in ("Tref", "V2W correct IM", "itmax=2000", "itmax=200 ", "itmax only",
                                        "first step with deviation", "period ", "UNRESOLVED", "Radau IIA",
                                        "Strang split", "YES", "NO")):
                    print(t, end="" if t.endswith("\n") else "\n")
sys.exit(0 if ok else 1)
