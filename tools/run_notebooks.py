"""Execute every notebook in notebooks/ in place (outputs saved), failing on the first error.

    pixi run notebooks            # all
    python tools/run_notebooks.py notebooks/02_crispr_benchmark.ipynb
"""
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

paths = [Path(p) for p in sys.argv[1:]] or sorted(Path("notebooks").glob("*.ipynb"))
for p in paths:
    t0 = time.time()
    nb = nbformat.read(p, as_version=4)
    NotebookClient(nb, timeout=1800, kernel_name="python3", resources={"metadata": {"path": str(p.parent)}}).execute()
    nbformat.write(nb, p)
    print(f"{p}: ok ({time.time() - t0:.0f} s)", flush=True)
