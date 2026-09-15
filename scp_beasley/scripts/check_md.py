
import nbformat
nb=nbformat.read('/mnt/d/exactTest/column-generation-solvers/scp_beasley/01_direct.ipynb', as_version=4)
for cell in nb.cells[:2]:
    print(repr(cell.source[:400]))
