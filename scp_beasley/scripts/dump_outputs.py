
import nbformat, glob, os
for path in sorted(glob.glob('/mnt/d/exactTest/column-generation-solvers/scp_beasley/*.ipynb')):
    nb=nbformat.read(path, as_version=4)
    print('====', os.path.basename(path))
    for idx,cell in enumerate(nb.cells):
        if cell.cell_type!='code': continue
        for out in cell.get('outputs',[]):
            if out.output_type=='stream':
                print(''.join(out.get('text',[])))
