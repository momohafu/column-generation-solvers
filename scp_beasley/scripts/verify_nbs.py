
import nbformat, glob, os
for path in sorted(glob.glob('/mnt/d/exactTest/column-generation-solvers/scp_beasley/*.ipynb')):
    nb=nbformat.read(path, as_version=4)
    print('====', os.path.basename(path), 'cells', len(nb.cells))
    for idx,cell in enumerate(nb.cells):
        if cell.cell_type=='code':
            errs=[]
            for out in cell.get('outputs',[]):
                if out.output_type=='error':
                    errs.append(out.ename+': '+out.evalue)
                elif out.output_type=='stream':
                    text=''.join(out.get('text',[]))
                    # print only selected lines containing key terms
                    for line in text.splitlines():
                        if any(k in line for k in ['python','m,n','termination','objective','best_bound','solve_time','num_selected','CG iterations','RMP LP objective','integer repair','LP-integer gap','iter 1','iter 2','iter 3','Benders wall','Lagrangian best','ratio greedy','reduced-cost','MIP repair','dual gap','LBBD wall','selected_columns']):
                            print('  cell',idx,line[:200])
            if errs:
                print('  ERROR', errs)
