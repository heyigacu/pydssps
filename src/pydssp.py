import pydssp
import os

def pydssp_from_pdb(pdb):
    os.system(f'pydssp {pdb} -o output.result > /dev/null 2>&1')
    with open('output.result', 'r') as f:
        line = f.readlines()[0].split()[0]
    os.remove('output.result')
    return line.strip()