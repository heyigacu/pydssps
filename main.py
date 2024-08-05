from src.chou_fas import ChouFasman
from src.s4pred import s4pred_dssp
from src.PSEA import PSEA_dssp_from_pdb
from src.dssp import dssp_from_pdb
from src.pydssp import pydssp_from_pdb

sequence = 'MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDN'
pdb = 'test.pdb'


print('Sequence:'+' '*(20-len('Sequence'))+sequence)
print('ChouFasman:'+' '*(20-len('ChouFasman'))+ChouFasman(sequence))
print('s4pred:'+' '*(20-len('s4pred'))+s4pred_dssp(sequence))
print('PSEA:'+' '*(20-len('PSEA'))+PSEA_dssp_from_pdb(pdb))
print('dssp:'+' '*(20-len('dssp'))+dssp_from_pdb(pdb))
print('pydssp:'+' '*(20-len('pydssp'))+pydssp_from_pdb(pdb))


