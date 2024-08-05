
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.Polypeptide import is_aa

from Bio.PDB.StructureBuilder import StructureBuilder
import biotite.structure.io as strucio
import os


def dssp_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    dssp = DSSP(structure[0], pdb_file)
    sequence = ""
    for chain in structure[0]:
        for residue in chain:
            if (chain.id, residue.id) in dssp:
                dssp_key = (chain.id, residue.id)
                dssp_data = dssp[dssp_key]
                if dssp_data[2] in ['H', 'G', 'I']:
                    sequence += 'H'
                elif dssp_data[2] in ['E', 'B']:
                    sequence += 'E'
                else:
                    sequence += '_'
    return sequence