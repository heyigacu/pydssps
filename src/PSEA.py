
import biotite.structure as struc
import biotite.structure.io as strucio


dssp_to_abc = {"I" : "c",
               "S" : "c",
               "H" : "a",
               "E" : "b",
               "G" : "c",
               "B" : "b",
               "T" : "c",
               "C" : "c"}

abc_to_EH_ = {'a':'H','b':'E','c':'_'}

def PSEA_dssp_from_pdb(pdb_file):
    array = strucio.load_structure(pdb_file)
    sse = struc.annotate_sse(array)
    sse = [abc_to_EH_[ss] for ss in sse]
    return "".join(sse)