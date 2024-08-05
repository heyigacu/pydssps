
import os 
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from Bio import SeqIO
PARENT_DIR = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
# model
# =============================================================================
class ResidueEmbedding(nn.Embedding):
    def __init__(self, vocab_size=21, embed_size=128, padding_idx=None):
        super().__init__(vocab_size, embed_size, padding_idx=padding_idx)

class GRUnet(nn.Module):
    def __init__(self,lstm_hdim=1024, embed_size=128, num_layers=3,bidirectional=True,lstm=False,outsize=3):
        super().__init__()
        """
            This version of the model has all the bells & whistles (e.g. 
            dropconnect) ripped out so its slimmed down for inference
            
        """
        self.lstm_hdim = lstm_hdim
        self.embed=ResidueEmbedding(vocab_size=22, embed_size=embed_size, padding_idx=21)
        self.lstm = nn.GRU(128, 1024, num_layers=3, bidirectional=True, batch_first=True,dropout=0.0)
        self.outlayer = nn.Linear(lstm_hdim*2, outsize)
        self.finalact=F.log_softmax
    def forward(self, x):
        """
            Assumes a batch size of one currently but can be changed
        """
        x=self.embed(x)
        x, _ = self.lstm(x)
        x=self.outlayer(x)
        x=self.finalact(x,dim=-1)
        return x.squeeze()     
       
class S4PRED(nn.Module):
    def __init__(self):
        super().__init__()
        """
            This loads the ensemble of models in a lazy way but its clear and 
            leaves the weight loading out of the run_model script. 
        """                        
        # Manually listing for clarity and hot swapping in future
        self.model_1=GRUnet()
        self.model_2=GRUnet()
        self.model_3=GRUnet()
        self.model_4=GRUnet()
        self.model_5=GRUnet()
    def forward(self, x):
        y_1=self.model_1(x)
        y_2=self.model_2(x)
        y_3=self.model_3(x)
        y_4=self.model_4(x)
        y_5=self.model_5(x)
        y_out=y_1*0.2+y_2*0.2+y_3*0.2+y_4*0.2+y_5*0.2
        return y_out        

def laod_s4pred():
    s4pred=S4PRED()
    s4pred.eval()
    s4pred.requires_grad=False
    weight_files=['/s4pred_weights/weights_1.pt',
                '/s4pred_weights/weights_2.pt',
                '/s4pred_weights/weights_3.pt',
                '/s4pred_weights/weights_4.pt',
                '/s4pred_weights/weights_5.pt']
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    s4pred.model_1.load_state_dict(torch.load(scriptdir + weight_files[0], map_location=lambda storage, loc: storage))
    s4pred.model_2.load_state_dict(torch.load(scriptdir + weight_files[1], map_location=lambda storage, loc: storage))
    s4pred.model_3.load_state_dict(torch.load(scriptdir + weight_files[2], map_location=lambda storage, loc: storage))
    s4pred.model_4.load_state_dict(torch.load(scriptdir + weight_files[3], map_location=lambda storage, loc: storage))
    s4pred.model_5.load_state_dict(torch.load(scriptdir + weight_files[4], map_location=lambda storage, loc: storage))
    return s4pred

# =============================================================================
# model
# =============================================================================
"""
C: coil 
H: helix 
E: strand 
"""
ind2char={0:"_", 1:"H", 2:"E"}


def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def format_horiz(sequence, ss, ss_conf):
    ''' Formats output for the PSIPRED HFORMAT .horiz files. 
        Care must be taken as there is a fixed column width of 60 char
    '''    
    lines=['# PSIPRED HFORMAT (S4PRED V1.2.4)']
    sub_seqs = list(chunkstring(sequence,60))
    sub_ss   = list(chunkstring("".join([ind2char[s.item()] for s in ss]),60))
    
    num_len  =  int(np.floor(len(sequence)/10))
    num_seq  = ''.join(f'{str((i+1)*10):>10}' for i in range(num_len+1))
    num_seq  = list(chunkstring(num_seq,60))
        
    # get confidences then floor them and convert to string 
    conf_idxs = ss_conf.argmax(-1)
    confs = ss_conf[np.arange(len(conf_idxs)),conf_idxs[:]]
    confs = "".join([str(x) for x in np.floor(confs*10).astype(np.int32)])
    confs = list(chunkstring(confs,60))

    for idx, subsq in enumerate(sub_seqs):
        lines.append(f'\nConf: {confs[idx]}')
        lines.append(f'Pred: {sub_ss[idx]}')
        lines.append(f'  AA: {subsq}')
        lines.append(f'      {num_seq[idx]}\n')
    return lines



def seq2int(sequence):
    aanumdict = {'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 
             'I':9, 'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16,
             'W':17, 'Y':18, 'V':19}
    ints = [aanumdict.get(res, 20) for res in sequence]
    return ints


def s4pred_forward_sequence(sequence):
    model = laod_s4pred()
    data = seq2int(sequence)
    with torch.no_grad():
        ss_conf=model(torch.tensor([data]))
        ss=ss_conf.argmax(-1).numpy()
        ss_conf=ss_conf.exp()
        tsum=ss_conf.sum(-1)
        tsum=torch.cat((tsum.unsqueeze(1),tsum.unsqueeze(1),tsum.unsqueeze(1)),1)
        ss_conf/=tsum
        ss_conf=ss_conf.numpy()
    return sequence, ss, ss_conf

def generate_horiz(path, sequence, ss, ss_conf):
    lines = format_horiz(sequence, ss, ss_conf)
    with open(path, 'w') as f:
        for line in lines:
            f.write(line+'\n')


def s4pred_dssp(sequence):
    sequence, ss, ss_conf = s4pred_forward_sequence(sequence)
    seq = ['_']*len(sequence)
    for i,ind in enumerate(ss): 
        seq[i] = ind2char[ind]
    return ''.join(seq)