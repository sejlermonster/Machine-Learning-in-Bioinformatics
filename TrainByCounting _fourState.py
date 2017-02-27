from Bio import SeqIO
import math
import string

observables = {'A':0, 'C':1, 'E':2, 'D':3, 'G':4, 'F':5, 'I':6, 'H':7, 'K':8, 'M':9, 'L':10, 'N':11, 'Q':12, 'P':13, 'S':14, 'R':15, 'T':16, 'W':17, 'V':18, 'Y':19 }

# Our initial propbabilities for the three different states
init_props = [0, 0, 0, 0]
#Our three different states
states = {'i':0, 'M':1, 'o':2, 'm':3}
num_of_states = 4
# An array that converts index to a state
index_to_states = {0:'i', 1: 'M', 2: 'o', 3:'M'}

# Transition probabiliy. This describes the probabaity of changing state given a state
# So in general we can see that if we observe a state there is a high probability we will stay in that state.
#Also we can see that we cannot go from state i(inside) to o(outside) we can only go through state M(membrane)
trans_probs = [[0.0, 0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0, 0.0]]

# Emition probabilities describes how our observations relates to our states For the different observations what is the probability of the different states.
#So we have 20 coloumns - one for each observation. And we have three rows on for each state.
emit_probs = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]

def CheckNextObservationChange(lzz, index):
    currVal = lzz[index]
    for i in range(index, len(lzz)): 
        if lzz[i] != currVal:
            if lzz[i] == states['o']:
                return states['M']
            else:
                return states['m']

    for i in range(len(lzz), -1,-1):
        if lzz[i] != currVal:
            if lzz[i] == states['i']:
                return states['M']
            else:
                return states['m']

def StateTranslater(lzz):
    result = []
    for i in range(len(lzz)):
        if lzz[i] != states['M']:
            result.append(lzz[i])
        else:
            result.append(CheckNextObservationChange(lzz,i))
    return result

def fasta(f):
    curr_val = ""
    curr_val_arr = ["", ""] 
    """
    Reads the fasta file f and returns a dictionary with the sequence names as keys and the
    sequences as the corresponding values. Lines starting with ';' in the fasta file are
    considered comments and ignored.
    """
    d = {}
    curr_key = ""
    lines = [string.strip(l) for l in open(f).readlines() if (l[0] != ';')]
    for l in lines:
        if l == '': 
            
            continue
        if l[0] == '>':
            if curr_key != "": 
                d[curr_key] = curr_val_arr
            curr_key = l[1:]
            curr_val_arr = ["", ""] 
            curr_val = ""
        else:
            if(l[0] == "#"):
                curr_val_arr[1] = l.strip("# ")
            else:
                curr_val_arr[0] = l
    d[curr_key] = curr_val_arr
    
    return d

def count(x, z):
    for i in range(0, len(x)):
        #Counting for initial propability
        init_props[z[i]] += 1
        #Counting emission probability
        emit_probs[z[i]][x[i]] += 1
        #Count transition probabilities
        if(i == len(x)-1): continue
        trans_probs[z[i]][z[i+1]] += 1 

def Normalize():
    #Divide counted initial propabilities
    TotalInitCount = 0.0
    for i in range(0, len(init_props)):
        TotalInitCount += init_props[i]
    for i in range(0, len(init_props)):
        init_props[i] = init_props[i]/TotalInitCount
    # Count emit_probs and divide by total
    for i in range(num_of_states):
        TotalEmitCount = 0.0
        for o in range(0, len(emit_probs[i])):
            TotalEmitCount += emit_probs[i][o]
        for j in range(0, len(emit_probs[i])):
            emit_probs[i][j] = emit_probs[i][j]/TotalEmitCount
    #Count transition  and divide by total
    for i in range(num_of_states):
        TotalTransCount = 0.0
        for o in range(0, len(trans_probs[i])):
            TotalTransCount += trans_probs[i][o]
        for j in range(0, len(trans_probs[i])):
            trans_probs[i][j] = trans_probs[i][j]/TotalTransCount

fastaData = []
fastaData.append(fasta("TrainingData/set160.0.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.1.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.2.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.3.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.4.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.5.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.6.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.7.labels.txt").values())
fastaData.append(fasta("TrainingData/set160.8.labels.txt").values())

Lxx = []
Lzz = []
for i in range(len(fastaData)):
    for j in range(len(fastaData[i])):
        Lxx.append(list(observables[c] for c in fastaData[i][j][0]))
        Lzz.append(list(states[c] for c in fastaData[i][j][1]))
    Lzz[i] = StateTranslater(Lzz[i])
    
for i in range(len(Lxx)):
    count(Lxx[i], Lzz[i])
Normalize()


print(init_props)
print(trans_probs)
print(emit_probs)
