from Bio import SeqIO
import math
import string

observables = {'A':0, 'C':1, 'E':2, 'D':3, 'G':4, 'F':5, 'I':6, 'H':7, 'K':8, 'M':9, 'L':10, 'N':11, 'Q':12, 'P':13, 'S':14, 'R':15, 'T':16, 'W':17, 'V':18, 'Y':19 }


# Our initial propbabilities for the three different states
#init_probs = [0.496359, 0.188095, 0.315546]
init_props = [0, 0, 0]
#Our three different states
states = {'i':0, 'M':1, 'o':2}
num_of_states = 3
# An array that converts index to a state
index_to_states = {0:'i', 1: 'M', 2: 'o'}

# Transition probabiliy. This describes the probabaity of changing state given a state
# So in general we can see that if we observe a state there is a high probability we will stay in that state.
#Also we can see that we cannot go from state i(inside) to o(outside) we can only go through state M(membrane)
trans_probs = [[0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0]]

# Emition probabilities describes how our observations relates to our states For the different observations what is the probability of the different states.
#So we have 20 coloumns - one for each observation. And we have three rows on for each state.
emit_probs = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]


#We load in our data. 
for recordX in SeqIO.parse ("TrainingData/set160.0.labelsX.txt", "fasta"):
    Lx = list(recordX)

#We load in our z for the test data, if we run it on the test data, or else ucomment this
for recordZ in SeqIO.parse ("TrainingData/set160.0.labelsZ.txt", "fasta"):
    Lz = list(recordZ)


Lxx = [observables[c] for c in Lx]
Lzz = [states[c] for c in Lz]


for i in range(0, len(Lxx)):
    #Counting for initial propability
    init_props[Lzz[i]] += 1
    #Counting emission probability
    emit_probs[Lzz[i]][Lxx[i]] += 1

#Divide counted initial propabilities
TotalCount = 0.0
for i in range(0, len(init_props)):
    TotalCount += init_props[i]
for i in range(0, len(init_props)):
    init_props[i] = init_props[i]/TotalCount

print "Intitial propbabilities:", init_props

# Count emit_probs and divide by total
for i in range(num_of_states):
    TotalEmitCount = 0.0
    for o in range(0, len(emit_probs[i])):
        TotalEmitCount += emit_probs[i][o]
    for j in range(0, len(emit_probs[i])):
        emit_probs[i][j] = emit_probs[i][j]/TotalEmitCount

print "Emission propbabilities:", emit_probs

#Count transition probabilities
for i in range(len(Lzz)-1):
    trans_probs[Lzz[i]][Lzz[i+1]] += 1

#Normalize
for i in range(num_of_states):
    TotalTransCount = 0.0
    for o in range(0, len(trans_probs[i])):
        TotalTransCount += trans_probs[i][o]
    for j in range(0, len(trans_probs[i])):
        trans_probs[i][j] = trans_probs[i][j]/TotalTransCount

print(trans_probs)
