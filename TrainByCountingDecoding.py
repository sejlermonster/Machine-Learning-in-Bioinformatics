from Bio import SeqIO
import math
import string
from compare_tm_pred import Compare

observables = {'A':0, 'C':1, 'E':2, 'D':3, 'G':4, 'F':5, 'I':6, 'H':7, 'K':8, 'M':9, 'L':10, 'N':11, 'Q':12, 'P':13, 'S':14, 'R':15, 'T':16, 'W':17, 'V':18, 'Y':19 }
index_to_observables = { 0:'A', 1:'C', 2:'E', 3:'D', 4:'G', 5:'F', 6:'I', 7:'H', 8:'K', 9:'M', 10:'L', 11:'N', 12:'Q', 13:'P', 14:'S', 15:'R', 16:'T', 17:'W', 18:'V', 19:'Y'}


# Our initial propbabilities for the three different states
#init_probs = [0.496359, 0.188095, 0.315546]
init_probs = [0, 0, 0]
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

# Return the index of the highet value in list
def GetMaxValueIndex(l):
    maxVal, maxIndex = l[0], 0
    for i in range(1, len(l)):
        if l[i] > maxVal:
            maxVal = l[i]
            maxIndex = i
    return maxIndex

# Returns the highest value in the list
def GetMaxValue(l):
    maxVal = l[0]
    for i in range(1, len(l)):
        if l[i] > maxVal:
            maxVal = l[i]
    return maxVal

def log(x):
    if x == 0:
        return float("-inf")
    else:
        return math.log(x)

# Compute most likely sequence given observations. 
def viterbi_decoding(obs):
        w = []
        #First(0th) column is calcuated based on initial propabilities and emit probabilities
        # So it is our initial probabilities + the emission probabilities of our 0th observation
        # All three rows are calculated as s will count up 3 times
        # Because we are working with log probabilities we use addition instead of multiplication
        w.append(["-inf"] * num_of_states) 
        for st in range(num_of_states):
            w[0][st] = (log(init_probs[st]) + log(emit_probs[st][obs[0]]))

        #We then initialize the values of the rest of the array
        for o in range(1, len(obs)):
            w.append(["-inf"] * num_of_states) # Each iteration we append a new column with num_of_states as amount of rows
            for s in range(num_of_states):
                #This step is described on slide 11  of HMM implementation
                # We find the most likely on for each state
                # We look at the three previous probabilities and add the transition probability for going to state s
                # in that way we find the highest probability of going to state s.
                highestVal = GetMaxValue([w[o-1][0] + log(trans_probs[0][s]),
                                          w[o-1][1] + log(trans_probs[1][s]),
                                          w[o-1][2] + log(trans_probs[2][s])])
                # We then take the highest val and add it to the emit_probs for the state s based on our observations
                w[o][s] = log(emit_probs[s][obs[o]]) + highestVal

        # We now backtrack through the table that we just created.
        # We create z which will contain the most likely sequence
        z = len(obs) * [None]
        z[len(z)-1] = GetMaxValueIndex(w[len(z)-1])
        
         # This step is desribed on slide 12 of HMM implementation
         # We have the nth row of w and find the probability for each. We add the tranistion probability of 
         # transitioning to the z nth+1 state
         # We find the max value of the previous
        for n in range(len(obs)-2, -1, -1):
            z[n] = GetMaxValueIndex([log(emit_probs[z[n+1]][obs[n+1]]) + w[n][0] + log(trans_probs[0][z[n+1]]), 
                                     log(emit_probs[z[n+1]][obs[n+1]]) + w[n][1] + log(trans_probs[1][z[n+1]]),
                                     log(emit_probs[z[n+1]][obs[n+1]]) + w[n][2] + log(trans_probs[2][z[n+1]])])
# We return the decoded(z) and we return the last column of w and the index of the last element in z
        return z, w[-1][z[-1]]
        
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
            #curr_val = curr_val + l
    
    return d

def count(x, z):
    #Counting for initial propability
    init_probs[z[0]] += 1
    for i in range(0, len(x)):
        #Counting emission probability
        emit_probs[z[i]][x[i]] += 1
        #Count transition probabilities
        if(i == len(x)-1): continue
        trans_probs[z[i]][z[i+1]] += 1 

def Normalize():
    #Divide counted initial propabilities
    TotalInitCount = 0.0
    for i in range(0, len(init_probs)):
        TotalInitCount += init_probs[i]
    for i in range(0, len(init_probs)):
        init_probs[i] = init_probs[i]/TotalInitCount
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
#fastaData.append(fasta("TrainingData/set160.9.labels.txt").values())

fastaData2 = []
fastaData2.append(fasta("TrainingData/set160.9.labels.txt"))

Lxx = []
Lzz = []
for i in range(len(fastaData)):
    for j in range(len(fastaData[i])):
        Lxx.append(list(observables[c] for c in fastaData[i][j][0]))
        Lzz.append(list(states[c] for c in fastaData[i][j][1]))
    
for i in range(len(Lxx)):
    count(Lxx[i], Lzz[i])
Normalize()

count = 0
fileCount = 0
resultName = "results/3state/result." + str(len(fastaData)-1)
f = open(resultName + ".txt", 'w')

for l in fastaData2:
    for key in l:
        z, prop = viterbi_decoding(list(observables[c] for c in l[key][0]))
        f.write(">" + key)
        f.write("\n")
        f.write("  ")
        f.write(l[key][0])
        f.write("\n")
        f.write('# ')
        f.write(string.join([index_to_states[c] for c in z]).replace(" ", "")) 
        f.write("\n")
        f.write("\n")
f.close()

Compare("TrainingData/set160.9.labels.txt", "results/3state/result." + str(len(fastaData)-1) + ".txt")
print "done"

#print(init_probs)
#print(trans_probs)
#print(emit_probs)




