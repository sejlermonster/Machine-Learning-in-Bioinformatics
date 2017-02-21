from Bio import SeqIO
import math
import string

observables = {'A':0, 'C':1, 'E':2, 'D':3, 'G':4, 'F':5, 'I':6, 'H':7, 'K':8, 'M':9, 'L':10, 'N':11, 'Q':12, 'P':13, 'S':14, 'R':15, 'T':16, 'W':17, 'V':18, 'Y':19 }

# Our initial propbabilities for the three different states
init_probs = [0.496359, 0.188095, 0.315546]
#Our three different states
states = {'i':0, 'M':1, 'o':2}
num_of_states = 3
# An array that converts index to a state
index_to_states = {0:'i', 1: 'M', 2: 'o'}

# Transition probabiliy. This describes the probabaity of changing state given a state
# So in general we can see that if we observe a state there is a high probability we will stay in that state.
trans_probs = [[0.990971, 0.009029, 0.000000],
               [0.023083, 0.953090, 0.023827],
               [0.000000, 0.013759, 0.986241]]

# Emition probabilities describes how our observations relates to our states For the different observations what is the probability of the different states.
#So we have 20 coloumns - one for each observation. And we have three rows on for each state.
emit_probs = [[0.043601, 0.011814, 0.053446, 0.065541, 0.049508, 0.049789, 0.054571, 0.024191, 0.055977, 0.035162, 0.103235, 0.045007, 0.029536, 0.048101, 0.075105, 0.059634, 0.068354, 0.016315, 0.067792, 0.043319],
              [0.102010, 0.019360, 0.009680, 0.011914, 0.033507, 0.103500, 0.118392, 0.003723, 0.000745, 0.039464, 0.138496, 0.014147, 0.011914, 0.026806, 0.067014, 0.012658, 0.073716, 0.037230, 0.119136, 0.056590],
              [0.082374, 0.008415, 0.059345, 0.059345, 0.069973, 0.031001, 0.049159, 0.019043, 0.081045, 0.025244, 0.068202, 0.047830, 0.032772, 0.052259, 0.073959, 0.086802, 0.056244, 0.007086, 0.062445, 0.027458]]


#Return the log joint probab
def log_joint_prob(x, z):
    logp = math.log(init_probs[z[0]]) + math.log(emit_probs[z[0]][x[0]])
    for i in range(1, len(x)):
        logp = logp + math.log(trans_probs[z[i - 1]][z[i]]) + math.log(emit_probs[z[i]][x[i]])
    return logp

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
        # So it is our initial probabilities + the emission probability of our 0th observation
        # All three rows are calculated as s will count of 3 times while o will stay at 0
        # Because we are working with log probabilities we use addition instead of multiplication
        w.append(["-inf"] * num_of_states) 
        for st in range(num_of_states):
            w[0][st] = (log(init_probs[st]) + log(emit_probs[st][obs[0]]))

        #We then initialize the values of the rest of the array
        for o in range(1, len(obs)):
            w.append(["-inf"] * num_of_states) # Each iteration we append a new column with num_of_states as amount of rows
            for s in range(num_of_states):
                #This step is desribed on slide 11  of HMM implementation
                #We find the most likely on for each state
                # This is done by finding the highest probability of each previous step and adding(because log prop).
                #  This gives su the probabilities  for each state added with the transition probability to go to state s
                # We then one which is most likely to by getting the max value.
                highestVal = GetMaxValue([w[o-1][0] + log(trans_probs[0][s]),
                                          w[o-1][1] + log(trans_probs[1][s]),
                                          w[o-1][2] + log(trans_probs[2][s])])
                # We then take the highest val and add it to the emit_probs for the state s based on our observations
                w[o][s] = log(emit_probs[s][obs[o]]) + highestVal

        # We now backtrack through the table that we just created.
        # We create z which will contain the most likely sequence
        z = len(obs) * [None]
        z[len(z)-1] = GetMaxValueIndex(w[len(z)-1])
        
         #This step is desribed on slide 12 of HMM implementation
         #We backtrack pretty much the same way as we created the table
         # We find the max value of the previous
        for n in range(len(obs)-2, -1, -1):
            z[n] = GetMaxValueIndex([log(emit_probs[z[n+1]][obs[n+1]]) + w[n][0] + log(trans_probs[0][z[n+1]]), 
                                     log(emit_probs[z[n+1]][obs[n+1]]) + w[n][1] + log(trans_probs[1][z[n+1]]),
                                     log(emit_probs[z[n+1]][obs[n+1]]) + w[n][2] + log(trans_probs[2][z[n+1]])])
# We return the decoded(z) and we return the last column of w and the index of the last element in z
        return z, w[-1][z[-1]]
            
for recordX in SeqIO.parse ("xval.txt", "fasta"):
    Lx = list(recordX)

for recordZ in SeqIO.parse ("zval.txt", "fasta"):
    Lz = list(recordZ)

Lxx = [observables[c] for c in Lx]
Lzz = [states[c] for c in Lz]

# Compute Viterbi decoding and its log probabilty
z, logpz = viterbi_decoding(Lxx)

# Print out the decoded sequence
print "Viterbi path z:", string.join([index_to_states[c] for c in z]) 

# Check if we decoded correctly
if (string.join([index_to_states[c] for c in z]) == string.join([index_to_states[c] for c in Lzz])):
    print "We decoded correctly! hurra", 
else:
    print "We did not Decode correctly! shit!", 

#Print out log probabilty calcualted by viterbi and by log_joint_prob function
print "The log probability of the decoded z:", logpz
print "Joint log probability of (x,z):", log_joint_prob(Lxx, z)
