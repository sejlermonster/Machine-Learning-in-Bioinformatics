from Bio import SeqIO
import math

observables = {'A':0, 'C':1, 'E':2, 'D':3, 'G':4, 'F':5, 'I':6, 'H':7, 'K':8, 'M':9, 'L':10, 'N':11, 'Q':12, 'P':13, 'S':14, 'R':15, 'T':16, 'W':17, 'V':18, 'Y':19 }

init_probs = [0.496359, 0.188095, 0.315546]
states = {'i':0, 'M':1, 'o':2}

trans_probs = [[0.990971, 0.009029, 0.000000],
               [0.023083, 0.953090, 0.023827],
               [0.000000, 0.013759, 0.986241]]

emit_probs = [[0.043601, 0.011814, 0.053446, 0.065541, 0.049508, 0.049789, 0.054571, 0.024191, 0.055977, 0.035162, 0.103235, 0.045007, 0.029536, 0.048101, 0.075105, 0.059634, 0.068354, 0.016315, 0.067792, 0.043319],
              [0.102010, 0.019360, 0.009680, 0.011914, 0.033507, 0.103500, 0.118392, 0.003723, 0.000745, 0.039464, 0.138496, 0.014147, 0.011914, 0.026806, 0.067014, 0.012658, 0.073716, 0.037230, 0.119136, 0.056590],
              [0.082374, 0.008415, 0.059345, 0.059345, 0.069973, 0.031001, 0.049159, 0.019043, 0.081045, 0.025244, 0.068202, 0.047830, 0.032772, 0.052259, 0.073959, 0.086802, 0.056244, 0.007086, 0.062445, 0.027458]]

# Returns the log transformed joint probability of x and z
def log_joint_prob(x, z):
    logp = math.log(init_probs[z[0]]) + math.log(emit_probs[z[0]][x[0]])
    for i in range(1, len(x)):
        logp = logp + math.log(trans_probs[z[i - 1]][z[i]]) + math.log(emit_probs[z[i]][x[i]])
    return logp

for recordX in SeqIO.parse ("xval.txt", "fasta"):
    Lx = list(recordX)

for recordZ in SeqIO.parse ("zval.txt", "fasta"):
    Lz = list(recordZ)

Lxx = [observables[c] for c in Lx]
Lzz = [states[c] for c in Lz]

print(log_joint_prob(Lxx, Lzz))