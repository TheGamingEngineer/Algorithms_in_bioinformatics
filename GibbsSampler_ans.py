import numpy as np
import matplotlib.pyplot as plt
from time import time

from argparse import ArgumentParser

parser = ArgumentParser(description="GibbsSampler")

parser.add_argument("-b", action="store", dest="beta", type=float, default=50.0,help="Weight on pseudo count (default: 50.0)")
parser.add_argument("-w", action="store_true", dest="sequence_weighting", help="Use Sequence weighting")
parser.add_argument("-f", action="store", dest="peptides_file", type=str, help="File with peptides")
parser.add_argument("-i", action="store", dest="iters_per_point", type=int,  default=6,help="Number of iteration per data point")
parser.add_argument("-s", action="store", dest="seed", type=int, default=1, help="Random number seed")
parser.add_argument("-Ts", action="store", dest="T_i", type=float, default=1.0, help="Start Temp")
parser.add_argument("-Te", action="store", dest="T_f", type=float, default=0.0001, help="End Temp")
parser.add_argument("-nT", action="store", dest="T_steps", type=int, default=5, help="Number of T steps")
args = parser.parse_args()
beta = args.beta
sequence_weighting = args.sequence_weighting
peptides_file = args.peptides_file
iters_per_point = args.iters_per_point
seed = args.seed
T_i = args.T_i
T_f = args.T_f
T_steps = args.T_steps


# ## DEFINE THE PATH TO YOUR COURSE DIRECTORY

# In[37]:

data_dir = "/Users/mniel/Courses/Algorithms_in_Bioinf/ipython/data/"

# ## Data imports

# ## Load data files

# In[39]:

alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)

bg_file = data_dir + "Matrices/bg.freq.fmt"

_bg = np.loadtxt(bg_file, dtype=float)

bg = {}
for i in range(0, len(alphabet)):
    bg[alphabet[i]] = _bg[i]

blosum_file = data_dir + "Matrices//blosum62.freq_rownorm"

_blosum62 = np.loadtxt(blosum_file, dtype=float).reshape((20, 20)).T

blosum62 = {}

for i, letter_1 in enumerate(alphabet):
    
    blosum62[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        
        blosum62[letter_1][letter_2] = _blosum62[i, j]


# ## Calculate log-odd matrix from a given peptide core alignment

# In[40]:

def initialize_matrix(core_len, alphabet):

    init_matrix = [0]*core_len

    for i in range(0, core_len):

        row = {}

        for letter in alphabet: 
            row[letter] = 0.0

        init_matrix[i] = row
        
    return init_matrix

def put_to_zero(matrix):
    
    for i in range(0, len(matrix)):

        for key in matrix[i].keys():
        
            matrix[i][key] = 0.0
    
    return matrix
            
def get_log_odds(peptides, alphabet, bg, scoring_scheme, core_len, c_matrix, f_matrix, g_matrix, p_matrix, w_matrix):

    # Amino Acid Count Matrix (c)

    c_matrix = put_to_zero(c_matrix)

    for position in range(0, core_len):

        # peptides has two elements; element[0] is the peptide sequence, element [1] is the core location
        for element in peptides:
            
            peptide = element[0]
            
            core_start = element[1]
            
            c_matrix[position][peptide[core_start+position]] += 1

            
    # Sequence Weighting
    weights = {}

    for element in peptides:

        peptide = element[0]
        core_start = element[1]
            
        # apply sequence weighting
        if sequence_weighting:

            w = 0.0
            neff = 0.0

            for position in range(0, core_len):

                r = 0

                for letter in alphabet:        

                    if c_matrix[position][letter] != 0:

                        r += 1

                s = c_matrix[position][peptide[core_start+position]]

                w += 1.0/(r * s)

                neff += r

            neff = neff / core_len

        # do not apply sequence weighting
        else:

            w = 1  

            neff = len(peptides)  


        weights[peptide] = w

    
    # Observed Frequencies Matrix (f)
    f_matrix = put_to_zero(f_matrix)

    for position in range(0, core_len):

        n = 0;

        for element in peptides:

            peptide = element[0]
            
            core_start = element[1]
              
            f_matrix[position][peptide[core_start+position]] += weights[peptide]

            n += weights[peptide]

        for letter in alphabet: 

            f_matrix[position][letter] = f_matrix[position][letter]/n

    
    # Pseudo Frequencies Matrix (g)
    g_matrix = put_to_zero(g_matrix)

    for position in range(0, core_len):

        for letter_1 in alphabet:
            
            for letter_2 in alphabet:

                 g_matrix[position][letter_1] += f_matrix[position][letter_2] * scoring_scheme[letter_1][letter_2]

                    
    # Combined Frequencies Matrix (p)

    alpha = neff - 1

    for position in range(0, core_len):

        for letter in alphabet:

            num = alpha*f_matrix[position][letter] + beta*g_matrix[position][letter]
            
            den = alpha + beta

            p_matrix[position][letter] = num / den
        

    # Log Odds Weight Matrix (w)
    for position in range(0, core_len):

        for letter in alphabet:

            if p_matrix[position][letter] != 0:
                
                w_matrix[position][letter] = np.log(p_matrix[position][letter]/bg[letter])/np.log(2)
    
    # Calculate the overall score of the peptides to the LO matrix
    _sum = 0
    for position in range(0, core_len):
        for letter in alphabet:
            _sum += f_matrix[position][letter] * w_matrix[position][letter]
    
    return w_matrix, _sum, p_matrix


# ## Score peptides to matrix

# In[41]:

def score_peptide(peptide, core_start, core_len, matrix):
    acum = 0
    for i in range(0, core_len):
        acum += matrix[i][peptide[i+core_start]]
    return acum


# ## Read peptides

# In[42]:

def load_peptide_data():
    
    # Remove peptides shorter than core_len
    raw_peptides = np.loadtxt(peptides_file, dtype=str).tolist()

    # only keep peptides with length equal to or longer than core_len
    peptides = []
    for i in range(0, len(raw_peptides)):
        if len(raw_peptides[i]) >= core_len:          
            peptides.append(raw_peptides[i])
        else:
            print "Peptide length too short discard", raw_peptides[i]
    
    peptides = sorted(peptides, key=len)
    min_pep_len = len(peptides[0])
    max_pep_len = len(peptides[-1])

    # random core start
    np.random.shuffle(peptides)
    cores_start = [0]*len(peptides)

    for i in range(0, len(cores_start)):

        if len(peptides[i]) != core_len:

            min_core_start = 0
            max_core_start = len(peptides[i]) - core_len

            cores_start[i] = np.random.randint(min_core_start, max_core_start)

    peptides = list(zip(peptides, cores_start))

    return peptides, min_pep_len, core_len


# ## Print out w-matrix in Psi_blast format

# In[43]:

def to_psi_blast(matrix):

    header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*header)) 

    letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    for i, row in enumerate(matrix):

        scores = []

        scores.append(str(i+1) + " A")

        for letter in letter_order:

            score = row[letter]

            scores.append(round(score, 4))

        print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*scores)) 


# ## Main loop

# In[ ]:

np.random.seed( seed )

core_len = 9

# get peptides
peptides, min_pep_len, core_len = load_peptide_data()

T_delta = (T_f - T_i) / T_steps
T = np.linspace(T_i,T_f,T_steps )
  

# Define number of iterations per temperature step
iters = len(peptides)*iters_per_point

print "# beta:", beta
print "# Sequence weighting:", sequence_weighting
print "# iters_per_point:", iters_per_point
print "# Seed:", seed
print "# Temperature:", T

# initialize matrices
c_matrix = initialize_matrix(core_len, alphabet)
f_matrix = initialize_matrix(core_len, alphabet)
g_matrix = initialize_matrix(core_len, alphabet)
p_matrix = initialize_matrix(core_len, alphabet)
w_matrix = initialize_matrix(core_len, alphabet)

# get initial log-odds matrix
log_odds_matrix, peptide_scores, _ = get_log_odds(peptides, alphabet, bg, blosum62, core_len, c_matrix, f_matrix, g_matrix, p_matrix, w_matrix)

# initial kld score
print( "Initial KLD score: " + str(peptide_scores))
kld = []
kld.append( peptide_scores )

# other stuff
t0 = time()
debug = False
#debug = True

for t in T:

    for i in range(0, iters):
            
        # extract peptide
        rand_index = np.random.randint(0,len(peptides)-1)
        peptide = peptides[rand_index][0]
        core_start_original = peptides[rand_index][1]

        if len(peptide) != core_len:
                
            max_core_start = len(peptide) - core_len

            core_start_shifted = np.random.randint(0, max_core_start)
                     
            # remove peptide from list
            peptides.remove(peptides[rand_index])

            # get base log_odds
            log_odds_matrix, peptide_scores, p_matrix = get_log_odds(peptides, alphabet, bg, blosum62, core_len, c_matrix, f_matrix, g_matrix, p_matrix, w_matrix)
                
            # score peptide against log_odds
            e_original = score_peptide(peptide, core_start_original, core_len, log_odds_matrix)
            if debug: print("Energy before shifting: " + str(e_original))

            # score shifted peptide against log_odds
	    # e_shift = XX
            e_shift = score_peptide(peptide, core_start_shifted, core_len, log_odds_matrix)
            if debug: print("Energy after shifting: " + str(e_shift))

            # energy differential
	    # de_shift = XX
            de = e_shift - e_original
    
            # probability of accepting move
            if ( de > 0):
                p = 1
            else:
		# e_shift = XX
                p = np.exp(de/t)
    
            # throw coin
            coin = np.random.uniform(0.0, 1.0, 1)[0]

	    # if coin XX p:
            if coin < p:
                if debug: print("RNG < P, Move accepted")
                peptides.append((peptide, core_start_shifted))
                kld.append(peptide_scores)

            else:
                if debug: print("RNG >= P, Move rejected")
                peptides.append((peptide, core_start_original))      
                      
        else:
            if debug: print("Can't shift peptide, it is a " + str(core_len) + "mer")
              
    print( "KLD score t: " + str(t) + " KLD: " + str(peptide_scores))
    
t1 = time()

print("Time elapsed (m):", (t1-t0)/60)

to_psi_blast(log_odds_matrix)
