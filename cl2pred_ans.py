import numpy as np
from pprint import pprint
from scipy.stats import pearsonr

from argparse import ArgumentParser

parser = ArgumentParser(description="Pep2score")

parser.add_argument("-mat", action="store", dest="mat_file", type=str, help="File with PSSM")
parser.add_argument("-f", action="store", dest="peptides_file", type=str, help="File with peptides (pep target=)")

args = parser.parse_args()

mat_file = args.mat_file
peptides_file = args.peptides_file

core_len = 9

# ## Initialize Matrix

def initialize_matrix(peptide_length, alphabet):

    init_matrix = [0]*peptide_length

    for i in range(0, peptide_length):

        row = {}

        for letter in alphabet: 
            row[letter] = 0.0

        #fancy way:  row = dict( zip( alphabet, [0.0]*len(alphabet) ) )

        init_matrix[i] = row
        
    return init_matrix


# ### Load Matrix from PSI-BLAST format

def from_psi_blast(file_name):

    f = open(file_name, "r")
    
    nline = 0
    for line in f:
    
        sline = str.split( line )
        
        if nline == 0:
        # recover alphabet
            alphabet = [str]*len(sline)
            for i in range(0, len(sline)):
                alphabet[i] = sline[i]
                
            matrix = initialize_matrix(core_len, alphabet)
        
        else:
            i = int(sline[0])
            
            for j in range(2,len(sline)):
                matrix[i-1][alphabet[j-2]] = float(sline[j])
                
        nline+= 1
            
    return matrix


# ### Score peptide to mat

def score_peptide(peptide, matrix):

    max_score = -99
    core_p1 = -9
    corelen = core_len
    for i in range(0, len(peptide)-corelen+1):
        score = 0
        for j in range(0, corelen):
            score += matrix[j][peptide[i+j]]
        if ( score > max_score):
            max_score = score
            core_p1 = i

    return max_score, core_p1


# ## Main

# Read evaluation data

evaluation = np.loadtxt(peptides_file, dtype=str).reshape(-1,3)
evaluation_peptides = evaluation[:, 0]
evaluation_targets = evaluation[:, 1].astype(float)

# Define which PSSM file to use (file save from pep2mat)

w_matrix = from_psi_blast(mat_file)

evaluation_predictions = []
for i in range(len(evaluation_peptides)):
    peptide = evaluation_peptides[i]
    score, p1 = score_peptide(peptide, w_matrix)	
    print peptide, p1, peptide[p1:p1+core_len], score, evaluation_targets[i]
