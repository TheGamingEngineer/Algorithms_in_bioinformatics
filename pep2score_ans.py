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
                
            matrix = initialize_matrix(peptide_length, alphabet)
        
        else:
            i = int(sline[0])
            
            for j in range(2,len(sline)):
                matrix[i-1][alphabet[j-2]] = float(sline[j])
                
        nline+= 1
            
    return matrix


# ### Score peptide to mat

def score_peptide(peptide, matrix):
    acum = 0
    for i in range(0, len(peptide)):
        # acum += XX
        acum += matrix[i][peptide[i]]
    return acum


# ## Main

# Read evaluation data

evaluation = np.loadtxt(peptides_file, dtype=str).reshape(-1,2)
evaluation_peptides = evaluation[:, 0]
evaluation_targets = evaluation[:, 1].astype(float)

peptide_length = len(evaluation_peptides[0])

# Define which PSSM file to use (file save from pep2mat)

w_matrix = from_psi_blast(mat_file)

evaluation_predictions = []
for i in range(len(evaluation_peptides)):
    score = score_peptide(evaluation_peptides[i], w_matrix)	
    evaluation_predictions.append(score)
    print evaluation_peptides[i], score, evaluation_targets[i]
    
pcc = pearsonr(evaluation_targets, evaluation_predictions)
print("PCC: ", pcc[0])
