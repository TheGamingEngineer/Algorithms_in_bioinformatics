import numpy as np

from argparse import ArgumentParser

parser = ArgumentParser(description="Pep2mat")

parser.add_argument("-q", action="store", dest="query_file", help="File with query sequence")
parser.add_argument("-db", action="store", dest="db_file", help="File with database sequence")
parser.add_argument("-go", action="store", dest="gap_open", type=float, default=-11.0, help="Value of gap open (-11.0)")
parser.add_argument("-ge", action="store", dest="gap_extension", type=float, default=-1.0, help="Value of gap extension (-1.0)")

args = parser.parse_args()

gap_open = args.gap_open
gap_extension = args.gap_extension
query_file = args.query_file
database_file = args.db_file


data_dir = "/Users/mniel/Courses/Algorithms_in_Bioinf/ipython/data/"

# ### Alphabet

# In[7]:

alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)


# ### Blosum Matrix

blosum_file = data_dir + "Matrices/BLOSUM50"

_blosum50 = np.loadtxt(blosum_file, dtype=float).reshape((24, -1)).T
blosum50 = {}

for i, letter_1 in enumerate(alphabet):
    
    blosum50[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        
        blosum50[letter_1][letter_2] = _blosum50[i, j]


# ## Alignment
# 
# This functions returns, apart from the final Alignment Matrix, all the intermedite Matrices (for plotting purposes).

# ### Alignment Matrix

# In[9]:

def smith_waterman_alignment(query="VLLP", database="VLILP", scoring_scheme={}, gap_open=-5, gap_extension=-1):
    
    # Matrix dimensions
    M = len(query)
    N = len(database)
    
    # E matrix (for backtracking)
    E_matrix = np.zeros((M+1, N+1), dtype=object)
    
    # D matrix (alignment matrix)
    D_matrix = np.zeros((M+1, N+1), np.int)

    # Initialize matrices (do we neede this?)
    for i in range(M, 0, -1):
        D_matrix[i-1, N] = 0
        E_matrix[i-1, N] = 0

    for j in range(N, 0, -1):
        D_matrix[M, j-1] = 0
        E_matrix[M, j-1] = 0
    
    
    D_matrix_max_score, D_matrix_i_max, D_matrix_j_max = -9, -9, -9
    for i in range(M-1, -1, -1): 
        for j in range(N-1, -1, -1):
                
            # digonal score
	    # diagonal_score = XX
            diagonal_score = D_matrix[i+1, j+1] + scoring_scheme[query[i]][database[j]]    
            
            # horizontal score
	    # max_horizontal_score = XX
            max_horizontal_score = D_matrix[i, j+1] + gap_open
            k_max_horizontal_score = 1
          
            for k in range(j+2, N):
                
		# score = XX
                score = D_matrix[i, k] + gap_open + (k-j-1)*gap_extension
                
                if score > max_horizontal_score: 
                    max_horizontal_score = score 
                    k_max_horizontal_score = k - j            
            
            
            # vertical score
            # max_vertical_score = XX	
            max_vertical_score = D_matrix[i+1, j] + gap_open
            k_max_vertical_score = 1
            
            for k in range(i+2, M):
    
		# score = XX
                score = D_matrix[k, j] + gap_open + (k-i-1)*gap_extension
               
                if score > max_vertical_score: 
                    max_vertical_score = score 
                    k_max_vertical_score = k - i
                  
                  
            ####################
            # E_matrix entries #
            ####################
            # E[i,j] = 0, negative number
            # E[i,j] = 1, match
            # E[i,j] = 2, gap opening in database
            # E[i,j] = 3, gap extension in database
            # E[i,j] = 4, gap opening in query
            # E[i,j] = 5, gap extension in query
            
            if diagonal_score >= max_vertical_score and diagonal_score >= max_horizontal_score:
                max_score = diagonal_score
                direction = "diagonal"
            elif max_horizontal_score > max_vertical_score:
                max_score = max_horizontal_score
                direction = "horizontal"
            else:
                max_score = max_vertical_score
                direction = "vertical"
                
            if max_score <= 0:
                max_score = 0
                direction = "none"

            # diagonal direction case
            if direction == "diagonal":
                E_matrix[i,j] = 1
                
            # vertical direction case
            elif direction == "vertical":

                # if k only moved one position, it means gap opening
                if k_max_vertical_score == 1: 
		    # E_matrix[i,j] = XX
                    E_matrix[i,j] = 2

                # else it is a gap extension
                else: 
		    # E_matrix[i,j] = XX
                    E_matrix[i,j] = 3
                        
            # horizontal direction case
            elif direction == "horizontal":

                # if k only moved one position, it means gap opening
                if k_max_horizontal_score == 1: 
                    E_matrix[i,j] = 4

                # else it is a gap extension
                else: 
                    E_matrix[i,j] = 5

            else:
                # max_score is negative, put E to zero
                E_matrix[i,j] = 0
                 
            # store max score
            D_matrix[i, j] = max_score
            
            # append partial alignment matrix to list
            #D_matrix_list.append(np.copy(D_matrix))
           
            # fetch global max score
            if max_score > D_matrix_max_score:
                D_matrix_max_score = max_score
                D_matrix_i_max = i
                D_matrix_j_max = j
            
    return D_matrix, E_matrix, D_matrix_i_max, D_matrix_j_max, D_matrix_max_score


# ## Alignment Matrix Traceback

def smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, query="VLLP", database="VLILP", gap_open=-5, gap_extension=-1):
    
    M = len(query)
    N = len(database)
    
    aligned_query = []
    aligned_database = []
    positions = []
    matches = 0
    
    # start from max_i, max_j
    i, j = i_max, j_max
    while i < M and j < N :

        positions.append([i,j])
        
        # E[i,j] = 0, stop back tracking
        if E_matrix[i, j] == 0:
            break
        
        # E[i,j] = 1, match
        if E_matrix[i, j] == 1:
            aligned_query.append(query[i])
            aligned_database.append(database[j])
            if (query[i] == database[j]):
                matches += 1
            i += 1
            j += 1
        
        
        # E[i,j] = 2, gap opening in database
        if E_matrix[i, j] == 2:
            aligned_database.append("-")
            aligned_query.append(query[i])
            i += 1

            
        # E[i,j] = 3, gap extension in database
        if E_matrix[i, j] == 3:
            
            count = i + 2
            score = D_matrix[count, j] + gap_open + gap_extension

            # Find length of gap (check if score == D_matrix[i, j])
            while((score - D_matrix[i, j])*(score - D_matrix[i, j]) >= 0.00001): 
                count += 1
		# score = XX
                score = D_matrix[count, j] + gap_open + (count-i-1)*gap_extension

            for k in range(i, count):
                aligned_database.append("-")
                aligned_query.append(query[i])
                i += 1
             
          
        # E[i,j] = 4, gap opening in query
        if E_matrix[i, j] == 4:
            aligned_query.append("-")
            aligned_database.append(database[j])
            j += 1
        
        
        # E[i,j] = 5, gap extension in query
        if E_matrix[i, j] == 5:
            
            count = j + 2
            score = D_matrix[i, count] + gap_open + gap_extension
            
            # Find length of gap (check if score == D_matrix[i, j])
            while((score - D_matrix[i, j])*(score - D_matrix[i, j]) >= 0.0001): 
                count += 1
                score = D_matrix[i, count] + gap_open + (count-j-1)*gap_extension

            for k in range(j, count):
                aligned_query.append("-")
                aligned_database.append(database[j])
                j += 1
                

    return aligned_query, aligned_database, matches


# ## Multiple alignments. Run an alignment of a sequence against a sequence database
# 

query_list = np.loadtxt(query_file, dtype=str).reshape(-1,2)
database_list = np.loadtxt(database_file, dtype=str).reshape(-1,2)

scoring_scheme = blosum50


# ### Align query against database. Might take a while. Go get some coffee 

from time import time

# this returns current timestamp in seconds
t0 = time()

for query in query_list:
    
    query_protein = query[0]
    query_sequence = query[1]
    
    for database in database_list:
    
        database_protein = database[0]
        database_sequence = database[1]
    
        D_matrix, E_matrix, i_max, j_max, max_score = smith_waterman_alignment(query_sequence, database_sequence, scoring_scheme, gap_open, gap_extension)
        aligned_query, aligned_database, matches = smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, query_sequence, database_sequence, gap_open, gap_extension)
        
        print "ALN", query_protein, len(query_sequence), database_protein, len(database_sequence), len(aligned_query), matches, max_score
        print "QAL", i_max, ''.join(aligned_query)
        print "DAL", j_max,''.join(aligned_database)
        print
        
t1 = time()

print( "Time (m):", (t1-t0)/60)


