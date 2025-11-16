import itertools
from pysat.card import CardEnc




def make_base_constraints(N:int, K:int, D:int, SymmetryBreaking = False, Decompose_Xors = False):

    # N, K, D refer to n,k,d parameters of the linear binary block code.
    # N is dimension = size of the codeword
    # K is the number of code wectors in the generator matrix
    # D is the minimum distance of the code, i.e. the minimal Hamming weight
    # We aim to minimize the first component of the weight spectrum, i.e. the number of codewords
    # that have the weight = D
    
    # We use the approach based on the generator matrix
    # It is of size n * k
    
    # XOR constraints: (left_part, right_part): left_part = XOR(right_part)    
    XOR_CONSTRAINTS = list()
    # Cardinality constraints: (combination, left_part, right_part): sum(left_part) >= right_part
    CARDINALITY_CONSTRAINTS = list()
    # Order constraints: (left_part, right_part): left_part > right_part
    ORDER_CONSTRAINTS = list()
    # Unit constraints: [a,b]: VAL(a) = b
    UNIT_CONSTRAINTS = dict()

    varcount = 1
    
    # Introduce the variables for the generator matrix
    # Note, that we imply that the generator matrix contains an identity submatrix, thus we ignore the first K columns
    # and thus only introduce N-K variables for each column
    GeneratorMatrix = list()
    for i in range(K):
        X = list()
        for j in range(N-K):
            X.append(varcount)
            varcount+=1
        GeneratorMatrix.append(X)
    
    
    # If we use Symmetry breaking, we i) fix some variables in the first row of the matrix, and 
    # ii) say that row 1 is > than row 2, etc. lexicographically, starting from the end of the row to the start of the row

    if SymmetryBreaking == True:        
        for u in range(N - K - D):
            UNIT_CONSTRAINTS[GeneratorMatrix[0][u]] = 0        
        for u in range(N-K - (D-1),N-K):
            UNIT_CONSTRAINTS[GeneratorMatrix[0][u]] = 1

        for u in range(N-K-1,0,-1):
            Column_1 = [GeneratorMatrix[j][u] for j in range(K)]
            Column_2 = [GeneratorMatrix[j][u-1] for j in range(K)]
            ORDER_CONSTRAINTS.append([Column_1,Column_2])


   
    # Next, we work with all possible combinations of K vectors of the generator matrix
    R = [u for u in range(K)]    
   
    # We process base vectors separately for simplicity
    for u in range(K):
        CARDINALITY_CONSTRAINTS.append([[u], GeneratorMatrix[u], D-1])

    # Auxiliary dictionary for the case when we decompose xors so that there are no xors with the number of inputs >=2
    XOR_DICT = dict()
    for u in range(K):
        XOR_DICT[tuple([u])] = GeneratorMatrix[u]

    for u in range(2,K+1):
        for COMB in itertools.combinations(R,u):
            # The value in the right part of the cardinality constraint is reduced by the number of vectors in a 
            # linear combination due to the implied identity matrix:
            # if we xor h (different) rows of an identity matrix - we will have a vector of Hamming weight h
            if D-len(COMB) > 0:
                # Introduce new variables that will encode this particular linear combination of base vectors                        
                XC = list()
                for j in range(N-K):
                    XC.append(varcount)
                    varcount+=1
                # Add XOR constraint connecting each component of this vector with components of base vectors involved
                if Decompose_Xors == False:
                    for j in range(N-K):                
                        XOR_CONSTRAINTS.append([XC[j],[GeneratorMatrix[i][j] for i in COMB]])
                
                else:              
                    assert(tuple(COMB[:-1]) in XOR_DICT)                     
                    RP_1 = XOR_DICT[tuple(COMB[:-1])]
                    RP_2 = GeneratorMatrix[COMB[-1]]
                    for j in range(N-K):            
                        XOR_CONSTRAINTS.append([XC[j],[RP_1[j], RP_2[j]]])

                    XOR_DICT[tuple(COMB)] = XC



                # Add cardinality constraint. We additionaly store the linear combination to be used later.
                
                CARDINALITY_CONSTRAINTS.append([COMB, XC, D-len(COMB)])
    
    return GeneratorMatrix, XOR_CONSTRAINTS, CARDINALITY_CONSTRAINTS, UNIT_CONSTRAINTS, ORDER_CONSTRAINTS, varcount

def optimization_formulation(N:int, K:int, D:int, SymmetryBreaking = False, Decompose_Xors = False, BLOCKS = []):
    # The concept of BLOCKS is as follows. We want to somehow restrict the 
    GM, XORS, CARDS, UNITS, ORDER, varcount = make_base_constraints(N,K,D,SymmetryBreaking, Decompose_Xors)
    USED_BLOCKS = BLOCKS
    if BLOCKS == []:
        USED_BLOCKS = [u for u in range(1,K+1)]
    
    Lvars = dict()
    
    R = [u for u in range(K)]    
    # Introduce variables only for the combinations corresponding to the considered blocks
    for u in USED_BLOCKS:
        for COMB in itertools.combinations(R,u):
            Lvars[tuple(COMB)] = varcount
            varcount += 1
    
    for u in range(len(CARDS)):
        # Update the cardinality constraint left part to account for newly introduced auxiliary vars
        if tuple(CARDS[u][0]) in Lvars:
            CARDS[u][1].append(Lvars[tuple(CARDS[u][0])])
        # Update the right part for all constraints
        CARDS[u][2]+=1
        
    OPTVARS = [Lvars[u] for u in Lvars]

    return GM, XORS, CARDS, UNITS, ORDER, OPTVARS, varcount
