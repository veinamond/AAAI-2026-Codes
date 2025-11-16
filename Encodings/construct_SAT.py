from constraint_encoding import *
import sys

def simple_and_to_cnf(lp,rp,clauses):
    lcl = list()
    for u in rp:
        lcl.append(-u)
        smcl = list()
        smcl.append(u)
        smcl.append(-lp)
        clauses.append(smcl)
    lcl.append(lp)
    clauses.append(lcl)

def simple_or_to_cnf(lp,rp,clauses):
    lcl = list()
    for u in rp:
        lcl.append(u)
        smcl = list()
        smcl.append(-u)
        smcl.append(lp)
        clauses.append(smcl)
    lcl.append(-lp)
    clauses.append(lcl)

def simple_xor_to_cnf(lp, rp, clauses):
    num = pow(2,len(rp)+1)
    vars = [u for u in rp]    
    vars.append(lp)
    for u in range(num):
        s = ('{0:0'+str(len(vars))+'b}').format(u)
        num1 = 0
        for v in s:
            if v=='1':
                num1 += 1
        if num1%2 == 1:
            cl = list()
            for v in range(len(vars)):
                if s[v]=='1':
                    cl.append(-vars[v])
                else:
                    cl.append(vars[v])
            clauses.append(cl)

def simple_geq_to_cnf(a, b, varcount, clauses):        
    assert(len(a) == len(b)) 
    equivalence_vars = list()

    for i in range(len(a)-1):
        equivalence_vars.append(varcount)
        varcount +=1
        simple_xor_to_cnf(equivalence_vars[-1],[a[i],b[i]], clauses)

    or_components = list()
    for i in range(len(a)):
        and_components = [a[i],-b[i]]
        for j in range(i-1,-1,-1):
            and_components.append(-equivalence_vars[j])
        
        and_variable = varcount
        varcount += 1
        simple_and_to_cnf(and_variable, and_components,clauses)
        or_components.append(and_variable)

    clauses.append(or_components) 
    return varcount   


def encode_to_SAT_noopt(N:int, K:int, D:int, SymmetryBreaking = True, Decompose_Xors = True, ENCTYPE = 1):
    # The main method to encode to SAT the problem of finding an (N,K,D)-code with error coefficient at most T

    GM, XORS, CARDS, UNITS, ORDER, nvars = make_base_constraints(N,K,D, SymmetryBreaking, Decompose_Xors)
    
    CLAUSES = list()
    for u in range(len(XORS)):
        simple_xor_to_cnf(XORS[u][0], XORS[u][1], CLAUSES)
    
    for u in range(len(ORDER)):
        nvars = simple_geq_to_cnf(ORDER[u][0], ORDER[u][1], nvars, CLAUSES)

    for u in UNITS:
        if UNITS[u] == 1:
            CLAUSES.append([u])
        else:
            CLAUSES.append([-u])
        
    
    for CC in CARDS:
        vars = CC[1]
        bound = CC[2]
        cl = CardEnc.atleast(vars, bound, nvars, None, encoding = ENCTYPE)            
        nvars = cl.nv + 1 
        for c in cl:
            CLAUSES.append(c)
    
    filename = "code_{}_{}_{}".format(N,K,D)
    if Decompose_Xors:
        filename+="_DX"
    if SymmetryBreaking:
        filename+="_SB"    
    
    filename+= "_noopt"
    
    filename+= "_e" + str(ENCTYPE)

    filename+= ".cnf"
    
    with open(filename,'w') as outfile:
        outfile.write("p cnf {} {}\n".format(nvars, len(CLAUSES)))
        for CL in CLAUSES:
            outfile.write(" ".join([str(u) for u in CL])+" 0\n")


def encode_to_SAT_opt(N:int, K:int, D:int, SymmetryBreaking = True, Decompose_Xors = True, BLOCKS = [], T = 0, ENCTYPE = 1):
    # The main method to encode to SAT the problem of finding an (N,K,D)-code with error coefficient at most T

    GM, XORS, CARDS, UNITS, ORDER, OPTVARS, nvars = optimization_formulation(N, K, D, SymmetryBreaking, Decompose_Xors, BLOCKS)
    
    CLAUSES = list()
    for u in range(len(XORS)):
        simple_xor_to_cnf(XORS[u][0], XORS[u][1], CLAUSES)
    
    for u in range(len(ORDER)):
        nvars = simple_geq_to_cnf(ORDER[u][0], ORDER[u][1], nvars, CLAUSES)

    for u in UNITS:
        if UNITS[u] == 1:
            CLAUSES.append([u])
        else:
            CLAUSES.append([-u])
        
    
    for CC in CARDS:
        vars = CC[1]
        bound = CC[2]
        cl = CardEnc.atleast(vars, bound, nvars, None, encoding = ENCTYPE)            
        nvars = cl.nv + 1 
        for c in cl:
            CLAUSES.append(c)
    
    # Add the constraint on the number of ones in OPTVARS
    cl = CardEnc.atmost(OPTVARS, T, nvars, None, encoding = ENCTYPE)            

    nvars = cl.nv + 1 
    for c in cl:
        CLAUSES.append(c)


    filename = "code_{}_{}_{}".format(N,K,D)
    if Decompose_Xors:
        filename+="_DX"
    if SymmetryBreaking:
        filename+="_SB"
    if BLOCKS == []:
        filename+="_all"
    else:
        filename+= "_" + "-".join([str(u) for u in BLOCKS])            
    
    filename+= "_T" + str(T)
    
    filename+= "_e" + str(ENCTYPE)

    filename+= ".cnf"
    
    with open(filename,'w') as outfile:
        outfile.write("p cnf {} {}\n".format(nvars, len(CLAUSES)))
        for CL in CLAUSES:
            outfile.write(" ".join([str(u) for u in CL])+" 0\n")
    

def encode_to_MaxSAT_opt(N:int, K:int, D:int, SymmetryBreaking = True, Decompose_Xors = True, BLOCKS = [], ENCTYPE = 1):
    # The main method to encode to MaxSAT the problem of finding an (N,K,D)-code with minimal error coefficient

    GM, XORS, CARDS, UNITS, ORDER, OPTVARS, nvars = optimization_formulation(N,K,D,SymmetryBreaking, Decompose_Xors, BLOCKS)
    
    CLAUSES = list()
    for u in range(len(XORS)):
        simple_xor_to_cnf(XORS[u][0], XORS[u][1], CLAUSES)
    
    nvars_old = nvars
    for u in range(len(ORDER)):
        nvars = simple_geq_to_cnf(ORDER[u][0], ORDER[u][1], nvars, CLAUSES)
        if nvars != nvars_old:
            print("well shit")

    for u in UNITS:
        if UNITS[u] == 1:
            CLAUSES.append([u])
        else:
            CLAUSES.append([-u])
        
    
    for CC in CARDS:
        vars = CC[1]
        bound = CC[2]
        cl = CardEnc.atleast(vars, bound, nvars, None, encoding = ENCTYPE)            
        nvars = cl.nv + 1 
        for c in cl:
            CLAUSES.append(c)
    
    filename = "code_{}_{}_{}".format(N,K,D)
    if Decompose_Xors:
        filename+="_DX"
    if SymmetryBreaking:
        filename+="_SB"
    if BLOCKS == []:
        filename+="_all"
    else:
        filename+= "_" + "-".join([str(u) for u in BLOCKS])            
    
    filename+= "_e" + str(ENCTYPE)

    filename+= ".wcnf"
    
    with open(filename,'w') as outfile:        
        for CL in CLAUSES:
            outfile.write("h " + " ".join([str(u) for u in CL])+" 0\n")
        
        for var in OPTVARS:
            outfile.write("1 {} 0\n".format(-var))
    
def print_help():
    print("To use the script invoke it as: construct_SAT.py [SAT, MAXSAT] [noopt, integer_value], N K D [DX,noDX] [SB,noSB] [1,2,3,6,7]")
    print("All arguments are mandatory")
    print("The second argument ([noopt,integer_value]) is used only if the first argument is SAT")
    print("noOPT means the basic formulation, without minimizing the error coefficient")
    print("If an integer is provided, then it will encode the problem of finding a code with T <= provided value")
    print("Otherwise, fill by any symbol")
    print("N, K, D are the characteristics of the code")
    print("DX means to decompose longer xors into xors with at most 2 inputs, strongly recommended to set it")
    print("SB means to use symmetry breaking. In the experiments it usually works better than without it for SAT / MaxSAT")
    print("The last number refers to the identifier of cardinality encoding in PySAT")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_help()
    else:    
        if len(sys.argv)  != 9:
            if len(sys.argv) == 2:
                if str(sys.argv[1].upper()) in ["-H","H","--HELP","-HELP","HELP"]:
                    print_help()
            else:
                print("Wrong number of arguments")
        else:
            noOPT = True
            targetT = 0
            SAT = True
            MAXSAT = True
            SB = True
            DX = True
            ET = 1
            BLOCKS = []            
            if str(sys.argv[1]).upper() == "SAT":
                MAXSAT = False
            elif str(sys.argv[1]).upper() == "MAXSAT":
                SAT = False
            else:
                print("Type of encoding is not specified")
                exit(0)
            
            if SAT:
                if str(sys.argv[2]).upper() != "NOOPT":
                    noOPT = False
                    targetT = int(sys.argv[2])                

            N = int(sys.argv[3])
            K = int(sys.argv[4])
            D = int(sys.argv[5])
            assert(N > 0)
            assert(K > 0)
            assert(D > 0)

            if str(sys.argv[6]).upper() in ["DX","NODX"]:
                if str(sys.argv[6]).upper() == "NODX":
                    DX = False
            else:
                print("DX is not set properly")
                exit(0)

            if str(sys.argv[7]).upper() in ["SB","NOSB"]:
                if str(sys.argv[7]).upper() == "NOSB":
                    SB = False
            else:
                print("SB is not set properly")
                exit(0)
            
            if int(sys.argv[8]) in [1,2,3,6,7]:
                ET = int(sys.argv[8])
            else:
                print("Cardinality encoding is specified incorrectly")
                exit(0)
            
            
            print("Symmetry breaking is set to {}, Decompose_XORs is set to {}".format(SB,DX))
            print("Cardinality encoding is set to {}".format(ET))
            if SAT:
                if noOPT:
                    print("Generating a SAT encoding for the problem of finding an ({},{},{})-code".format(N,K,D))                    
                    encode_to_SAT_noopt(N,K,D,SB,DX,ET)
                else:                    
                    print("Generating a SAT encoding for the problem of finding an ({},{},{})-code with target error coefficient at most {}".format(N,K,D,targetT))
                    encode_to_SAT_opt(N,K,D,SB,DX,BLOCKS, targetT, ET)
            else:
                    print("Generating a MaxSAT encoding for the problem of finding an ({},{},{})-code with the goal at minimizing error coefficient".format(N,K,D))
                    encode_to_MaxSAT_opt(N, K, D, SB, DX,BLOCKS, ET)
            

            

    

