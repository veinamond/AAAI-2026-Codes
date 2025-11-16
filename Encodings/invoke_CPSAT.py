from ortools.sat.python import cp_model
from constraint_encoding import *
import os
import sys
Nglobal = 0
Kglobal = 0
Dglobal = 0

class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    def __init__(self, variables: list[cp_model.IntVar]):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables = variables
        self.__solution_count = 0

    def on_solution_callback(self) -> None:
        assert(Nglobal > 0 and Kglobal > 0 and Dglobal > 0) 
        BV = [[] for _ in range(Kglobal)]
        cnt = 1
        for i in range (Kglobal):
            for j in range(Nglobal-Kglobal):
                BV[i].append(cnt)
                cnt += 1

        self.__solution_count += 1
        for v in self.__variables:
            ind = v.index
            t = v.name
            t = t.strip("x[]")            
            t_int = int(t) - 1            
            if t_int != ind:
                print("?FUCK?")
            i = t_int //(Nglobal-Kglobal)
            j = t_int - i*(Nglobal-Kglobal)
            BV[i][j]=self.value(v)            

        check_solution_CPSAT(Nglobal,Kglobal,Dglobal,BV)    
        
        print("Current bound: {} @ {}".format(self.objective_value, self.WallTime()))

    @property
    def solution_count(self) -> int:
        return self.__solution_count


def count_ones(L:list):
    r = 0
    for u in L:
        if u == 1:
            r+= 1
    return r

def check_solution_CPSAT(N:int, K:int, D: int, GM:list, verbosity = True, outfile=""):
    # The main difference is in how we obtain the values of variables
    valid_solution = True
    Tval = 0
    # assume that M is a list of ints containing a model in DIMACS format, e.g. -1 -2 3 4 5 -6, etc
    
    R = [u for u in range(K)]
    
    spectrum = dict()
    for u in range (N+1):
        spectrum [u] = 0
    
    for u in range(1,K+1):
        for COMB in itertools.combinations(R,u):
            if D - len(COMB) > 0:
                vector_sum = [0 for _ in range(N-K)]
                for i in range(N-K):
                    for j in COMB:
                        vector_sum[i] ^= GM[j][i]
                
                HW = count_ones(vector_sum) + len(COMB)

                spectrum[HW] += 1
                
                if HW == D:
                    Tval += 1
                if HW < D:
                    valid_solution = False
    
    for i in range(K):
        BVi = GM[i]
        unit_i = [0 for _ in range(K)]
        unit_i[i] = 1
        fff = unit_i+BVi
        print("".join([str(t) for t in fff]))
    
    if verbosity:
        if valid_solution == True:
            print("Solution is valid")
            print("T = {}".format(spectrum[D]))
        else:
            print("Solution is invalid")
        

    print("Weight spectrum: ")
    for u in spectrum:
        if spectrum[u]!=0:
            print("Codewords of weight {} : {}".format(u, spectrum[u]))
    #print("Statistics for combinations: "+ str(sizes_comb))
    #print("All sizes: "+ str(sizes_all))
    
    outfile = "{}_{}_{}_{}.gen".format(N,K,D,Tval)
    cnt_file = 0
    while os.path.exists(outfile):
        outfile = "{}_{}_{}_{}_{}.gen".format(N,K,D,Tval,cnt_file)
        cnt_file+=1    

    
    with open(outfile, 'w') as wr:    
        wr.write("{} {}\n".format(N,K))
        for i in range(K):
            xva = GM[i]
            unit_i = [0 for _ in range(K)]
            unit_i[i] = 1
            fff = unit_i+xva             
            wr.write("".join([str(t) for t in fff])+"\n")
    return Tval, valid_solution

def simple_geq_to_CPSAT(a, b, varcount, m):        
    assert(len(a) == len(b)) 
    equivalence_vars = list()   
    
    for i in range(len(a)-1):
        equivalence_vars.append(m.new_bool_var("x(%s)" % varcount))        
        varcount +=1
        m.AddBoolXOr([a[i],b[i], equivalence_vars[-1], True])        

    or_components = list()    
    for i in range(len(a)):
        and_components = [a[i],b[i].Not()]
        for j in range(i-1,-1,-1):
            and_components.append(equivalence_vars[j].Not())
        
        and_variable = m.new_bool_var("x(%s)" % varcount)
        varcount += 1

        or_components.append(and_variable)

        ff = [u.Not() for u in and_components]
        ff.append(and_variable)

        m.add_bool_or(ff)

        for u in and_components:
            m.AddImplication(or_components[-1],u)

    m.add_bool_or(or_components)
        

    return varcount   

    
    

def solve_using_ortools(N:int, K:int, D:int, SB : bool, DX : bool, bl:list,  limsec : int, nthreads : int):
    GM, XORS, CARDS, UNITS, ORDER, OPTVARS, nvars = optimization_formulation(N,K,D,SymmetryBreaking=SB, Decompose_Xors=DX, BLOCKS = bl)

    global Nglobal, Kglobal, Dglobal
    Nglobal = N
    Kglobal = K
    Dglobal = D
  
    m = cp_model.CpModel()
  
    # Make variables for Generator Matrix
    
    x = {}
    
    for i in range(K):
        for j in range(N-K):            
            x[GM[i][j]] = m.new_bool_var("x[%s]" % (GM[i][j]))            
        
  

    for u in range(len(XORS)):        
        x[XORS[u][0]] = m.new_bool_var("x(%s)" % XORS[u][0])        
        rp = [x[i] for i in XORS[u][1]] + [x[XORS[u][0]]] + [True]        
        m.AddBoolXOr(rp)

    for u in OPTVARS:
        x[u] = m.new_bool_var("x[%s]"%u)
    
    for u in UNITS:        
        m.add(x[u]==UNITS[u])

    
    for u in range(len(ORDER)):        
        rp1 = [x[i] for i in ORDER[u][0]]
        rp2 = [x[i] for i in ORDER[u][1]]
        nvars = simple_geq_to_CPSAT(rp1, rp2, nvars, m)
        
    for CC in CARDS:
        vars = [x[u] for u in CC[1]]
        bound = CC[2]
        m.add(sum(vars) >= bound)
        
    m.minimize(sum([x[i] for i in OPTVARS]))
   
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = nthreads
    solver.parameters.max_time_in_seconds = limsec

    solution_printer = VarArraySolutionPrinter([x[GM[i][j]] for i in range(K) for j in range(N-K)])

    status = solver.solve(m, solution_printer)
        
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f"Maximum of objective function: {solver.objective_value}\n")
        BV = [[] for _ in range(K)]
        for i in range(K):
            for j in range(N-K):
                
                BV[i].append(int(solver.value(x[GM[i][j]])))
        
        check_solution_CPSAT(N,K,D,BV)    
    else:
        print("No solution found.")
    
    print("Done!")

def print_help():
    print("To use the script invoke it as: python3 invoke_CPSAT.py N K D [DX|noDX] [SB|noSB] limsec nthreads")
    print("All arguments are mandatory")
    print("N, K, D are the characteristics of the code")
    print("DX means to decompose longer xors into xors with at most 2 inputs, strongly recommended to unset it, e.g. to use noDX")
    print("SB means to use symmetry breaking")
    print("The limsec is the limit in seconds given to the CPSAT")
    print("The nthreads is the number of threads given to the CPSAT. Note, that in the parallel mode CPSAT is not deterministic.")

N = 0
K = 0
D = 0


Decompose_XORs = False

SymmetryBreaking = False
limsec = 600 
nthreads = 1
BLOCKS = []


if len(sys.argv) < 7:
    print_help()
    exit(0)
N = int(sys.argv[1])
K = int(sys.argv[2])
D = int(sys.argv[3])
assert(N > 0 and K > 0 and D > 0)

if str(sys.argv[4]).upper() in ["DX","NODX"]:
    if str(sys.argv[4]).upper() == "NODX":
            Decompose_XORs = False
else:
    print("DX is not set properly")
    exit(0)

if str(sys.argv[5]).upper() in ["SB","NOSB"]:
    if str(sys.argv[5]).upper() == "NOSB":
        SymmetryBreaking = False
else:
    print("SB is not set properly")
    exit(0)


limsec = int(sys.argv[6])
nthreads = int(sys.argv[7])

assert(limsec >= 1)
assert(nthreads >= 1)
print('N : ' + str(N))
print('K : ' + str(K))
print('D : ' + str(D))
print('DecXors : ' + str(Decompose_XORs))
print('SymBr : ' + str(SymmetryBreaking))
print('nthreads : ' + str(nthreads))
print('limsec : ' + str(limsec))

solve_using_ortools(N, K, D, SymmetryBreaking,Decompose_XORs,BLOCKS,limsec,nthreads)

