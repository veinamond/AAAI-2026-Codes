import sys
import itertools
import os

def Read_gen(fn):
    N = 0
    K = 0
    D = 0
    GM = list()
    with open(fn,'r') as infile:
        for line in infile:
            if N==0:
                p = line.strip("\r\n").split(" ")
                N = int(p[0])
                K = int(p[1])                
            else:
                p = line.strip("\r\n")
                h = [int(u) for u in p]
                GM.append(h)
    
    assert(len(GM) == K)
    for RM in GM:
        assert(len(RM)==N)

    spectrum = dict()
    for u in range(N+1):
        spectrum[u] = 0
        
    R = [u for u in range(K)]
    
    for u in range(2,K+1):
        for COMB in itertools.combinations(R,u):    
            if len(COMB) > 0:
                vector_sum = [0 for _ in range(N)]
                for i in range(N):
                    for j in COMB:
                        vector_sum[i]^=GM[j][i]                
                ccnt = count_ones(vector_sum)                
                spectrum[ccnt]+=1
    
    for u in range(N+1):
        if spectrum[u]!=0:
            D = u
            break
    Tval = spectrum[D]

    for i in range(K):
        BVi = GM[i]
        unit_i = [0 for _ in range(K)]
        unit_i[i] = 1
        fff = unit_i+BVi
        print("".join([str(t) for t in fff]))

    print()
    print("D = {}".format(D))    
    print("T = {}".format(spectrum[D]))
    
    print("Weight spectrum: ")
    for u in spectrum:
        if spectrum[u]!=0:
            print("Codewords of weight {} : {}".format(u, spectrum[u]))
    return Tval

def read_kissat_output(filename):
    model = list()
    reading_ss = False
    is_SAT = False
    ss = ""
    SS_read = False
    with open(filename,'r') as infile:        
        for line in infile:            
            #print(line)
            if line=="c\n" or  line.startswith("c"):
                #if reading_ss == True:
                    #print(line)
                reading_ss = False
            if line.startswith("s SATISFIABLE"):
                is_SAT = True
            elif line.startswith("v"):
                if not SS_read:
                    reading_ss = True
                    SS_read = True
                if reading_ss:
                    ss += line.strip(" \r\n")
            elif reading_ss == True and len(line) > 1:
                ss += " " + line.strip(" \r\n")
                
    
    if is_SAT == False:
        return False, {} 

    p = ss.split(" ")
    g = list()
    model.append(0)
    #176
    #-1
    #77
    incomplete_token = False

    for i in range(len(p)):
        #print(i, p[i])
        p[i] = p[i].replace("v","")
        t = p[i]
        if t == "-":
            p[i+1] = t+p[i+1]
        else:
            if t != "" and t !="v":                            
                if (abs(int(t))) == abs(model[-1]) + 1:                
                    model.append(int(t))
                else:                                
                    if i != len(p) - 1:                    
                        #print("last token {} current token {}".format(str(model[-1]),u))
                        assert(i > 1)
                        if incomplete_token:
                            f = p[i-1] + p[i]
                            if abs(int(f)) != abs(model[-1])+1:
                                p[i] = f
                            else:
                                model.append(int(f))
                                incomplete_token = False
                        else:
                            incomplete_token = True
                    
    
    return True, model

def count_ones(L:list):
    r = 0
    for u in L:
        if u == 1:
            r+= 1
    return r

def check_solution(N:int, K:int, D: int, M:list, verbosity = True, outfile=""):
    valid_solution = True
    Tval = 0
    # assume that M is a list of ints containing a model in DIMACS format, e.g. -1 -2 3 4 5 -6, etc
    cnt = 1
    mm = dict()
    for u in M:
        if u < 0:
            mm[abs(u)] = 0
        else:
            mm[abs(u)] = 1

    GMvars = list()
    for i in range(K):
        bvvars = [cnt + t for t in range(N-K)]        
        GMvars.append(bvvars)
        cnt += N-K
    
    GM = list()
    
    for bvvars in GMvars:
        BVvals = [mm[u] for u in bvvars]
        GM.append(BVvals)

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

def print_help():
    print("To use the script invoke it as: check_solution.py kissat_output N K D")
    print("Alternatively, it can be used to check .gen files and outputting their info")
    print("In this case just use it as check_solution.py /path_to_file/file.gen")
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_help()
    else:
        if len(sys.argv) == 2:
            if sys.argv[1].endswith(".gen"):
                Read_gen(sys.argv[1])
            elif str(sys.argv[1].upper()) in ["-H","H","--HELP","-HELP","HELP"]:
                print_help()
        else:
            if len(sys.argv) < 5:
                print("Not enough arguments")
            else:
                fn = sys.argv[1]
                N = int(sys.argv[2])
                K = int(sys.argv[3])
                D = int(sys.argv[4])
                R, m = read_kissat_output(fn)
                if R == False:
                    print("CNF is unsatisfiable")
                else:
                    check_solution(N,K,D, m)







