The repository contains supplementary material for the paper **Using Constraint Solvers to Construct Binary Codes with Good Error Correction Performance** accepted to the AAAI 2026 conference.

It contains:

- Codes found by our methods (for which we measure the error correcting performance) as well as the example codes from Figure 1 : in the Codes folder.
	- The *.gen files contain generator matrices for the respective codes. The filenames usually show the parameters of the code, the name of the method used to find it (e.g. kissat), the type of cardinality encoding (e.g. e7), and the value of the error coefficient (last number before .gen).
	- The *.fer file with the same name as the *.gen file contains the results of modeling FER for the respective code. There are usually other auxiliary values in these logs, however in all cases SNR is the value in the first column, and FER is the value in the third column.

- Scripts used to construct SAT / MaxSAT encodings and invoke CP-SAT: in the Encodings folder.
	- The constraint_encoding.py script details the construction of a relatively general set of constraints for the problem of finding an (n,k,d)-code.
	- The construct_SAT.py script is used to generate SAT and MaxSAT encodings. It has many parameters which can be shown by running script without parameters or with -h.
	- The invoke_CPSAT.py script uses the constraints constructed by constraint_encoding.py script to build a constraint model for CPSAT and invoke it. We could have created the model in CP-SAT straightaway, but decided on this variant for consistency.
		- Note, that we do not put blocks into command line parameters of the scripts, but the corresponding encodings (with blocks) can be constructed when importing scripts as source files.
	- The check_solution.py extracts the generator matrix of a code from the output of a SAT solver, in particular it is geared to use Kissat, but since all SAT solvers output satisfying assignments in more or less the same format, we expect it to work well for other solvers, possibly with minor modifications. When given a *.gen file, the script computes minimal distance D and the weight spectrum.
		- Note, that check_solution and invoke_CPSAT.py will put the *.gen file with the found solutions into the current working folder.
- Decoder implementing the Viterbi algorithm that can be used to compute FER for a code. This is a computationally intensive procedure which takes several hours per code.
	- To build the decode invoke make 
	- To construct the FER values that can be used for the plot invoke the compute_FER.sh script in the same folder followed by the .gen file with the generator matrix for a code. It will compute a range of FER values for different SNR.
- The PDF file with additional information on experiments, including the howto for Magma, QextNewEdition and Cube-and-Conquer, and the detatiled table listing the runtimes of SAT/MaxSAT/CP-SAT on code-finding problems with the corresponding best found codes.

We are not aware of any special requirements to run the code, apart from python-sat. We tested it on different Linux distributions (e.g. Ubuntu) and they should also work in Windows and MacOS. It is possible that for old Python versions something will break, but python 3.10 and newer should be fine. 
The decoder should also work in Windows without a problem, and can be built using e.g. MVS.


