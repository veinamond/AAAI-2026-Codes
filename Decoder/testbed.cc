#include <iostream>
#include <fstream>
#include <random>

#include "FastTTLViterbi.h"
using namespace std;

inline double EbN0StdDev(double SNR, double codeRate)
{
	return sqrt(0.5 * pow(10, -SNR / 10) / codeRate);
}

void PrintUsage() {
	cout << "Usage: SimViterbi.exe gen_matrix_file SNR max_err max_iterations\n";
	cout << "Generator matrix format:\nLength Dimension\nGenerator matrix with no spaces\nExample:\n8 4\n11110000\n11001100\n10101010\n11111111";
}

int main(int argc, char* argv[])
{
	unsigned I = 1;
	if (argc != 5) {
		PrintUsage();
		return 1;
	}
	try {
		ifstream gen_matrix(argv[I++]);
		if (!gen_matrix.good()) {
			string message = "Cannot open input file " + string(argv[1]);
			throw exception(message.c_str());
		}
		double SNR = atof(argv[I++]);
		size_t err = atoi(argv[I++]);
		size_t max_it = atoi(argv[I++]);

		unsigned N, K;
		gen_matrix >> N;
		gen_matrix >> K;

		if (K > N) throw exception("Dimension cannot exceed Length");

		tBit* pGenMatrix = new tBit[N * K];
		for (unsigned i = 0; i < K; ++i) {
			string row;
			gen_matrix >> row;
			if(row.length() != N) throw exception("Insufficient number of elements in generator matrix");
			for (unsigned j = 0; j < N; ++j) {
				if (!(row[j] == '0' || row[j] == '1')) throw exception("Wrong file content");
				pGenMatrix[i * N + j] = (row[j] == '1');
			}
		}

		FastTTLViterbi decoder(N, K, pGenMatrix, 1);
		cout << "Viterbi decoder for (" << N << ", " << K << ") code is initialized" << endl;
		double StdDev = EbN0StdDev(SNR, double(K) / N);

		size_t NumErrors = 0;
		size_t NumMLErrors = 0;

		std::mt19937_64 generator;
		generator.seed(8);
		std::uniform_int_distribution<unsigned> distribution(0, 1);
		std::normal_distribution<KType> noise(0, StdDev);

		tBit* pInfoBits = new tBit[K];
		tBit* pCodeword = new tBit[N];
		tBit* pOutputCodeword = new tBit[N];
		KType* pModulated = new KType[N];
		KType* pNoisy = new KType[N];
		size_t it = 0;
		for (it; it < max_it; ++it) {
			if (NumErrors >= err) break;

			//generate info bits
			for (unsigned i = 0; i < K; ++i)
				pInfoBits[i] = distribution(generator);

			fill_n(pCodeword, N, 0);

			//encode by matrix multiplication
			for (unsigned i = 0; i < K; ++i) {
				if (pInfoBits[i] != 0)
					XOR(pCodeword, pGenMatrix + i * N, N);
			}

			// modulate, 0 -> 1, 1 -> -1
			for (unsigned i = 0; i < N; ++i) {
				pModulated[i] = (1 - 2 * (bool)pCodeword[i]);
			}

			//add channel noise
			for (unsigned i = 0; i < N; ++i) {
				double r = noise(generator);
				pNoisy[i] = pModulated[i] + r;
			}

			KType score;
			//run decoder
			decoder.forwardPass(pNoisy, pOutputCodeword, score);

			//check error and score
			bool hasError = false;
			for (unsigned i = 0; i < N; ++i) {
				if (pOutputCodeword[i] != pCodeword[i]) {
					hasError = true;
					break;
				}
			}

			if (hasError) {
				NumErrors++;

				KType origScore = 0;

				for (unsigned i = 0; i < N; ++i) {
					bool sgn = pNoisy[i] > 0;
					if ((bool)pCodeword[i] == 0 && !sgn) {
						origScore += pNoisy[i];
					}
					else if ((bool)pCodeword[i] == 1 && sgn) {
						origScore -= pNoisy[i];
					}
				}

				if (origScore < score) NumMLErrors++;
			}

			if (it % 10000 == 0 && it != 0) {
				cout << "sim iteration = " << it << ", FER = " << (NumErrors / (double)(it)) << ", Num of errors = " << NumErrors << "\n";
			}
		}
		cerr << "#SNR StdDev FER MLFER NumErrors\n";
		cerr << SNR << "\t" << StdDev << "\t" << (NumErrors / (double)(it)) << "\t" << (NumMLErrors / (double)(it))<<"\t"<<NumErrors << endl;

		delete[] pGenMatrix;
		delete[] pInfoBits;
		delete[] pModulated;
		delete[] pNoisy;
		delete[] pCodeword;
		delete[] pOutputCodeword;
	}
	catch (std::exception& ex) {
		cerr << ex.what() << endl;
		return 1;
	}
	catch (...) {
		cerr << "Undefined exception!" << endl;
		return 2;
	}
	return 0;
}
