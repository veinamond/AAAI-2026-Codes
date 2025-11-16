#include "Trellis.h"

Trellis::Trellis(int n, int k, tBit* generatorMatrix){

	//cout << "Trellis for (" << n << ", " << k << ") code initialized\n";
	m_Length = n;
	m_Dimension = k;
	m_pGeneratorMatrix = generatorMatrix;

	unsigned* pPermutation = new unsigned[m_Length];
	for (unsigned i = 0; i < m_Length; i++) pPermutation[i] = i;
	Gauss(m_pGeneratorMatrix, m_Dimension, m_Length, false, pPermutation);

	MinSpan(m_pGeneratorMatrix, m_Dimension, m_Length);

	cout << "Generator matrix in min span form:\n";
	for (unsigned i = 0; i < m_Dimension; ++i) {
		for (unsigned j = 0; j < m_Length; ++j) {
			cout << (bool)m_pGeneratorMatrix[i * m_Length + j];
		}
		cout << endl;
	}

	m_pActiveNum = new int[m_Length + 1];
	m_pNode = new Vertex*[m_Length + 1];

	//find trellis width
	int* varNumBeginCol = new int[m_Length];
	int* varNumEndCol = new int[m_Length];
	bool* state = new bool[m_Dimension];
	bool* isVarActiveNow = new bool[m_Dimension];
	bool* isVarActiveNext = new bool[m_Dimension];

	memset(isVarActiveNow, false, sizeof(bool)* m_Dimension);
	memset(isVarActiveNext, false, sizeof(bool)* m_Dimension);
	memset(varNumBeginCol, -1, sizeof(int)* m_Length);
	memset(varNumEndCol, -1, sizeof(int)* m_Length);
	memset(state, false, sizeof(bool)* m_Dimension);

	int counter = 0;
	for (unsigned i = 0; i < m_Dimension; ++i){
		for (unsigned j = 0; j < m_Length; ++j){
			if ((m_pGeneratorMatrix[i * m_Length + j] != 0)){
				varNumBeginCol[j] = counter;
				counter++;
				break;
			}
		}
	}

	counter = 0;
	for (unsigned i = 0; i < m_Dimension; ++i){
		for (unsigned j = m_Length - 1; j >= 0; j--){
			if ((m_pGeneratorMatrix[i * m_Length + j] != 0)){
				varNumEndCol[j] = counter;
				counter++;
				break;
			}
		}
	}
	counter = 0;

	int max_width = 0;
	m_pActiveNum[0] = 0;
	for (unsigned i = 1; i < m_Length + 1; ++i){
		m_pActiveNum[i] = m_pActiveNum[i - 1] + ((varNumBeginCol[i - 1] != -1) ? 1 : 0) - ((varNumEndCol[i - 1] != -1) ? 1 : 0);
		if (m_pActiveNum[i] > max_width){
			max_width = m_pActiveNum[i];
		}
	}

	for (unsigned i = 0; i < m_Length + 1; ++i){
		int t = 1 << m_pActiveNum[i];
		m_pNode[i] = new Vertex[t];
	}

	bool* bit_view_now = new bool[max_width];
	bool* bit_view_next = new bool[max_width];

	if (m_pActiveNum[1] == 0){
		m_pNode[0][0] = Vertex(0, 0, 0, 1);
	}
	else{
		m_pNode[0][0] = Vertex(0, 1, 0, 1);
	}


	isVarActiveNow[varNumBeginCol[0]] = true;
	if (varNumEndCol[1] == -1 || varNumEndCol[1] != varNumBeginCol[0]){
		isVarActiveNext[varNumBeginCol[0]] = true;
	}

	for (unsigned i = 1; i < m_Length; ++i){
		memset(state, false, sizeof(bool)* m_Dimension);
		memset(bit_view_now, false, sizeof(bool) * max_width);
		memset(bit_view_next, false, sizeof(bool) * max_width);

		//choose the active variables
		int t = varNumBeginCol[i - 1];
		if (t != -1){
			isVarActiveNow[t] = true;
		}

		int add = varNumBeginCol[i];
		if (add != -1){
			isVarActiveNext[add] = true;
		}

		t = varNumEndCol[i];
		if (t != -1){
			isVarActiveNext[t] = false;
		}

		int active_num = m_pActiveNum[i];
		int active_num_next = m_pActiveNum[i + 1];
		int loop_cond = 1 << active_num;
		//loop on vertices
		for (int j = 0; j < loop_cond; ++j){
			memset(state, false, sizeof(bool)* m_Dimension);
			int num = j;
			for (int k = 0; k < active_num; ++k){
				bit_view_now[k] = (num % 2 == 1);
				num /= 2;
			}

			int counter = 0;
			for (unsigned k = 0; k < m_Dimension; ++k){
				if (isVarActiveNow[k]){
					state[k] = bit_view_now[counter];
					counter++;
				}
			}

			int now_index = 0;
			int degree_two = 1;
			for (int k = active_num - 1; k >= 0; --k){
				now_index += degree_two * bit_view_now[k];
				degree_two *= 2;
			}
			counter = 0;
			for (unsigned k = 0; k < m_Dimension; ++k){
				if (isVarActiveNext[k]){
					bit_view_next[counter] = state[k];
					counter++;
				}
			}

			int index = 0;
			degree_two = 1;
			for (int k = active_num_next - 1; k >= 0; --k){
				index += degree_two * bit_view_next[k];
				degree_two *= 2;
			}

			bool result = false;
			for (unsigned k = 0; k < m_Dimension; ++k){
				result ^= ((m_pGeneratorMatrix[k * m_Length + i] != 0) && state[k]);
			}
			m_pNode[i][now_index] = Vertex(index, -1, result, false);

			//copy-paste
			if (add != -1){
				memset(bit_view_next, false, sizeof(bool) * max_width);
				state[add] = true;

				counter = 0;
				for (unsigned k = 0; k < m_Dimension; ++k){
					if (isVarActiveNext[k]){
						bit_view_next[counter] = state[k];
						counter++;
					}
				}
				int index = 0;
				degree_two = 1;
				for (int k = active_num_next - 1; k >= 0; --k){
					index += degree_two * bit_view_next[k];
					degree_two *= 2;
				}

				bool result = false;
				for (unsigned k = 0; k < m_Dimension; ++k){
					result ^= ((m_pGeneratorMatrix[k* m_Length + i] != 0) && state[k]);
				}
				m_pNode[i][now_index].right = index;
				m_pNode[i][now_index].sr = result;
			}
		}

		for (unsigned j = 0; j < m_Dimension; ++j){
			isVarActiveNow[j] = isVarActiveNext[j];
		}
	}

	delete[] state;
	delete[] bit_view_now;
	delete[] bit_view_next;
	delete[] isVarActiveNow;
	delete[] isVarActiveNext;
	delete[] varNumBeginCol;
	delete[] varNumEndCol;
	delete[] pPermutation;

	reverse_edges();
	return;
}

void Trellis::reverse_edges(){
	for (int i = m_Length - 1; i >= 0; --i){
		//cleaning next level
		int loop_cond = (1 << m_pActiveNum[i + 1]);
		for (int j = 0; j < loop_cond; ++j){
			m_pNode[i + 1][j] = Vertex();
		}
		loop_cond = (1 << m_pActiveNum[i]);
		for (int j = 0; j < loop_cond; ++j){
			Vertex v(m_pNode[i][j]);
			if (m_pNode[i + 1][v.left].left == -1){
				m_pNode[i + 1][v.left].left = j;
				m_pNode[i + 1][v.left].sl = v.sl;
			}
			else{
				m_pNode[i + 1][v.left].right = j;
				m_pNode[i + 1][v.left].sr = v.sl;
			}

			if (v.right != -1){
				if (m_pNode[i + 1][v.right].left == -1){
					m_pNode[i + 1][v.right].left = j;
					m_pNode[i + 1][v.right].sl = v.sr;
				}
				else{
					m_pNode[i + 1][v.right].right = j;
					m_pNode[i + 1][v.right].sr = v.sr;
				}
			}
		}
	}
	m_pNode[0][0] = Vertex();
}

Trellis::~Trellis(){
#ifdef TRELLIS_LOGS
	cout << "Trellis object destructor" << endl;
#endif
	for (unsigned i = 0; i <= m_Length; ++i){
		delete[] m_pNode[i];
	}
	delete[] m_pNode;
	delete[] m_pActiveNum;
}