#include "FastTTLViterbi.h"

FastTTLViterbi::FastTTLViterbi(int n, int k, tBit* generatorMatrix, unsigned listSize):Trellis(n,k,generatorMatrix)
{
	m_ListSize = listSize;
	m_pPath = new unsigned[m_ListSize * (m_Length + 1)];
	m_pM1 = new KType*[m_Length + 1];
	m_pM2 = new KType*[m_Length + 1];
	m_pV = new unsigned*[m_Length + 1];
	m_pVS = new unsigned*[m_Length + 1];
	for (unsigned i = 0; i < m_Length + 1; ++i){
		if (m_pActiveNum[i] > 31){
			std::cerr << "Error, trellis.m_pActiveNum" << endl;
		}
		unsigned t = 1 << m_pActiveNum[i];
		m_pM1[i] = new KType[t];
		m_pM2[i] = new KType[t];
		m_pV[i] = new unsigned[t];
		m_pVS[i] = new unsigned[t];
	}
	m_pHardDecisionPoints = new bool[m_Length];
	m_pAbsDecisionPoints = new KType[m_Length];
}

FastTTLViterbi::~FastTTLViterbi(){
#ifdef TRELLIS_LOGS
	cout << "OuterDecoderTreeTrellis Destructor" << endl;
#endif
	for (unsigned i = 0; i <= m_Length; ++i){
		delete[] m_pM1[i];
		delete[] m_pM2[i];
		delete[] m_pV[i];
		delete[] m_pVS[i];
	}
	delete[] m_pPath;
	delete[] m_pM1;
	delete[] m_pM2;
	delete[] m_pV;
	delete[] m_pVS;
	delete[] m_pAbsDecisionPoints;
	delete[] m_pHardDecisionPoints;
}

bool FastTTLViterbi::forwardPass(const KType* pLLR, tBit* pOutCodeword, KType& pMetrics){
	/*
	for trellis representation, where node know predecessors
	*/
	m_NextCodewordIndex = 0;
	B.clear();
	for (unsigned i = 0; i < m_Length; ++i){
		if (pLLR[i] > 0){
			m_pHardDecisionPoints[i] = false;
			m_pAbsDecisionPoints[i] = pLLR[i];
		}
		else{
			m_pHardDecisionPoints[i] = true;
			m_pAbsDecisionPoints[i] = -((KType)pLLR[i]);
		}
	}

	m_pM1[0][0] = 0;
	m_pM2[0][0] = 0;

	for (unsigned i = 1; i <= m_Length; ++i){
		int loop_cond = (1 << m_pActiveNum[i]);
		for (int j = 0; j < loop_cond; ++j){
			Vertex* node = &m_pNode[i][j];
			KType left_branch_score = m_pM1[i - 1][node->left];
			KType m = m_pAbsDecisionPoints[i - 1];
			if (node->sl){
				if (!m_pHardDecisionPoints[i-1]){
					left_branch_score -= m;
				}
			}
			else{
				if (m_pHardDecisionPoints[i-1]){
					left_branch_score -= m;
				}
			}
			if (node->right == -1){
				m_pM1[i][j] = left_branch_score;
				m_pM2[i][j] = left_branch_score;
				m_pV[i][j] = node->left;
				m_pVS[i][j] = node->left;
			}
			else{
				KType right_branch_score = m_pM1[i - 1][node->right];
				if (node->sr){
					if (!m_pHardDecisionPoints[i-1]){
						right_branch_score -= m;
					}
				}
				else{
					if (m_pHardDecisionPoints[i-1]){
						right_branch_score -= m;
					}
				}

				if (left_branch_score > right_branch_score){
					m_pM1[i][j] = left_branch_score;
					m_pV[i][j] = node->left;
					m_pM2[i][j] = right_branch_score;
					m_pVS[i][j] = node->right;
				}
				else{
					m_pM1[i][j] = right_branch_score;
					m_pV[i][j] = node->right;
					m_pM2[i][j] = left_branch_score;
					m_pVS[i][j] = node->left;
				}
			}
		}
	}

	return backwardPass(pOutCodeword, pMetrics);
}

KType FastTTLViterbi::backwardPass(tBit* pOutCodeword, KType& pMetrics){

	//_ASSERT(m_NextCodewordIndex < m_ListSize);

	//step 1
	unsigned i, j, t;
	KType m;

	if (m_NextCodewordIndex == 0){
		i = 0;
		t = m_Length;
		m = 0;
		m_pPath[m_NextCodewordIndex * (m_Length + 1) + t] = 0;
	}
	else{
		auto b_end = B.end();
		auto&& itMax = B.begin();
		for (auto&& it = B.begin(); it != B.end(); ++it){
			if (it->m > itMax->m){
				itMax = it;
			}
		}
		Element S(*itMax);
		//step 2
		t = S.t;
		i = m_pPath[S.q * (m_Length + 1) + t];
		m = S.u;
		for (unsigned p = t; p <= m_Length; ++p){//memcpy
			m_pPath[m_NextCodewordIndex * (m_Length + 1) + p] = m_pPath[S.q * (m_Length + 1) + p];
		}
		//step 3
		B.erase(itMax);
		//step 4
		j = m_pVS[t][i];
		//step 5
		m += m_pM2[t][i] - m_pM1[t - 1][j];
		m_pPath[m_NextCodewordIndex * (m_Length + 1) + t - 1] = j;
		t--;
		i = j;
	}
	//step 9
	while (t != 0){
		//step 6
		KType tmp = m_pM2[t][i] + m;
		//comp_count++;
		if (m_pV[t][i] != m_pVS[t][i]){
			B.push_back(Element(tmp, m_NextCodewordIndex, t, m));
		}

		//step 7
		j = m_pV[t][i];
		//step 8
		m += m_pM1[t][i] - m_pM1[t - 1][j];
		m_pPath[m_NextCodewordIndex * (m_Length + 1) + t - 1] = j;
		t--;
		i = j;
	}

	pMetrics = m;

	int current_node = 0;
	for (int i = m_Length; i > 0; --i){
		unsigned next_node = m_pPath[m_NextCodewordIndex * (m_Length + 1) + i - 1];
		Vertex* node = &m_pNode[i][current_node];
		if (node->left == next_node){
			pOutCodeword[i - 1] = node->sl;
		}
		else{
			pOutCodeword[i - 1] = node->sr;
		}
		current_node = next_node;
	}

	m_NextCodewordIndex++;

	return m;
}