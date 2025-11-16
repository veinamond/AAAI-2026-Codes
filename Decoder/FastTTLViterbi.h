#pragma once

#include "Trellis.h"
#include <vector>

using std::vector;
/// Element of the stack
struct Element {
	KType m;
	KType u;
	unsigned q;
	unsigned t;
	Element(){ m = 0; q = 0; t = 0; u = 0; }
	Element(KType _m, unsigned _q, unsigned _t, KType _u){
		m = _m;
		q = _q;
		t = _t;
		u = _u;
	}
};

/**
 * Fast Tree Trellis List Viterbi decoder. Can be used for arbitrary code besides (n,n) code.
 */
class FastTTLViterbi: public Trellis
{
public:
	bool* m_pHardDecisionPoints;
	KType* m_pAbsDecisionPoints;

	FastTTLViterbi(int n, int k, tBit* generatorMatrix, unsigned listSize);
	virtual ~FastTTLViterbi();

	//preprocessing stage for TreeTrellis algorithm. Return latency
	virtual bool forwardPass(const KType* pLLR, tBit* pOutCodeword, KType& pMetrics);

	//returns k-th best codeword. Return latency and metrcis.
	virtual KType backwardPass(tBit* pOutCodeword, KType& pMetrics);

private:
	unsigned m_NextCodewordIndex;
	unsigned m_ListSize;
	//state variables for trellis vertices
	KType** m_pM1;
	KType** m_pM2;
	//pointer to previous vertices
	unsigned** m_pV;
	unsigned** m_pVS;
	//store path throught trellis corresponds to codeword
	unsigned* m_pPath;

	vector<Element> B;

	FastTTLViterbi(const FastTTLViterbi&);
};