#pragma once
#include "algo.h"
#include <cstring>
//Describes vertice for trellis diagramm
struct Vertex{
	int left;
	int right;
	bool sl;
	bool sr;

	Vertex(){ left = -1; right = -1;  sl = false; sr = false; }
	Vertex(int _left, int _right, bool _sl, bool _sr){
		left = _left;
		right = _right;
		sl = _sl;
		sr = _sr;
	}
	Vertex(const Vertex& node){
		left = node.left;
		right = node.right;
		sl = node.sl;
		sr = node.sr;
	}
};

//Represent trellis diagramm
class Trellis{
public:
	unsigned m_Length;
	unsigned m_Dimension;
	Vertex** m_pNode;
	int* m_pActiveNum;
	Trellis(){
		m_Length = 0;
		m_Dimension = 0;
		m_pGeneratorMatrix = nullptr;
		m_pNode = nullptr;
		m_pActiveNum = nullptr;
	}
	Trellis(int n, int k, tBit* generatorMatrix);
	~Trellis();
	void reverse_edges();
private:
	tBit* m_pGeneratorMatrix;
	Trellis(const Trellis&);
};