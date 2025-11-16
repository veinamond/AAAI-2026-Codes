#pragma once

#include "algo.h"

void XOR(tBit* pDest, const tBit* pSrc, size_t Size)
{
	for (size_t i = 0; i < Size; i++) pDest[i] ^= pSrc[i];
}

void XOR(tBit* pDest, const tBit* pSrc1, const tBit* pSrc2, size_t Size)
{
	for (size_t i = 0; i < Size; i++) pDest[i] = pSrc1[i] ^ pSrc2[i];
}

void Gauss(tBit* pMatrix, unsigned& NumOfRows, unsigned NumOfColumns, bool ReversePass, unsigned* pPermutation)
{
	for (unsigned i = 0; i < NumOfRows; i++) {
		// identify the leading column
		unsigned c = i;
		bool Success = false;
		for (; c < NumOfColumns; c++) {
			unsigned C = (pPermutation) ? pPermutation[c] : c;
			for (unsigned j = i; j < NumOfRows; j++) {
				if (pMatrix[j * NumOfColumns + C]) {
					Success = true;
					if (j > i)
						XOR(pMatrix + i * NumOfColumns, pMatrix + j * NumOfColumns, NumOfColumns);
					break;
				}
			}
			if (Success) {
				if ((c != i) && pPermutation) std::swap(pPermutation[c], pPermutation[i]);
				break;
			}
		}
		if (!Success) {
			NumOfRows = i;
			break;
		}
		unsigned LoopStart = (ReversePass) ? 0 : (i + 1);
		unsigned C = (pPermutation) ? pPermutation[i] : c;
		for (unsigned j = LoopStart; j < NumOfRows; j++) {
			if (j == i) continue;
			if (pMatrix[j * NumOfColumns + C])
				XOR(pMatrix + j * NumOfColumns, pMatrix + i * NumOfColumns, NumOfColumns);
		}
	}
}

void MinSpan(uint8_t* pMatrix, unsigned NumOfRows, unsigned NumOfColumns)
{
	for (int i = NumOfRows - 1; i > 0; --i) {
		// identify the rightmost 1
		int c =NumOfColumns-1;
		while (c >= 0 && pMatrix[i * NumOfColumns + c] == 0) c -= 1;
		if (c < 0) continue;

		for (int j = 0; j < i; ++j) {
			if (pMatrix[j * NumOfColumns + c] != 0) {
				XOR(pMatrix + j * NumOfColumns, pMatrix + i * NumOfColumns, NumOfColumns);
			}
		}
	}
}