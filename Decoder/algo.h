#pragma once

#include <iostream>
#include <cstdint>

using std::cout;
using std::endl;

typedef unsigned char tBit;
typedef float KType;
typedef uint64_t Word;

void XOR(tBit* pDest, const tBit* pSrc, size_t Size);

void XOR(tBit* pDest, const tBit* pSrc1, const tBit* pSrc2, size_t Size);

void Gauss(tBit* pMatrix, unsigned& NumOfRows, unsigned NumOfColumns, bool ReversePass, unsigned* pPermutation);

void MinSpan(uint8_t* pMatrix, unsigned NumOfRows, unsigned NumOfColumns);