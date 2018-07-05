#pragma once
#include <Arduino.h>
#include "memory.h"


void SerialPrintLeadingZeros(uint8_t Number);
void SerialPrintPaddedHex(uint8_t Number);
void SerialPrintPaddedHexArray(void *Number, size_t Length);
void SerialPrintPaddedBin(uint8_t Number);
void SerialPrintFormattedHex(uint16_t Value, size_t Index);
void SerialPrintFormattedHex(uint8_t Value, size_t Index);

void printLLNumber(uint64_t n, uint8_t base);

void Print2DMatrix(float **X, const size_t RowStart, const size_t ColStart, const size_t RowEnd, const size_t ColEnd);
void PrintArray(float *a, const size_t Start, const size_t End);
void Print1DMatrix(float *X, const size_t Ncol, const size_t RowStart, const size_t ColStart, const size_t RowEnd, const size_t ColEnd);
void Print1DMatrix(SPISRAMFLOATARRAY &X, const size_t Ncol, const size_t RowStart, const size_t ColStart, const size_t RowEnd, const size_t ColEnd);
void Print1DMatrix(SPISRAMFLOATARRAY &X, const size_t Ncol);
void Print1DMatrix(SPISRAMFLOATARRAY &X);
void Print1DMatrixMATLAB(SPISRAMFLOATARRAY &X, size_t Ncol, String VarName);
void Print1DMatrixArrayMATLAB(std::vector <SPISRAMFLOATARRAY> &X, std::vector <size_t> &Ncol, String VarName);
void Print1DMatrixArrayMATLAB(std::vector <SPISRAMFLOATARRAY> &X, size_t Ncol, String VarName);
