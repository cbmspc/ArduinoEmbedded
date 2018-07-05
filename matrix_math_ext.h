#pragma once
#include "memory.h"

void CholeskyDecomp(SPISRAMFLOATARRAY &B, size_t n);
void Determinant(SPISRAMFLOATARRAY &A, float &d, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList);
void EigSolve(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &eig, float &l, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList);
void GeneralizedSelfAdjointEigSolver(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &B, SPISRAMFLOATARRAY &w, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList);
void LDA(SPISRAMFLOATARRAY &X1, SPISRAMFLOATARRAY &X2, size_t NumDim, SPISRAMFLOATARRAY &w, std::vector <SPISRAMClass *> &SRAMPtrList);
void LDA(SPISRAMFLOATARRAY &mu1, SPISRAMFLOATARRAY &mu2, SPISRAMFLOATARRAY &S1, SPISRAMFLOATARRAY &S2, SPISRAMFLOATARRAY &w, size_t n1, size_t n2, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList);
void LDLDecomp(SPISRAMFLOATARRAY &B, SPISRAMFLOATARRAY &D, size_t n);
void LowerTriangularInverse(SPISRAMFLOATARRAY &L, size_t n);
void MatrixMultiply(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &B, SPISRAMFLOATARRAY &C, size_t n, size_t m, size_t p);
void MatrixTranspose(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &AT, size_t n, size_t m);
void MatrixTransposeSquare(SPISRAMFLOATARRAY &A, size_t n);
void MatrixNormalize(SPISRAMFLOATARRAY &A, size_t Ncol);
void Normalize(SPISRAMFLOATARRAY &w);
float Norm(SPISRAMFLOATARRAY &w);
void PCA_SmallSampleSize(SPISRAMFLOATARRAY &X, SPISRAMFLOATARRAY &PA, size_t &nk, size_t NumSamples, size_t NumVariables, float RetainThres, std::vector <SPISRAMClass *> &SRAMPtrList);
void PCA(SPISRAMFLOATARRAY &S, size_t &nk, size_t n, float RetainThres, std::vector <SPISRAMClass *> &SRAMPtrList, size_t nkmax = SIZE_MAX);
void PCA(SPISRAMFLOATARRAY &X, SPISRAMFLOATARRAY &PA, size_t &nk, size_t NumSamples, size_t NumVariables, float RetainThres, std::vector <SPISRAMClass *> &SRAMPtrList);
void MatrixCovariance(SPISRAMFLOATARRAY &X, const size_t Nrow, const size_t Ncol, SPISRAMFLOATARRAY &S, std::vector <SPISRAMClass *> &SRAMPtrList, uint8_t SpecialOps = 0);
void MatrixDemean(SPISRAMFLOATARRAY &X, const size_t Nrow, const size_t Ncol, SPISRAMFLOATARRAY &Xdm, std::vector <SPISRAMClass *> &SRAMPtrList, bool DeMeanRowsInstead);
bool ResizeMatrix(SPISRAMFLOATARRAY &X, size_t Nrow, size_t Ncol, const size_t NewNrow, const size_t NewNcol);
void ResizeMatrix(SPISRAMFLOATARRAY &X, size_t Nrow, size_t Ncol, const size_t NewNrow, const size_t NewNcol, size_t &FinalNrow, size_t &FinalNcol);
void ConcatenateMatrixToRight(SPISRAMFLOATARRAY &A, size_t &NcolA, SPISRAMFLOATARRAY &B, size_t NcolB);
void ModifiedGramSchmidt(std::vector <SPISRAMFLOATARRAY> &V, size_t NumVectors, size_t VectorDim, std::vector <SPISRAMFLOATARRAY> &U);
void Orth(SPISRAMFLOATARRAY &A, size_t VectorDim, size_t &NumVectors, std::vector <SPISRAMClass *> &SRAMPtrList);
void ClasswisePCA(std::vector <SPISRAMFLOATARRAY> &Data, uint16_t NumClass, std::vector <size_t> &NumKept, size_t NumDim, float RetainThres, std::vector <SPISRAMFLOATARRAY> &CPA, std::vector <SPISRAMClass *> &SRAMPtrList);
void DisplaySRAMFloatUsage(std::vector <SPISRAMClass *> &SRAMPtrList);
void DisplaySRAMFloatContents(std::vector <SPISRAMClass *> &SRAMPtrList);
void MatrixSVD(SPISRAMFLOATARRAY &A, size_t m, size_t n, SPISRAMFLOATARRAY &w, SPISRAMFLOATARRAY &V, std::vector <SPISRAMClass *> &SRAMPtrList);
void MatrixMean(SPISRAMFLOATARRAY &X, const size_t Nrow, const size_t Ncol, SPISRAMFLOATARRAY &M, std::vector <SPISRAMClass *> &SRAMPtrList);
