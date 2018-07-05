// Arduino Matrix Math using external memory

#include <math.h>
#include "matrix_math_ext.h"
#include "formatted_display.h"


void Determinant(SPISRAMFLOATARRAY &A, float &d, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	SPISRAMFLOATARRAY B(n*n, SRAMPtrList);
	SPISRAMFLOATARRAY D(n, SRAMPtrList);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			B[i*n + j] = A[i*n + j];
		}
	}

	LDLDecomp(B, D, n);

	d = B[0] * D[0] * B[0];
	for (size_t i = 1; i < n; i++) {
		d *= B[i*n + i] * D[i] * B[i*n + i];
	}

}

void CholeskyDecomp(SPISRAMFLOATARRAY &B, size_t n)
{
	//Cholesky Decomp
	//where B = L * L'
	//replace B with L

	for (size_t i = 0; i < n; i++) {
		B[i*n + i] = B[i*n + i];
		for (size_t j = 0; j < i; j++) {
			B[i*n + j] = B[i*n + j];
			for (size_t k = 0; k < j; k++) {
				B[i*n + j] -= B[i*n + k] * B[j*n + k];
			}
			B[i*n + j] /= B[j*n + j];
			B[j*n + i] = 0.0;
			B[i*n + i] -= pow(B[i*n + j],2.0);
		}
		B[i*n + i] = sqrt(B[i*n + i]);
	}
}

void LDLDecomp(SPISRAMFLOATARRAY &B, SPISRAMFLOATARRAY &D, size_t n)
{
	//LDL Decomp
	//where B = L * D * L'
	//replace B with L, D (vector) with D

	for (size_t i = 0; i < n; i++) {
		D[i] = B[i*n + i];
		for (size_t j = 0; j < i; j++) {
			B[i*n + j] = B[i*n + j];
			for (size_t k = 0; k < j; k++) {
				B[i*n + j] -= B[i*n + k] * D[k] * B[j*n + k];
			}
			B[i*n + j] /= D[j];
			B[j*n + i] = 0.0;
			D[i] -= B[i*n + j] * D[j] * B[i*n + j];
		}
		B[i*n + i] = 1.0;
	}
}

void LowerTriangularInverse(SPISRAMFLOATARRAY &L, size_t n)
{
	//Solve L * L^-1 = I
	//for L^-1 by
	//replacing L with L^-1

	//float Li[n][n];

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			L[i*n + j] = -L[i*n + j] * L[j*n + j];
			for (size_t k = j + 1; k < i; k++) {
				L[i*n + j] -= L[i*n + k] * L[k*n + j];
			}
			L[i*n + j] /= L[i*n + i];
			L[j*n + i] = 0.0;
		}
		L[i*n + i] = 1.0 / L[i*n + i];
	}
}

void MatrixMultiply(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &B, SPISRAMFLOATARRAY &C, size_t n, size_t m, size_t p)
{
	// C is n x p
	// A is n x m
	// B is m x p

	for (size_t i = 0; i<n; i++) {
		//single row from A
		for (size_t k = 0; k<p; k++) {
			float sum = 0;
			for (size_t j = 0; j < m; j++) {
				sum += A[i*m + j] * B[j*p + k];
			}
			C[i*p + k] = sum;
			//C[i*p + k] = A[i*m] * B[k];
			//for (size_t j = 1; j<m; j++) {
			//	C[i*p + k] += A[i*m + j] * B[j*p + k];
			//}
		}
	}
}

void MatrixTransposeSquare(SPISRAMFLOATARRAY &A, size_t n)
{
	// A is n x n

	float el;

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			el = A[i*n + j];
			A[i*n + j] = A[j*n + i];
			A[j*n + i] = el;
		}
	}
}

void MatrixTranspose(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &AT, size_t n, size_t m)
{
	// A = input matrix (n x m)
	// AT = output matrix = the transpose of AT (m x n)
	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < m; j++)
			AT[j*n + i] = A[i*m + j];
}

void EigSolve(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &eig, float &l, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList)
{

	// A is nxn, eig is nx1, l is 1x1
	int it;
	bool repeat;

	SPISRAMFLOATARRAY eig_tmp(n, SRAMPtrList);
	SPISRAMFLOATARRAY l_n(1, SRAMPtrList);
	SPISRAMFLOATARRAY l_d(1, SRAMPtrList);

	// initialize eigenvector
	float tmp = 1.0 / sqrt((float)n);
	for (size_t i = 0; i < n; i++) {
		eig[i] = tmp;
	}

	it = 0;
	do {
		it++;

		// calculate the matrix-by-vector product A*eig
		MatrixMultiply(A, eig, eig_tmp, n, n, 1);

		// calculate the length of the resultant vector
		Normalize(eig_tmp);

		// normalize b to unit vector for next iteration
		repeat = false;
		for (size_t i = 0; i < n; i++) {
			if (fabs(eig[i] - eig_tmp[i]) > 0.0001 && it<100) {
				repeat = true;
			}
			eig[i] = eig_tmp[i];
		}
	} while (repeat);

	MatrixMultiply(A, eig, eig_tmp, n, n, 1);

	MatrixMultiply(eig, eig_tmp, l_n, 1, n, 1);

	MatrixMultiply(eig, eig, l_d, 1, n, 1);
	
	l = l_n[0] / l_d[0];

}

void Normalize(SPISRAMFLOATARRAY &w)
{
	// Vector normalization of w
	float norm = Norm(w);
	size_t n = w.size();
	for (size_t i = 0; i < n; i++) {
		w[i] /= norm;
	}
}


void MatrixNormalize(SPISRAMFLOATARRAY &A, size_t Ncol)
{
	// Normalize each column in the matrix
	size_t Nrow = A.size() / Ncol;
	for (size_t j = 0; j < Ncol; j++) {
		float s = 0.0;
		for (size_t i = 0; i < Nrow; i++) {
			s += A[i*Ncol + j] * A[i*Ncol + j];
		}
		s = sqrt(s);
		for (size_t i = 0; i < Nrow; i++) {
			A[i*Ncol + j] /= s;
		}
	}
}


float Norm(SPISRAMFLOATARRAY &w)
{
	// Vector norm of w
	size_t n = w.size();
	float norm = 0.0;
	for (size_t i = 0; i < n; i++) {
		norm += w[i] * w[i];
	}
	norm = sqrt(norm);
	return norm;
}


void GeneralizedSelfAdjointEigSolver(SPISRAMFLOATARRAY &A, SPISRAMFLOATARRAY &B, SPISRAMFLOATARRAY &w, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	SPISRAMFLOATARRAY C(n*n, SRAMPtrList);
	SPISRAMFLOATARRAY w_tmp(n, SRAMPtrList);
	float l;

	CholeskyDecomp(B, n); // L where L * L^T = B

	LowerTriangularInverse(B, n); // L^-1

	MatrixMultiply(B, A, C, n, n, n); //L^-1 * A, output C

	MatrixTransposeSquare(B, n); // (L^-1)^T

	MatrixMultiply(C, B, A, n, n, n); //L^-1 * A * (L^-1)^T, output A

	EigSolve(A, w_tmp, l, n, SRAMPtrList);

	MatrixMultiply(B, w_tmp, w, n, n, 1); // eigvecs = (L^-1)^T * eigvecs

	Normalize(w);
}



void PCA(SPISRAMFLOATARRAY &X, SPISRAMFLOATARRAY &PA, size_t &nk, size_t NumSamples, size_t NumVariables, float RetainThres, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	/* Automatically chooses between small sample size vs. large sample size versions of PCA. 
	Note: Small-sample-size calculation currently doesn't work, since Colin's PCA requires a true covariance. Inputting dmdata*dmdata' will produce NaNs
	*/
	//if (NumSamples <= NumVariables) {
	//	PCA_SmallSampleSize(X, PA, nk, NumSamples, NumVariables, RetainThres, SRAMPtrList);
	//}
	//else {
		MatrixCovariance(X, NumSamples, NumVariables, PA, SRAMPtrList);
		PCA(PA, nk, NumVariables, RetainThres, SRAMPtrList, NumSamples - 1);
		ResizeMatrix(PA, NumVariables, NumVariables, NumVariables, nk);

		SerialUSB.print("NumVariables=");
		SerialUSB.println(NumVariables);
		SerialUSB.print("New PA size=");
		SerialUSB.println(PA.size());

	//}
}


void PCA_SmallSampleSize(SPISRAMFLOATARRAY &X, SPISRAMFLOATARRAY &PA, size_t &nk, size_t NumSamples, size_t NumVariables, float RetainThres, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	/*	Performs traditional PCA when the number of samples is equal to or less than the number of variables
	*/

	SPISRAMFLOATARRAY XT(NumSamples*NumVariables, SRAMPtrList);
	MatrixTranspose(X, XT, NumSamples, NumVariables);

	SPISRAMFLOATARRAY S;
	MatrixCovariance(XT, NumVariables, NumSamples, S, SRAMPtrList, 0);

	PCA(S, nk, NumSamples, RetainThres, SRAMPtrList);

	if (nk >= NumSamples) {
		nk = NumSamples - 1;
	}

	ResizeMatrix(S, NumSamples, NumSamples, NumSamples, nk);

	PA.resize(NumVariables*nk, SRAMPtrList);

	MatrixMultiply(XT, S, PA, NumVariables, NumSamples, nk);
	MatrixNormalize(PA, nk);
	

}


void PCA(SPISRAMFLOATARRAY &S, size_t &nk, size_t n, float RetainThres, std::vector <SPISRAMClass *> &SRAMPtrList, size_t nkmax)
{
	/*	Dimension Reduction using traditional Principal Component Analysis
		Given covariance matrix S = X'*X with dimension n, produce first nk principal axes that explain > RetainThres of variance
		S is overwritten by the matrix of principal axes
	*/

	SerialUSB.print(" Performing a PCA with dimension ");
	SerialUSB.print(n);
	SerialUSB.print(" RetainThres=");
	SerialUSB.print(RetainThres);
	SerialUSB.print(" nkmax=");
	SerialUSB.println(nkmax);
	
	float l, s;
	SPISRAMFLOATARRAY S_tmp(n*n, SRAMPtrList);
	SPISRAMFLOATARRAY w(n, SRAMPtrList);
	SPISRAMFLOATARRAY W(n*n, SRAMPtrList);
	SPISRAMFLOATARRAY WW(n*n, SRAMPtrList);
	//DisplaySRAMFloatUsage(SRAMPtrList);

	// initialize S_tmp
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			S_tmp[i*n + j] = S[i*n + j];
			S[i*n + j] = 0.0;
		}
	}

	nk = 0;
	s = 0.0;

	EigSolve(S_tmp, w, l, n, SRAMPtrList); // get first eigenvec/eigenval

	do
	{
		// fill nkth column of S with nkth principal axis

		//SerialUSB.println("Filling nkth column of S with nkth principal axis");
		//SerialUSB.print(" nk=");
		//SerialUSB.println(nk);
		//SerialUSB.print(" w=");
		//Print1DMatrix(w, n);

		for (size_t i = 0; i < n; i++) {
			S[i*n + nk] = w[i];
		}

		nk++;

		if (nk >= nkmax) {
			break;
		}

		// prepare S_tmp for next principal component

		MatrixMultiply(w, w, W, n, 1, n);
		MatrixMultiply(W, W, WW, n, n, n);
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				S_tmp[i*n + j] += l*WW[i*n + j] - 2*l*W[i*n + j];
			}
		}

		// check if enough variance is already explained by nk components, by getting nk+1 eigenvalue
		s += l;
		EigSolve(S_tmp, w, l, n, SRAMPtrList);

	} while ((nk < n) && (s / (s + (n-nk)*l) <= RetainThres));
}



// LDA wrapper for Colin-LDA.
void LDA(SPISRAMFLOATARRAY &X1, SPISRAMFLOATARRAY &X2, size_t NumDim, SPISRAMFLOATARRAY &w, std::vector <SPISRAMClass *> &SRAMPtrList) {
	size_t n1 = X1.size() / NumDim;
	size_t n2 = X2.size() / NumDim;
	SPISRAMFLOATARRAY mu1;
	SPISRAMFLOATARRAY mu2;
	SPISRAMFLOATARRAY S1;
	SPISRAMFLOATARRAY S2;
	w.resize(NumDim, SRAMPtrList);
	MatrixMean(X1, n1, NumDim, mu1, SRAMPtrList);
	MatrixMean(X2, n2, NumDim, mu2, SRAMPtrList);
	MatrixCovariance(X1, n1, NumDim, S1, SRAMPtrList);
	MatrixCovariance(X2, n2, NumDim, S2, SRAMPtrList);
	LDA(mu1, mu2, S1, S2, w, n1, n2, NumDim, SRAMPtrList);
}




void LDA(SPISRAMFLOATARRAY &mu1, SPISRAMFLOATARRAY &mu2, SPISRAMFLOATARRAY &S1, SPISRAMFLOATARRAY &S2, SPISRAMFLOATARRAY &w, size_t n1, size_t n2, size_t n, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	// DOES NOT change S1 and S2 (covariance matrices of class 1 and 2)
	// useful output is w
	// utilizes a PCA step to ensure nondegenerate matrix

	int ntot = n1 + n2;
	float p1 = (float)n1 / (float)ntot;
	float p2 = (float)n2 / (float)ntot;
	float d;
	SPISRAMFLOATARRAY mu(n, SRAMPtrList);
	SPISRAMFLOATARRAY Sw(n*n, SRAMPtrList);
	SPISRAMFLOATARRAY Sb(n*n, SRAMPtrList);

	// calculate overall mean
	for (size_t j = 0; j < n; j++) {
		mu[j] = p1*mu1[j] + p2*mu2[j];
	}

	// calculate within-class scattering matrix
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			Sw[i*n + j] = p1*S1[i*n + j] + p2*S2[i*n + j];
		}
	}

	// calculate between-class scattering matrix
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			Sb[i*n + j] = p1*mu1[i] * mu1[j] + p2*mu2[i] * mu2[j] - mu[i] * mu[j];
		}
	}


	// check if Sw is rank deficient
	// if it's not use lda
	// if it is run PCA
	Determinant(Sw, d, n, SRAMPtrList);
	if (d > .0000000001) {  // 10^-10
		// run lda
		GeneralizedSelfAdjointEigSolver(Sb, Sw, w, n, SRAMPtrList);
	}
	else {
		// run PCA
		SPISRAMFLOATARRAY S(n*n, SRAMPtrList);

		// calculate total covariance matrix (S)
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				S[i*n + j] = Sw[i*n + j] + Sb[i*n + j];
			}
		}

		size_t nk = 0;
		PCA(S, nk, n, 0.99, SRAMPtrList);

		SPISRAMFLOATARRAY K(n*nk, SRAMPtrList);
		SPISRAMFLOATARRAY KT(nk*n, SRAMPtrList);

		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < nk; j++) {
				K[i*nk + j] = S[i*n + j];
			}
		}

		MatrixTranspose(K, KT, n, nk);

		// project data to principal components before running lda
		SPISRAMFLOATARRAY Sk_tmp(n*nk, SRAMPtrList);
		SPISRAMFLOATARRAY Swk(nk*nk, SRAMPtrList);
		SPISRAMFLOATARRAY Sbk(nk*nk, SRAMPtrList);
		SPISRAMFLOATARRAY wk(nk, SRAMPtrList);

		MatrixMultiply(Sw, K, Sk_tmp, n, n, nk); //Swk
		MatrixMultiply(KT, Sk_tmp, Swk, nk, n, nk); //Swk

		MatrixMultiply(Sb, K, Sk_tmp, n, n, nk); //Sbk
		MatrixMultiply(KT, Sk_tmp, Sbk, nk, n, nk); //Sbk

													// run lda
		GeneralizedSelfAdjointEigSolver(Sbk, Swk, wk, nk, SRAMPtrList);

		// find projection from original space that incorporates 
		// principal components and lda vector
		MatrixMultiply(K, wk, w, n, nk, 1);

	}

	// scale w so that mu * w = 1
	// find s where mu * w = s s.t. mu * (w/s) = 1
	SPISRAMFLOATARRAY s(1, SRAMPtrList);
	MatrixMultiply(mu, w, s, 1, n, 1);
	for (size_t i = 0; i < n; i++) {
		w[i] /= s[0];
	}

}


// Calculate the mean value of each column in the matrix
// X is Nrow x Ncol
// M will be automatically resized to 1 x Ncol
void MatrixMean(SPISRAMFLOATARRAY &X, const size_t Nrow, const size_t Ncol, SPISRAMFLOATARRAY &M, std::vector <SPISRAMClass *> &SRAMPtrList) {
	M.resize(Ncol, SRAMPtrList);
	for (size_t j = 0; j < Ncol; j++) {
		M[j] = 0;
		for (size_t i = 0; i < Nrow; i++) {
			M[j] += X[i*Ncol + j];
		}
		M[j] /= Nrow;
	}
}


void MatrixDemean(SPISRAMFLOATARRAY &X, const size_t Nrow, const size_t Ncol, SPISRAMFLOATARRAY &Xdm, std::vector <SPISRAMClass *> &SRAMPtrList, bool DeMeanRowsInstead)
{
	// De-mean. Calculate mean for each column (or row, if DeMeanRowsInstead = true), subtract it from X, store in Xdm
	Xdm.resize(Nrow*Ncol, SRAMPtrList);
	if (DeMeanRowsInstead) {
		for (size_t i = 0; i < Nrow; i++) {
			float m = 0.0;
			for (size_t j = 0; j < Ncol; j++) {
				m += X[i*Ncol + j];
			}
			m /= Ncol;
			for (size_t j = 0; j < Ncol; j++) {
				Xdm[i*Ncol + j] = X[i*Ncol + j] - m;
			}
		}
	}
	else {
		for (size_t j = 0; j < Ncol; j++) {
			float m = 0.0;
			for (size_t i = 0; i < Nrow; i++) {
				m += X[i*Ncol + j];
			}
			m /= Nrow;
			for (size_t i = 0; i < Nrow; i++) {
				Xdm[i*Ncol + j] = X[i*Ncol + j] - m;
			}
		}
	}
}

void MatrixCovariance(SPISRAMFLOATARRAY &X, const size_t Nrow, const size_t Ncol, SPISRAMFLOATARRAY &S, std::vector <SPISRAMClass *> &SRAMPtrList, uint8_t SpecialOps)
/* Calculate the covariance of matrix X (Nrow by Ncol). S = transpose(X) * X
   Each row is a sample. Each column is a variable.
   Note: For real matrix only. Only calculates the upper triangle.

   S is automatically allocated/resized.

   SpecialOps is a set of bit switches
   0 = default covariance
   1 = de-mean in the other dimension (useful for small-sample-size PCA)
   2 = disable scaling by 1/(n-1)
   3 = both 1 and 2

*/
{

	// De-mean
	SPISRAMFLOATARRAY Xdm;
	bool DeMeanRowsInstead = false;
	if (SpecialOps & 1) {
		DeMeanRowsInstead = true;
	}
	MatrixDemean(X, Nrow, Ncol, Xdm, SRAMPtrList, DeMeanRowsInstead);

	// Calculate covariance, scaled to 1/(n-1)
	float nrm1 = Nrow - 1;
	if (SpecialOps & 2) {
		nrm1 = 1;
	}

	S.resize(Ncol*Ncol, SRAMPtrList);

	for (size_t i = 0; i < Ncol; i++) {
		for (size_t j = 0; j < Ncol; j++) {
			if (j >= i) {
				float sum = 0;
				for (size_t k = 0; k < Nrow; k++) {
					sum += Xdm[k*Ncol + i] * Xdm[k*Ncol + j] / nrm1;
				}
				S[i*Ncol + j] = sum;
			}
			else {
				S[i*Ncol + j] = S[j*Ncol + i];
			}
		}
	}
}



// Resize the matrix, anchor at the first element
// Note: Will not change Nrow and Ncol, but will instead return true or false
bool ResizeMatrix(SPISRAMFLOATARRAY &X, size_t Nrow, size_t Ncol, const size_t NewNrow, const size_t NewNcol)
{
	size_t FinalNrow;
	size_t FinalNcol;
	ResizeMatrix(X, Nrow, Ncol, NewNrow, NewNcol, FinalNrow, FinalNcol);
	if (NewNrow == FinalNrow && NewNcol == FinalNcol) { return true; }
	return false;
}

// Resize the matrix, anchor at the first element
void ResizeMatrix(SPISRAMFLOATARRAY &X, size_t Nrow, size_t Ncol, const size_t NewNrow, const size_t NewNcol, size_t &FinalNrow, size_t &FinalNcol)
{
	// Resize the rows first
	if (X.resize(NewNrow*Ncol)) {
		// Clean up by filling zeros to the expanded area
		if (NewNrow > Nrow) {
			for (size_t i = Nrow; i < NewNrow; i++) {
				for (size_t j = 0; j < Ncol; j++) {
					X[i*Ncol + j] = 0.0;
				}
			}
		}
		Nrow = NewNrow;
		FinalNrow = Nrow;
		//SerialUSB.println("Resized row");
	}
	else {
		//SerialUSB.println("Cannot resize row");
	}

	// Then resize the columns
	if (NewNcol > Ncol) {
		if (X.resize(Nrow*NewNcol)) {
			// For an increased size, we need to reshuffle data from the end first
			// First row never needs to be shuffled
			for (size_t i = Nrow - 1; i > 0; i--) {
				for (size_t j = Ncol - 1; j > 0; j--) {
					X[i*NewNcol + j] = X[i*Ncol + j];
					X[i*Ncol + j] = 0;
				}
				X[i*NewNcol] = X[i*Ncol];
				X[i*Ncol] = 0;
			}

			// Zero the expanded area
			for (size_t i = 0; i < Nrow; i++) {
				for (size_t j = Ncol; j < NewNcol; j++) {
					X[i*NewNcol + j] = 0;
				}
			}

			Ncol = NewNcol;
			FinalNcol = Ncol;
		}
	}
	else if (NewNcol < Ncol) {
		// For a decreased size, we need to reshuffle data from the start first
		// First row never needs to be shuffled
		// Then reduce the array allocation
		for (size_t i = 1; i < Nrow; i++) {
			for (size_t j = 0; j < Ncol; j++) {
				X[i*NewNcol + j] = X[i*Ncol + j];
				X[i*Ncol + j] = 0;
			}
		}
		Ncol = NewNcol;
		X.resize(Nrow*NewNcol);
		FinalNcol = Ncol;
	}
	//SerialUSB.println("Resized col");
}



// Classwise PCA
// CPA should have NumClass number of arrays, just like Data
// Data[0] should be (trial x dim) orientation (each row is a sample)
// Elements of CPA will be resized based on NumKept. For example, if CPCA reduces NumDim=80 to NumKept[0]=5 and NumKept[1]=3, then CPA[0] is 80 by 5 and CPA[1] is 80 by 3.
void ClasswisePCA(std::vector <SPISRAMFLOATARRAY> &Data, uint16_t NumClass, std::vector <size_t> &NumKept, size_t NumDim, float RetainThres, std::vector <SPISRAMFLOATARRAY> &CPA, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	SPISRAMFLOATARRAY Data_b(NumClass * NumDim, SRAMPtrList);
	std::vector <uint16_t> NumTrial(NumClass);
	uint16_t NumTrialTotal = 0;
	size_t nk;
	std::vector <SPISRAMFLOATARRAY> ClassMean;
	ClassMean.resize(2);
	SPISRAMFLOATARRAY TotalMean(NumDim, SRAMPtrList);
	CPA.resize(NumClass);
	NumKept.resize(NumClass);

	for (uint16_t j = 0; j < NumDim; j++) {
		TotalMean[j] = 0;
	}

	SerialUSB.println("CPCA: Calculating total mean and class mean.");
	for (uint16_t c = 0; c < NumClass; c++) {
		uint16_t s = Data[c].size() / NumDim;
		NumTrial[c] = s;
		NumTrialTotal += s;
		ClassMean[c].resize(NumDim, SRAMPtrList);
		
		// Calculate class mean
		for (uint16_t j = 0; j < NumDim; j++) {
			float sum = 0.0;
			for (uint16_t t = 0; t < NumTrial[c]; t++) {
				sum += Data[c][t*NumDim + j];
			}
			ClassMean[c][j] = sum / NumTrial[c];
			TotalMean[j] += sum;
		}
	}
	for (uint16_t j = 0; j < NumDim; j++) {
		TotalMean[j] /= NumTrialTotal;
	}
	

	// Calculate the between class principal axis

	for (uint16_t c = 0; c < NumClass; c++) {
		//Data_b
		SerialUSB.print("CPCA: Calculating between-class PA for class ");
		SerialUSB.println(c);

		float p = sqrt((float)NumTrial[c] / (float)NumTrialTotal);
		//SerialUSB.print(" p=");
		//SerialUSB.println(p);
		for (uint16_t j = 0; j < NumDim; j++) {
			Data_b[c*NumDim + j] = p * (ClassMean[c][j] - TotalMean[j]);
		}
	}

	// TotalMean and ClassMean are no longer needed
	//DisplaySRAMFloatUsage(SRAMPtrList);
	SerialUSB.println("Releasing TotalMean and ClassMean");
	TotalMean.resize(0);
	for (uint16_t c = 0; c < NumClass; c++) {
		ClassMean[c].resize(0);
	}

	//SerialUSB.print(" Data_b = ");
	//Print1DMatrix(Data_b, NumDim);

	// Data_b itself is correct. pca(Data_b) is a problem if want to use shortcut
	

	//DisplaySRAMFloatUsage(SRAMPtrList);

	SerialUSB.println("CPCA: Calculating total pca(Data_b)");
	SPISRAMFLOATARRAY S_b; // S_b will be resized by PCA
	PCA(Data_b, S_b, nk, NumClass, NumDim, 999.999, SRAMPtrList);

	SerialUSB.print(" nk=");
	SerialUSB.println(nk);


	SerialUSB.print(" S_b = ");
	Print1DMatrix(S_b, nk);

	//DisplaySRAMFloatUsage(SRAMPtrList);

	// Within class covariance and PCA
	//std::vector <SPISRAMFLOATARRAY> S_w; S_w.resize(2);
	// Use CPA as both working storage and final output
	for (uint16_t c = 0; c < NumClass; c++) {
		SerialUSB.print("CPCA: Calculating within-class PA for class ");
		SerialUSB.println(c);
		//DisplaySRAMFloatUsage(SRAMPtrList);

		MatrixCovariance(Data[c], NumTrial[c], NumDim, CPA[c], SRAMPtrList);
		PCA(CPA[c], NumKept[c], NumDim, RetainThres, SRAMPtrList);
		ResizeMatrix(CPA[c], NumDim, NumDim, NumDim, NumKept[c]);
		//DisplaySRAMFloatUsage(SRAMPtrList);

		SerialUSB.print(" NumKept{");
		SerialUSB.print(c);
		SerialUSB.print("}=");
		SerialUSB.println(NumKept[c]);
	}

	Print1DMatrixArrayMATLAB(CPA, NumKept, "S_w");

	
	// Add the between-class PA to each within-class PA

	SerialUSB.println("CPCA: Combining between-class PA and within-class PA.");

	for (uint16_t c = 0; c < NumClass; c++) {
		ConcatenateMatrixToRight(CPA[c], NumKept[c], S_b, nk);
	}

	// Orthonormalize the concatenated matrices

	SerialUSB.println("CPCA: Orth. ");
	for (uint16_t c = 0; c < NumClass; c++) {

		SerialUSB.print("Orthonormalizing CPA[");
		SerialUSB.print(c);
		SerialUSB.print("] with NumDim=");
		SerialUSB.print(NumDim);
		SerialUSB.print(", NumKept[");
		SerialUSB.print(c);
		SerialUSB.print("]=");
		SerialUSB.println(NumKept[c]);
		//TODO: Orth is broken

		Orth(CPA[c], NumDim, NumKept[c], SRAMPtrList);
	}

	//delete S_w;
	//delete ClassMean;
	//delete NumTrial;
}


// Resize and concatenate matrix A with matrix B, saving result to A.
// NcolA will be updated to reflect the new number of columns if the concatenation is successful
void ConcatenateMatrixToRight(SPISRAMFLOATARRAY &A, size_t &NcolA, SPISRAMFLOATARRAY &B, size_t NcolB)
{
	size_t NrowA = A.size() / NcolA;
	size_t NrowB = B.size() / NcolB;
	if (NrowA != NrowB) {
		return;
	}
	size_t OldNcolA = NcolA;
	ResizeMatrix(A, NrowA, NcolA, NrowA, NcolA + NcolB, NrowA, NcolA);
	if (NcolA > OldNcolA) {
		// resize successful
		// start copying
		for (size_t i = 0; i < NrowA; i++) {
			for (size_t j = 0; j < NcolB; j++) {
				A[i*NcolA + OldNcolA + j] = B[i*NcolB + j];
			}
		}
	}
	else {
		SerialUSB.println("Resize failure in Concatenation");
	}
}



// Compatibility wrapper to do orthonormalization on a matrix of tall vectors
// In this case, the dimension is Nrow, and the number of basis vectors is Ncol
// Input array will be overwritten with the matrix formed by the tall orthonormal basis vectors
// Null vectors will not be included
// Possible bug: During testing, Orth lost two vectors.
void Orth(SPISRAMFLOATARRAY &A, size_t VectorDim, size_t &NumVectors, std::vector <SPISRAMClass *> &SRAMPtrList)
{
	// Using SVD:
	SPISRAMFLOATARRAY w;
	SPISRAMFLOATARRAY V;
	MatrixSVD(A, VectorDim, NumVectors, w, V, SRAMPtrList);
	size_t r;
	for (size_t i = w.size()-1; i-- > 0; ) {
		SerialUSB.print("w[");
		SerialUSB.print(i);
		SerialUSB.print("]=");
		SerialUSB.println(w[i]);
		if (w[i] > 1e-6) {
			r = i;
			break;
		}
	}
	ResizeMatrix(A, VectorDim, NumVectors, VectorDim, r);
	NumVectors = r;

	// Using Gram Schmidt:
	//// Create temp arrays
	//std::vector <SPISRAMFLOATARRAY> V(NumVectors);
	//std::vector <SPISRAMFLOATARRAY> U(NumVectors);
	//for (size_t j = 0; j < NumVectors; j++) {
	//	V[j].resize(VectorDim, SRAMPtrList);
	//	U[j].resize(VectorDim, SRAMPtrList);
	//	for (size_t i = 0; i < VectorDim; i++) {
	//		V[j][i] = A[i*NumVectors + j];
	//	}
	//}

	//ModifiedGramSchmidt(V, NumVectors, VectorDim, U);


	//size_t n = 0;
	//size_t *keep = new size_t[NumVectors];

	//

	//for (size_t j = 0; j < NumVectors; j++) {
	//	//SerialUSB.print(" Uj=");
	//	//Print1DMatrix(U[j], VectorDim);
	//	//SerialUSB.print(" Norm(U[");
	//	//SerialUSB.print(j);
	//	//SerialUSB.print("])=");
	//	//SerialUSB.println(Norm(U[j]));
	//	if (Norm(U[j]) > 0.99) {
	//		keep[j] = n;
	//		n++;
	//	}
	//	else {
	//		keep[j] = n;
	//	}
	//}

	////SerialUSB.print(" n=");
	////SerialUSB.println(n);

	//// Update A
	//for (size_t i = 0; i < VectorDim; i++) {
	//	for (size_t j = 0; j < NumVectors; j++) {
	//		A[i*n + keep[j]] = U[j][i];
	//	}
	//}
	//A.resize(n*VectorDim);
	//NumVectors = n;
	//delete keep;
	//delete U;
	//delete V;
}




// Perform the modified Gram-Schmidt orthonormalization process
// based on http://web.mst.edu/~hilgers/classes/CS328/notes/modgs/node2.html
// Requirement: V must have rank == VectorDim, otherwise unexpected or incorrect result results
void ModifiedGramSchmidt(std::vector <SPISRAMFLOATARRAY> &V, size_t NumVectors, size_t VectorDim, std::vector <SPISRAMFLOATARRAY> &U)
{
	U[0] = V[0];
	Normalize(U[0]);

	for (size_t j = 1; j < NumVectors; j++) {
		float vu = V[j].dotprod(U[0]);
		for (size_t i = 0; i < VectorDim; i++) {
			U[j][i] = V[j][i] - vu * U[0][i];
		}
	}

	for (size_t k = 1; k < NumVectors - 1; k++) {
		Normalize(U[k]);
		for (size_t j = k + 1; j < NumVectors; j++) {
			float uu = U[j].dotprod(U[k]);
			for (size_t i = 0; i < VectorDim; i++) {
				U[j][i] = U[j][i] - uu * U[k][i];
			}
		}
	}

}


void DisplaySRAMFloatUsage(std::vector <SPISRAMClass *> &SRAMPtrList) {
	size_t b, c;
	for (uint8_t s = 0; s < SRAMPtrList.size(); s++) {
		b = SPISRAMFLOATARRAY::GetTotalBytesAllocated(SRAMPtrList[s]);
		c = SPISRAMFLOATARRAY::CountAllocatedArrays(SRAMPtrList[s]);
		SerialUSB.print("SRAM ");
		SerialUSB.print((uint32_t)SRAMPtrList[s]);
		SerialUSB.print(": ");
		SerialUSB.print(b);
		SerialUSB.print(" bytes allocated in ");
		SerialUSB.print(c);
		SerialUSB.println(" arrays.");
	}
}


void DisplaySRAMFloatContents(std::vector <SPISRAMClass *> &SRAMPtrList) {
	for (uint8_t s = 0; s < SRAMPtrList.size(); s++) {
		std::vector <SPISRAMFLOATARRAY*> List;
		SPISRAMFLOATARRAY::GetArrayObjPtrs(SRAMPtrList[s], List);
		SerialUSB.print("SRAM ");
		SerialUSB.print((uint32_t)SRAMPtrList[s]);
		SerialUSB.print(" has ");
		SerialUSB.print(List.size());
		SerialUSB.println(" arrays.");

		for (std::vector<SPISRAMFLOATARRAY*>::iterator it = List.begin(); it < List.end(); it++) {
			SerialUSB.print("Array: Address=");
			SerialUSB.print((*it)->address());
			SerialUSB.print(" Size=");
			SerialUSB.print((*it)->size());
			SerialUSB.print(" Content=");
			for (size_t i = 0; i < (*it)->size(); i++) {
				SerialUSB.print((**it)[i], 6);
				SerialUSB.print(", ");
			}
			SerialUSB.println("");
		}
	}
}






static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))



void MatrixSVD(SPISRAMFLOATARRAY &A, size_t m, size_t n, SPISRAMFLOATARRAY &w, SPISRAMFLOATARRAY &V, std::vector <SPISRAMClass *> &SRAMPtrList)
/* Singular Value Decomposition
Adapted from Numerical Recipes in C: The Art of Scientific Computing, chapter 2.6, pages 59--70
Given a matrix A[0..m-1][0..n-1], this routine computes its singular value decomposition, A = U·W·V'.
The matrix U replaces A on output.
The diagonal matrix of singular values w is output as a vector w[0..n-1].
The matrix V (not the transpose V') is output as V[0..n-1][0..n-1].
*/
{
	float pythag(float a, float b);
	int flag, i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z;
	//SPISRAMFLOATARRAY rv1(n, SRAMPtrList);
	std::vector <float> rv1(n);  //float *rv1 = new float[n];
	w.resize(n, SRAMPtrList);
	V.resize(n*n, SRAMPtrList);
	g = scale = anorm = 0.0; //Householder reduction to bidiagonal form.
	for (i = 1;i <= n;i++) {
		l = i + 1;
		rv1[i - 1] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i;k <= m;k++) scale += fabs(A[(k - 1)*n + i - 1]);
			if (scale) {
				for (k = i;k <= m;k++) {
					A[(k - 1)*n + i - 1] /= scale;
					s += A[(k - 1)*n + i - 1] * A[(k - 1)*n + i - 1];
				}
				f = A[(i - 1)*n + i - 1];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				A[(i - 1)*n + i - 1] = f - g;
				for (j = l;j <= n;j++) {
					for (s = 0.0, k = i;k <= m;k++) s += A[(k - 1)*n + i - 1] * A[(k - 1)*n + j - 1];
					f = s / h;
					for (k = i;k <= m;k++) A[(k - 1)*n + j - 1] += f * A[(k - 1)*n + i - 1];
				}
				for (k = i;k <= m;k++) A[(k - 1)*n + i - 1] *= scale;
			}
		}
		w[i - 1] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l;k <= n;k++) scale += fabs(A[(i - 1)*n + k - 1]);
			if (scale) {
				for (k = l;k <= n;k++) {
					A[(i - 1)*n + k - 1] /= scale;
					s += A[(i - 1)*n + k - 1] * A[(i - 1)*n + k - 1];
				}
				f = A[(i - 1)*n + l - 1];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				A[(i - 1)*n + l - 1] = f - g;
				for (k = l;k <= n;k++) rv1[k - 1] = A[(i - 1)*n + k - 1] / h;
				for (j = l;j <= m;j++) {
					for (s = 0.0, k = l;k <= n;k++) s += A[(j - 1)*n + k - 1] * A[(i - 1)*n + k - 1];
					for (k = l;k <= n;k++) A[(j - 1)*n + k - 1] += s * rv1[k - 1];
				}
				for (k = l;k <= n;k++) A[(i - 1)*n + k - 1] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i - 1]) + fabs(rv1[i - 1])));
	}
	for (i = n;i >= 1;i--) {
		//Accumulation of right-hand transformations.
		if (i < n) {
			if (g) {
				for (j = l;j <= n;j++) //Double division to avoid possible underflow.
					V[(j - 1)*n + i - 1] = (A[(i - 1)*n + j - 1] / A[(i - 1)*n + l - 1]) / g;
				for (j = l;j <= n;j++) {
					for (s = 0.0, k = l;k <= n;k++) s += A[(i - 1)*n + k - 1] * V[(k - 1)*n + j - 1];
					for (k = l;k <= n;k++) V[(k - 1)*n + j - 1] += s * V[(k - 1)*n + i - 1];
				}
			}
			for (j = l;j <= n;j++) {
				V[(j - 1)*n + i - 1] = 0.0;
				V[(i - 1)*n + j - 1] = 0.0;
			}
		}
		V[(i - 1)*n + i - 1] = 1.0;
		g = rv1[i - 1];
		l = i;
	}
	for (i = IMIN(m, n);i >= 1;i--) {
		//Accumulation of left-hand transformations.
		l = i + 1;
		g = w[i - 1];
		for (j = l;j <= n;j++) A[(i - 1)*n + j - 1] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l;j <= n;j++) {
				for (s = 0.0, k = l;k <= m;k++) s += A[(k - 1)*n + i - 1] * A[(k - 1)*n + j - 1];
				f = (s / A[(i - 1)*n + i - 1])*g;
				for (k = i;k <= m;k++) A[(k - 1)*n + j - 1] += f * A[(k - 1)*n + i - 1];
			}
			for (j = i;j <= m;j++) A[(j - 1)*n + i - 1] *= g;
		}
		else for (j = i;j <= m;j++) A[(j - 1)*n + i - 1] = 0.0;
		A[(i - 1)*n + i - 1] += 1.0;
	}
	for (k = n;k >= 1;k--) {
		//Diagonalization of the bidiagonal form : Loop over singular values, and over allowed iterations.
		for (its = 1;its <= 30;its++) {
			flag = 1;
			for (l = k;l >= 1;l--) {
				//Test for splitting.
				nm = l - 1; //Note that rv1[1-1] is always zero.
				if ((float)(fabs(rv1[l - 1]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((float)(fabs(w[nm - 1]) + anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0; //Cancellation of rv1[l-1], if l > 1.
				s = 1.0;
				for (i = l;i <= k;i++) {
					f = s * rv1[i - 1];
					rv1[i - 1] = c * rv1[i - 1];
					if ((float)(fabs(f) + anorm) == anorm) break;
					g = w[i - 1];
					h = pythag(f, g);
					w[i - 1] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1;j <= m;j++) {
						y = A[(j - 1)*n + nm - 1];
						z = A[(j - 1)*n + i - 1];
						A[(j - 1)*n + nm - 1] = y * c + z * s;
						A[(j - 1)*n + i - 1] = z * c - y * s;
					}
				}
			}
			z = w[k - 1];
			if (l == k) {
				//Convergence.
				if (z < 0.0) {
					//Singular value is made nonnegative.
					w[k - 1] = -z;
					for (j = 1;j <= n;j++) V[(j - 1)*n + k - 1] = -V[(j - 1)*n + k - 1];
				}
				break;
			}
			if (its == 30) SerialUSB.println("no convergence in 30 svdcmp iterations");
			x = w[l - 1]; //Shift from bottom 2 - by - 2 minor.
			nm = k - 1;
			y = w[nm - 1];
			g = rv1[nm - 1];
			h = rv1[k - 1];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, 1.0);
			f = ((x - z)*(x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0; //Next QR transformation :
			for (j = l;j <= nm;j++) {
				i = j + 1;
				g = rv1[i - 1];
				y = w[i - 1];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j - 1] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1;jj <= n;jj++) {
					x = V[(jj - 1)*n + j - 1];
					z = V[(jj - 1)*n + i - 1];
					V[(jj - 1)*n + j - 1] = x * c + z * s;
					V[(jj - 1)*n + i - 1] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j - 1] = z; //Rotation can be arbitrary if z = 0.
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1;jj <= m;jj++) {
					y = A[(jj - 1)*n + j - 1];
					z = A[(jj - 1)*n + i - 1];
					A[(jj - 1)*n + j - 1] = y * c + z * s;
					A[(jj - 1)*n + i - 1] = z * c - y * s;
				}
			}
			rv1[l - 1] = 0.0;
			rv1[k - 1] = f;
			w[k - 1] = x;
		}
	}
	//delete rv1;
}

float pythag(float a, float b)
/* Computes(a2 + b2)^(1/2) without destructive underflow or overflow. */
{
	float absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

