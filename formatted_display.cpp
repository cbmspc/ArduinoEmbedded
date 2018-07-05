#include <Arduino.h>
#include "formatted_display.h"
#define Serial SerialUSB



void SerialPrintLeadingZeros(uint8_t Number) {
	if (Number < 128) Serial.print("0");
	if (Number < 64) Serial.print("0");
	if (Number < 32) Serial.print("0");
	if (Number < 16) Serial.print("0");
	if (Number < 8) Serial.print("0");
	if (Number < 4) Serial.print("0");
	if (Number < 2) Serial.print("0");
}

void SerialPrintPaddedHex(uint8_t Number) {
	if (Number < 16) Serial.print("0");
	Serial.print(Number, HEX);
}


void SerialPrintPaddedHexArray(void *Number, size_t Length) {
	uint8_t *Num = reinterpret_cast<uint8_t*>(Number);
	for (size_t i = 0; i < Length; i++) {
		SerialPrintPaddedHex(Num[i]);
		if (i != Length - 1) {
			Serial.print(" ");
		}
	}
}

void SerialPrintPaddedBin(uint8_t Number) {
	SerialPrintLeadingZeros(Number);
	Serial.print(Number, BIN);
}


void SerialPrintFormattedHex(uint16_t Value, size_t Index) {
	SerialPrintFormattedHex((uint8_t)(Value >> 8), Index);
	SerialPrintFormattedHex((uint8_t)Value, Index + 1);
}

void SerialPrintFormattedHex(uint8_t Value, size_t Index) {
	if (Index % 32 == 0) {
		if (Index < 0x10000000) Serial.print("0");
		if (Index < 0x1000000) Serial.print("0");
		if (Index < 0x100000) Serial.print("0");
		if (Index < 0x10000) Serial.print("0");
		if (Index < 0x1000) Serial.print("0");
		if (Index < 0x100) Serial.print("0");
		if (Index < 0x10) Serial.print("0");
		Serial.print(Index, HEX);
		Serial.print("h");
		Serial.print("  ");
	}
	if (Value < 16) Serial.print("0");
	Serial.print(Value, HEX);
	Serial.print(" ");
	if ((Index + 1) % 32 == 0) {
		Serial.print("  ");
		if (Index < 0x10000000) Serial.print("0");
		if (Index < 0x1000000) Serial.print("0");
		if (Index < 0x100000) Serial.print("0");
		if (Index < 0x10000) Serial.print("0");
		if (Index < 0x1000) Serial.print("0");
		if (Index < 0x100) Serial.print("0");
		if (Index < 0x10) Serial.print("0");
		Serial.print(Index, HEX);
		Serial.println("h");
	}
	else if ((Index + 1) % 16 == 0) {
		Serial.print(" ");
	}
}



void printLLNumber(uint64_t n, uint8_t base)
{
	/* Written by robtillaart at http://forum.arduino.cc/index.php?topic=143584.0 */
	unsigned char buf[16 * sizeof(long)];
	unsigned int i = 0;

	if (n == 0)
	{
		Serial.print((char)'0');
		return;
	}

	while (n > 0)
	{
		buf[i++] = n % base;
		n /= base;
	}

	for (; i > 0; i--)
		Serial.print((char)(buf[i - 1] < 10 ?
			'0' + buf[i - 1] :
			'A' + buf[i - 1] - 10));
}




void Print2DMatrix(float **X, const size_t RowStart, const size_t ColStart, const size_t RowEnd, const size_t ColEnd)
{
	int k = 0;
	Serial.print("[");
	for (int i = RowStart; i <= RowEnd; i++) {
		for (int j = ColStart; j <= ColEnd; j++) {
			Serial.print(X[i][j], 6);
			if (j < ColEnd) {
				// Not yet the last number in a row
				Serial.print(", ");
			}
			k++;
			if (k % 20 == 0 && j < ColEnd) {
				// 20 numbers have been printed. If this is not the last number in a row, break a line
				Serial.println("... ");
			}
		}
		if (i < RowEnd) {
			// End of a row but not the end of matrix
			k = 0;
			Serial.println("; ");
		}
	}
	// End of matrix
	Serial.println("]; ");
}

void PrintArray(float *a, const size_t Start, const size_t End)
{
	Print2DMatrix(&a, 0, Start, 0, End);
}

void Print1DMatrix(float *X, const size_t  Ncol, const size_t RowStart, const size_t ColStart, const size_t RowEnd, const size_t ColEnd)
/* Print a flattened matrix (a matrix that's reshaped into an array) */
{
	int k = 0;
	Serial.print("[");
	for (int i = RowStart; i <= RowEnd; i++) {
		for (int j = ColStart; j <= ColEnd; j++) {
			Serial.print(X[i*Ncol + j], 6);
			if (j < ColEnd) {
				// Not yet the last number in a row
				Serial.print(", ");
			}
			k++;
			if (k % 20 == 0 && j < ColEnd) {
				// 20 numbers have been printed. If this is not the last number in a row, break a line
				Serial.println("... ");
			}
		}
		if (i < RowEnd) {
			// End of a row but not the end of matrix
			k = 0;
			Serial.println("; ");
		}
	}
	// End of matrix
	Serial.println("]; ");
}


void Print1DMatrix(SPISRAMFLOATARRAY &X, const size_t Ncol, const size_t RowStart, const size_t ColStart, const size_t RowEnd, const size_t ColEnd)
/* Print a flattened matrix (a matrix that's reshaped into an array) */
{
	if (X.size() == 0) {
		Serial.println("[Empty matrix]");
		return;
	}
	size_t k = 0;
	size_t q = 0;
	size_t qsize = X.size();
	Serial.print("[");
	for (size_t i = RowStart; i <= RowEnd; i++) {
		for (size_t j = ColStart; j <= ColEnd; j++) {
			if (q++ >= qsize) {
				Serial.println("Index out of bound]; ");
				return;
			}
			Serial.print(X[i*Ncol + j], 6);
			if (j < ColEnd) {
				// Not yet the last number in a row
				Serial.print(", ");
			}
			k++;
			if (k % 20 == 0 && j < ColEnd) {
				// 20 numbers have been printed. If this is not the last number in a row, break a line
				Serial.println("... ");
			}
		}
		if (i < RowEnd) {
			// End of a row but not the end of matrix
			k = 0;
			Serial.println("; ");
		}
	}
	// End of matrix
	Serial.println("]; ");
}

void Print1DMatrix(SPISRAMFLOATARRAY &X, const size_t Ncol) {
	// X knows its own size.
	size_t Nrow = X.size() / Ncol;
	Print1DMatrix(X, Ncol, 0, 0, Nrow - 1, Ncol - 1);
}

void Print1DMatrix(SPISRAMFLOATARRAY &X) {
	Print1DMatrix(X, X.size());
}

void Print1DMatrixMATLAB(SPISRAMFLOATARRAY &X, size_t Ncol, String VarName) {
	Serial.print(" ");
	Serial.print(VarName);
	Serial.print(" = ");
	Print1DMatrix(X, Ncol);
}

void Print1DMatrixArrayMATLAB(std::vector <SPISRAMFLOATARRAY> &X, std::vector <size_t> &Ncol, String VarName) {
	size_t Nrow;
	for (size_t i = 0; i < X.size(); i++) {
		Serial.print(" ");
		Serial.print(VarName);
		Serial.print("{");
		Serial.print(i + 1);
		Serial.print("} = ");
		Print1DMatrix(X[i], Ncol[i]);
	}
}

void Print1DMatrixArrayMATLAB(std::vector <SPISRAMFLOATARRAY> &X, size_t Ncol, String VarName) {
	size_t Nrow;
	for (size_t i = 0; i < X.size(); i++) {
		Serial.print(" ");
		Serial.print(VarName);
		Serial.print("{");
		Serial.print(i + 1);
		Serial.print("} = ");
		Print1DMatrix(X[i], Ncol);
	}
}