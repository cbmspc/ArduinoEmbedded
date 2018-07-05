/*
Handles the SPISRAM chip, including writing/reading.
*/

#pragma once

#if defined(ARDUINO) && ARDUINO >= 100
	#include "arduino.h"
#else
	#include "WProgram.h"
#endif
#undef max
#undef min

#include <vector>
#include <SPI.h>
#include "pins_arduino.h"
#include "wiring_private.h"

namespace std {
	void __throw_bad_alloc();
	void __throw_length_error(char const*e);
}

class SPISRAMClass {
private:
	SPIClass * SPIPTR = 0;
	uint8_t _CS_PIN;
	uint8_t _MISO_PIN;
	uint8_t _MOSI_PIN;
	uint8_t _SCK_PIN;
	uint8_t _SPI_DATAMODE;
	SercomSpiTXPad _SERCOM_TXPAD;
	SercomRXPad _SERCOM_RXPAD;
	_EPioType _MISO_PIOTYPE;
	_EPioType _MOSI_PIOTYPE;
	_EPioType _SCK_PIOTYPE;
	bool _HasPinControl;
	size_t _SRAMSIZEBYTES;

	void SelectDevice();
	void UnselectDevice();
	void EnableWriteSRAM();
	void ResetRegister();
	void Reverse_Endian(void* A, void* B, size_t Length);

	//void DefaultWrite(uint8_t *Buffer, size_t DataLength, uint32_t Address);
	//void DefaultRead(uint8_t *Buffer, size_t DataLength, uint32_t Address);
	bool FastWrite(void *Buffer, size_t DataLength, uint32_t Address);
	bool FastRead(void *Buffer, size_t DataLength, uint32_t Address);

	bool _SequentialWriteOpened;
	uint32_t _SequentialWriteBytesWritten;
	uint32_t _SequentialWriteNextAddress;

	bool _SequentialReadOpened;
	uint32_t _SequentialReadBytesRead;
	uint32_t _SequentialReadNextAddress;


	// This proxy is used by the bracket operator access in public
	class proxy {
	private:
		size_t _index;
		SPISRAMClass *_parentPtr;
	public:
		proxy(size_t i, SPISRAMClass *parentPtr);
		operator uint8_t() const;
		void operator=(uint8_t rhs);
		void operator=(proxy rhs);
	};






public:
	SPISRAMClass();
	~SPISRAMClass();

	// This merely sets up the class, but does not assume pin control.
	void Initialize(SERCOM *SPISRAM_SERCOM,
		uint8_t SPISRAM_CS_PIN, uint8_t SPISRAM_MISO_PIN, uint8_t SPISRAM_SCK_PIN, uint8_t SPISRAM_MOSI_PIN,
		uint8_t SPI_DATAMODE,
		SercomSpiTXPad SPISRAM_SERCOM_TXPAD, SercomRXPad SPISRAM_SERCOM_RXPAD,
		_EPioType SPISRAM_MISO_PIOTYPE, _EPioType SPISRAM_SCK_PIOTYPE, _EPioType SPISRAM_MOSI_PIOTYPE,
		size_t SRAMSIZEBYTES);


/*	
	The proper procedure for handing off pin control is:
	1. First core: ReleasePinControl -- set pins to input or input_pullup
	2. First core sends a core-to-core command to second core
	3. Second core: AssumePinControl -- reconfigures mux to SPI mode and begins SPI
*/
	//void BeginPinControlTransferProcedure();
	void AssumePinControl();
	void ReleasePinControl();

	bool get_JEDEC_ID(uint8_t *Response, size_t Length);
	bool get_SRAM_Mode_Register(uint8_t &Response);

	bool FillSRAM(uint8_t Fill);

	inline uint32_t get_SRAM_SIZE() {
		return _SRAMSIZEBYTES;
	}

	/* Conduct a comprehensive self test. */
	uint8_t SelfTest();

	bool isReady();

	bool SeqWriteOpen(uint32_t StartAddress);
	uint32_t SeqWrite(void *Buf, size_t BufferLengthBytes);
	uint32_t get_SeqBytesWritten();
	void SeqWriteClose();

	bool SeqReadOpen(uint32_t StartAddress);
	uint32_t SeqRead(void *Buf, size_t BufferLengthBytes);
	uint32_t get_SeqBytesRead();
	void SeqReadClose();


	/* Bracket operator access */
	proxy operator[](uint32_t i);   // Indexed read and write

	/* Calculate a simple 16-bit hash value for the SRAM content. This function will not reset Sequential read position. */
	uint16_t CalculateHash(uint32_t StartAddress, uint32_t DataLength);

	/* Support for SPISRAMFLOATARRAY */
	float ReadFloat(uint32_t StartAddress);
	void WriteFloat(float Data, uint32_t StartAddress);

	// Error codes for the SelfTest
	enum ERRORCODES {
		OK,
		NOT_INITIALIZED,
		NO_PIN_CONTROL,
		INCORRECT_JEDEC_ID,
		INCORRECT_STATUS_REGISTER,
		WRITE_CHANNEL_IS_OPEN,
		READ_CHANNEL_IS_OPEN,
		INCORRECT_READ_TOTAL,
		INCORRECT_WRITE_TOTAL,
		READBACK_ERROR
	};

};


/* Allows bracket access that directly accepts or returns float 
Warning:	Arrays like this do not have the concepts of rows and columns like a matrix. 
			Instead, everything is 1D, and rows/columns are left to the caller to handle.
*/
class SPISRAMFLOATARRAY {
public:
	struct ATracker {
		SPISRAMFLOATARRAY *ObjPtr;
		SPISRAMClass *SRAMPtr;
		uint32_t StartAddress;
		uint32_t Bytes;
		bool operator<(const ATracker &rhs) const;   // Comparison operator for ascend sorting by StartAddress
		//bool operator>(const ATracker &rhs) const;   // Comparison operator for descend sorting by array size
	};
	struct SRAMRange {
		SPISRAMClass *SRAMPtr;
		uint32_t AllowedStartAddress;
		uint32_t AllowedEndAddress;
	};
private:
	static std::vector <ATracker> ATrack; // Tracks the SRAM Pointers, start address, and length of each class object
	static std::vector <SRAMRange> RTrack; // Tracks the allowable address ranges for this class
	uint32_t _StartAddress;
	uint32_t _Length;  // Length is for float. Each float is 4 bytes. ATracker tracks the RAM usage, and Length tracks bracket index boundary
	SPISRAMClass *_SRAMPtr;
	bool AttemptAllocate(uint32_t ArrayLength, SPISRAMClass *SRAMPtr);
	void Allocate(uint32_t ArrayLength, uint32_t ArrayBytes, uint32_t ArrayStartAddress, SPISRAMClass *SRAMPtr);
	void Deallocate(uint32_t ArrayStartAddress, SPISRAMClass *SRAMPtr);
	bool DecreaseAllocation(uint32_t NewArrayLength);
	bool IncreaseAllocationInPlace(uint32_t NewArrayLength);
	bool Relocate(uint32_t NewArrayLength, SPISRAMClass *NewSRAMPtr);
	//static bool TestPinControl(SPISRAMClass *SRAMPtr, uint32_t TestAddress);

	class proxy {
	private:
		size_t _index;
		SPISRAMFLOATARRAY *_parentPtr;
	public:
		proxy(size_t i, SPISRAMFLOATARRAY *parentPtr);
		operator float() const; // This allows the proxy to be casted to float

		// Float math operators
		float operator=(float rhs);
		float operator+=(float rhs);
		float operator-=(float rhs);
		float operator*=(float rhs);
		float operator/=(float rhs);
		proxy operator=(proxy rhs);
		proxy operator+=(proxy rhs);
		proxy operator-=(proxy rhs);
		proxy operator*=(proxy rhs);
		proxy operator/=(proxy rhs);

		// Increment/decrement operators
		proxy& operator++(); // Prefix increment
		proxy operator++(int); // Postfix increment
		proxy& operator--(); // Prefix decrement
		proxy operator--(int); // Postfix decrement

	};


public:
	SPISRAMFLOATARRAY();
	SPISRAMFLOATARRAY(uint32_t ArrayLength, SPISRAMClass *SRAMPtr);
	//SPISRAMFLOATARRAY(uint32_t ArrayLength, SPISRAMClass **SRAMPtrList, const uint8_t NumSRAMPtrs);
	SPISRAMFLOATARRAY(uint32_t ArrayLength, std::vector <SPISRAMClass *> &SRAMPtrList);
	SPISRAMFLOATARRAY(uint32_t ArrayLength, SPISRAMClass *SRAMPtr, uint32_t ManualStartAddress);
	SPISRAMFLOATARRAY(const SPISRAMFLOATARRAY &old_obj); // Copy constructor
	~SPISRAMFLOATARRAY();

	void manual_allocate(uint32_t ArrayLength, SPISRAMClass *SRAMPtr, uint32_t ManualStartAddress);
	static void GetSRAMRange(uint32_t &AllowedStartAddress, uint32_t &AllowedEndAddress, SPISRAMClass *SRAMPtr);
	static void SetSRAMRange(const uint32_t &AllowedStartAddress, const uint32_t &AllowedEndAddress, SPISRAMClass *SRAMPtr);
	//bool resize(uint32_t NewArrayLength, SPISRAMClass **SRAMPtrList, const uint8_t NumSRAMPtrs); // Resize the array, trying from a list of SRAM pointers
	bool resize(uint32_t NewArrayLength, std::vector <SPISRAMClass *> &SRAMPtrList); // Resize the array, trying from a list of SRAM pointers
	bool resize(uint32_t NewArrayLength, SPISRAMClass *NewSRAMPtr); // Resize the array, using the specific SRAM pointer if needed to relocate
	bool resize(uint32_t NewArrayLength); // Resize the array, using the same SRAM if needed to relocate
	uint32_t size();
	uint32_t address();
	SPISRAMClass* SRAMPtr();

	static const size_t GetTotalBytesAllocated(SPISRAMClass *SRAMPtr);
	static const size_t CountAllocatedArrays(SPISRAMClass *SRAMPtr);
	static void GetArrayObjPtrs(SPISRAMClass *SRAMPtr, std::vector <SPISRAMFLOATARRAY*> &List);
	static size_t Defrag(SPISRAMClass *SRAMPtr); // Consolidate within one chip
	//static size_t Defrag(SPISRAMClass **SRAMPtrList, const size_t NumSRAMPtrs); // Consolidate all arrays towards the first SRAM chip
	static size_t Defrag(std::vector <SPISRAMClass *> &SRAMPtrList); // Consolidate all arrays towards the first SRAM chip

	static bool TestPinControl(SPISRAMClass *SRAMPtr, uint32_t TestAddress);


	// Bracket access
	proxy operator[](uint32_t i);

	// Array operations
	SPISRAMFLOATARRAY& operator=(SPISRAMFLOATARRAY &rhs);  // LHS array becomes a copy of the RHS array
	float dotprod(SPISRAMFLOATARRAY &rhs);   // Dot product with array
	float mean(); // Arithmetic mean of all array elements
	float var(); // Variance of all array elements
	SPISRAMFLOATARRAY& operator*=(const float rhs); // Multiply by a scalar number
	SPISRAMFLOATARRAY& operator/=(const float rhs); // Divide by a scalar number


};

