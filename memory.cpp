// 
// 
// 

#include <Arduino.h>
#undef max
#undef min
#include "memory.h"
#include <algorithm>
#include <functional>

#define OOB_RETURN_VALUE NAN

#define SRAM_CLOCK_HZ 12000000
#define SELFTEST_SIZE 512    // Must be at least 2 and also divides _SRAMSIZEBYTES. There will be 3 variables using this many bytes

#define BUFFERSIZE 2048

// These op codes should be standard among all compatible SPI RAM modules
#define OpCodeWriteStatusRegister	0x01
#define OpCodeWrite					0x02
#define OpCodeRead					0x03
#define OpCodeReadStatusRegister	0x05
#define OpCodeWriteEnable			0x06
#define OpCodeGetDeviceID			0x9f



// This fix is needed for Arduino to use std::vector with SPISRAMFLOATARRAY
namespace std {
	void __throw_bad_alloc() {
		SerialUSB.println("Unable to allocate memory");
	}
	void __throw_length_error(char const*e) {
		SerialUSB.print("Length Error :"); SerialUSB.println(e);
	}
}



SPISRAMClass::SPISRAMClass() {
	/* Dummy initializer. Must call Initialize before doing anything else to the SRAM. */
	_SRAMSIZEBYTES = 0;
	_HasPinControl = false;
}

void SPISRAMClass::Initialize(SERCOM *SPISRAM_SERCOM,
	uint8_t SPISRAM_CS_PIN, uint8_t SPISRAM_MISO_PIN, uint8_t SPISRAM_SCK_PIN, uint8_t SPISRAM_MOSI_PIN,
	uint8_t SPI_DATAMODE,
	SercomSpiTXPad SPISRAM_SERCOM_TXPAD, SercomRXPad SPISRAM_SERCOM_RXPAD,
	_EPioType SPISRAM_MISO_PIOTYPE, _EPioType SPISRAM_SCK_PIOTYPE, _EPioType SPISRAM_MOSI_PIOTYPE,
	size_t SRAMSIZEBYTES) {

	if (_SRAMSIZEBYTES || !SRAMSIZEBYTES) { return; } // Prevents repetitively initializing

													  // MISO, SCK, MOSI
	_CS_PIN = SPISRAM_CS_PIN;
	_MISO_PIN = SPISRAM_MISO_PIN;
	_MOSI_PIN = SPISRAM_MOSI_PIN;
	_SCK_PIN = SPISRAM_SCK_PIN;
	_SPI_DATAMODE = SPI_DATAMODE;
	_SERCOM_TXPAD = SPISRAM_SERCOM_TXPAD;
	_SERCOM_RXPAD = SPISRAM_SERCOM_RXPAD;
	_MISO_PIOTYPE = SPISRAM_MISO_PIOTYPE;
	_MOSI_PIOTYPE = SPISRAM_MOSI_PIOTYPE;
	_SCK_PIOTYPE = SPISRAM_SCK_PIOTYPE;
	_SRAMSIZEBYTES = SRAMSIZEBYTES;

	delay(10); //Delay is required for proper SRAM initialization (empirically tested)

	SPIPTR = new SPIClass(SPISRAM_SERCOM, _MISO_PIN, _SCK_PIN, _MOSI_PIN, _SERCOM_TXPAD, _SERCOM_RXPAD);

	// Moved the pin modes to AssumePinControl



	ReleasePinControl(); // By default, set the pins to input mode until we decide to gain pin control

}

SPISRAMClass::~SPISRAMClass() {
	if (_HasPinControl) {
		UnselectDevice();
	}
	_SRAMSIZEBYTES = 0;
	_HasPinControl = false;
	(*SPIPTR).end();
	delete SPIPTR;
}


bool SPISRAMClass::isReady() {
	if (!_SRAMSIZEBYTES) { return false; } // Not initialized
	if (!_HasPinControl) { return false; } // No control
	return true;
}

uint8_t SPISRAMClass::SelfTest() {
	if (!_SRAMSIZEBYTES) { return ERRORCODES::NOT_INITIALIZED; } // Not initialized
	if (!_HasPinControl) { return ERRORCODES::NO_PIN_CONTROL; } // No control
	if (_SequentialWriteOpened) {
		return ERRORCODES::WRITE_CHANNEL_IS_OPEN;
	}
	if (_SequentialReadOpened) {
		return ERRORCODES::READ_CHANNEL_IS_OPEN;
	}
	uint8_t Buffer1[SELFTEST_SIZE];
	uint8_t Buffer2[SELFTEST_SIZE];
	uint8_t Buffer3[SELFTEST_SIZE];
	get_JEDEC_ID(Buffer1, 2);
	if (Buffer1[0] != 0 || Buffer1[1] != 0) {
		// JEDEC ID not a useful test
		//SerialUSB.print("JEDEC ID = ");
		//SerialUSB.print(Buffer1[0], HEX);
		//SerialUSB.print(" ");
		//SerialUSB.print(Buffer1[1], HEX);
		//SerialUSB.println("");
		//return ERRORCODES::INCORRECT_JEDEC_ID; // Incorrect JEDEC ID response
	}
	ResetRegister();
	get_SRAM_Mode_Register(Buffer1[0]);
	if (Buffer1[0] != 0) {
		SerialUSB.print("Unexpected SRAM register: ");
		SerialUSB.println(Buffer1[0], DEC);
		return ERRORCODES::INCORRECT_STATUS_REGISTER; // Wrong values in the status register
	}

	// Read original -> write random -> read verify -> restore original

	//SerialUSB.print("Now testing ");
	//SerialUSB.print(BytesTested);
	//SerialUSB.print(" -- ");
	//SerialUSB.print(BytesTested + SELFTEST_SIZE - 1);
	//SerialUSB.println(" bytes");


	for (uint32_t TestOffset = 0; TestOffset < _SRAMSIZEBYTES; TestOffset+=SELFTEST_SIZE) {

		//SerialUSB.print("Now testing SRAM range starting at ");
		//SerialUSB.println(TestOffset);

		uint32_t ExpectedTotalBytesRead = 0;
		uint32_t ExpectedTotalBytesWritten = 0;

		// First, read the original content into Buffer1 and keep it there
		SeqReadOpen(TestOffset);
		SeqRead(Buffer1, SELFTEST_SIZE);
		ExpectedTotalBytesRead += SELFTEST_SIZE;
		if (get_SeqBytesRead() != ExpectedTotalBytesRead) {
			SeqReadClose();
			SeqWriteClose();
			return ERRORCODES::INCORRECT_READ_TOTAL;
		}

		// Overwrite with a random content
		for (size_t j = 0; j < SELFTEST_SIZE; j++) {
			Buffer2[j] = rand();
		}
		SeqWriteOpen(TestOffset);
		SeqWrite(Buffer2, SELFTEST_SIZE);
		ExpectedTotalBytesWritten += SELFTEST_SIZE;
		if (get_SeqBytesWritten() != ExpectedTotalBytesWritten) {
			SeqReadClose();
			SeqWriteClose();
			return ERRORCODES::INCORRECT_WRITE_TOTAL;
		}

		// Re-read what was written and verify
		SeqReadOpen(TestOffset);
		ExpectedTotalBytesRead = 0;
		SeqRead(Buffer3, SELFTEST_SIZE);
		for (size_t j = 0; j < SELFTEST_SIZE; j++) {
			//SerialUSB.print(j);
			//SerialUSB.print(" ");
			//SerialUSB.print(Buffer2[j], HEX);
			//SerialUSB.print(" ");
			//SerialUSB.print(Buffer3[j], HEX);
			//SerialUSB.println("");

			if (Buffer2[j] != Buffer3[j]) {
				SeqReadClose();
				SeqWriteClose();
				SerialUSB.println("SRAM read back error. Either this core does not really have control, the SRAM is defective, or interference is too strong.");
				return ERRORCODES::READBACK_ERROR; // Content of the re-read does not match.
			}
		}

		// Write back the original content
		SeqWriteOpen(TestOffset);
		SeqWrite(Buffer1, SELFTEST_SIZE);

		SeqReadClose();
		SeqWriteClose();

	}


	return ERRORCODES::OK; // Self test passed
}


//void SPISRAMClass::BeginPinControlTransferProcedure() {
//	/* This begins the pin control handoff. Both the releasing core and assuming core run this function first. */
//	_HasPinControl = false;
//	SeqWriteClose();
//	SeqReadClose();
//	pinMode(_CS_PIN, OUTPUT);
//	digitalWrite(_CS_PIN, HIGH);
//}


void SPISRAMClass::AssumePinControl() {
	/* This completes the step to gain pin control. */
	
	_HasPinControl = false;
	SeqWriteClose();
	SeqReadClose();

	pinMode(_CS_PIN, OUTPUT);
	digitalWrite(_CS_PIN, HIGH);
	(*SPIPTR).begin();

	// Required because they are different than those in variant.cpp
	pinPeripheral(_MISO_PIN, _MISO_PIOTYPE);
	pinPeripheral(_SCK_PIN, _SCK_PIOTYPE);
	pinPeripheral(_MOSI_PIN, _MOSI_PIOTYPE);
	UnselectDevice();
	ResetRegister();
	//EraseSRAM();
	_HasPinControl = true;
}


void SPISRAMClass::ReleasePinControl() {
	/* This releases the pin control. */
	_HasPinControl = false;

	//(*SPIPTR).end();
	pinMode(_CS_PIN, INPUT_PULLUP);
	pinMode(_MISO_PIN, INPUT);
	pinMode(_SCK_PIN, INPUT);
	pinMode(_MOSI_PIN, INPUT);

	SeqWriteClose();
	SeqReadClose();

}



void SPISRAMClass::Reverse_Endian(void* A, void* B, size_t Length) {
	uint8_t *AA = reinterpret_cast<uint8_t *>(A);
	uint8_t *BB = reinterpret_cast<uint8_t *>(B);
	for (size_t i = 0; i < Length; i++) {
		BB[i] = AA[Length - i - 1];
	}
}


void SPISRAMClass::SelectDevice() {
	(*SPIPTR).beginTransaction(SPISettings(SRAM_CLOCK_HZ, MSBFIRST, _SPI_DATAMODE));
	digitalWrite(_CS_PIN, LOW);
}

void SPISRAMClass::UnselectDevice() {
	digitalWrite(_CS_PIN, HIGH);
	(*SPIPTR).endTransaction();
}

void SPISRAMClass::EnableWriteSRAM() {
	SelectDevice();
	(*SPIPTR).transfer(OpCodeWriteEnable); // Enable write
	UnselectDevice();
}

void SPISRAMClass::ResetRegister() {

	/* Cypress SRAM requires enabling write first. */
	EnableWriteSRAM();
	SelectDevice();
	(*SPIPTR).transfer(OpCodeWriteStatusRegister); // Write mode register

	/* For Cypress SRAM: */
	(*SPIPTR).transfer(0b00000000);

	/* For Microchip SRAM:
	(*SPIPTR).transfer(0b01000000); // Write 01000000 to enter Sequential mode
	*/

	UnselectDevice();

	//SelectDevice();
	//(*SPIPTR).transfer(0xFF); // Reset to SPI bus mode
	//UnselectDevice();

}

bool SPISRAMClass::get_JEDEC_ID(uint8_t *Response, size_t Length) {
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization
	SelectDevice();
	(*SPIPTR).transfer(OpCodeGetDeviceID); // Get JEDEC ID
	(*SPIPTR).transfer(0x00);
	(*SPIPTR).transfer(Response, Length);
	UnselectDevice();
	return true;
}

bool SPISRAMClass::get_SRAM_Mode_Register(uint8_t &Response) {
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization
	SelectDevice();
	(*SPIPTR).transfer(OpCodeReadStatusRegister);
	Response = (*SPIPTR).transfer(0x00);
	UnselectDevice();
	return true;
}


bool SPISRAMClass::FillSRAM(uint8_t Fill) {
	/* Slow way to fill SRAM with the same byte. */
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization
	EnableWriteSRAM();
	SelectDevice();
	(*SPIPTR).transfer(OpCodeWrite);
	(*SPIPTR).transfer(0x00); // 24-bit address
	(*SPIPTR).transfer(0x00);
	(*SPIPTR).transfer(0x00);
	for (size_t i = 0; i < _SRAMSIZEBYTES; i++) {
		(*SPIPTR).transfer(Fill);
	}
	UnselectDevice();
	return true;
}

//void SPISRAMClass::DefaultWrite(uint8_t *Buffer, size_t DataLength, uint32_t Address) {
//	uint32_t AddressToSend;
//	Reverse_Endian(&Address, &AddressToSend, 3);
//	EnableWriteSRAM();
//	SelectDevice();
//	(*SPIPTR).transfer(0x02);
//	(*SPIPTR).transfer(&AddressToSend, 3); // 24-bit address out of a 32-bit uint. Only the 3 LSB will be in use
//	for (size_t i = 0; i < DataLength; i++) {
//		(*SPIPTR).transfer(Buffer[i]);
//	}
//	UnselectDevice();
//}
//
//
//void SPISRAMClass::DefaultRead(uint8_t *Buffer, size_t DataLength, uint32_t Address) {
//	uint32_t AddressToSend;
//	Reverse_Endian(&Address, &AddressToSend, 3);
//	SelectDevice();
//	(*SPIPTR).transfer(0x03);
//	(*SPIPTR).transfer(&AddressToSend, 3); // 24-bit address out of a 32-bit uint. Only the 3 LSB will be in use
//	(*SPIPTR).transfer(Buffer, DataLength);
//	UnselectDevice();
//}


bool SPISRAMClass::FastWrite(void *Buffer, size_t DataLength, uint32_t Address) {
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization

	//SerialUSB.print("F[");
	//SerialUSB.print((uint32_t) Buffer);
	//SerialUSB.print("]; ");



	uint8_t cmd = OpCodeWrite;
	uint32_t AddressToSend;
	Reverse_Endian(&Address, &AddressToSend, 3);
	EnableWriteSRAM();
	SelectDevice();
	(*SPIPTR).transferfastwriteonly(&cmd, 1);
	(*SPIPTR).transferfastwriteonly(&AddressToSend, 3);
	(*SPIPTR).transferfastwriteonly(Buffer, DataLength);
	UnselectDevice();
	return true;
}


bool SPISRAMClass::FastRead(void *Buffer, size_t DataLength, uint32_t Address) {
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization
	uint8_t cmd = OpCodeRead;
	uint32_t AddressToSend;
	Reverse_Endian(&Address, &AddressToSend, 3);
	SelectDevice();
	(*SPIPTR).transferfastwriteonly(&cmd, 1);
	(*SPIPTR).transferfastwriteonly(&AddressToSend, 3);
	(*SPIPTR).transferfastreadonly(Buffer, DataLength);
	UnselectDevice();
	return true;
}



bool SPISRAMClass::SeqWriteOpen(uint32_t StartAddress) {
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization
	_SequentialWriteBytesWritten = 0;
	_SequentialWriteNextAddress = StartAddress;
	_SequentialWriteOpened = true;
	return true;
}


uint32_t SPISRAMClass::SeqWrite(void *Buf, size_t BufferLength) {
	if (!_SequentialWriteOpened) { return 0; }
	uint8_t *Buffer = reinterpret_cast<uint8_t *>(Buf);
	FastWrite(Buffer, BufferLength, _SequentialWriteNextAddress);
	_SequentialWriteBytesWritten += BufferLength;
	_SequentialWriteNextAddress += BufferLength;
	_SequentialWriteNextAddress %= _SRAMSIZEBYTES;
	return BufferLength;
}


uint32_t SPISRAMClass::get_SeqBytesWritten() {
	return _SequentialWriteBytesWritten;
}


void SPISRAMClass::SeqWriteClose() {
	_SequentialWriteOpened = false;
}


bool SPISRAMClass::SeqReadOpen(uint32_t StartAddress) {
	if (!_SRAMSIZEBYTES || !_HasPinControl) { return false; } // Checks initialization
	_SequentialReadBytesRead = 0;
	_SequentialReadNextAddress = StartAddress;
	_SequentialReadOpened = true;
	return true;
}


uint32_t SPISRAMClass::SeqRead(void *Buf, size_t BufferLength) {
	if (!_SequentialReadOpened) { return 0; }
	uint8_t *Buffer = reinterpret_cast<uint8_t *>(Buf);
	FastRead(Buffer, BufferLength, _SequentialReadNextAddress);
	_SequentialReadBytesRead += BufferLength;
	_SequentialReadNextAddress += BufferLength;
	_SequentialReadNextAddress %= _SRAMSIZEBYTES;
	return BufferLength;
}


uint32_t SPISRAMClass::get_SeqBytesRead() {
	return _SequentialReadBytesRead;
}


void SPISRAMClass::SeqReadClose() {
	_SequentialReadOpened = false;
}



/* Bracket overloads */

SPISRAMClass::proxy SPISRAMClass::operator[](uint32_t i) {
	// Read
	return proxy(i, this);
}

SPISRAMClass::proxy::proxy(size_t i, SPISRAMClass *parentPtr) {
	_index = i;
	_parentPtr = parentPtr;
}

SPISRAMClass::proxy::operator uint8_t() const {
	uint8_t Buf[1];
	(*_parentPtr).FastRead(Buf, 1, _index);
	return Buf[0];
}

void SPISRAMClass::proxy::operator=(uint8_t rhs) {
	// Write
	(*_parentPtr).FastWrite(&rhs, 1, _index);
}

void SPISRAMClass::proxy::operator=(proxy rhs) {
	// Read and then write
	uint8_t temp = rhs;  // A temporary variable is required to evaluate the referenced operator function
	(*_parentPtr).FastWrite(&temp, 1, _index);
}





uint16_t SPISRAMClass::CalculateHash(uint32_t StartAddress, uint32_t DataLength) {
	// modified from http://www.sal.wisc.edu/st5000/documents/tables/crc16.c
	uint16_t crc = 0xffff;

	if (!_SRAMSIZEBYTES || !_HasPinControl) { return crc; } // Checks initialization


	size_t BufSize = SELFTEST_SIZE;
	unsigned char *Buffer = new unsigned char[BufSize];
	size_t n = 0;
	unsigned char *p;

	while (BufSize) {

		if (BufSize > DataLength) {
			BufSize = DataLength;
		}

		FastRead(Buffer, BufSize, StartAddress);
		n = BufSize;
		p = Buffer;
		while (n-- > 0) {
			crc = (unsigned char)(crc >> 8) | (crc << 8);
			crc ^= *p++;
			crc ^= (unsigned char)(crc & 0xff) >> 4;
			crc ^= (crc << 8) << 4;
			crc ^= ((crc & 0xff) << 4) << 1;
		}
		StartAddress += BufSize;
		DataLength -= BufSize;

	}

	delete Buffer;
	return(crc);
}


/* Support for SPISRAMFLOATARRAY in SPISRAMClass */
float SPISRAMClass::ReadFloat(uint32_t StartAddress) {
	float Data;
	FastRead(&Data, 4, StartAddress);
	return Data;
}

void SPISRAMClass::WriteFloat(float Data, uint32_t StartAddress) {
	FastWrite(&Data, 4, StartAddress);
}


/* SPISRAMFLOATARRAY */

// Define the static members once, and only once
std::vector <SPISRAMFLOATARRAY::ATracker> SPISRAMFLOATARRAY::ATrack(0); 
std::vector <SPISRAMFLOATARRAY::SRAMRange> SPISRAMFLOATARRAY::RTrack(0);


SPISRAMFLOATARRAY::SPISRAMFLOATARRAY() {
	// Default constructor does not do anything. 
	// If you created an object with the default constructor, do one of the following:
	// * Use resize(ArrayLength, SRAMPtr) to allocate an unused area automatically.
	// * Use manual_allocate(ArrayLength, SRAMPtr, ManualStartAddress) to manually allocate to a specific address.

	_Length = 0;
	_SRAMPtr = nullptr;
}


SPISRAMFLOATARRAY::SPISRAMFLOATARRAY(uint32_t ArrayLength, SPISRAMClass *SRAMPtr) {
	// Constructor that automatically allocates a free range to a new array to a specific SRAM
	AttemptAllocate(ArrayLength, SRAMPtr);
}


SPISRAMFLOATARRAY::SPISRAMFLOATARRAY(uint32_t ArrayLength, std::vector <SPISRAMClass *> &SRAMPtrList) {
	// Constructor that automatically allocates a free range to a new array to the first SRAM in the list that is capable of holding this.
	for (uint8_t s = 0; s < SRAMPtrList.size(); s++) {
		if (AttemptAllocate(ArrayLength, SRAMPtrList[s])) {
			return;
		}
	}

	// Defrag all SRAMs and try again
	if (Defrag(SRAMPtrList)) {
		for (uint8_t s = 0; s < SRAMPtrList.size(); s++) {
			if (AttemptAllocate(ArrayLength, SRAMPtrList[s])) {
				return;
			}
		}
	}

	SerialUSB.print("Cannot allocate an array with length ");
	SerialUSB.print(ArrayLength);
	SerialUSB.println(" in any available SRAM even after defrag.");
}


SPISRAMFLOATARRAY::SPISRAMFLOATARRAY(uint32_t ArrayLength, SPISRAMClass *SRAMPtr, uint32_t ManualStartAddress) {
	// Constructor that specifies the start address. 
	// This is useful if you already have data in the SRAM and wants to read it back easily
	manual_allocate(ArrayLength, SRAMPtr, ManualStartAddress);
}

// Manually specify the start address of the array. Useful if there is already data on the SRAM
// Cannot bypass SRAMRange
// Cannot bypass allocation check
void SPISRAMFLOATARRAY::manual_allocate(uint32_t ArrayLength, SPISRAMClass *SRAMPtr, uint32_t ManualStartAddress) {
	if (_Length > 0) {
		// Array already exists. Releasing current one.
		_Length = 0;
		Deallocate(_StartAddress, _SRAMPtr);
	}

	uint32_t ArrayBytes = ArrayLength * sizeof(float);
	bool conflict = false;

	/* First check the SRAMRange */
	uint32_t AllowedStartAddress;
	uint32_t AllowedEndAddress;
	GetSRAMRange(AllowedStartAddress, AllowedEndAddress, SRAMPtr);
	if (ManualStartAddress < AllowedStartAddress || ManualStartAddress + ArrayBytes - 1 > AllowedEndAddress) {
		conflict = true;
	}

	if (!conflict && !ATrack.empty()) {
		for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
			if (it->SRAMPtr == SRAMPtr) {
				if (ManualStartAddress >= it->StartAddress + it->Bytes) {
					// no conflict
					continue;
				}
				else if (ManualStartAddress + ArrayBytes <= it->StartAddress) {
					// no conflict
					continue;
				}
				else {
					// conflict
					conflict = true;
					break;
				}
			}
		}
	}
	if (!conflict) {
		Allocate(ArrayLength, ArrayBytes, ManualStartAddress, SRAMPtr);
	}
	else {
		_Length = 0;
		SerialUSB.print("Cannot allocate an array with length ");
		SerialUSB.print(ArrayLength);
		SerialUSB.print(" to the specified address ");
		SerialUSB.print(ManualStartAddress);
		SerialUSB.print(" on SRAM ");
		SerialUSB.println((uint32_t)SRAMPtr);
	}
}

SPISRAMFLOATARRAY::SPISRAMFLOATARRAY(const SPISRAMFLOATARRAY &old_obj) {
	// Copy constructor. Used when object is passed as value
	if (AttemptAllocate(old_obj._Length, old_obj._SRAMPtr)) {
		//SerialUSB.println("Copy allocation successful. Copying content to new location.");
		for (size_t i = 0; i < _Length; i++) {
			_SRAMPtr->WriteFloat(old_obj._SRAMPtr->ReadFloat(old_obj._StartAddress + i * sizeof(float)), _StartAddress + i * sizeof(float));
		}
	}
	else {
		_Length = 0;
	}
}

SPISRAMFLOATARRAY::~SPISRAMFLOATARRAY() {
	// Destructor
	//SerialUSB.print("Destructing ");
	//SerialUSB.println((uint32_t)this);
	Deallocate(_StartAddress, _SRAMPtr);
	_Length = 0;
	_SRAMPtr = nullptr;
}


/* Sets the allowable assignment address range for SRAM chips. Already allocated arrays are not affected
*/
void SPISRAMFLOATARRAY::SetSRAMRange(const uint32_t &AllowedStartAddress, const uint32_t &AllowedEndAddress, SPISRAMClass *SRAMPtr) {
	if (RTrack.size() == 0) {
		// Add it in
		SRAMRange a;
		a.SRAMPtr = SRAMPtr;
		a.AllowedStartAddress = AllowedStartAddress;
		a.AllowedEndAddress = AllowedEndAddress;
		RTrack.push_back(a);
		return;
	}
	else {
		// Look for the stored info
		bool found = false;
		for (std::vector<SRAMRange>::iterator it = RTrack.begin(); it < RTrack.end(); it++) {
			if (it->SRAMPtr == SRAMPtr) {
				found = true;
				it->AllowedStartAddress = AllowedStartAddress;
				it->AllowedEndAddress = AllowedEndAddress;
				return;
			}
		}
		if (!found) {
			// If not found, add it in
			SRAMRange a;
			a.SRAMPtr = SRAMPtr;
			a.AllowedStartAddress = AllowedStartAddress;
			a.AllowedEndAddress = AllowedEndAddress;
			RTrack.push_back(a);
			return;
		}

	}
}


void SPISRAMFLOATARRAY::GetSRAMRange(uint32_t &AllowedStartAddress, uint32_t &AllowedEndAddress, SPISRAMClass *SRAMPtr) {
	if (RTrack.size() == 0) {
		// Use default and add it in
		AllowedStartAddress = 0;
		AllowedEndAddress = SRAMPtr->get_SRAM_SIZE() - 1;
		SRAMRange a;
		a.SRAMPtr = SRAMPtr;
		a.AllowedStartAddress = AllowedStartAddress;
		a.AllowedEndAddress = AllowedEndAddress;
		RTrack.push_back(a);
		return;
	}
	else {
		// Look for the stored info
		bool found = false;
		for (std::vector<SRAMRange>::iterator it = RTrack.begin(); it < RTrack.end(); it++) {
			if (it->SRAMPtr == SRAMPtr) {
				found = true;
				AllowedStartAddress = it->AllowedStartAddress;
				AllowedEndAddress = it->AllowedEndAddress;
				return;
			}
		}
		if (!found) {
			// If not found, use default and add it in
			AllowedStartAddress = 0;
			AllowedEndAddress = SRAMPtr->get_SRAM_SIZE() - 1;
			SRAMRange a;
			a.SRAMPtr = SRAMPtr;
			a.AllowedStartAddress = AllowedStartAddress;
			a.AllowedEndAddress = AllowedEndAddress;
			RTrack.push_back(a);
			return;
		}

	}
}


bool SPISRAMFLOATARRAY::AttemptAllocate(uint32_t ArrayLength, SPISRAMClass *SRAMPtr) {
	//SerialUSB.print("Object attempts allocation: ");
	//SerialUSB.println((uint32_t)this);
	_Length = 0; // Default fail state. If Length is still zero by the end of constructor, allocation has failed.
	uint32_t ProposedStartAddress;
	uint32_t ArrayBytes = ArrayLength * sizeof(float);

	/* First check the SRAMRange */
	uint32_t AllowedStartAddress;
	uint32_t AllowedEndAddress;
	GetSRAMRange(AllowedStartAddress, AllowedEndAddress, SRAMPtr);

	if (ATrack.size() == 0) {
		// ATrack is empty, hence we can allocate anywhere. Choose ArrayStartAddress = 0.
		ProposedStartAddress = AllowedStartAddress;
		if (ArrayBytes <= AllowedEndAddress - AllowedStartAddress + 1) {
			//SerialUSB.print("Attempt empty. ");
			if (TestPinControl(SRAMPtr, ProposedStartAddress)) {
				Allocate(ArrayLength, ArrayBytes, ProposedStartAddress, SRAMPtr);
				return true;
			}
			else {
				return false;
			}
		}
		else {
			// If requested array size > RAM size, there is no way to fit no matter what
			return false;
		}
	}
	// Otherwise we have to find room before, between, or after existing tenants
	// Remember the list is not ordered by StartAddress so we have to look through all of them


	// First, try if we can allocate to the start of SRAM. This means that the start address of all tenants must be >= (ProposedStartAddress + ArrayBytes)
	// Example: If we propose to allocate this new array to SRAM1[0] thru SRAM1[999], then all existing tenants must have start address >= 1000.

	{ // Scope guard
		ProposedStartAddress = AllowedStartAddress;
		uint32_t MinimumRequiredStartAddressForExistingTenants = ProposedStartAddress + ArrayBytes;
		bool successful = true;
		for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
			if (it->SRAMPtr == SRAMPtr) {
				// Only care if SRAMPtr is the same

				if (it->StartAddress < MinimumRequiredStartAddressForExistingTenants) {
					// Found a tenant with start address conflicting the proposed range
					successful = false;
					break;
				}
			}
		}
		if (successful) {
			//SerialUSB.print("Attempt begin. ");
			if (ArrayBytes <= AllowedEndAddress - AllowedStartAddress + 1) {
				if (TestPinControl(SRAMPtr, ProposedStartAddress)) {
					Allocate(ArrayLength, ArrayBytes, ProposedStartAddress, SRAMPtr);
					return true;
				}
				else {
					return false;
				}
			}
			else {
				return false;
			}
		}
	}

	// Next, try if we can allocate to the end of SRAM. 
	// For each tenant, we'll move the ProposedStartAddress to the EndAddress+1 if it is bigger. 
	// After the loop, ProposedStartAddress should be MaximumEndAddressOfAllTenants + 1.
	{
		ProposedStartAddress = AllowedStartAddress;
		for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
			if (it->SRAMPtr == SRAMPtr) {
				// Only care if SRAMPtr is the same

				if (it->StartAddress + it->Bytes > ProposedStartAddress) {
					ProposedStartAddress = it->StartAddress + it->Bytes;
				}
			}
		}

		// Now check if ProposedStartAddress + ArrayBytes - 1 <= AllowedEndAddress
		if (ProposedStartAddress + ArrayBytes - 1 <= AllowedEndAddress) {
			//SerialUSB.print("Attempt end. ");
			if (TestPinControl(SRAMPtr, ProposedStartAddress)) {
				Allocate(ArrayLength, ArrayBytes, ProposedStartAddress, SRAMPtr);
				return true;
			}
			else {
				return false;
			}
		}
	}

	// Next, try to find a gap between existing tenants large enough for the new array
	{
		for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
			if (it->SRAMPtr == SRAMPtr) {
				// Only care if SRAMPtr is the same

				// For each tenant, propose the new array to reside immediately after this tenant
				// Then we check if this is possible
				ProposedStartAddress = it->StartAddress + it->Bytes;

				// Since tenants are not sorted, we have to loop again to see if another tenant is using this space
				// A tenant does not collide with the proposed range if StartAddress > ProposedEndAddress or EndAddress < ProposedStartAddress
				// We need all tenants to not collide

				

				//SerialUSB.print("Proposing newcomer Addr=");
				//SerialUSB.print(ProposedStartAddress);
				//SerialUSB.print(" Bytes=");
				//SerialUSB.println(ArrayBytes);

				bool successful = false;
				if (ProposedStartAddress >= AllowedStartAddress && ProposedStartAddress + ArrayBytes - 1 <= AllowedEndAddress) {
					// First of all, the proposed range must fit in SRAM. Otherwise it's automatic fail

					successful = true;
					for (std::vector<ATracker>::iterator jt = ATrack.begin(); jt < ATrack.end(); jt++) {
						if (jt->SRAMPtr == SRAMPtr) {

							if (jt->StartAddress >= ProposedStartAddress + ArrayBytes) {
								// not collide

								//SerialUSB.print("A tenant Addr=");
								//SerialUSB.print(jt->StartAddress);
								//SerialUSB.print(" Bytes=");
								//SerialUSB.print(jt->Bytes);
								//SerialUSB.print(" does not collide 1 with newcomer Addr=");
								//SerialUSB.print(ProposedStartAddress);
								//SerialUSB.print(" Bytes=");
								//SerialUSB.println(ArrayBytes);

								continue;
							}
							else if (jt->StartAddress + jt->Bytes <= ProposedStartAddress) {
								// not collide
								//SerialUSB.print("A tenant Addr=");
								//SerialUSB.print(jt->StartAddress);
								//SerialUSB.print(" Bytes=");
								//SerialUSB.print(jt->Bytes);
								//SerialUSB.print(" does not collide 2 with newcomer Addr=");
								//SerialUSB.print(ProposedStartAddress);
								//SerialUSB.print(" Bytes=");
								//SerialUSB.println(ArrayBytes);
								continue;
							}
							else {
								// collide
								//SerialUSB.print("A tenant Addr=");
								//SerialUSB.print(jt->StartAddress);
								//SerialUSB.print(" Bytes=");
								//SerialUSB.print(jt->Bytes);
								//SerialUSB.print(" collides with newcomer Addr=");
								//SerialUSB.print(ProposedStartAddress);
								//SerialUSB.print(" Bytes=");
								//SerialUSB.println(ArrayBytes);

								successful = false;
								break;
							}
						}
					}
				}
				else {
					successful = false;
				}

				if (successful) {
					//SerialUSB.print("Attempt gap. ");
					if (TestPinControl(SRAMPtr, ProposedStartAddress)) {
						Allocate(ArrayLength, ArrayBytes, ProposedStartAddress, SRAMPtr);
						return true;
					}
					else {
						return false;
					}
				}
			}
		}
	}

	//SerialUSB.print("Cannot allocate an array with length ");
	//SerialUSB.print(ArrayLength);
	//SerialUSB.print(" in SRAM pointed from ");
	//SerialUSB.println((uint32_t)SRAMPtr);
	return false;
}

void SPISRAMFLOATARRAY::Allocate(uint32_t ArrayLength, uint32_t ArrayBytes, uint32_t ArrayStartAddress, SPISRAMClass *SRAMPtr) {
	_StartAddress = ArrayStartAddress;
	_Length = ArrayLength;
	_SRAMPtr = SRAMPtr;
	ATracker a;
	a.ObjPtr = this;
	a.Bytes = ArrayBytes;
	a.StartAddress = ArrayStartAddress;
	a.SRAMPtr = SRAMPtr;

	//SerialUSB.print("ALLOCATE ObjPtr=");
	//SerialUSB.print((uint32_t)a.ObjPtr);
	//SerialUSB.print(" Bytes=");
	//SerialUSB.print((uint32_t)a.Bytes);
	//SerialUSB.print(" StartAddress=");
	//SerialUSB.print((uint32_t)a.StartAddress);
	//SerialUSB.print(" SRAMPtr=");
	//SerialUSB.println((uint32_t)a.SRAMPtr);

	ATrack.push_back(a);

	//SerialUSB.print("Allocated ");
	//SerialUSB.println(ATrack.size());

	//SerialUSB.print("Object ");
	//SerialUSB.print((uint32_t)this);
	//SerialUSB.print(" allocated an array Len=");
	//SerialUSB.print(ArrayLength);
	//SerialUSB.print(", Bytes=");
	//SerialUSB.print(ArrayBytes);
	//SerialUSB.print(" StartAddress=");
	//SerialUSB.print(ArrayStartAddress);
	//SerialUSB.print(" SRAMPtr=");
	//SerialUSB.print((uint32_t)SRAMPtr);
	//SerialUSB.println("");

}

// Locate this object in the ATrack list and delete it
// Does not delete/deallocate the object itself
void SPISRAMFLOATARRAY::Deallocate(uint32_t ArrayStartAddress, SPISRAMClass *SRAMPtr) {
	// Locate this object in the ATrack list and delete it

	//SerialUSB.print("Object ");
	//SerialUSB.print((uint32_t)this);
	//SerialUSB.println(" attempts deallocation.");
	//SerialUSB.print("ATrack size=");
	//SerialUSB.println(ATrack.size());


	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
		if (it->ObjPtr == this && it->SRAMPtr == SRAMPtr && it->StartAddress == ArrayStartAddress) {
			// Found it
			//SerialUSB.print("Object ");
			//SerialUSB.print((uint32_t)this);
			//SerialUSB.print(" deallocated the array StartAddress=");
			//SerialUSB.print(ArrayStartAddress);
			//SerialUSB.print(" SRAMPtr=");
			//SerialUSB.print((uint32_t)SRAMPtr);
			//SerialUSB.println("");
			ATrack.erase(it);
		}
	}
}

bool SPISRAMFLOATARRAY::DecreaseAllocation(uint32_t NewArrayLength) {
	// Reduce the allocated length of the array
	if (NewArrayLength >= _Length) {
		return false;
	}
	uint32_t NewArrayBytes = NewArrayLength * sizeof(float);
	// Locate this object in the ATrack list and modify the allocation bytes
	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
		if (it->SRAMPtr == _SRAMPtr && it->StartAddress == _StartAddress) {
			// Found it
			if (it->Bytes <= NewArrayBytes) {
				return false;
			}
			it->Bytes = NewArrayBytes;
			_Length = NewArrayLength;
			return true;
		}
	}
}

bool SPISRAMFLOATARRAY::IncreaseAllocationInPlace(uint32_t NewArrayLength) {

	//SerialUSB.println("Inside the InPlace");

	// Increase the allocated length of the array
	if (NewArrayLength <= _Length) {
		return false;
	}
	uint32_t NewArrayBytes = NewArrayLength * sizeof(float);

	//SerialUSB.println("Inplace2");

	// Locate this object in the ATrack list

	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {

		//SerialUSB.println("Inplace3");

		if (it->SRAMPtr == _SRAMPtr && it->StartAddress == _StartAddress) {
			// Found it

			//SerialUSB.println("Inplace4");

			if (it->Bytes >= NewArrayBytes) {
				//SerialUSB.println("Bad1");
				return false;
			}

			//SerialUSB.println("Inplace5");

			// Now make sure the new bytes don't exceed SRAM size limit
			uint32_t AllowedStartAddress;
			uint32_t AllowedEndAddress;
			GetSRAMRange(AllowedStartAddress, AllowedEndAddress, _SRAMPtr);

			//SerialUSB.println("Inplace6");

			if (_StartAddress + NewArrayBytes - 1 > AllowedEndAddress) {
				//SerialUSB.println("Bad2");
				return false;
			}

			//SerialUSB.println("Inplace7");

			// Now look at every object's start address
			bool noconflict = true;

			//SerialUSB.println("Inplace8");

			for (std::vector<ATracker>::iterator jt = ATrack.begin(); jt < ATrack.end(); jt++) {

				//SerialUSB.println("Inplace9");

				if (jt->SRAMPtr == _SRAMPtr) {

					//SerialUSB.println("Inplace10");

					if (_StartAddress != jt->StartAddress && _StartAddress + NewArrayBytes > jt->StartAddress) {
						noconflict = false;
						//SerialUSB.println("Bad3");
						return false;
					}
				}
			}
			if (noconflict) {
				// No other arrays have a start address in the proposed expansion range
				it->Bytes = NewArrayBytes;
				_Length = NewArrayLength;
				//SerialUSB.print("Good4. Array resized to ");
				//SerialUSB.print(it->Bytes);
				//SerialUSB.print(" bytes and ");
				//SerialUSB.print(_Length);
				//SerialUSB.println(" length");
				return true;
			}
		}
	}
	return false;

}




bool SPISRAMFLOATARRAY::Relocate(uint32_t NewArrayLength, SPISRAMClass *NewSRAMPtr) {
	// Relocate an array to NewSRAMPtr
	// This is meant to reallocate the entire array to a new StartAddress on the same or different SRAM.
	// It can be used to increase the size of the array.

	uint32_t NewArrayBytes = NewArrayLength * sizeof(float);
	uint32_t _OldLength = _Length;
	uint32_t _OldStartAddress = _StartAddress;
	SPISRAMClass *_OldSRAMPtr = _SRAMPtr;

	if (AttemptAllocate(NewArrayLength, NewSRAMPtr)) {
		// Allocation successful. Now test pin control at the new location on the new chip
		if (!TestPinControl(_SRAMPtr, _StartAddress)) {
			// No pin control. Abort relocation. Restore old referrers.
			Deallocate(_StartAddress, _SRAMPtr);
			_Length = _OldLength;
			_StartAddress = _OldStartAddress;
			_SRAMPtr = _OldSRAMPtr;
			SerialUSB.println("Relocate failed due to lack of pin control.");
			return false;
		}
		// Allocation successful. Object now refers to new location. Need to copy data from old location.
		for (size_t i = 0; i < _OldLength; i++) {
			_SRAMPtr->WriteFloat(_OldSRAMPtr->ReadFloat(_OldStartAddress + i * sizeof(float)), _StartAddress + i * sizeof(float));
		}
		// User function is responsible to zero out the content of the expanded area.

		// Once copied, deallocate the old location in ATrack.
		Deallocate(_OldStartAddress, _OldSRAMPtr);

		return true;
	}
	else {
		// Allocation failure. Restore old referrers.
		//SerialUSB.println("Cannot relocate.");
		_Length = _OldLength;
		_StartAddress = _OldStartAddress;
		_SRAMPtr = _OldSRAMPtr;
		return false;
	}


	return false;

}



// Resize using a list of SRAM modules (i.e. try first SRAM, then next SRAM, etc.)
// Warning: This function does not change the contents, but only change the allocation.
bool SPISRAMFLOATARRAY::resize(uint32_t NewArrayLength, std::vector <SPISRAMClass *> &SRAMPtrList) {
	for (uint8_t s = 0; s < SRAMPtrList.size(); s++) {
		if (resize(NewArrayLength, SRAMPtrList[s])) {
			return true;
		}
	}

	// Defrag all SRAMs and try again
	if (Defrag(SRAMPtrList)) {
		for (uint8_t s = 0; s < SRAMPtrList.size(); s++) {
			if (resize(NewArrayLength, SRAMPtrList[s])) {
				return true;
			}
		}
	}

	SerialUSB.print("Cannot allocate an array with length ");
	SerialUSB.print(NewArrayLength);
	SerialUSB.println(" in any available SRAM even after defrag.");

	return false;
}



// Resize using the specific SRAM if needed, i.e. this won't move content to new SRAM if in-place resize is successful
bool SPISRAMFLOATARRAY::resize(uint32_t NewArrayLength, SPISRAMClass *NewSRAMPtr) {
	if (NewArrayLength == 0 && _Length > 0) {
		// This is the same as deallocating.
		_Length = 0;
		Deallocate(_StartAddress, _SRAMPtr);
		return true;
	}
	if (NewArrayLength < _Length) {
		//SerialUSB.println("AlloDec");
		return DecreaseAllocation(NewArrayLength);
	}
	if (NewArrayLength > _Length) {
		//SerialUSB.println("AlloInc");
		// First attempt to expand in place
		if (_Length && IncreaseAllocationInPlace(NewArrayLength)) {
			//SerialUSB.println("InPlaceTrue");
			return true;
		}
		else {
			//SerialUSB.println("Attempt relocate.");
			// Then try to completely relocate to a new address range.
			return Relocate(NewArrayLength, NewSRAMPtr);
		}
	}
	if (NewArrayLength == _Length) {
		// No need to resize
		return true;
	}
	return false;
}


// Resizes within the same SRAM (will fail if it has never been allocated)
bool SPISRAMFLOATARRAY::resize(uint32_t NewArrayLength) {
	if (_SRAMPtr == nullptr) {
		return false;
	}
	SPISRAMClass *NewSRAMPtr = _SRAMPtr;
	return resize(NewArrayLength, NewSRAMPtr);
}


// Size of the array (not bytes used)
uint32_t SPISRAMFLOATARRAY::size() {
	uint32_t t = _Length;
	return t;
}

uint32_t SPISRAMFLOATARRAY::address() {
	uint32_t t = _StartAddress;
	return t;
}

SPISRAMClass* SPISRAMFLOATARRAY::SRAMPtr() {
	SPISRAMClass *t = _SRAMPtr;
	return t;
}

const size_t SPISRAMFLOATARRAY::GetTotalBytesAllocated(SPISRAMClass *SRAMPtr) {
	if (ATrack.empty()) {
		return 0;
	}
	size_t bytes = 0;
	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
		if (it->SRAMPtr == SRAMPtr) {
			bytes += it->Bytes;
		}
	}
	return bytes;
}

const size_t SPISRAMFLOATARRAY::CountAllocatedArrays(SPISRAMClass *SRAMPtr) {
	if (ATrack.empty()) {
		return 0;
	}
	size_t count = 0;
	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
		if (it->SRAMPtr == SRAMPtr) {
			count++;
		}
	}
	return count;
}

void SPISRAMFLOATARRAY::GetArrayObjPtrs(SPISRAMClass *SRAMPtr, std::vector <SPISRAMFLOATARRAY*> &List) {

	if (ATrack.empty()) {
		//SerialUSB.println("ATrack is empty.");
		return;
	}
	size_t count = 0;
	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
		//SerialUSB.println("debug1.");
		if (it->SRAMPtr == SRAMPtr) {
			//SerialUSB.println("debug2.");
			List.push_back(it->ObjPtr);
		}
	}
	return;
}


// Defrag all arrays in all available SRAMs, including moving arrays from one chip to another
// The goal is to create the largest contiguous free space
size_t SPISRAMFLOATARRAY::Defrag(std::vector <SPISRAMClass *> &SRAMPtrList) {

	if (SRAMPtrList.empty()) {
		return 0;
	}

	size_t nm = 0;

	// Try to move as many arrays as possible to the first SRAM (in the SRAMPtrList), then second SRAM, etc.
	for (size_t s = 0; s < SRAMPtrList.size() - 1; s++) {
		nm += Defrag(SRAMPtrList[s]);
		for (size_t t = s + 1; t < SRAMPtrList.size(); t++) {
			for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
				if (it->SRAMPtr == SRAMPtrList[t]) {
					//SerialUSB.print(" Trying to move [&SRAM=");
					//SerialUSB.print((uint32_t)it->SRAMPtr);
					//SerialUSB.print(" StartAddress=");
					//SerialUSB.print(it->StartAddress);
					//SerialUSB.print(" Bytes=");
					//SerialUSB.print(it->Bytes);
					//SerialUSB.print("] to &SRAM ");
					//SerialUSB.println((uint32_t)SRAMPtrList[s]);
					if (it->ObjPtr->Relocate(it->ObjPtr->_Length, SRAMPtrList[s])) {
						nm++;
					}
				}
			}
		}
	}
	return nm;

}



// Defrag all arrays in the entire SRAM. Returns the number of arrays moved.
size_t SPISRAMFLOATARRAY::Defrag(SPISRAMClass *SRAMPtr) {
	size_t nm = 0;
	uint32_t AllowedStartAddress;
	uint32_t AllowedEndAddress;
	GetSRAMRange(AllowedStartAddress, AllowedEndAddress, SRAMPtr);
	uint32_t WorkingAddress = AllowedStartAddress; // Address before this is defragged. Addresses at and after this are not.


	SerialUSB.print("Defragging SRAM ");
	SerialUSB.print((uint32_t)SRAMPtr);
	SerialUSB.print(" ..");

	// Sort ATrack in ascending StartAddress order
	std::sort(ATrack.begin(), ATrack.end());

	for (std::vector<ATracker>::iterator it = ATrack.begin(); it < ATrack.end(); it++) {
		if (it->SRAMPtr == SRAMPtr) {

			if (it->StartAddress > WorkingAddress && WorkingAddress + it->Bytes - 1 <= AllowedEndAddress) {
				// There is a gap before this array, and the new range is permitted

				if (!TestPinControl(SRAMPtr, WorkingAddress)) {
					// No pin control
					return nm;
				}

				// SRAM copy
				//SerialUSB.print("Moving array ");
				//SerialUSB.print((uint32_t)it->ObjPtr);
				//SerialUSB.print(" from ");
				//SerialUSB.print(it->StartAddress);
				//SerialUSB.print(":");
				//SerialUSB.print(it->StartAddress + it->Bytes - 1);
				//SerialUSB.print(" to ");
				//SerialUSB.print(WorkingAddress);
				//SerialUSB.print(":");
				//SerialUSB.println(WorkingAddress + it->Bytes - 1);
				for (size_t i = 0; i < it->Bytes; i++) {
					SRAMPtr->WriteFloat(SRAMPtr->ReadFloat(it->StartAddress + i), WorkingAddress + i);
				}
				// Set offset for *it to the new start address
				// Change both the address in ATrack and in the class object
				it->StartAddress = WorkingAddress;
				it->ObjPtr->_StartAddress = WorkingAddress;
				WorkingAddress = it->StartAddress + it->Bytes;
				nm++;
			}
			else {
				//SerialUSB.print("Array ");
				//SerialUSB.print((uint32_t)it->ObjPtr);
				//SerialUSB.print(" does not need to move from ");
				//SerialUSB.print(it->StartAddress);
				//SerialUSB.print(":");
				//SerialUSB.println(it->StartAddress + it->Bytes - 1);
				WorkingAddress = it->StartAddress + it->Bytes;
			}

		}
	}
	SerialUSB.println(" Done.");
	return nm;
}



bool SPISRAMFLOATARRAY::ATracker::operator<(const ATracker &rhs) const {
	if (SRAMPtr < rhs.SRAMPtr) {
		return true;
	}
	else if (SRAMPtr > rhs.SRAMPtr) {
		return false;
	}
	else {
		if (StartAddress < rhs.StartAddress) {
			return true;
		}
		else {
			return false;
		}
	}
}


//bool SPISRAMFLOATARRAY::ATracker::operator>(const ATracker &rhs) const {
//	if (Bytes > rhs.Bytes) {
//		return true;
//	}
//	else {
//		return false;
//	}
//}



// Test if this processor has pin control by doing a test read/write
bool SPISRAMFLOATARRAY::TestPinControl(SPISRAMClass *SRAMPtr, uint32_t TestAddress) {
	bool status = false;
	float original = SRAMPtr->ReadFloat(TestAddress);
	float testdata = original + PI;
	if (!std::isfinite(original)) {
		testdata = PI;
	}
	SRAMPtr->WriteFloat(testdata, TestAddress);
	float readback = SRAMPtr->ReadFloat(TestAddress);
	if (readback == testdata) {
		status = true;
	}
	// Write back the original content no matter the test result
	SRAMPtr->WriteFloat(original, TestAddress);

	if (!status) {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::TestPinControl(SPISRAMClass *SRAMPtr, uint32_t TestAddress): No pin control. ");
		SerialUSB.print("original=");
		SerialUSB.print(original, 3);
		SerialUSB.print(", testdata=");
		SerialUSB.print(testdata, 3);
		SerialUSB.print(", readback=");
		SerialUSB.println(readback, 3);
	}

	return status;
}



SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::operator[](uint32_t i) {
	// Read
	return proxy(i, this);
}

SPISRAMFLOATARRAY::proxy::proxy(size_t i, SPISRAMFLOATARRAY *parentPtr) {
	_index = i;
	_parentPtr = parentPtr;
}

SPISRAMFLOATARRAY::proxy::operator float() const {
	// Read
	if (_index < _parentPtr->_Length) {
		return _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator float() const: Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return OOB_RETURN_VALUE;
	}
}

float SPISRAMFLOATARRAY::proxy::operator=(float rhs) {
	// Write
	if (_index < _parentPtr->_Length) {
		_parentPtr->_SRAMPtr->WriteFloat(rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator=(float rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}
}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator=(proxy rhs) {
	// Read and then write
	if (_index < _parentPtr->_Length) {
		float temp = rhs;
		_parentPtr->_SRAMPtr->WriteFloat(temp, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator=(proxy rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}
}

float SPISRAMFLOATARRAY::proxy::operator+=(float rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x + rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator+=(float rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}
}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator+=(proxy rhs) {
	//SerialUSB.println("proxy+=");
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x + rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator+=(proxy rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}
}

float SPISRAMFLOATARRAY::proxy::operator-=(float rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x - rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator-=(float rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}
}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator-=(proxy rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x - rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator-=(proxy rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}
}

float SPISRAMFLOATARRAY::proxy::operator*=(float rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x * rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator*=(float rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}

}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator*=(proxy rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x * rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator*=(proxy rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}

}

float SPISRAMFLOATARRAY::proxy::operator/=(float rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x / rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator/=(float rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}

}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator/=(proxy rhs) {
	if (_index < _parentPtr->_Length) {
		// Read existing
		float x = _parentPtr->_SRAMPtr->ReadFloat(_parentPtr->_StartAddress + _index * sizeof(float));
		// Write
		_parentPtr->_SRAMPtr->WriteFloat(x / rhs, _parentPtr->_StartAddress + _index * sizeof(float));
		return *this;
	}
	else {
		SerialUSB.print("Warning in SPISRAMFLOATARRAY::proxy::operator/=(proxy rhs): Index out of bound: i=");
		SerialUSB.print(_index);
		SerialUSB.print(" len=");
		SerialUSB.println(_parentPtr->_Length);
		return *this;
	}

}

SPISRAMFLOATARRAY::proxy& SPISRAMFLOATARRAY::proxy::operator++() {
	this->operator+=(1.0);
	return *this;
}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator++(int) {
	proxy tmp(*this);
	operator++();
	return tmp;
}

SPISRAMFLOATARRAY::proxy& SPISRAMFLOATARRAY::proxy::operator--() {
	this->operator-=(1.0);
	return *this;
}

SPISRAMFLOATARRAY::proxy SPISRAMFLOATARRAY::proxy::operator--(int) {
	proxy tmp(*this);
	operator--();
	return tmp;
}


// Operators that act on the entire array.

// Array equal
SPISRAMFLOATARRAY& SPISRAMFLOATARRAY::operator=(SPISRAMFLOATARRAY &rhs) {
	size_t len = this->size();
	if (rhs.size() != len) {
		if (!(this->resize(len))) {
			// Can't resize
			SerialUSB.print("Warning in SPISRAMFLOATARRAY::operator=(SPISRAMFLOATARRAY &rhs): Cannot resize.");
			return *this;
		}
	}
	//SerialUSB.print(" copying from ");
	//SerialUSB.print((uint32_t)&rhs);
	//SerialUSB.print(" to ");
	//SerialUSB.println((uint32_t)this);
	for (size_t i = 0; i < len; i++) {
		(*this)[i] = rhs[i];
	}
	return *this;
}

float SPISRAMFLOATARRAY::dotprod(SPISRAMFLOATARRAY &rhs) {
	float out = 0.0;
	size_t len = this->size();
	if (rhs.size() < len) {
		len = rhs.size();
	}
	for (size_t i = 0; i < len; i++) {
		out += (*this)[i] * rhs[i];
	}
	return out;
}

float SPISRAMFLOATARRAY::mean() {
	float out = 0.0;
	size_t len = this->size();
	for (size_t i = 0; i < len; i++) {
		out += (*this)[i];
	}
	out /= len;
	return out;
}

float SPISRAMFLOATARRAY::var() {
	float m = this->mean();
	float out = 0.0;
	float t;
	size_t len = this->size();
	for (size_t i = 0; i < len; i++) {
		t = (*this)[i] - m;
		out += t * t;
	}
	out /= (len - 1);
	return out;
}

SPISRAMFLOATARRAY& SPISRAMFLOATARRAY::operator*=(const float rhs) {
	for (size_t i = 0; i < this->size(); i++) {
		(*this)[i] *= rhs;
	}
	return *this;
}

SPISRAMFLOATARRAY& SPISRAMFLOATARRAY::operator/=(const float rhs) {
	for (size_t i = 0; i < this->size(); i++) {
		(*this)[i] /= rhs;
	}
	return *this;
}

