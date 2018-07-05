// 
// 
// 

#include "storage.h"

#define NAND_CLOCK_HZ 12000000
#define REGULARFILE_FILETYPE 82

SPINANDClass::SPINANDClass() {
	/* Dummy initializer. Must call Initialize before doing anything else to the NAND. */
	_NAND_BYTES_PER_PAGE = 0;
	_HasPinControl = false;
}


void SPINANDClass::Initialize(SERCOM *SPINAND_SERCOM,
	uint8_t SPINAND_CS_PIN, uint8_t SPINAND_MISO_PIN, uint8_t SPINAND_SCK_PIN, uint8_t SPINAND_MOSI_PIN,
	uint8_t SPI_DATAMODE,
	SercomSpiTXPad SPINAND_SERCOM_TXPAD, SercomRXPad SPINAND_SERCOM_RXPAD,
	_EPioType SPINAND_MISO_PIOTYPE, _EPioType SPINAND_SCK_PIOTYPE, _EPioType SPINAND_MOSI_PIOTYPE,
	size_t NAND_BYTES_PER_PAGE, size_t NAND_PAGES_PER_BLOCK, size_t NAND_BLOCKS_PER_DIE, size_t NAND_DIES_PER_DEVICE) {

	if (_NAND_BYTES_PER_PAGE || !NAND_BYTES_PER_PAGE) { return; } // Prevents initializing again

																  // MISO, SCK, MOSI
	_CS_PIN = SPINAND_CS_PIN;
	_MISO_PIN = SPINAND_MISO_PIN;
	_MOSI_PIN = SPINAND_MOSI_PIN;
	_SCK_PIN = SPINAND_SCK_PIN;
	_SPI_DATAMODE = SPI_DATAMODE;
	_SERCOM_TXPAD = SPINAND_SERCOM_TXPAD;
	_SERCOM_RXPAD = SPINAND_SERCOM_RXPAD;
	_MISO_PIOTYPE = SPINAND_MISO_PIOTYPE;
	_MOSI_PIOTYPE = SPINAND_MOSI_PIOTYPE;
	_SCK_PIOTYPE = SPINAND_SCK_PIOTYPE;

	delay(10);

	SPIPTR = new SPIClass(SPINAND_SERCOM, _MISO_PIN, _SCK_PIN, _MOSI_PIN, _SERCOM_TXPAD, _SERCOM_RXPAD);

	_NAND_BYTES_PER_PAGE = NAND_BYTES_PER_PAGE;
	_NAND_PAGES_PER_BLOCK = NAND_PAGES_PER_BLOCK;
	_NAND_BLOCKS_PER_DIE = NAND_BLOCKS_PER_DIE;
	_NAND_DIES_PER_DEVICE = NAND_DIES_PER_DEVICE;
	_NAND_PAGES_PER_DIE = _NAND_PAGES_PER_BLOCK * _NAND_BLOCKS_PER_DIE;
	_NAND_PAGES_PER_DEVICE = _NAND_PAGES_PER_DIE * _NAND_DIES_PER_DEVICE;
	_NAND_BLOCKS_PER_DEVICE = ((_NAND_BLOCKS_PER_DIE * _NAND_DIES_PER_DEVICE) / 8) * 8; // Round down to multiple of 8 so DeviceBlockMap won't overflow

	_HighestGlobalFileSerial = UINT32_MAX;
	_HighestGlobalPageSerial = UINT64_MAX;
	_HighestGlobalPageIndex = UINT32_MAX;
	_HighestGlobalPageScanned = false;
	_DeviceBlockMapBuilt = false;
	_FileWriteChannelOpen = false;
	_FileReadChannelOpen = false;


	// Do not allocate until we AssumePinControl
	DeviceBlockMap = nullptr;
	//DeviceBlockMap = new uint8_t[_NAND_BLOCKS_PER_DEVICE / 8]; // A bit map to mark for used/unused blocks

	ReleasePinControl(); // By default, set the pins to input mode until we decide to gain pin control

}

SPINANDClass::~SPINANDClass() {
	if (_HasPinControl) {
		UnselectDevice();
	}
	ReleasePinControl();
	_NAND_BYTES_PER_PAGE = 0;
	(*SPIPTR).end();
	delete DeviceBlockMap;
	delete SPIPTR;
}

//void SPINANDClass::BeginPinControl() {
//	/* This begins the pin control. The proper procedure for handing off pin control is:
//	1. First core sets CS pin to high
//	2. First core informs second core
//	3. Second core sets CS pin to output mode and to high (this function)
//	4. First core releases pin control (sets all its pins to input mode)
//	5. Second core reconfigures mux to SPI mode and begins SPI (AssumePinControl function)
//	*/
//	pinMode(_CS_PIN, OUTPUT);
//	digitalWrite(_CS_PIN, HIGH);
//}


void SPINANDClass::AssumePinControl() {
	/* This completes the step to gain pin control. */
	pinMode(_CS_PIN, OUTPUT);
	digitalWrite(_CS_PIN, HIGH);
	(*SPIPTR).begin();
	pinPeripheral(_MISO_PIN, _MISO_PIOTYPE);
	pinPeripheral(_SCK_PIN, _SCK_PIOTYPE);
	pinPeripheral(_MOSI_PIN, _MOSI_PIOTYPE);
	UnselectDevice();
	delayMicroseconds(1250);
	_SelectedDie = 1;
	SelectDie(_SelectedDie);
	DisableAllBlockLocks();
	EnableECC();
	_SelectedDie = 0;
	SelectDie(_SelectedDie);
	DisableAllBlockLocks();
	EnableECC();

	if (DeviceBlockMap == nullptr) {
		DeviceBlockMap = new uint8_t[_NAND_BLOCKS_PER_DEVICE / 8];
	}

	_HasPinControl = true;
}


void SPINANDClass::ReleasePinControl() {
	/* This releases the pin control. */
	_HasPinControl = false;

	FileWriteClose();
	FileReadClose();
	_DeviceBlockMapBuilt = false;
	_HighestGlobalPageScanned = false;


	//(*SPIPTR).end();
	pinMode(_CS_PIN, INPUT_PULLUP);
	pinMode(_MISO_PIN, INPUT);
	pinMode(_SCK_PIN, INPUT);
	pinMode(_MOSI_PIN, INPUT);

	delete DeviceBlockMap;
	DeviceBlockMap = nullptr;

}


void SPINANDClass::Reverse_Endian(void *A, void *B, size_t Length) {
	uint8_t *AA = reinterpret_cast<uint8_t *>(A);
	uint8_t *BB = reinterpret_cast<uint8_t *>(B);
	for (size_t i = 0; i < Length; i++) {
		BB[i] = AA[Length - i - 1];
	}
}


void SPINANDClass::SelectDevice() {
	(*SPIPTR).beginTransaction(SPISettings(NAND_CLOCK_HZ, MSBFIRST, _SPI_DATAMODE));
	digitalWrite(_CS_PIN, LOW);
}

void SPINANDClass::UnselectDevice() {
	digitalWrite(_CS_PIN, HIGH);
	(*SPIPTR).endTransaction();
}

void SPINANDClass::EnableWrite() {
	SelectDevice();
	(*SPIPTR).transfer(0x06); // Enable write
	UnselectDevice();
}

bool SPINANDClass::WakeupDevice() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	SelectDevice();
	(*SPIPTR).transfer(0xAB);
	(*SPIPTR).transfer(0);
	UnselectDevice();
	return true;
}

bool SPINANDClass::WaitForReady() {
	uint8_t DataBack;
	size_t q = 9999;
	SelectDevice();
	(*SPIPTR).transfer(0x0f); // get feature
	(*SPIPTR).transfer(0xc0); // status register address
	while (q--) {
		DataBack = (*SPIPTR).transfer(0x00); // status register data out
		if (DataBack == 0) {
			break;
		}
	}
	UnselectDevice();
	if (q == 0) {
		return 1;
	}
	else {
		return 0;
	}
}

uint8_t SPINANDClass::get_blocklock() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return UINT8_MAX; } // Checks initialization
	uint8_t Response;
	SelectDevice();
	(*SPIPTR).transfer(0x0f);
	(*SPIPTR).transfer(0xa0);
	Response = (*SPIPTR).transfer(0x00);
	UnselectDevice();
	return Response;
}

uint8_t SPINANDClass::get_config() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return UINT8_MAX; } // Checks initialization
	uint8_t Response;
	SelectDevice();
	(*SPIPTR).transfer(0x0f);
	(*SPIPTR).transfer(0xb0);
	Response = (*SPIPTR).transfer(0x00);
	UnselectDevice();
	return Response;
}

uint8_t SPINANDClass::get_status() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return UINT8_MAX; } // Checks initialization
	uint8_t Response;
	SelectDevice();
	(*SPIPTR).transfer(0x0f);
	(*SPIPTR).transfer(0xc0);
	Response = (*SPIPTR).transfer(0x00);
	UnselectDevice();
	return Response;
}

uint8_t SPINANDClass::get_dieselect() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return UINT8_MAX; } // Checks initialization
	uint8_t Response;
	SelectDevice();
	(*SPIPTR).transfer(0x0f);
	(*SPIPTR).transfer(0xd0);
	Response = (*SPIPTR).transfer(0x00);
	UnselectDevice();
	return Response;
}


void SPINANDClass::EnableECC() {
	SelectDevice();
	(*SPIPTR).transfer(0x0f);
	(*SPIPTR).transfer(0xb0);
	(*SPIPTR).transfer(0b00010000);
	UnselectDevice();
}


void SPINANDClass::DisableAllBlockLocks() {
	SelectDevice();
	(*SPIPTR).transfer(0x1f);
	(*SPIPTR).transfer(0xa0);
	(*SPIPTR).transfer(0b00000000);
	UnselectDevice();
}


bool SPINANDClass::get_JEDEC_ID(uint8_t *Response, size_t Length) {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	SelectDevice();
	(*SPIPTR).transfer(0x9f); // Get JEDEC ID
	(*SPIPTR).transfer(0x00);
	(*SPIPTR).transfer(Response, Length);
	UnselectDevice();
	return true;
}

uint8_t SPINANDClass::SelectDie(size_t die) {
	SelectDevice();
	(*SPIPTR).transfer(0x1f);
	(*SPIPTR).transfer(0xd0);
	switch (die) {
	case 0:
		(*SPIPTR).transfer(0x00);
		_SelectedDie = 0;
		break;
	case 1:
		(*SPIPTR).transfer(0x40);
		_SelectedDie = 1;
		break;
	}
	UnselectDevice();
	return _SelectedDie;
}


bool SPINANDClass::EraseDeviceFull() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	for (size_t i = 0; i < _NAND_DIES_PER_DEVICE; i++) {
		EraseDie(i);
	}
	_HighestGlobalPageIndex = UINT32_MAX;
	_HighestGlobalPageSerial = UINT64_MAX;
	_HighestGlobalFileSerial = UINT32_MAX;
	_DeviceBlockMapBuilt = false;
	_HighestGlobalPageScanned = false;

	return true;
}


void SPINANDClass::EraseDie(size_t die) {
	uint32_t Address;
	SelectDie(die);
	for (size_t block = 0; block < _NAND_BLOCKS_PER_DIE; block++) {
		Address = ((uint32_t)block) << 6;
		EraseBlock(Address);
	}
}


void SPINANDClass::EraseBlock(uint32_t RowAddress) {
	uint32_t AddressToSend;
	Reverse_Endian(&RowAddress, &AddressToSend, 3);
	EnableWrite();
	SelectDevice();
	(*SPIPTR).transfer(0xd8); // Erase block
	(*SPIPTR).transfer(&AddressToSend, 3);
	UnselectDevice();
	WaitForReady();
}


uint16_t SPINANDClass::EraseDevice() {
	/* Erase only the used blocks and return number of blocks erased. */
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return 0; } // Checks initialization
	if (!_DeviceBlockMapBuilt) {
		BuildDeviceBlockMap();
	}
	if (!_HighestGlobalPageScanned) {
		ScanForHighestPage();
	}
	uint32_t RA;
	uint16_t Erased = 0;
	uint8_t Die;
	for (uint16_t GBI = 0; GBI < _NAND_BLOCKS_PER_DEVICE; GBI++) {
		if (IsBlockUsed(GBI)) {
			GlobalBlockIndex_to_RA_and_Die(GBI, RA, Die);
			if (Die != _SelectedDie) {
				SelectDie(Die);
			}
			EraseBlock(RA);
			MarkBlockUnused(GBI);
			Erased++;
		}
	}

	_HighestGlobalPageIndex = InitializeHighestGPI(_HighestGlobalPageIndex);
	_HighestGlobalPageSerial = UINT64_MAX;
	_HighestGlobalFileSerial = UINT32_MAX;

	return Erased;
}


void SPINANDClass::EraseData(uint32_t GPI) {
	/* Uses Global Page Index to make it more accessible to public.
	Warning: This still erases the entire block that the page is in.
	GPI will be set to the first page of the block during the function.
	*/

	uint64_t PageSerial;
	uint32_t RA, FileSerial, FilePageIndex;
	uint16_t GBI, BytesUsedInThisPage;
	uint8_t Die, FileType;
	char FileDescription[NANDFILEDESCRIPTIONLENGTH];
	bool PlaneBit;

	GlobalPageIndex_to_GlobalBlockIndex(GPI, GBI);
	GlobalBlockIndex_to_GlobalPageIndex(GBI, GPI);


	GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
	FastRetrieveFromNAND(RA);
	ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription, PlaneBit);
	if (PageSerial == UINT64_MAX) {
		/* This block is already empty. */
		return;
	}
	if (PageSerial + _NAND_PAGES_PER_BLOCK >= _HighestGlobalPageSerial) {
		/* If the block being erased contains a page with the highest serial, the highest serial is no longer valid. */
		_HighestGlobalPageIndex = UINT32_MAX;
		_HighestGlobalPageSerial = UINT64_MAX;
		_HighestGlobalFileSerial = UINT32_MAX;
		_HighestGlobalPageScanned = false;
	}

	EraseBlock(RA);
	/* Set the block map to unused again */
	MarkBlockUnused(GBI);
}

uint32_t SPINANDClass::InitializeHighestGPI() {
	randomSeed(micros());
	uint32_t GPI = random(_NAND_BLOCKS_PER_DEVICE)*_NAND_PAGES_PER_BLOCK;
	GPI--;
	return GPI;
}


uint32_t SPINANDClass::InitializeHighestGPI(uint32_t GPI) {
	uint16_t GBI;
	GlobalPageIndex_to_GlobalBlockIndex(GPI, GBI);
	GBI++; // Don't reuse the block we last used.
	GlobalBlockIndex_to_GlobalPageIndex(GBI, GPI);
	GPI--;
	return GPI;
}


void SPINANDClass::MarkBlockUnused(uint16_t GBI) {
	if (DeviceBlockMap != nullptr) {
		DeviceBlockMap[GBI / 8] |= (1 << (GBI % 8));
	}
}

void SPINANDClass::MarkBlockUsed(uint16_t GBI) {
	if (DeviceBlockMap != nullptr) {
		DeviceBlockMap[GBI / 8] &= (UINT8_MAX - (1 << (GBI % 8)));
	}
}

bool SPINANDClass::IsBlockUsed(uint16_t GBI) {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	if (!_DeviceBlockMapBuilt) {
		BuildDeviceBlockMap();
	}
	if (!_HighestGlobalPageScanned) {
		ScanForHighestPage();
	}

	if (DeviceBlockMap != nullptr) {
		if (DeviceBlockMap[GBI / 8] & (1 << (GBI % 8))) {
			// 1 bit is empty
			// 0 bit is used (partially or full)
			return false;
		}
		return true;
	}
	else {
		// Device block map is not available due to insufficient CPU RAM. We'll do it the slow way.
		// NOT TESTED
		uint64_t PageSerial;
		uint32_t GPI, RA, FileSerial, FilePageIndex;
		uint16_t BytesUsedInThisPage;
		uint8_t Die, FileType;
		char FileDescription[NANDFILEDESCRIPTIONLENGTH];
		bool PlaneBit;
		GlobalBlockIndex_to_GlobalPageIndex(GBI, GPI);
		GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
		FastRetrieveFromNAND(RA);
		ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription, PlaneBit);
		if (PageSerial == UINT64_MAX) {
			return false;
		}
		else {
			return true;
		}
	}

	// Fail-safe returns true
	return true;
}

uint16_t SPINANDClass::GetNumFreeBlocks() {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return 0; } // Checks initialization
	uint16_t Free = 0;
	for (uint16_t GBI = 0; GBI < _NAND_BLOCKS_PER_DEVICE; GBI++) {
		if (!IsBlockUsed(GBI)) {
			Free++;
		}
	}
	return Free;
}

void SPINANDClass::FastWriteToCache(void *Buf, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit, bool RetainExistingDataInCache) {
	/*
	It's important to select the correct plane and the starting byte address (column address) to write
	Do not pre-transpose the ColumnAddress.
	[cmd code]  [3 dummy bits]  [1 plane select bit]  [12 bits column address]
	*/
	uint8_t cmd = 0x02;
	if (RetainExistingDataInCache) {
		cmd = 0x84;
	}
	uint8_t AddressToSend[2];
	uint8_t *Buffer = reinterpret_cast<uint8_t *>(Buf);
	Reverse_Endian(&ColumnAddress, AddressToSend, 2);
	AddressToSend[0] |= (PlaneBit ? 0x10 : 0x00);
	EnableWrite();
	SelectDevice();
	(*SPIPTR).transferfastwriteonly(&cmd, 1);
	(*SPIPTR).transferfastwriteonly(AddressToSend, 2);
	(*SPIPTR).transferfastwriteonly(Buffer, DataLength);
	UnselectDevice();
}


void SPINANDClass::FastReadFromCache(void *Buf, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit, bool RetainExistingDataInCache) {
	/*
	[cmd code]  [3 dummy bits]  [1 plane select bit]  [12 bits column address]  [8 dummy bits]
	*/
	uint8_t cmd = 0x03;
	if (RetainExistingDataInCache) {
		cmd = 0x30;
	}
	uint8_t dummybyte = 0x00;
	uint8_t AddressToSend[2];
	uint8_t *Buffer = reinterpret_cast<uint8_t *>(Buf);
	Reverse_Endian(&ColumnAddress, AddressToSend, 2);
	AddressToSend[0] |= (PlaneBit ? 0x10 : 0x00);

	SelectDevice();
	(*SPIPTR).transferfastwriteonly(&cmd, 1);
	(*SPIPTR).transferfastwriteonly(AddressToSend, 2);
	(*SPIPTR).transferfastwriteonly(&dummybyte, 1);
	(*SPIPTR).transferfast(Buffer, DataLength);
	UnselectDevice();
}


void SPINANDClass::FastCommitToNAND(uint32_t RowAddress) {
	uint8_t cmd = 0x10;
	uint32_t AddressToSend;
	Reverse_Endian(&RowAddress, &AddressToSend, 3);
	EnableWrite();
	SelectDevice();
	(*SPIPTR).transferfastwriteonly(&cmd, 1);
	(*SPIPTR).transferfastwriteonly(&AddressToSend, 3);
	UnselectDevice();
	WaitForReady();
}


void SPINANDClass::FastRetrieveFromNAND(uint32_t RowAddress) {
	uint8_t cmd = 0x13;
	uint32_t AddressToSend;
	Reverse_Endian(&RowAddress, &AddressToSend, 3);
	EnableWrite();
	SelectDevice();
	(*SPIPTR).transferfastwriteonly(&cmd, 1);
	(*SPIPTR).transferfastwriteonly(&AddressToSend, 3);
	UnselectDevice();
	WaitForReady();
}


bool SPINANDClass::FileWriteOpen(char *FileDescription, uint32_t &FirstPageGPI) {
	/* Find the first free page to allow writing.
	You can only have one FileWrite channel open at any time. Trying to open another channel will abandon the existing channel
	*/

	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization

	if (!_DeviceBlockMapBuilt) {
		BuildDeviceBlockMap();
	}

	if (!_HighestGlobalPageScanned) {
		ScanForHighestPage();
	}

	if (_FileWriteChannelOpen) {
		delete _FileWriteChannelFileDescription;
	}

	if (_HighestGlobalPageSerial == UINT64_MAX - 1) {
		/* Reached maximum page serial */
		return false;
	}

	if (_HighestGlobalFileSerial == UINT32_MAX - 1) {
		/* Reached maximum file serial. */
		return false;
	}


	uint32_t GPI = _HighestGlobalPageIndex;

	/* Find a free page. */
	size_t i = 0;
	while (i++ < _NAND_PAGES_PER_DEVICE) {
		GPI++;
		GPI %= _NAND_PAGES_PER_DEVICE;

		if (!IsPageUsed(GPI)) {
			/* Found an unused page. Will write to it on the next FileWrite call. */
			if (!_FileWriteChannelOpen) {
				_FileWriteChannelFirstPageGPI = GPI;
				_FileWriteChannelFileDescription = new char[NANDFILEDESCRIPTIONLENGTH];
				for (int j = 0; j < NANDFILEDESCRIPTIONLENGTH; j++) {
					_FileWriteChannelFileDescription[j] = FileDescription[j];
				}
				_FileWriteChannelTotalBytesWritten = 0;
				_FileWriteChannelFilePageIndex = 0;
				_HighestGlobalFileSerial++;
				_FileWriteChannelOpen = true;
				FirstPageGPI = _FileWriteChannelFirstPageGPI;
				return true;
			}
		}
	}

	/* If we cannot find an empty page after checking all pages on the device, return false. */
	FirstPageGPI = UINT32_MAX;
	_FileWriteChannelOpen = false;
	return false;

}


void SPINANDClass::FileWriteClose() {
	/* This closes the currently open file channel

	*/
	if (_FileWriteChannelOpen) {
		delete _FileWriteChannelFileDescription;
		_FileWriteChannelOpen = false;
	}

}



uint32_t SPINANDClass::FileWrite(void *Buf, uint32_t BufferLength) {
	/* Automatically find enough free pages and write data into them.
	At this time, sequential write is required.
	Writing to a new unused block must also start on the first page of that block.

	You must first open a file channel with FileOpen before you can write to it.

	Warning: You will still use up one whole page if you write less than 2048 bytes for each call.

	*/

	if (!_FileWriteChannelOpen) {
		/* You MUST have a file channel open to write anything. */
		return 0;
	}

	uint64_t PageSerial;
	uint32_t GPI = _FileWriteChannelFirstPageGPI;
	uint32_t BytesWritten = 0; // Bytes written for this call only. _FileWriteChannelTotalBytesWritten is for all calls while the channel is open
	uint32_t RA;
	uint16_t GBI, GBIPREV = UINT16_MAX;
	uint8_t Die;
	bool PlaneBit;

	/* Create a local buffer */
	uint8_t *Buffer = reinterpret_cast<uint8_t *>(Buf);
	uint16_t CurrentDataLength;

	/* Find a free page. */
	size_t i = 0;
	while (i++ < _NAND_PAGES_PER_DEVICE) {
		CurrentDataLength = BufferLength - BytesWritten;
		if (CurrentDataLength > _NAND_BYTES_PER_PAGE) {
			CurrentDataLength = _NAND_BYTES_PER_PAGE;
		}

		if (!CurrentDataLength) {
			/* Finished writing everything. */
			break;
		}

		if (_HighestGlobalPageSerial == UINT64_MAX - 1) {
			/* Reached maximum page serial */
			FileWriteClose();
			break;
		}


		if (IsPageUsed(GPI)) {
			GPI++;
			GPI %= _NAND_PAGES_PER_DEVICE;
		}
		else {
			/* Found an unused page. Will write to it. */
			GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
			PageSerial = _HighestGlobalPageSerial + 1;

			//SerialUSB.print("NANDWriting ");
			//SerialUSB.println(CurrentDataLength);
			FastWriteToCache(Buffer, CurrentDataLength, 0, PlaneBit, false);

			WriteMetadata(PageSerial, _HighestGlobalFileSerial, _FileWriteChannelFilePageIndex, CurrentDataLength, REGULARFILE_FILETYPE, _FileWriteChannelFileDescription, PlaneBit);

			//char TempBuf[2132];
			//FastReadFromCache(TempBuf, 2132, 0, PlaneBit);

			FastCommitToNAND(RA);

			/* Update indices after writing */
			_HighestGlobalPageSerial = PageSerial;
			_HighestGlobalPageIndex = GPI;
			BytesWritten += CurrentDataLength;
			_FileWriteChannelTotalBytesWritten += CurrentDataLength;
			Buffer += CurrentDataLength;
			_FileWriteChannelFilePageIndex++;


			/* Mark block as used if necessary */
			GlobalPageIndex_to_GlobalBlockIndex(GPI, GBI);
			if (GBI != GBIPREV) {
				GBIPREV = GBI;
				MarkBlockUsed(GBI);
			}
		}
	}
	return BytesWritten;
}


bool SPINANDClass::FileReadOpen(uint32_t FirstPageGPI) {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization

	if (!_DeviceBlockMapBuilt) {
		BuildDeviceBlockMap();
	}

	if (!_HighestGlobalPageScanned) {
		ScanForHighestPage();
	}
	uint16_t BytesUsedInThisPage;
	bool PlaneBit;
	RetrieveMetadata(FirstPageGPI, _FileReadChannelFileSerial, _FileReadChannelFilePageIndex, BytesUsedInThisPage, PlaneBit);

	if (_FileReadChannelFileSerial == UINT32_MAX) {
		/* There is no file at this GPI */
		_FileReadChannelOpen = false;
		return false;
	}
	if (_FileReadChannelFilePageIndex != 0) {
		/* The GPI indicated is not the start of file. */
		_FileReadChannelOpen = false;
		return false;
	}

	_FileReadChannelFirstPageGPI = FirstPageGPI;
	_FileReadChannelTotalBytesRead = 0;
	_FileReadChannelCurrentPageStartPos = 0;
	_FileReadChannelOpen = true;
	return true;
}


void SPINANDClass::FileReadClose() {
	/* This closes the currently open FileRead channel. */
	if (_FileReadChannelOpen) {
		_FileReadChannelOpen = false;
	}
}


uint32_t SPINANDClass::FileRead(void *Buf, uint32_t BufferLength) {
	/* Read data from NAND starting for BufferLength and place into Buf
	Seek data based on page serial
	Must call FileReadOpen first
	*/

	if (!_FileReadChannelOpen) {
		return 0;
	}

	uint32_t FileSerial;
	uint32_t GPI;
	uint32_t FilePageIndex;
	uint16_t BytesUsedInThisPage;
	uint32_t RA;
	uint8_t Die;
	bool PlaneBit;

	/* Create a byte counter and a buffer */
	size_t BytesRead = 0;
	uint8_t *Buffer = reinterpret_cast<uint8_t *>(Buf);
	size_t CurrentDataLength;

	/* Read */
	size_t PagesTraversed = 0;
	GPI = (_FileReadChannelFirstPageGPI + _FileReadChannelFilePageIndex) % _NAND_PAGES_PER_DEVICE;


	while (PagesTraversed < _NAND_PAGES_PER_DEVICE) {
		CurrentDataLength = BufferLength - BytesRead;

		if (!CurrentDataLength) {
			/* Finished reading everything. */
			break;
		}

		GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
		RetrieveMetadata(GPI, FileSerial, FilePageIndex, BytesUsedInThisPage, PlaneBit);


		if (BytesUsedInThisPage > _NAND_BYTES_PER_PAGE) {
			BytesUsedInThisPage = _NAND_BYTES_PER_PAGE;
		}


		if (FileSerial == _FileReadChannelFileSerial && FilePageIndex == _FileReadChannelFilePageIndex) {
			/* The page content was already loaded into cache by RetrieveMetadata

			Given:
			This page contains at least 0 bytes of the file we requested

			Scenarios:
			Fully reading all used bytes in this page
			Partially reading from start
			Partially reading from somewhere not the start
			Reading across pages

			*/

			if (_FileReadChannelCurrentPageStartPos + CurrentDataLength > BytesUsedInThisPage) {
				/* Fulfilling the current request would exceed the available data in this page.
				Solution: Read partially.
				*/
				CurrentDataLength = BytesUsedInThisPage - _FileReadChannelCurrentPageStartPos;
			}

			if (CurrentDataLength > 0) {
				FastReadFromCache(Buffer, CurrentDataLength, _FileReadChannelCurrentPageStartPos, PlaneBit, false);

				BytesRead += CurrentDataLength;
				Buffer += CurrentDataLength; // Move the Buffer pointer
				_FileReadChannelTotalBytesRead += CurrentDataLength;
				_FileReadChannelCurrentPageStartPos += CurrentDataLength;
			}

			if (_FileReadChannelCurrentPageStartPos >= BytesUsedInThisPage) {
				/* Exhausted all bytes in the current page.
				Move onto the next page
				*/
				GPI++;
				GPI %= _NAND_PAGES_PER_DEVICE;
				_FileReadChannelFilePageIndex++;
				PagesTraversed++;
				_FileReadChannelCurrentPageStartPos = 0;
			}


		}
		else {
			/* That wasn't the page with the correct file serial or file page index. This could mean
			a) A part of the file has been erased by an erase block command
			b) The caller is trying to retrieve more than the file
			Either way, we should give up since it is corrupted.
			*/
			FileReadClose();
			return BytesRead;
		}
	}

	if (PagesTraversed >= _NAND_PAGES_PER_DEVICE) {
		/* There is no possible way there can be more pages to read for this file. Close the channel. */
		FileReadClose();
	}

	return BytesRead;
}



uint32_t SPINANDClass::FindFiles(uint32_t *FileGPICache, uint32_t CacheSizeSlots, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength, uint32_t &FirstMatchingFileFPGPI, uint32_t &LastMatchingFileFPGPI) {
	/* Returns the number of slots used in FileGPICache */

	if (!_DeviceBlockMapBuilt) {
		BuildDeviceBlockMap();
	}

	if (!_HighestGlobalPageScanned) {
		ScanForHighestPage();
	}


	uint32_t GPISTART;
	uint32_t CacheSlot = 0;
	for (uint16_t GBI = 0; GBI < _NAND_BLOCKS_PER_DEVICE; GBI++) {
		if (!IsBlockUsed(GBI)) {
			continue;
		}

		GlobalBlockIndex_to_GlobalPageIndex(GBI, GPISTART);

		if (GPISTART < MinimumFileGPI && GPISTART + _NAND_PAGES_PER_BLOCK < MinimumFileGPI) {
			continue;
		}

		for (uint32_t GPI = GPISTART; GPI < GPISTART + _NAND_PAGES_PER_BLOCK; GPI++) {

			if (GPI < MinimumFileGPI) {
				continue;
			}

			uint64_t PageSerial;
			uint32_t FileSerial;
			uint32_t FilePageIndex;
			uint16_t BytesUsedInThisPage;
			uint8_t FileType;
			char FileDescription[NANDFILEDESCRIPTIONLENGTH];
			RetrieveMetadata(GPI, PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription);

			if (PageSerial == UINT64_MAX || FileSerial == UINT32_MAX) {
				break;
			}

			uint64_t PageSerialMax = 0;
			uint64_t PageSerialMin = UINT64_MAX;

			if (FilePageIndex == 0 && memcmp(FileDescription, MatchingFileDescription, MatchLength) == 0) {
				// This file is a match

				if (CacheSizeSlots && CacheSlot < CacheSizeSlots) {
					FileGPICache[CacheSlot] = GPI;
				}
				if (!CacheSizeSlots) {
					if (PageSerial > PageSerialMax) {
						PageSerialMax = PageSerial;
						LastMatchingFileFPGPI = GPI;
						//SerialUSB.print("Debug: This file ");
						//SerialUSB.print(GPI);
						//SerialUSB.print(" ");
						//for (int cc = 0; cc < NANDFILEDESCRIPTIONLENGTH; cc++) {
						//	SerialUSB.print((uint8_t)FileDescription[cc], HEX);
						//	SerialUSB.print(" ");
						//}
						//SerialUSB.print(" matches string ");
						//for (int cc = 0; cc < NANDFILEDESCRIPTIONLENGTH; cc++) {
						//	SerialUSB.print((uint8_t)MatchingFileDescription[cc], HEX);
						//	SerialUSB.print(" ");
						//}
						//SerialUSB.println("");
						//SerialUSB.print("MatchLength = ");
						//SerialUSB.println(MatchLength);
					}
					if (PageSerial < PageSerialMin) {
						PageSerialMin = PageSerial;
						FirstMatchingFileFPGPI = GPI;
					}
				}
				CacheSlot++;
			}
			if (CacheSizeSlots && CacheSlot == CacheSizeSlots) { break; }

		}
		if (CacheSizeSlots && CacheSlot == CacheSizeSlots) { break; }

	}

	return CacheSlot;
}


uint32_t SPINANDClass::FindFiles(uint32_t *FileGPICache, uint32_t CacheSizeSlots, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength) {
	uint32_t FirstMatchingFileFPGPI;
	uint32_t LastMatchingFileFPGPI;
	return FindFiles(FileGPICache, CacheSizeSlots, MinimumFileGPI, MatchingFileDescription, MatchLength, FirstMatchingFileFPGPI, LastMatchingFileFPGPI);
}

uint32_t SPINANDClass::FindFiles(uint32_t *FileGPICache, uint32_t CacheSizeSlots, uint32_t MinimumFileGPI) {
	char MatchingFileDescription[1];
	size_t MatchLength = 0;
	uint32_t FirstFileGPI;
	uint32_t LastFileGPI;
	return FindFiles(FileGPICache, CacheSizeSlots, MinimumFileGPI, MatchingFileDescription, MatchLength, FirstFileGPI, LastFileGPI);
}



uint32_t SPINANDClass::FindFiles(std::vector <uint32_t> &FileGPICache, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength, uint32_t &FirstMatchingFileFPGPI, uint32_t &LastMatchingFileFPGPI) {
	uint32_t n = CountFiles(0, MatchingFileDescription, MatchLength);
	FileGPICache.resize(n);
	while (n && FileGPICache.size() < n) {
		SerialUSB.print("Warning in SPINANDClass::FindFiles: Failure to allocate ");
		SerialUSB.print(n);
		SerialUSB.print(" slots. Trying to allocate ");
		n--;
		SerialUSB.print(n);
		SerialUSB.println(" slots instead.");
		FileGPICache.resize(n);
	}
	if (n) {
		return FindFiles(FileGPICache.data(), n, MinimumFileGPI, MatchingFileDescription, MatchLength);
	}
	else {
		return 0;
	}
}

uint32_t SPINANDClass::FindFiles(std::vector <uint32_t> &FileGPICache, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength) {
	uint32_t FirstMatchingFileFPGPI;
	uint32_t LastMatchingFileFPGPI;
	return FindFiles(FileGPICache, MinimumFileGPI, MatchingFileDescription, MatchLength, FirstMatchingFileFPGPI, LastMatchingFileFPGPI);
}

uint32_t SPINANDClass::FindFiles(std::vector <uint32_t> &FileGPICache, uint32_t MinimumFileGPI) {
	char MatchingFileDescription[1];
	size_t MatchLength = 0;
	uint32_t FirstFileGPI;
	uint32_t LastFileGPI;
	FindFiles(FileGPICache, MinimumFileGPI, MatchingFileDescription, MatchLength, FirstFileGPI, LastFileGPI);
}

bool SPINANDClass::FindFirstAndLastFile(uint32_t &FirstFileGPI, uint32_t &LastFileGPI, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength) {
	/* Find the first and last file that matches the search criteria. First and last refer to the files with the lowest and highest PageSerial. */
	uint32_t *FileGPICache = nullptr;
	if (FindFiles(FileGPICache, 0, MinimumFileGPI, MatchingFileDescription, MatchLength, FirstFileGPI, LastFileGPI) > 0) {
		return true;
	}
	return false;
}


uint32_t SPINANDClass::CountFiles(uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength) {
	/* Returns the number of files matching the search criteria, without outputting the file GPIs */
	uint32_t *FileGPICache = nullptr;
	uint32_t FirstFileGPI;
	uint32_t LastFileGPI;
	return FindFiles(FileGPICache, 0, MinimumFileGPI, MatchingFileDescription, MatchLength, FirstFileGPI, LastFileGPI);
}


uint32_t SPINANDClass::CountFiles(uint32_t MinimumFileGPI) {
	/* Returns the number of files matching the search criteria, without outputting the file GPIs */
	uint32_t *FileGPICache = nullptr;
	return FindFiles(FileGPICache, 0, MinimumFileGPI);
}


bool SPINANDClass::FileGetInfo(uint32_t FirstPageGPI, char *FileDescription) {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	uint64_t PageSerial;
	uint32_t FileSerial;
	uint32_t FilePageIndex;
	uint16_t BytesUsedInThisPage;
	uint8_t FileType;
	RetrieveMetadata(FirstPageGPI, PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription);
	return true;
}


bool SPINANDClass::FileGetInfo(uint32_t FirstPageGPI, char *FileDescription, uint32_t &FileSize, uint32_t &PagesUsed) {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	uint64_t PageSerial;
	uint32_t FileSerial;
	return FileGetInfo(FirstPageGPI, PageSerial, FileSerial, FileDescription, FileSize, PagesUsed);
}

bool SPINANDClass::FileGetInfo(uint32_t FirstPageGPI, uint64_t &PageSerial, uint32_t &FileSerial, char *FileDescription, uint32_t &FileSize, uint32_t &PagesUsed) {
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	uint64_t PageSerial_Test;
	uint32_t FileSerial_Test;
	uint32_t FilePageIndex, FilePageIndex_Test;
	uint32_t GPI;
	uint16_t BytesUsedInThisPage;
	uint8_t FileType;
	char FileDesc[NANDFILEDESCRIPTIONLENGTH];
	FileSize = 0;
	PagesUsed = 0;
	RetrieveMetadata(FirstPageGPI, PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription);

	if (PageSerial == UINT64_MAX || FileSerial == UINT32_MAX) {
		return false;
	}

	if (FilePageIndex != 0) {
		// The GPI provided was not of the first page. Backtrack by FPI to find it.
		FirstPageGPI -= FilePageIndex;
		FirstPageGPI %= _NAND_PAGES_PER_DEVICE;
		RetrieveMetadata(FirstPageGPI, PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription);
	}

	FilePageIndex_Test = 0;
	for (uint32_t i = 0; i < _NAND_PAGES_PER_DEVICE; i++) {
		// Calculate file size
		GPI = (FirstPageGPI + i) % _NAND_PAGES_PER_DEVICE;
		RetrieveMetadata(GPI, PageSerial_Test, FileSerial_Test, FilePageIndex, BytesUsedInThisPage, FileType, FileDesc);
		if (FileSerial_Test == FileSerial && FilePageIndex_Test++ == FilePageIndex) {
			// All pages used by the same file will have the same FileSerial and increasing FilePageIndex
			FileSize += BytesUsedInThisPage;
			PagesUsed++;
		}
		else {
			break;
		}
	}

	return true;
}



size_t SPINANDClass::GetFileDescriptionLength() {
	return NANDFILEDESCRIPTIONLENGTH;
}


bool SPINANDClass::IsPageUsed(uint32_t GPI) {
	/* Checks if a page is used.
	Warning: Overwrites the cache.
	Side effect: Can also use this to read into cache.
	*/
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return false; } // Checks initialization
	uint64_t PageSerial;
	uint32_t FileSerial;
	RetrieveMetadata(GPI, PageSerial, FileSerial);
	if (PageSerial == UINT64_MAX) {
		return false;
	}
	else {
		return true;
	}
}


bool SPINANDClass::GetPlaneBit(uint32_t RA) {
	return ((RA & 0b1000000) >> 6);
}

//void SPINANDClass::WritePageSerial(uint64_t PageSerial, bool PlaneBit) {
//	/* Writes the page serial to the cache */
//	FastWriteToCache(&PageSerial, 8, 0x0838, PlaneBit);
//}


void SPINANDClass::WriteMetadata(uint64_t PageSerial, uint32_t FileSerial, uint32_t FilePageIndex, uint16_t BytesUsedInThisPage, uint8_t FileType, char *FileDescription, bool PlaneBit) {
	/* This writes the metadata to the cache. Do this after writing content to cache

	offset	address				allocation		purpose
	0		0x0820 -- 0x082c	(13 bytes)		File Description
	13		0x082d -- 0x082d	(1 byte)		File Type
	14		0x082e -- 0x082f	(2 bytes)		Number of bytes written in this page
	16		0x0830 -- 0x0833	(4 bytes)		File Page Index
	20		0x0834 -- 0x0837	(4 bytes)		Global File Serial
	24		0x0838 -- 0x083f	(8 bytes)		Global Page Serial


	Indicate unspecified input arguments with 0xff's
	*/

	uint8_t *filedesc = reinterpret_cast<uint8_t *>(FileDescription);
	uint8_t *bytesused = reinterpret_cast<uint8_t *>(&BytesUsedInThisPage);
	uint8_t *filepageindex = reinterpret_cast<uint8_t *>(&FilePageIndex);
	uint8_t *fileserial = reinterpret_cast<uint8_t *>(&FileSerial);
	uint8_t *pageserial = reinterpret_cast<uint8_t *>(&PageSerial);
	uint8_t Metadata[32];
	int i;
	for (i = 0; i < NANDFILEDESCRIPTIONLENGTH; i++) { Metadata[i] = filedesc[i]; }
	Metadata[13] = FileType;
	for (i = 0; i < 2; i++) { Metadata[14 + i] = bytesused[i]; }
	for (i = 0; i < 4; i++) { Metadata[16 + i] = filepageindex[i]; }
	for (i = 0; i < 4; i++) { Metadata[20 + i] = fileserial[i]; }
	for (i = 0; i < 8; i++) { Metadata[24 + i] = pageserial[i]; }

	FastWriteToCache(Metadata, 32, 0x0820, PlaneBit, true);
}


void SPINANDClass::ReadMetadata(uint64_t &PageSerial, uint32_t &FileSerial, uint32_t &FilePageIndex, uint16_t &BytesUsedInThisPage, uint8_t &FileType, char *FileDescription, bool PlaneBit) {
	/* Read the page serial, file serial, file page index, bytesusedinthispage, file type, and file description from the cache.
	Data must already be loaded from NAND to cache.
	*/
	uint8_t Metadata[32];
	uint8_t *filedesc = reinterpret_cast<uint8_t *>(FileDescription);
	uint8_t *bytesused = reinterpret_cast<uint8_t *>(&BytesUsedInThisPage);
	uint8_t *filepageindex = reinterpret_cast<uint8_t *>(&FilePageIndex);
	uint8_t *fileserial = reinterpret_cast<uint8_t *>(&FileSerial);
	uint8_t *pageserial = reinterpret_cast<uint8_t *>(&PageSerial);
	FastReadFromCache(Metadata, 32, 0x0820, PlaneBit, false);
	int i;
	for (i = 0; i < NANDFILEDESCRIPTIONLENGTH; i++) { filedesc[i] = Metadata[i]; }
	FileType = Metadata[13];
	for (i = 0; i < 2; i++) { bytesused[i] = Metadata[14 + i]; }
	for (i = 0; i < 4; i++) { filepageindex[i] = Metadata[16 + i]; }
	for (i = 0; i < 4; i++) { fileserial[i] = Metadata[20 + i]; }
	for (i = 0; i < 8; i++) { pageserial[i] = Metadata[24 + i]; }
}


void SPINANDClass::RetrieveMetadata(uint32_t GPI, uint64_t &PageSerial, uint32_t &FileSerial, uint32_t &FilePageIndex, uint16_t &BytesUsedInThisPage, uint8_t &FileType, char *FileDescription) {
	uint32_t RA;
	uint8_t Die;
	bool PlaneBit;
	GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
	FastRetrieveFromNAND(RA);
	ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription, PlaneBit);
}


void SPINANDClass::RetrieveMetadata(uint32_t GPI, uint32_t &FileSerial, uint32_t &FilePageIndex, uint16_t &BytesUsedInThisPage, bool &PlaneBit) {
	uint64_t PageSerial;
	uint32_t RA;
	uint8_t Die;
	char FileDescription[NANDFILEDESCRIPTIONLENGTH];
	GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
	FastRetrieveFromNAND(RA);
	ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, Die, FileDescription, PlaneBit);
}


void SPINANDClass::RetrieveMetadata(uint32_t GPI, uint64_t &PageSerial, uint32_t &FileSerial) {
	uint32_t FilePageIndex, RA;
	uint16_t BytesUsedInThisPage;
	uint8_t Die, FileType;
	bool PlaneBit;
	char FileDescription[NANDFILEDESCRIPTIONLENGTH];
	GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
	FastRetrieveFromNAND(RA);
	ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription, PlaneBit);
}



void SPINANDClass::GlobalPageIndex_to_PlaneBit(uint32_t GlobalPageIndex, bool &PlaneBit) {
	PlaneBit = ((GlobalPageIndex & 0b1000000) >> 6);
}


void SPINANDClass::GlobalPageIndex_to_RA_and_Die(uint32_t GlobalPageIndex, uint32_t &RowAddress, uint8_t &DieIndex) {
	DieIndex = GlobalPageIndex / _NAND_PAGES_PER_DIE;
	RowAddress = GlobalPageIndex % _NAND_PAGES_PER_DIE;
}


void SPINANDClass::GlobalBlockIndex_to_RA_and_Die(uint16_t GlobalBlockIndex, uint32_t &RowAddress, uint8_t &DieIndex) {
	DieIndex = GlobalBlockIndex / _NAND_BLOCKS_PER_DIE;
	RowAddress = (GlobalBlockIndex % _NAND_BLOCKS_PER_DIE) * _NAND_PAGES_PER_BLOCK;
}


void SPINANDClass::GlobalBlockIndex_to_GlobalPageIndex(uint16_t GlobalBlockIndex, uint32_t &GlobalPageIndex) {
	GlobalPageIndex = _NAND_PAGES_PER_BLOCK * GlobalBlockIndex;
}

void SPINANDClass::GlobalPageIndex_to_GlobalBlockIndex(uint32_t GlobalPageIndex, uint16_t &GlobalBlockIndex) {
	GlobalBlockIndex = GlobalPageIndex / _NAND_PAGES_PER_BLOCK;
	/* Sanitize GlobalBlockIndex since it cannot exceed _NAND_BLOCKS_PER_DEVICE */
	//GlobalBlockIndex &= (_NAND_BLOCKS_PER_DEVICE - 1); 
	GlobalBlockIndex %= _NAND_BLOCKS_PER_DEVICE;
}


void SPINANDClass::GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(uint32_t GlobalPageIndex, uint32_t &RowAddress, uint8_t &DieIndex, bool &PlaneBit) {
	GlobalPageIndex_to_RA_and_Die(GlobalPageIndex, RowAddress, DieIndex);
	if (DieIndex != _SelectedDie) {
		SelectDie(DieIndex);
	}
	PlaneBit = GetPlaneBit(RowAddress);
}

void SPINANDClass::BuildDeviceBlockMap() {
	/* Builds a bitmap of used/unused blocks. Only the first page in each block is checked.
	Also checks the page serial of the first page in each block.
	*/
	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return; } // Checks initialization

	uint64_t PageSerial;
	uint32_t GPI, RA, FileSerial, FilePageIndex;
	uint16_t BytesUsedInThisPage;
	uint8_t Die, FileType;
	char FileDescription[NANDFILEDESCRIPTIONLENGTH];
	bool PlaneBit, Found = 0;

	for (size_t GBI = 0; GBI < _NAND_BLOCKS_PER_DEVICE; GBI++) {
		GlobalBlockIndex_to_GlobalPageIndex(GBI, GPI);
		GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
		FastRetrieveFromNAND(RA);
		ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription, PlaneBit);
		if (PageSerial == UINT64_MAX) {
			MarkBlockUnused(GBI);
		}
		else {
			MarkBlockUsed(GBI);
		}
		if (PageSerial != UINT64_MAX && (PageSerial > _HighestGlobalPageSerial || _HighestGlobalPageSerial == UINT64_MAX)) {
			_HighestGlobalPageSerial = PageSerial;
			_HighestGlobalPageIndex = GPI;
			Found = 1;
		}
	}
	if (!Found) {
		_HighestGlobalPageSerial = UINT64_MAX;
		//_HighestGlobalPageIndex = UINT32_MAX;
		_HighestGlobalPageIndex = InitializeHighestGPI();
	}
	_DeviceBlockMapBuilt = true;
}


void SPINANDClass::ScanForHighestPage() {
	/*	Scan for the used page with the highest serial number in the entire device and record its serial number and index location.
	Also scan for the file with the highest serial number in the entire device and record its serial number.
	Requires DeviceBlockMap to be built first. Otherwise it will be built.
	This function will only scan the pages in the block with the highest starting serial.
	*/

	if (!_NAND_BYTES_PER_PAGE || !_HasPinControl) { return; } // Checks initialization

	if (!_DeviceBlockMapBuilt) {
		BuildDeviceBlockMap();
	}

	bool PlaneBit, Found = 0;
	uint64_t PageSerial;
	uint32_t FileSerial;
	uint32_t RA, FilePageIndex;
	uint16_t BytesUsedInThisPage;
	uint8_t FileType, Die;
	uint32_t GPIMAXX = _HighestGlobalPageIndex + _NAND_PAGES_PER_BLOCK;
	char FileDescription[NANDFILEDESCRIPTIONLENGTH];

	for (uint32_t GPI = _HighestGlobalPageIndex; GPI < GPIMAXX; GPI++) {
		GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(GPI, RA, Die, PlaneBit);
		FastRetrieveFromNAND(RA);
		ReadMetadata(PageSerial, FileSerial, FilePageIndex, BytesUsedInThisPage, FileType, FileDescription, PlaneBit);
		if (PageSerial != UINT64_MAX && (PageSerial > _HighestGlobalPageSerial || _HighestGlobalPageSerial == UINT64_MAX)) {
			_HighestGlobalPageSerial = PageSerial;
			_HighestGlobalPageIndex = GPI;
		}
		if (FileSerial != UINT32_MAX && (FileSerial > _HighestGlobalFileSerial || _HighestGlobalFileSerial == UINT32_MAX)) {
			_HighestGlobalFileSerial = FileSerial;
		}
	}
	_HighestGlobalPageScanned = true;
}




















//void SPINANDClass::DefaultWriteToCache(uint8_t *Buffer, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit) {
//	uint8_t AddressToSend[2];
//	Reverse_Endian(&ColumnAddress, AddressToSend, 2);
//	AddressToSend[0] |= (PlaneBit ? 0x10 : 0x00);
//
//	EnableWrite();
//	SelectDevice();
//	(*SPIPTR).transfer(0x02);
//	(*SPIPTR).transfer(AddressToSend, 2);
//	for (size_t i = 0; i < DataLength; i++) {
//		(*SPIPTR).transfer(Buffer[i]);
//	}
//	UnselectDevice();
//}
//
//
//void SPINANDClass::DefaultCommitToNAND(uint32_t RowAddress) {
//	uint32_t AddressToSend;
//	Reverse_Endian(&RowAddress, &AddressToSend, 3);
//	EnableWrite();
//	SelectDevice();
//	(*SPIPTR).transfer(0x10);
//	(*SPIPTR).transfer(&AddressToSend, 3);
//	UnselectDevice();
//	WaitForReady();
//}
//
//
//void SPINANDClass::DefaultRetrieveFromNAND(uint32_t RowAddress) {
//	uint32_t AddressToSend;
//	Reverse_Endian(&RowAddress, &AddressToSend, 3);
//	EnableWrite();
//	SelectDevice();
//	(*SPIPTR).transfer(0x13);
//	(*SPIPTR).transfer(&AddressToSend, 3);
//	UnselectDevice();
//	WaitForReady();
//}
//
//
//void SPINANDClass::DefaultReadFromCache(uint8_t *Buffer, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit) {
//	uint8_t dummybyte = 0x00;
//	uint8_t AddressToSend[2];
//	Reverse_Endian(&ColumnAddress, AddressToSend, 2);
//	AddressToSend[0] |= (PlaneBit ? 0x10 : 0x00);
//
//	SelectDevice();
//	(*SPIPTR).transfer(0x03);
//	(*SPIPTR).transfer(AddressToSend, 2); 
//	(*SPIPTR).transfer(&dummybyte, 1);
//	(*SPIPTR).transfer(Buffer, DataLength);
//	UnselectDevice();
//}




uint8_t SPINANDClass::SelfTest() {
	if (!_NAND_BYTES_PER_PAGE) { return ERRORCODES::NOT_INITIALIZED; } // Not initialized
	if (!_HasPinControl) { return ERRORCODES::NO_PIN_CONTROL; } // No control
	if (_FileWriteChannelOpen) {
		return ERRORCODES::WRITE_CHANNEL_IS_OPEN;
	}
	if (_FileReadChannelOpen) {
		return ERRORCODES::READ_CHANNEL_IS_OPEN;
	}
	uint8_t *Buffer = new uint8_t[_NAND_BYTES_PER_PAGE];
	if (Buffer == nullptr) {
		delete Buffer;
		return ERRORCODES::TEST_BUFFER_ALLOCATION_ERROR;
	}
	get_JEDEC_ID(Buffer, 2);
	if (Buffer[0] != 0x2C || Buffer[1] != 0x36) {
		SerialUSB.print("JEDEC ID = ");
		SerialUSB.print(Buffer[0], HEX);
		SerialUSB.print(" ");
		SerialUSB.print(Buffer[1], HEX);
		SerialUSB.println("");
		delete Buffer;
		return ERRORCODES::INCORRECT_JEDEC_ID; // Incorrect JEDEC ID response
	}
	//ResetRegister();
	//get_SRAM_Mode_Register(Buffer1[0]);
	//if (Buffer1[0] != 0) {
	//	SerialUSB.print("Unexpected SRAM register: ");
	//	SerialUSB.println(Buffer1[0], DEC);
	//	return ERRORCODES::INCORRECT_STATUS_REGISTER; // Wrong values in the status register
	//}

	// Test the cache's one-shot read/write capability
	for (int j = 0; j < 4; j++) {
		size_t match = 0;
		bool PlaneBit = (bool)(j%2);
		SelectDie(j / 2);

		// Create a pseudorandom binary sequence
		uint16_t seed = random(0xffff);
		build_test_sequence(Buffer, _NAND_BYTES_PER_PAGE, seed, true);

		//SerialUSB.print("Writing to cache: ");
		//for (size_t i = 0; i < 10; i++) {
		//	SerialUSB.print(Buffer[i]);
		//	SerialUSB.print(" ");
		//}
		//SerialUSB.println("");

		FastWriteToCache(Buffer, _NAND_BYTES_PER_PAGE, 0, PlaneBit, false);

		for (size_t i = 0; i < _NAND_BYTES_PER_PAGE; i++) {
			Buffer[i] = 0;
		}

		FastReadFromCache(Buffer, _NAND_BYTES_PER_PAGE, 0, PlaneBit, false);
		//SerialUSB.print("Getting from cache: ");
		//for (size_t i = 0; i < 10; i++) {
		//	SerialUSB.print(Buffer[i]);
		//	SerialUSB.print(" ");
		//}
		//SerialUSB.println("");

		// Compare
		match = build_test_sequence(Buffer, _NAND_BYTES_PER_PAGE, seed, false);

		if (match != _NAND_BYTES_PER_PAGE) {
			delete Buffer;
			switch (j) {
			case 0: 
				return ERRORCODES::CACHE_ERROR_DIE_0_PLANE_0;
			case 1: 
				return ERRORCODES::CACHE_ERROR_DIE_0_PLANE_1;
			case 2: 
				return ERRORCODES::CACHE_ERROR_DIE_1_PLANE_0;
			case 3: 
				return ERRORCODES::CACHE_ERROR_DIE_1_PLANE_1;
			}
		}
	}


	delete Buffer;
	return ERRORCODES::OK; // Self test passed
}



size_t SPINANDClass::build_test_sequence(uint8_t *Buffer, size_t Length, uint16_t seed, bool write) {
	// Adapted from http://www.wattnotions.com/prbs-array-c-code-generator-pretty-nasty/
	uint8_t *p = Buffer;
	uint16_t a = seed;
	size_t match = 0;
	bool bit;
	if (Length == 0) {
		return 0;
	}
	unsigned period = 0;
	do {
		bit = ((a >> 0) ^ (a >> 2)) & 1;
		a = (a >> 1) | (bit << 10);
		if (Length > 2) {
			if (write) {
				memcpy(p, &a, 2);
				match += 2;
			}
			else {
				if (!memcmp(p, &a, 2)) {
					match += 2;
				}
			}
			p += 2;
			Length -= 2;
		}
		else {
			if (write) {
				memcpy(p, &a, 1);
				match += 2;
			}
			else {
				if (!memcmp(p, &a, 1)) {
					match++;
				}
			}
			p++;
			Length--;
		}
	} while (Length);
	return match;

}


uint8_t SPINANDClass::Benchmark(uint32_t WriteBytes, BenchmarkResult &result) {
	if (!_NAND_BYTES_PER_PAGE) { return ERRORCODES::NOT_INITIALIZED; } // Not initialized
	if (!_HasPinControl) { return ERRORCODES::NO_PIN_CONTROL; } // No control
	if (_FileWriteChannelOpen) {
		return ERRORCODES::WRITE_CHANNEL_IS_OPEN;
	}
	if (_FileReadChannelOpen) {
		return ERRORCODES::READ_CHANNEL_IS_OPEN;
	}
	if (WriteBytes < 2) {
		return ERRORCODES::WRITEBYTES_TOO_SMALL;
	}
	if (WriteBytes > ((uint32_t)GetNumFreeBlocks() * (uint32_t)_NAND_PAGES_PER_BLOCK * (uint32_t)_NAND_BYTES_PER_PAGE)) {
		return ERRORCODES::NOT_ENOUGH_FREE_BLOCKS;
	}

	uint8_t *Buffer = new uint8_t[WriteBytes];
	if (Buffer == nullptr) {
		delete Buffer;
		return ERRORCODES::TEST_BUFFER_ALLOCATION_ERROR;
	}
	get_JEDEC_ID(Buffer, 2);
	if (Buffer[0] != 0x2C || Buffer[1] != 0x36) {
		SerialUSB.print("JEDEC ID = ");
		SerialUSB.print(Buffer[0], HEX);
		SerialUSB.print(" ");
		SerialUSB.print(Buffer[1], HEX);
		SerialUSB.println("");
		delete Buffer;
		return ERRORCODES::INCORRECT_JEDEC_ID; // Incorrect JEDEC ID response
	}



	uint32_t FPGPI, BytesWritten, BytesRead;
	uint32_t t_start;
	//uint32_t t_cmd;
	//uint32_t t_writeopen, t_write, t_writeclose;
	//uint32_t t_readopen, t_read, t_readclose;


	// Benchmark delay between two commands
	t_start = micros();
	result.t_cmd = micros() - t_start;

	// Benchmark "IsBlockUsed"
	t_start = micros();
	IsBlockUsed(0);
	result.t_isblockused = micros() - t_start;

	// Benchmark opening a write channel
	t_start = micros();
	if (!FileWriteOpen("BENCHMARK    ", FPGPI)) {
		delete Buffer;
		return ERRORCODES::FILEWRITEOPEN_FAIL;
	}
	result.t_writeopen = micros() - t_start;

	// Create a pseudorandom binary sequence
	uint16_t seed = random(0xffff);
	build_test_sequence(Buffer, WriteBytes, seed, true);

	// Benchmark writing
	t_start = micros();
	BytesWritten = FileWrite(Buffer, WriteBytes);
	result.t_write = micros() - t_start;

	// Benchmark closing the write channel
	t_start = micros();
	FileWriteClose();
	result.t_writeclose = micros() - t_start;


	// Verify correct number of bytes written
	if (BytesWritten != WriteBytes) {
		delete Buffer;
		return ERRORCODES::FILEWRITE_WRONG_BYTES_WRITTEN;
	}



	// Benchmark finding a file (not implemented yet)


	// Benchmark opening a read channel
	t_start = micros();
	if (!FileReadOpen(FPGPI)) {
		delete Buffer;
		return ERRORCODES::FILEREADOPEN_FAIL;
	}
	result.t_readopen = micros() - t_start;

	// Benchmark reading
	t_start = micros();
	BytesRead = FileRead(Buffer, WriteBytes);
	result.t_read = micros() - t_start;

	// Benchmark closing the read channel
	t_start = micros();
	FileReadClose();
	result.t_readclose = micros() - t_start;

	
	// Verify correct number of bytes read
	if (BytesRead != WriteBytes) {
		delete Buffer;
		return ERRORCODES::FILEREAD_WRONG_BYTES_READ;
	}


	// Verify read-back content
	size_t match = build_test_sequence(Buffer, WriteBytes, seed, false);
	if (match != WriteBytes) {
		delete Buffer;
		return ERRORCODES::INCORRECT_READ_BACK_CONTENT;
	}


	delete Buffer;
	return ERRORCODES::OK;
}



//

//bool strequal(const char *A, const char *B, size_t L) {
//	for (size_t i = 0; i < L; i++) {
//		if (A[i] != B[i]) {
//			return false;
//		}
//	}
//	return true;
//}
//
