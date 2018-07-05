/*
Handles the SPINAND chip, including writing/reading.
*/
#pragma once
#if defined(ARDUINO) && ARDUINO >= 100
	#include "arduino.h"
#else
	#include "WProgram.h"
#endif
#include <SPI.h>
#include "pins_arduino.h"
#include "wiring_private.h"
#undef max
#undef min
#include <vector>


#define NANDFILEDESCRIPTIONLENGTH 13


class SPINANDClass {
public:
	SPINANDClass();
	void Initialize(SERCOM *SPINAND_SERCOM,
		uint8_t SPINAND_CS_PIN, uint8_t SPINAND_MISO_PIN, uint8_t SPINAND_SCK_PIN, uint8_t SPINAND_MOSI_PIN,
		uint8_t SPI_DATAMODE,
		SercomSpiTXPad SPINAND_SERCOM_TXPAD, SercomRXPad SPINAND_SERCOM_RXPAD,
		_EPioType SPINAND_MISO_PIOTYPE, _EPioType SPINAND_SCK_PIOTYPE, _EPioType SPINAND_MOSI_PIOTYPE,
		size_t NAND_BYTES_PER_PAGE, size_t NAND_PAGES_PER_BLOCK, size_t NAND_BLOCKS_PER_DIE, size_t NAND_DIES_PER_DEVICE);
	~SPINANDClass();

	bool WakeupDevice();
	bool get_JEDEC_ID(uint8_t *Response, size_t Length);
	uint8_t get_blocklock();
	uint8_t get_config();
	uint8_t get_status();
	uint8_t get_dieselect();

	uint16_t EraseDevice();
	bool EraseDeviceFull();

	bool FileWriteOpen(char *FileDescription, uint32_t &FirstPageGPI);
	void FileWriteClose();
	uint32_t FileWrite(void *Buf, uint32_t BufferLength);

	bool FileReadOpen(uint32_t FirstPageGPI);
	void FileReadClose();
	uint32_t FileRead(void *Buf, uint32_t BufferLength);
	bool FileGetInfo(uint32_t FirstPageGPI, char *FileDescription);
	bool FileGetInfo(uint32_t FirstPageGPI, char *FileDescription, uint32_t &FileSize, uint32_t &PagesUsed);
	bool FileGetInfo(uint32_t FirstPageGPI, uint64_t &PageSerial, uint32_t &FileSerial, char *FileDescription, uint32_t &FileSize, uint32_t &PagesUsed);

	uint32_t FindFiles(uint32_t *FileGPICache, uint32_t CacheSizeSlots, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength, uint32_t &FirstMatchingFileFPGPI, uint32_t &LastMatchingFileFPGPI);
	uint32_t FindFiles(uint32_t *FileGPICache, uint32_t CacheSizeSlots, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength);
	uint32_t FindFiles(uint32_t *FileGPICache, uint32_t CacheSizeSlots, uint32_t MinimumFileGPI);
	uint32_t FindFiles(std::vector <uint32_t> &FileGPICache, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength, uint32_t &FirstMatchingFileFPGPI, uint32_t &LastMatchingFileFPGPI);
	uint32_t FindFiles(std::vector <uint32_t> &FileGPICache, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength);
	uint32_t FindFiles(std::vector <uint32_t> &FileGPICache, uint32_t MinimumFileGPI);
	uint32_t CountFiles(uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength);
	uint32_t CountFiles(uint32_t MinimumFileGPI);
	bool FindFirstAndLastFile(uint32_t &FirstFileGPI, uint32_t &LastFileGPI, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength);
	//bool FindLastFile(uint32_t &FileGPI, uint32_t MinimumFileGPI, const char *MatchingFileDescription, size_t MatchLength);
	size_t GetFileDescriptionLength();

	bool IsBlockUsed(uint16_t GBI);
	bool IsPageUsed(uint32_t GPI);

	uint16_t GetNumFreeBlocks();



	// Note: Self test will not actually write to the NAND itself. But it will still overwrite the cache.
	uint8_t SelfTest();

	struct BenchmarkResult {
		uint32_t t_cmd;
		uint32_t t_isblockused;
		uint32_t t_writeopen;
		uint32_t t_write;
		uint32_t t_writeclose;
		uint32_t t_readopen;
		uint32_t t_read;
		uint32_t t_readclose;
	};

	// Benchmark
	uint8_t Benchmark(uint32_t WriteBytes, BenchmarkResult &result);

	/*
	The proper procedure for handing off pin control is:
	1. First core: ReleasePinControl -- set pins to input or input_pullup
	2. First core sends a core-to-core command to second core
	3. Second core: AssumePinControl -- reconfigures mux to SPI mode and begins SPI
	*/
	void AssumePinControl();
	void ReleasePinControl();


	inline bool HighestGlobalPageScanned() {
		return _HighestGlobalPageScanned;
	}
	inline uint64_t HighestGlobalPageSerial() {
		return _HighestGlobalPageSerial;
	}
	inline uint32_t HighestGlobalPageIndex() {
		return _HighestGlobalPageIndex;
	}
	inline uint32_t HighestGlobalFileSerial() {
		return _HighestGlobalFileSerial;
	}
	inline uint32_t GBI_MAX() {
		return _NAND_BLOCKS_PER_DEVICE;
	}
	inline uint32_t GPI_MAX() {
		return _NAND_PAGES_PER_DEVICE;
	}


	// PRBS-11 generation/verification code
	size_t build_test_sequence(uint8_t *Buffer, size_t Length, uint16_t seed, bool write);



	enum ERRORCODES {
		OK,
		NOT_INITIALIZED,
		NO_PIN_CONTROL,
		TEST_BUFFER_ALLOCATION_ERROR,
		INCORRECT_JEDEC_ID,
		INCORRECT_STATUS_REGISTER,
		WRITE_CHANNEL_IS_OPEN,
		READ_CHANNEL_IS_OPEN,
		CACHE_ERROR_DIE_0_PLANE_0,
		CACHE_ERROR_DIE_0_PLANE_1,
		CACHE_ERROR_DIE_1_PLANE_0,
		CACHE_ERROR_DIE_1_PLANE_1,
		FILEWRITEOPEN_FAIL,
		FILEREADOPEN_FAIL,
		WRITEBYTES_TOO_SMALL,
		NOT_ENOUGH_FREE_BLOCKS,
		FILEWRITE_WRONG_BYTES_WRITTEN,
		FILEREAD_WRONG_BYTES_READ,
		INCORRECT_READ_BACK_CONTENT
	};

private:
	SPIClass *SPIPTR = 0;
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
	uint16_t _NAND_BYTES_PER_PAGE;
	uint16_t _NAND_PAGES_PER_BLOCK;
	size_t _NAND_BLOCKS_PER_DIE;
	size_t _NAND_DIES_PER_DEVICE;
	size_t _NAND_BLOCKS_PER_DEVICE;
	uint32_t _NAND_PAGES_PER_DIE;
	uint32_t _NAND_PAGES_PER_DEVICE;
	uint64_t _HighestGlobalPageSerial;
	uint32_t _HighestGlobalPageIndex;
	uint32_t _HighestGlobalFileSerial;
	uint8_t _SelectedDie;
	uint8_t *DeviceBlockMap;
	bool _HasPinControl;
	bool _HighestGlobalPageScanned;
	bool _DeviceBlockMapBuilt;

	bool _FileWriteChannelOpen;
	uint32_t _FileWriteChannelFirstPageGPI;
	char *_FileWriteChannelFileDescription;
	//uint64_t _FileWriteChannelFileDate;
	uint32_t _FileWriteChannelTotalBytesWritten;
	uint32_t _FileWriteChannelFilePageIndex;
	//uint32_t _FileWriteChannelFileType;

	bool _FileReadChannelOpen;
	uint32_t _FileReadChannelFirstPageGPI;
	uint16_t _FileReadChannelCurrentPageStartPos;
	uint32_t _FileReadChannelTotalBytesRead;
	uint32_t _FileReadChannelFileSerial;
	uint32_t _FileReadChannelFilePageIndex;

	void SelectDevice();
	void UnselectDevice();
	void EnableWrite();
	void Reverse_Endian(void *A, void *B, size_t Length);
	bool WaitForReady();

	void FastWriteToCache(void *Buf, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit, bool RetainExistingDataInCache);
	void FastReadFromCache(void *Buf, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit, bool RetainExistingDataInCache);
	void FastCommitToNAND(uint32_t RowAddress);
	void FastRetrieveFromNAND(uint32_t RowAddress);

	uint8_t SelectDie(size_t die);
	void EnableECC();
	void DisableAllBlockLocks();
	void EraseData(uint32_t GPI);
	void EraseDie(size_t die);
	void EraseBlock(uint32_t Address);

	uint32_t InitializeHighestGPI();
	uint32_t InitializeHighestGPI(uint32_t GPI);


	void WriteMetadata(uint64_t PageSerial, uint32_t FileSerial, uint32_t FilePageIndex, uint16_t BytesUsedInThisPage, uint8_t FileType, char *FileDescription, bool PlaneBit);
	void ReadMetadata(uint64_t &PageSerial, uint32_t &FileSerial, uint32_t &FilePageIndex, uint16_t &BytesUsedInThisPage, uint8_t &FileType, char *FileDescription, bool PlaneBit);
	void RetrieveMetadata(uint32_t GPI, uint64_t &PageSerial, uint32_t &FileSerial, uint32_t &FilePageIndex, uint16_t &BytesUsedInThisPage, uint8_t &FileType, char *FileDescription);
	void RetrieveMetadata(uint32_t GPI, uint32_t &FileSerial, uint32_t &FilePageIndex, uint16_t &BytesUsedInThisPage, bool &PlaneBit);
	void RetrieveMetadata(uint32_t GPI, uint64_t &PageSerial, uint32_t &FileSerial);


	void BuildDeviceBlockMap();
	void ScanForHighestPage();

	void MarkBlockUnused(uint16_t GBI);
	void MarkBlockUsed(uint16_t GBI);
	bool GetPlaneBit(uint32_t RA);

	void GlobalPageIndex_to_RA_and_Die(uint32_t GlobalPageIndex, uint32_t &RowAddress, uint8_t &DieIndex);
	void GlobalBlockIndex_to_RA_and_Die(uint16_t GlobalBlockIndex, uint32_t &RowAddress, uint8_t &DieIndex);
	void GlobalBlockIndex_to_GlobalPageIndex(uint16_t GlobalBlockIndex, uint32_t &GlobalPageIndex);
	void GlobalPageIndex_to_GlobalBlockIndex(uint32_t GlobalPageIndex, uint16_t &GlobalBlockIndex);
	void GlobalPageIndex_to_PlaneBit(uint32_t GlobalPageIndex, bool &PlaneBit);
	void GlobalPageIndex_to_RA_and_selDie_and_PlaneBit(uint32_t GlobalPageIndex, uint32_t &RowAddress, uint8_t &DieIndex, bool &PlaneBit);

	//void DefaultWriteToCache(uint8_t *Buffer, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit);
	//void DefaultReadFromCache(uint8_t *Buffer, size_t DataLength, uint16_t ColumnAddress, bool PlaneBit);
	//void DefaultCommitToNAND(uint32_t RowAddress);
	//void DefaultRetrieveFromNAND(uint32_t RowAddress);



};


//bool strequal(const char *A, const char *B, size_t L);
