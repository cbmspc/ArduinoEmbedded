#pragma once

#if defined(ARDUINO) && ARDUINO >= 100
	#include "arduino.h"
#else
	#include "WProgram.h"
#endif


void mcpy(void *dest, const void *src, size_t len);
int mcmp(void *str1, void *str2, size_t count);
