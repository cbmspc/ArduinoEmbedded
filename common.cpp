#include "common.h"

// memcpy and memcmp give a lot of fake errors in Visual Studio IDE. Use these to avoid the fake errors.
void mcpy(void *dest, const void *src, size_t len) {
	unsigned char *d = (unsigned char*)dest;
	const unsigned char *s = (const unsigned char*)src;
	while (len--) {
		*d++ = *s++;
	}
}

// memcpy and memcmp give a lot of fake errors in Visual Studio IDE. Use these to avoid the fake errors.
int mcmp(void *str1, void *str2, size_t count) {
	register const unsigned char *s1 = (const unsigned char*)str1;
	register const unsigned char *s2 = (const unsigned char*)str2;
	while (count-- > 0)	{
		if (*s1++ != *s2++) {
			return s1[-1] < s2[-1] ? -1 : 1;
		}
	}
	return 0;
}

