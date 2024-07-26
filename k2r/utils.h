#include <stdint.h>
#include <math.h>

int32_t makehist_32(int32_t *S, int32_t *hist, uint32_t len);
double  entropy(int32_t *hist, uint32_t histlen, uint32_t len);
int32_t makehist_16(int16_t *S, int32_t *hist, uint32_t len);
int32_t makehist_u16(uint16_t *S, int32_t *hist, uint32_t len);