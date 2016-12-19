#include <stdint.h>

typedef struct {
  uint32_t high;
  uint8_t low;
} uint40_t;

static inline uint64_t read_uint40(uint40_t i) {
  return ((uint64_t) i.high <<8 | (uint64_t) low);
}

static inline uint40_t write_uint40(uint64_t i) {
  uint40_t j = {
    .high = (uint32_t) (i>>8),
    .low = (uint8_t) (i&0xff),
  };
  return j;
}

static inline uint64_t *get_primary_CT(char *seq) {
  uint132_t k, k4;
  for (i=0; i<15; ++i) {
    k4 = (k4<<2) | char_to_nt4_table[seq[i]];
    if (i%5==4) {
      k = ((k<<8) | fivebase43CT_switch[k4]);
      k4 = 0;
    }
  }
}
