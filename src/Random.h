/*
 * Random.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <stdlib.h>
#include <stdint.h>
#define HAS_RAND48 1

class Random {
private:
        static int bits_num;
        static uint bits_data;
        static uint64_t state48;	// drand48 state, advanced inline (see fraction())
public:
        static int time_seed();

        static void reset(int seed = -1);
        static uint bits();
        static uint bits(int num);

        static bool boolean();

        static float fraction();			// returns [0,1]
        static float fraction_truncated();  // returns [0,1)
        static float peek_fraction(int k);	// what fraction() returns k calls from now; no state change
};

#if HAS_RAND48
// mrand48() uses the libc rand48 state, which fraction() does not advance (the shuffler does not use bits())
inline uint Random::bits() {
        uint raw = mrand48();
        return(raw ^ (raw >> 16));
}
// Same stream as glibc drand48() after seed48(): X = (0x5DEECE66D * X + 0xB) mod 2^48,
// returned as X / 2^48 (exact in a double).
inline float Random::fraction() {
        state48 = (state48 * 0x5DEECE66DULL + 0xBULL) & 0xFFFFFFFFFFFFULL;
        return(float(double(state48) * 0x1p-48));
}
// k steps of the LCG at once: X_k = A_k * X + C_k (mod 2^48)
struct Lcg48Jump { uint64_t a; uint64_t c; };
constexpr Lcg48Jump lcg48_jump(int k) {
        Lcg48Jump j = {1, 0};
        for (int i = 0; i < k; i++) {
                j.a = (j.a * 0x5DEECE66DULL) & 0xFFFFFFFFFFFFULL;
                j.c = (j.c * 0x5DEECE66DULL + 0xBULL) & 0xFFFFFFFFFFFFULL;
        }
        return(j);
}
constexpr Lcg48Jump lcg48_jumps[8] = {lcg48_jump(0), lcg48_jump(1), lcg48_jump(2), lcg48_jump(3),
        lcg48_jump(4), lcg48_jump(5), lcg48_jump(6), lcg48_jump(7)};
inline float Random::peek_fraction(int k) {	// 0 < k < 8
        uint64_t state = (lcg48_jumps[k].a * state48 + lcg48_jumps[k].c) & 0xFFFFFFFFFFFFULL;
        return(float(double(state) * 0x1p-48));
}
inline float Random::fraction_truncated() {
	float ret = Random::fraction();
	while (ret == 1) {
		ret = Random::fraction();
	}
	return ret;
}

#else // HAS_RAND48
inline float Random::fraction() {
        return(float(bits() / double(UINT_MAX)));
}
#endif // HAS_RAND48


#endif /* RANDOM_H_ */
