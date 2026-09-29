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
};

#if HAS_RAND48
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
