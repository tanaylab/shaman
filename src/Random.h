/*
 * Random.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <stdlib.h>
#define HAS_RAND48 1

class Random {
private:
        static int bits_num;
        static uint bits_data;
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
inline float Random::fraction() {
        return(float(drand48()));
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
