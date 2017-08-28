#include "Random.h"
#include "macro.h"
#include <time.h>
int Random::bits_num = 0;
uint Random::bits_data = 0;
#if !HAS_RAND48
static const int IM1 = 2147483563;
static const int IM2 = 2147483399;
static const int IMM1 = (IM1 - 1);
static const int IA1 = 40014;
static const int IA2 = 40692;
static const int IQ1 = 53668;
static const int IQ2 = 52774;
static const int IR1 = 12211;
static const int IR2 = 3791;
static const int NTAB = 32;
static const int NDIV = (1 + IMM1 / NTAB);
static int prev;
static int table[NTAB];
static int state1;
static int state2;
#endif // HAS_RAND48

#if HAS_RAND48
void Random::reset(int seed) {

        if(seed == -1)
                seed = time_seed();

        static uint2 seeds[3];

        seeds[0] = (seed >> 16);
        seeds[1] = (seed & 0xFFFF);
        seeds[2] = seeds[0] ^ seeds[1];

        seed48(seeds);
}
#else // HAS_RAND48
void Random::reset(int seed) {

        if(seed == -1)
                seed = time_seed();

        state2 = state1 = abs(seed) + !seed;

        for(int i = int(NTAB + 7); i >= 0; i--) {
                int k = state1 / IQ1;
                state1 = IA1 * (state1 - k * IQ1) - k * IR1;
                if(state1 < 0)
                        state1 += IM1;
                if(i < NTAB)
                        table[i] = state1;
        }

        prev = table[0];
}
#if INT_MAX != 2147483647
# error The following code assumes 32 bit integers
#endif // INT_MAX

uint Random::bits() {
        DBG_ASSERT(state1 || state2, "You forgot to reset the Random seed");
        int k = state1 / IQ1;
        state1 = IA1 * (state1 - k * IQ1) - k * IR1;
        if(state1 < 0)
                state1 += IM1;

        k = state2 / IQ2;
        state2 = IA2 * (state2 - k * IQ2) - k * IR2;
        if(state2 < 0)
                state2 += IM2;

        int i = int(prev / NDIV);
        prev = table[i] - state2;
        table[i] = state1;
        if(prev < 1)
                prev += IMM1;

        return(uint(prev ^ (prev << 16)));
}
#endif // SCM_HAS_RAND48
int Random::time_seed() {

        time_t now;
        time(&now);
        return(int(now));
}
uint Random::bits(int num) {

        ASSERT(0 <= num && num <= 32,
                "asked for invalid number of bits (" << num << ')');

        if(bits_num < num) {
                bits_num = 32;
                bits_data = bits();
        }

        int value = (num == 32 ? bits_data : bits_data & (uint(1 << num) - 1));
        bits_data >>= num;
        bits_num -= num;

        return(value);
}
bool Random::boolean() {
        if(!bits_num) {
                bits_num = 32;
                bits_data = bits();
        }

        bool value = bits_data & 1;
        bits_data >>= 1;
        bits_num--;

        return(value);
}
