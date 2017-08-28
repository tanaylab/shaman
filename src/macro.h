/*
 * macro.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef MACRO_H_
#define MACRO_H_

#include <iostream>


#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif


typedef unsigned int uint;
typedef unsigned short uint2;

#endif /* MACRO_H_ */
