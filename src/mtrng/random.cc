#include <cassert>
#include <cstdlib>

// These are factorized out, because for C++ TR1 the header names might be
// <tr1/random> and the namespace is std::tr1.
#ifndef HAVE_CXX11RANDOM
    #define HAVE_CXX11RANDOM (__cplusplus >= 201103L)
    #define CXX11RANDOM_HEADER <random>
    #define CXX11RANDOM_NAMESPACE std
#endif

#if HAVE_CXX11RANDOM

#include CXX11RANDOM_HEADER
using CXX11RANDOM_NAMESPACE::mt19937;
using CXX11RANDOM_NAMESPACE::uniform_real_distribution;

/**
 * Wrapper around C++11 standard Mersenne twister random source.
 *
 * This class wraps around the C++11 standard random source std::mt19973, as
 * defined in the `random` header.  A wrapper is used as we (1) want to
 * substitute the source with any other random source in the case that the
 * compiler does not support C++11 and (2) provide a stripped-down interface
 * to random numbers in the interval [0,1), which can then be
 * easily implemented.
 */
class Cxx11Random
{
public:
    static unsigned entropy() { return mt19937::word_size; }

    Cxx11Random(unsigned seed=0) : engine_(seed), distrib_(0., 1.) {}

    void seed(unsigned seed) { engine_.seed(seed); }

    double operator() () { return distrib_(engine_); }

    unsigned raw() { return engine_(); }

private:
    mt19937 engine_;
    uniform_real_distribution<double> distrib_;
};

#endif /* HAVE_CXX11RANDOM */

class CStdlibRandom;

/** Return log2(x+1), thus avoiding overflow for x == -1 */
static unsigned log2p1(unsigned x) {
    ++x;
    if (x == 0)
        return log2p1(-2) + 1;
    unsigned res = 0;
    while (x > 1) {
        x >>= 1;
        ++res;
    }
    return res;
}

class CStdlibRandom;
static CStdlibRandom *current = NULL;

/**
 * Wrapper around the C standard random number generator `stdlib.h`
 *
 * This class wraps around the C standard PRNG functions `rand()` and
 * `srand()`.  This PRNG is almost always of the shift-modulo type, which
 * produces random numbers of poor quality, low entropy and short period.
 * However, since it is the only PRNG for old C++ versions, we provide a
 * wrapper here.
 *
 * This class is a singleton (it may only be instantiated *once*), since the
 * underlying C library function store the current state in a global.
 */
class CStdlibRandom
{
public:
    static unsigned entropy() { return log2p1(RAND_MAX); }

    CStdlibRandom(unsigned seed=0) {
        assert(current == NULL && "class is a singleton");
        current = this;
        this->seed(seed);
    }

    ~CStdlibRandom() {
        if (current == this)
            current = NULL;
    }

    void seed(unsigned seed) { srand(seed); }

    double operator() () {
        assert(current != NULL && "class is a singleton");
        return rand()/(RAND_MAX + 1.0);
    }
};


// Make sure we always have a RNG
#if HAVE_CXX11RANDOM
    typedef Cxx11Random DefaultRandom;
#else
    #pragma message "WARNING: Falling back on poor Stdlib RNG. Check results."
    typedef CStdlibRandom DefaultRandom;
#endif


extern "C" {

void *new_random() {
    return new DefaultRandom();
}

void delete_random(void *rng) {
    delete (DefaultRandom *)rng;
}

void seed_random(void *rng, unsigned seed) {
    ((DefaultRandom *)rng)->seed(seed);
}

double get_random(void *rng) {
    return ((DefaultRandom *)rng)->operator()();
}

}
