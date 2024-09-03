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
using CXX11RANDOM_NAMESPACE::uniform_int_distribution;

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

    double rand () { return distrib_(engine_); }

    int randint (int min, int max) { uniform_int_distribution<int> dist(min, max); return dist(engine_); }

    unsigned raw() { return engine_(); }

private:
    mt19937 engine_;
    uniform_real_distribution<double> distrib_;
};

#endif /* HAVE_CXX11RANDOM */

#if HAVE_CXX11RANDOM
    typedef Cxx11Random DefaultRandom;
#else
    #error "Support for the RNGs provided by the C++ standard library in C++11 and higher is required"
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
    return ((DefaultRandom *)rng)->rand();
}

int get_randint(void *rng, int min, int max) {
    return ((DefaultRandom *)rng)->randint(min, max);
}

}
