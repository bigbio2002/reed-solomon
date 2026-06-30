#ifndef _GALOIS_H_
#define _GALOIS_H_

/* https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders */

#include <vector>

extern unsigned int GF_exponents[512];
extern unsigned int GF_logarithms[256];

typedef std::vector<unsigned int> Polynomial;

void init_tables(const unsigned int prim=0x11d);
unsigned int GF_multiply(const unsigned int x, const unsigned int y);
unsigned int GF_divide(const unsigned int x, const unsigned int y);
unsigned int GF_pow(const unsigned int x, const unsigned int power);
unsigned int GF_inverse(const unsigned int x);

Polynomial GF_polynomial_scale(const Polynomial &p, const unsigned int x);
Polynomial GF_polynomial_add(const Polynomial &p, const Polynomial &q);
Polynomial GF_polynomial_multiply(const Polynomial &p, const Polynomial &q);
unsigned int GF_polynomial_eval(const Polynomial &poly, const unsigned int x);
std::pair<Polynomial,Polynomial> GF_polynomial_divide(const Polynomial &dividend, const Polynomial &divisor);
Polynomial GF_polynomial_divide_remainder(const Polynomial &dividend, const Polynomial &divisor);

unsigned int qr_check_format(unsigned int fmt);
unsigned int hamming_weight(unsigned int x);
int qr_decode_format(unsigned int fmt);

unsigned int GF_add(const unsigned int x, const unsigned int y);
unsigned int GF_subtract(const unsigned int x, const unsigned int y);
unsigned int carryless_multiply(const unsigned int x, const unsigned int y);
unsigned int GF_multiply_noLUT(unsigned int x, unsigned int y, const unsigned int prim=0, const unsigned int field_char_full=256, const bool carryless=true);

#endif
