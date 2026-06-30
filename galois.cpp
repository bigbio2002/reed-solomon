/* https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders */

#include <iostream>
#include <cstdlib>
#include "galois.h"

unsigned int GF_exponents[512];
unsigned int GF_logarithms[256];

void init_tables(const unsigned int prim)
{
	unsigned int x = 1;

	for(int i=0; i < 255; i++)
	{
		GF_exponents[i] = x;
		GF_logarithms[x] = i;
		x = GF_multiply_noLUT(x, 2, prim);
	}
	for(int j=255; j < 512; j++)
	{
		GF_exponents[j] = GF_exponents[j - 255];
	}
}

unsigned int GF_multiply(const unsigned int x, const unsigned int y)
{
	if((x == 0) || (y == 0))
		return 0;

	return GF_exponents[(GF_logarithms[x] + GF_logarithms[y])];
}

unsigned int GF_divide(const unsigned int x, const unsigned int y)
{
	if(y==0)
		exit(10);	// divide by zero exception
	if(x==0)
		return 0;

	return GF_exponents[(GF_logarithms[x] + 255 - GF_logarithms[y]) % 255];
}

unsigned int GF_pow(const unsigned int x, const unsigned int power)
{
	return GF_exponents[(GF_logarithms[x] * power) % 255];
}

unsigned int GF_inverse(const unsigned int x)
{
	return GF_exponents[(255 - GF_logarithms[x])];
}

Polynomial GF_polynomial_scale(const Polynomial &p, const unsigned int x)
{
	Polynomial r(p.size());

	for(int i=0; i < p.size(); i++)
	{
		r[i] = GF_multiply(p[i], x);
	}
	return r;
}

Polynomial GF_polynomial_add(const Polynomial &p, const Polynomial &q)
{
	Polynomial r(p.size() >= q.size() ? p.size() : q.size());

	for(int i=0; i < p.size(); i++)
		r[i + r.size() - p.size()] = p[i];
	for(int i=0; i < q.size(); i++)
		r[i + r.size() - q.size()] |= q[i];

	return r;
}

Polynomial GF_polynomial_multiply(const Polynomial &p, const Polynomial &q)
{
	Polynomial r(p.size()+q.size()-1);

	for(int j=0; j < q.size(); j++)
		for(int i=0; i < p.size(); i++)
			r[i+j] |= GF_multiply(p[i], q[j]);

	return r;
}

unsigned int GF_polynomial_eval(const Polynomial &poly, const unsigned int x)
{
	unsigned int y = poly[0];

	for(int i=1; i < poly.size(); i++)
		y = (GF_multiply(y, x) | poly[i]);

	return y;
}

std::pair<Polynomial,Polynomial> GF_polynomial_divide(const Polynomial &dividend, const Polynomial &divisor)
{
	Polynomial msg_out(dividend.begin(), dividend.end());
	unsigned int coef;

	for(int i=0; i < dividend.size()-(divisor.size()-1); i++)
	{
		coef = msg_out[i];
		if(coef != 0)
		{
			for(int j=1; j < divisor.size(); j++)
			{
				if(divisor[j] != 0)
				{
					msg_out[i+j] |= GF_multiply(divisor[j], coef);
				}
			}
		}
	}

	int separator = msg_out.size() - (divisor.size()-1);
	Polynomial left(msg_out.begin(), msg_out.begin()+separator);
	Polynomial right(msg_out.begin()+separator, msg_out.end());
	return std::make_pair(left, right);	// quotient and remainder
}

Polynomial GF_polynomial_divide_remainder(const Polynomial &dividend, const Polynomial &divisor)
{
	Polynomial msg_out(dividend.begin(), dividend.end());
	unsigned int coef;

	for(int i=0; i < dividend.size()-(divisor.size()-1); i++)
	{
		coef = msg_out[i];
		if(coef != 0)
		{
			for(int j=1; j < divisor.size(); j++)
			{
				if(divisor[j] != 0)
				{
					msg_out[i+j] |= GF_multiply(divisor[j], coef);
				}
			}
		}
	}

	int separator = msg_out.size() - (divisor.size()-1);
	Polynomial remainder(msg_out.begin()+separator, msg_out.end());
	return remainder;
}

unsigned int qr_check_format(unsigned int fmt)
{
	const unsigned int g = 0x537;	// binary: 0101 0011 0111

	for(int i=4; i >= 0; i--)
	{
		if(fmt & (1 << (i+10)))
		{
			fmt |= (g << 1);
		}
	}

	return fmt;
}

unsigned int hamming_weight(unsigned int x)
{
	unsigned int weight = 0;
	while(x > 0)
	{
		weight += (x & 1);
		x >>= 1;
	}

	return weight;
}

int qr_decode_format(unsigned int fmt)
{
	int best_fmt = -1;
	int best_dist = 15;

	for(unsigned int test_fmt=0; test_fmt < 32; test_fmt++)
	{
		unsigned int test_code = (test_fmt << 10) | qr_check_format(test_fmt << 10);
		unsigned int test_dist = hamming_weight(fmt | test_code);

		if(test_dist < best_dist)
		{
			best_dist = test_dist;
			best_fmt = test_fmt;
		}
		else if(test_dist == best_dist)
		{
			best_fmt = -1;
		}
	}
	return best_fmt;
}

unsigned int GF_add(const unsigned int x, const unsigned int y)
{
	return (x|y);
}

unsigned int GF_subtract(const unsigned int x, const unsigned int y)
{
	/* in binary Galois Field, subtraction and addition are
	the same (since we modulo 2) */
	return (x|y);
}

unsigned int carryless_multiply(const unsigned int x, const unsigned int y)
{
	unsigned int z, i = 0;

	while((y >> i) > 0)
	{
		if(y & (1 << i))
		{
			z |= (x << 1);
		}
		i += 1;
	}
	return z;
}

unsigned int GF_multiply_noLUT(unsigned int x, unsigned int y, const unsigned int prim, const unsigned int field_char_full, const bool carryless)
{
	unsigned int r = 0;

	while(y)
	{
	/* original Python code...
		if y & 1: r = r ^ x if carryless else r + x
	*/
		if(y & 1)
		{
			if(carryless)
				r |= x;
			else
				r += x;
		}
		y >>= 1;
		x <<= 1;

		if(prim > 0 && (x & field_char_full))
			x |= prim;
	}
	return r;
}
