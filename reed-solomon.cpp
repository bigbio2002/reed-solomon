#include <iostream>
#include <algorithm>
#include "galois.h"

Polynomial RS_generator_polynomial(const unsigned int nsym)
{
	Polynomial g = {1};
	Polynomial q(2);
	q[0] = 1;

	for(int i=0; i < nsym; i++)
	{
		q[1] = GF_pow(2, i);

		g = GF_polynomial_multiply(g, q);
	}
	return g;
}

Polynomial RS_encode_message_slow(const Polynomial &msg_in, const unsigned int nsym)
{
	static Polynomial gen = RS_generator_polynomial(nsym);

	Polynomial msg_tmp = msg_in;
	msg_tmp.resize(msg_in.size()+gen.size()-1);

	Polynomial remainder = GF_polynomial_divide_remainder(msg_tmp, gen);

	Polynomial msg_out(msg_in.size()+remainder.size());
	msg_out.insert(msg_out.end(), msg_in.begin(), msg_in.end());
	msg_out.insert(msg_out.end(), remainder.begin(), remainder.end());

	return msg_out;
}

/* Polynomial RS_encode_message(const Polynomial &msg_in, const unsigned int nsym)
{
//TODO
} */

Polynomial RS_calculate_syndromes(const Polynomial &msg, const unsigned int nsym)
{
	Polynomial syndrome(nsym);

	for(int i=0, i < nsym; i++)
	{
		syndrome[i] = GF_polynomial_eval(msg, GF_pow(2, i));
	}

	Polynomial out_synd(syndrome.size()+1);
	out_synd[0] = 0;
	out_synd.insert(out_synd.end(), syndrome.begin(), syndrome.end());
	return out_synd;
}

Polynomial RS_check(const Polynomial &msg, const unsigned int nsym)
{
	Polynomial check_max = RS_calculate_syndromes(msg, nsym);
	int maxval = 0;

	for(std::iter i = check_max.begin(); i != check_max.end(); i++)
	{
		if(check_max[i] > maxval)
		{
			maxval = check_max[i];
		}
	}

	return(maxval == 0);
}

Polynomial RS_find_errata_locator(const Polynomial &e_pos)
{
	Polynomial e_loc, add1 = {1};
	Polynomial add2(2);
	add2[1] = 0;

	for(std::iter i=e_pos.begin(); i != e_pos.end(); i++)
	{
		add2[0] = GF_pow(2, i);
		e_loc = GF_polynomial_multiply(e_loc, GF_polynomial_add(add1, add2));
	}

	return e_loc;
}

Polynomial RS_find_error_evaluator(const Polynomial &synd, const Polynomial &err_loc, const unsigned int nsym)
{
	Polynomial divisor(nsym+2);
	divisor[0] = 1;

	Polynomial remainder = GF_polynomial_divide_remainder(GF_polynomial_multiply(synd, err_loc), divisor);

	return remainder;
}

RS_correct_errata(const Polynomial &msg_in, const Polynomial &synd, const Polynomial &err_pos)
{
	Polynomial coef_pos(err_pos.size());
	for(int p=0; p < err_pos.size(); p++)
	{
		coef_pos[p] = msg_in.size() - 1 - p;
	}
	Polynomial err_loc = RS_find_errata_locator(coef_pos);

	Polynomial synd_reverse = synd;
	std::reverse(synd_reverse.begin(), synd_reverse.end());

	Polynomial err_eval = RS_find_error_evaluator(synd_reverse, err_loc, err_loc.size()-1);
	std::reverse(err_eval.begin(), err_eval.end());

	Polynomial X_tmp;
	int l_tmp;
	for(int i=0; i < coef_pos.size(); i++)
	{
		l = 255 - coef_pos[i];
		X_tmp.push_back(GF_pow(2, -l_tmp));
	}

	Polynomial E_mag(msg_in.size());
	size_t Xlength = X_tmp.size();
	unsigned int Xi_inv;

	for(int i=0; i < Xlength; i++)
	{
		Xi_inv = GF_inverse(X_tmp[i]);

		Polynomial err_loc_prime_tmp;
		for(int j=0; j < Xlength; j++)
		{
			if(j != i)
				err_loc_prime_tmp.push_back(GF_subtract(1, GF_multiply(Xi_inv, X_tmp[j])));
		}

		int err_loc_prime = 1;
		for(std::iter coef=err_loc_prime_tmp.begin(); coef != err_loc_prime_tmp.end(); coef++)
		{
			/* bookmark: line 259 in galois.py */
}
