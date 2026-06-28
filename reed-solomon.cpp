#include <iostream>
#include <iterator>
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

	for(int i=0; i < nsym; i++)
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

	for(std::iter it = check_max.begin(); it != check_max.end(); it++)
	{
		if(check_max[it] > maxval)
		{
			maxval = check_max[it];
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

Polynomial RS_correct_errata(const Polynomial &msg_in, const Polynomial &synd, const Polynomial &err_pos)
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
		l_tmp = 255 - coef_pos[i];
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

		unsigned int err_loc_prime = 1;
		for(int coef=0; coef < err_loc_prime_tmp.size(); i++)
			err_loc_prime = GF_multiply(err_loc_prime, err_loc_prime_tmp[coef]);

		std::reverse(err_eval.begin(), err_eval.end());

		unsigned int y = GF_polynomial_eval(err_eval, Xi_inv);
		y = GF_multiply(GF_pow(X_tmp[i], 1), y);

		if(err_loc_prime == 0)
			exit(1);	//could not find error magnitude

		unsigned int magnitude = GF_divide(y, error_loc_prime);
		E_mag[err_pos[i]] = magnitude;
	}

	Polynomial msg_out = GF_polynomial_add(msg_in, E_mag);

	return msg_out;
}

Polynomial RS_find_error_locator(const Polynomial &synd, const unsigned int nsym, const unsigned int erase_count=0)
{
	Polynomial err_loc, old_loc = { 1 };

	unsigned int synd_shift = synd.size()-nsym;
	unsigned int K;

	for(int i=0; i < nsym - erase_count; i++)
	{
		K = i+synd_shift;

		unsigned int delta = synd[K];

		for(int j=1; j < err_loc.size(); j++)
			delta |= GF_multiply(err_loc[-(j+1)], synd[K-j]);

		old_loc.push_back(0);

		if(delta != 0)
		{
			if(old_loc.size() > err_loc.size())
			{
				Polynomial new_loc = GF_polynomial_scale(old_loc, delta);
				old_loc = GF_polynomial_scale(err_loc, GF_inverse(delta));
				err_loc = new_loc;
			}

			err_loc = GF_polynomial_add(err_loc, GF_polynomial_scale(old_loc, delta));
		}
	}

	while(err_loc.size() && err_loc[0] == 0)
		err_loc[0].erase();

	unsigned int errs = err_loc.size()-1;
	if((errs - erase_count)*2 + erase_count > nsym)
		exit(99);

	return err_loc;
}

Polynomial RS_find_error_locator(const Polynomial &synd, const unsigned int nsym, const unsigned int erase_count=0, const Polynomial &erase_loc)
{
	Polynomial err_loc, old_loc = erase_loc;

	unsigned int synd_shift = synd.size()-nsym;
	unsigned int K;

	for(int i=0; i < nsym - erase_count; i++)
	{
		K = erase_count+i+synd_shift;

		unsigned int delta = synd[K];

		for(int j=1; j < err_loc.size(); j++)
			delta |= GF_multiply(err_loc[-(j+1)], synd[K-j]);

		old_loc.push_back(0);

		if(delta != 0)
		{
			if(old_loc.size() > err_loc.size())
			{
				Polynomial new_loc = GF_polynomial_scale(old_loc, delta);
				old_loc = GF_polynomial_scale(err_loc, GF_inverse(delta));
				err_loc = new_loc;
			}

			err_loc = GF_polynomial_add(err_loc, GF_polynomial_scale(old_loc, delta));
		}
	}

	while(err_loc.size() && err_loc[0] == 0)
		err_loc[0].erase();

	unsigned int errs = err_loc.size()-1;
	if((errs - erase_count)*2 + erase_count > nsym)
		exit(99);

	return err_loc;
}

Polynomial RS_find_errors(const Polynomial &err_loc, const unsigned int nmess)
{
	unsigned int errs = err_loc.size()-1;
	Polynomial err_pos;

	for(int i=0; i < nmess; i++)
	{
		if(GF_polynomial_eval(err_loc, GF_pow(2, i)) == 0)
			err_pos.push_back(nmess-1-i);
	}

	if(err_pos.size() != errs)
		exit(10);

	return err_pos;
}

Polynomial RS_forney_syndromes(const Polynomial &synd, const Polynomial &pos, const unsigned int nmess)
{
	Polynomial erase_pos_reversed;
	for(int p=0; p < pos.size(); p++)
		erase_pos_reversed.push_back(pos[nmess-1-p]);

	Polynomial fsynd(synd.size()-1);
	for(int l=1; l < synd.size(); l++)
		fsynd[l-1] = synd[l];

	unsigned int x;
	for(int i=0; i < pos.size(); i++)
	{
		x = GF_pow(2, erase_pos_reversed[i]);
		for(int j=0; j < fsynd.size()-1; j++)
		{
			fsynd[j] = GF_multiply(fsynd[j], x) | fsynd[j+1];
		}
	}

	return fsynd;
}

Polynomial RS_correct_message(const Polynomial &msg_in, const unsigned int nsym, const Polynomial &erase_pos)
{
	if(msg_in.size() > 255)
		exit(56);

	Polynomial msg_out = msg_in;

	for(int e_pos=0; e_pos < erase_pos.size(); e_pos++)
		msg_out[e_pos] = 0;

	if(erase_pos.size() > nsym)
		exit(3);

	Polynomial synd = RS_calculate_syndromes(msg_out, nsym);
	for(int i=0; i < synd.size(); i++)
	{
		if(synd[i])	/* scan through synd. if any of the coefficients arent zero, we can't return the message immediately, so abort the scan */
			break;

		/* are we on the last iteration of the for loop? we've scanned through synd, so if we haven't broken out yet, that means all coefficients are zero and we can return early */
		if(i >= synd.size()-1)
			return /* TODO */;
	}

	Polynomial fsynd = RS_forney_syndromes(synd, erase_pos, msg_out.size());
	Polynomial err_loc = RS_find_error_locator(fsynd, nsym, erase_pos.size());
	std::reverse(err_loc.begin(), err_loc.end());
	Polynomial err_pos = RS_find_errors(err_loc, msg_out.size());
	if(err_pos.empty())
		exit(34);

	Polynomial positions = erase_pos;
	positions.insert(positions.end(), err_pos.begin(), err_pos.end());
	msg_out = RS_correct_errata(msg_out, synd, positions);

	synd = RS_calculate_syndromes(msg_out, nsym);
	for(int i=0; i < synd.size(); i++)
	{
		if(synd[i] > 0)
			exit(9);
	}

	return /* TODO */;
}
