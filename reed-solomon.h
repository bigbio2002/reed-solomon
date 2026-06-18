#ifndef _REED_SOLOMON_H_
#define _REED_SOLOMON_H_

#include "galois.h"

Polynomial RS_generator_polynomial(const unsigned int nsym);
Polynomial RS_encode_message_slow(const Polynomial &msg_in, const unsigned int nsym);
/* Polynomial RS_encode_message(const Polynomial &msg_in, const unsigned int nsym); */
Polynomial RS_calculate_syndromes(const Polynomial &msg, const unsigned int nsym);
Polynomial RS_check(const Polynomial &msg, const unsigned int nsym);
Polynomial RS_find_errata_locator(const Polynomial &e_pos);
Polynomial RS_find_error_evaluator(const Polynomial &synd, const Polynomial &err_loc, const unsigned int nsym);

#endif
