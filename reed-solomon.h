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

Polynomial RS_correct_errata(const Polynomial &msg_in, const Polynomial &synd, const Polynomial &err_pos);

Polynomial RS_find_error_locator(const Polynomial &synd, const unsigned int nsym, const unsigned int erase_count);
Polynomial RS_find_error_locator(const Polynomial &synd, const unsigned int nsym, const unsigned int erase_count, const Polynomial &erase_loc);

Polynomial RS_find_errors(const Polynomial &err_loc, const unsigned int nmess);
Polynomial RS_forney_syndromes(const Polynomial &synd, const Polynomial &pos, const unsigned int nmess);

Polynomial RS_correct_message(const Polynomial &msg_in, const unsigned int nsym, const Polynomial &erase_pos);

#endif
