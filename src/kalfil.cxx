// ----------------------------------------------------------------------------
// kalfil.c -- bayesian morse code decoder 
//
// Copyright (C) 2012-2014
//		     (C) Mauri Niininen, AG1LE
//
// This file is part of Bayesian Morse code decoder   

// bmorse is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bmorse is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bmorse.  If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------------

#include "bmorse.h"
#include <math.h> 
#include <stdio.h>


int morse::kalfil_(const real z, const integer ip, const real rn,
	const integer ixs, const integer kelem, const integer jnode, const real
	dur, const integer ilrate, real *pin, real *lkhdj)
{
/*   THIS SUBROUTINE COMPUTES THE ARRAY OF KALMAN FILTER */
/*   RECURSIONS USED TO DETERMINE THE LIKELIHOODS. */

/*   VARIABLES: */
/*       Z -	INPUT MEASUREMENT */
/*       IP -	INPUT PATH IDENTITY */
/*       RN -	INPUT NOISE POWER ESTIMATE */
/*       ILX -	INPUT SAVED KEYSTATE ON PATH IP */
/*       IXS -	INPUT KEYSTAT OF NEW NODE */
/*       KELEM -	INPUT ELEM STATE OF NEW NODE */
/*       ISRATE 	INPUT SPEED STATE OF NEW NODE */
/*       DUR - 	INPUT CURRENT DURATION OF ELEMENT ON IP */
/*       ILRATE 	INPUT SPEED STATE ON PATH IP */
/*       PIN - 	TRANS PROB FROM PATH IP TO NODE N */
/*       LKHDJ - OUTPUT CALCULATED LIKELIHOOD VALUE */

/*   SUBROUTINES USED */
/*       MODEL - OBTAINS THE SIGNAL-STATE-DEPENDENT LINEAR */
/*       	     MODEL FOR THE KALMAN FILTER RECURSIONS */

/*   IF TRANSITION PROBABILITY IS VERY SMALL, DON'T */
/*   BOTHER WITH LIKELIHOOD CALCULATION: */

    const float TRANS_PROB_THRESHOLD = 1e-4f;
    if (*pin <= TRANS_PROB_THRESHOLD) {
		*lkhdj = 0.f;
		return 0;
    }

/*   OBTAIN STATE-DEPENDENT MODEL PARAMETERS: */
    real qa, hz, phi;
    model_(dur, kelem, ilrate, ixs, &phi, &qa, &hz);

/* 	GET PREVIOUS ESTIMATES FOR PATH IP */

    real ykk = ykkip[ip - 1];
    real pkk = pkkip[ip - 1];

/*  IMPLEMENT KALMAN FILTER FOR THIS TRANSITION */

    real ypred = phi * ykk;
    real ppred = phi * pkk * phi + qa;
    real pz = hz * ppred + rn;
    real pzinv = 1.f / pz;
    real g = ppred * hz * pzinv;
    real pest = (1.f - g * hz) * ppred;
    real zr = z - hz * ypred;

    ykksv[jnode - 1] = ypred + g * zr;
    pkksv[jnode - 1] = pest;
    if (ykksv[jnode - 1] <= .01f) {
		ykksv[jnode - 1] = .01f;
    }
/* Computing 2nd power */
    real a = .5f*pzinv*(zr * zr);
    if (a > 1e3f) {
		*lkhdj = 0.;
		return 0;
    }
    *lkhdj = (1.f / sqrt(pz)) * exp(-a);
 //   printf("\nz:%f a:%f lkhdj:%f israte:%d ilrate:%d dur:%f",*z,a,*lkhdj,*israte,*ilrate,*dur);
    return 0;
} /* kalfil_ */

