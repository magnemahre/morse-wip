// ----------------------------------------------------------------------------
// likhd.c --  bayesian morse code decoder 
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
#include <stdio.h>

int morse::likhd_(const real z, const real rn, const integer ip, const integer lambda,
	 const real dur, const integer ilrate, real *p, real *lkhd)
{
 
    /* Local variables */
    static real lkhdj;
 


/* 	THIS SUBROUTINE CALCULATES,FOR EACH PATH */
/* 	EXTENSION TO STATE N, THE LIKELIHOOD OF THAT */
/* 	TRANSITION GIVEN THE MEASUREMENTZ. IT USES */
/* 	AN ARRAY OF LINEAR (KALMAN) FILTERS TO DO SO. */

/* 	VARIABLES: */
/* 	Z- 	INPUT MEASUREMENT */
/* 	RN-	INPUT NOISE POWER ESTIMATE */
/* 	IP-	INPUT SAVED PATH IDENTITY */
/* 	LAMBDA-	INPUT SAVED LTR STATE IDENTITY */
/* 	DUR-	INPUT SAVED DURATION OF ELEMENT ON PATH IP */
/* 	ILRATE-	INPUT SAVED DATA RATE (SPEED) */
/* 	P-	INPUT TRANSITION PROBABILITIES */
/* 	LKHD-	OUTPUT COMPUTED LIKELIHOODS FOR EACH TRANS */

/*  SUBROUTINES USED: */
/* 	KALFIL-KALMAN FILTER FOR EACH NEW PATH */

/*   OBTAIN SAVED KEYSTATE: */
    /* Parameter adjustments */
    --lkhd;
    p -= 26;

    /* Function Body */
    if (lambda != 0) {
        integer kelem = ilami[ielmst[lambda - 1] - 1];
        integer ilx = ilamx[kelem - 1];
/* 	FOR EACH STATE: */
        for (int k = 1; k <= 6; ++k) {
		for (int i = 1; i <= 5; ++i) {
/* 	OBTAIN KEYSTATE, RATE STATE, STATE N, NEW NODE: */
			integer ixs = isx[k - 1];
			integer israte = i;
			integer n = (i - 1) * 6 + k;
			integer j = (ip - 1) * 30 + n;
			real pin = p[ip + n * 25];
/* 	COMPUTE AND STORE LIKELIHOOD: */
			kalfil_(z, ip, rn, ixs, kelem, j, dur, ilrate,&pin, &lkhdj);
			lkhd[j] = lkhdj;
//			goto L100;
			if (pin > 1e-6f) {

//printf("\nz:%f ip:%3d rn:%f ilx:%d ixs:%d kelem:%d j:%4d israte:%d dur:%4.1f ilrate:%d pin:%f lkhd:%f",(double)*z,(int)*ip,(double)*rn,(int)ilx,(int)ixs,(int)kelem,(int)j,(int)israte,(double)*dur,(int)*ilrate,(double)pin,(double)lkhdj);
			}
			
		}
        }
    }

    return 0;
} /* likhd_ */

