// ----------------------------------------------------------------------------
// model.c --  bayesian morse code decoder 
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

#include <math.h>
#include "bmorse.h"




/* Table of constant values */

//static doublereal c_b5 = 10.;

int morse::model_(const real dur, const integer ielm, const integer ilr, const integer ixs, real *phi, real *qa, real *hz)
{

/* 	THIS SUBROUTINE COMPUTES THE PARAMETERS OF THE */
/* 	OBSERVATION STATE TRANSITION MATRIX PHI, THE */
/* 	MEASUREMENT MATRIX, AND THE COVARIANCES. */

/* 	VARIABLES: */
/* 		DUR-	INPUT ELEMENT DURATION */
/* 		IELM-	INPUT ELEMENT TYPE */
/* 		ILR-	INPUT SAVED RATE */
/* 		ISR-	INPUT RATE OF NEW STATE */
/* 		IXS-	INPUT KEYSTATE OF NEW STATE */
/* 		PHI-	OUTPUT STATE TRANSITION MATRIX ENTRY FOR SIGNAL AMPLITUDE STATE */
/* 		QA-	OUTPUT COVARIANCE FOR AMPLITUDE STATE */
/* 		HZ-	OUTPUT MEASUREMENT MATRIX VALUE */

/* 	COMPUTE MEASUREMENT COEFFICIENT: */
    *hz = (real) ixs;
    
/* 	COMPUTE PHI AND AMPLITUDE STATE VARIANCE (Q): */
    real r1 = 1200.f / ilr;
    real bauds = dur / r1;
    if (bauds >= 14.f) {
		bauds = 14.f;
    }

    if (ielm < 3) {
		*qa = 1e-4f;
		*phi = 1.f;
		return 0;
    }

    if (ixs != 0) {
		*phi = 1.f;
		*qa = exp((bauds - 14.f) * .6f) * .15f;
		*qa += bauds * .01f * exp((1.f - bauds) * .2f);
		return 0;
    }

    real xsamp = r1 * 22.4f;
    doublereal d1 = (doublereal) (-2 / xsamp);
    *phi = pow(10.0, d1);

    if (bauds >= 14.f) {
		*phi = 1.f;
    }
    *qa = 0.f;
    return 0;
} /* model_ */

