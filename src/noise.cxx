// ----------------------------------------------------------------------------
// noise.c --  bayesian morse code decoder 
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
// ----------------------------------------------------------------------------

#include "bmorse.h"

int morse::noise_(double zin, real *rn, real *z)
{
    /* Initialized data */

    static real ylong[200] = { 1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f };
    static real yshort[50] = { 1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	    1.f,1.f,1.f,1.f,1.f,1.f,1.f };
    static integer kl = 1;
    static integer kkl = 1;
    static integer ks = 1;
    static integer kks = 1;
    static real ymin1 = 1.f;
    static real ymin2 = 1.f;
    static real ymavg = .05f;

    /* System generated locals */
    integer i1;


/* 	THIS SUBROUTINE ESTIMATES THE NOISE POWER IN THE */
/* 	ENVELOPE DETECTED OUTPUT FOR USE BY THE KALMAN */
/* 	FILTERS. AN ESTIMATE OF THE NOISE POWER IS */
/* 	ALSO SUBTRACTED FROM THE ENVELOPE DETECTED SIGNAL */
/* 	IN ORDER TO PRODUCE A ZERO-MEAN NOISE PROCESS */
/* 	AT THE INPUT TO THE MORSE PROCESSOR. */

    ++kl;
    if (kl == 201) {
	kl = 1;
    }
    ++ks;
    if (ks == 51) {
	ks = 1;
    }
    ++kkl;
    if (kkl >= 200) {
	kkl = 200;
    }
    ++kks;
    if (kks >= 50) {
	kks = 50;
    }
    if (kks > 2) {
        ylong[kl - 1] = zin;
        yshort[ks - 1] = zin;
        ymin1 = zin;
        ymin2 = zin;
    }

    i1 = kkl;
    for (int i = 1; i <= i1; ++i) {
	if (ylong[i - 1] > ymin1) {
	    continue;
	}
	ymin1 = ylong[i - 1];
    }
    i1 = kks;
    for (int i = 1; i <= i1; ++i) {
	if (yshort[i - 1] > ymin2) {
	    continue;
	}
	ymin2 = yshort[i - 1];
    }
    real ymin = ymin1;
    if (ymin2 < ymin1) {
	ymin = ymin2;
    }
    ymavg += (ymin - ymavg) * .004f;
    *rn = ymavg * .3f;
    if (*rn < .005f) {
	*rn = .005f;
    }
    *z = (zin - ymavg * 2.4f - .05f) * 1.1f;
    return 0;
} /* noise_ */

