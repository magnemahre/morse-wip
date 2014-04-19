// ----------------------------------------------------------------------------
// proces.c --  bayesian morse code decoder 
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

int morse::proces_(real *z, real *rn, integer *xhat, real *px, integer *elmhat,  real *spdhat, integer *imax, real *pmax, int spd)
{
    /* Initialized data */

    static integer isave = 25;
    static integer lambda[25]; 
	static integer ilrate[25]; 
    static real dur[25];
    static integer pathsv[25]; 
    static integer sort[25];
    
    static real p[750];
    static integer lamsav[750];
    static real dursav[750];
    static integer ilrsav[750];
    static real pin[750]	/* was [25][30] */, lkhd[750];


    /* Local variables */
//    static real pin[30][25];		// N: 5 spd x 6 morse element states I: 25 paths  - transition probability from path I to state N 
//    static real lkhd[30][25];		// 5 speeds x 6 morse element states x 25 paths = 750 likelyhoods 
    
    static real pelm = 0;
    static integer ipmax;
    static int init =1;
 

/* 	THIS SUBROUTINE IMPLEMENTS THE PROCESSING ALGORITHM */
/* 	FOR JOINT DEMODULATION, DECODING, AND TRANSLATION OF */
/* 	THE RECEIVED MORSE PROCESS. IT TAKES IN A NEW MEASURE- */
/* 	MENT, Z, OF THE DETECTED SIGNAL EVERY 5 MSEC AND PRO- */
/* 	DUCES AN ESTIMATE OF THE CURRENT KEYSTATE, ELEMENT */
/* 	STATE, AND LETTER OF THE RECEIVED SIGNAL. */

/* 	VIS */
/* 	DEFINITIONS OF VARIABLE NAMES: */
/* 	Z- INPUT SAMPLE OF DETECTED SIGNAL */
/* 	RN- INPUT NOISE POWER ESTIMATE */
/* 	XHAT- OUTPUT ESTIMATE OF KEYSTATE */
/* 	ELMHAT- OUTPUT ESTIMATE OF ELEMENT STATE */
/* 	LTRHAT- OUTPUT ESTIMATE OF LETTER STATE */


/* 	ISAVE- NO. OF PREVIOUS PATHS SAVED */
/* 	IPATH- IDENTITY OF SAVED PATH */
/* 	LAMBDA (I) - IDENTITY OF LTR STATE OF SAVED PATH I */
/* 	DUR (I)- DURATION OF ELEMENT ON PATH I */
/* 	ILRATE(I)- IDENTITY OF DATA RATE ON PATH I */
/* 	PIN (I,N)- COMPUTED TRANS PROB FROM PATH I TO STATE N */
/* 	LAMSAV(J)- IDENTITY OF LTR STATE AT NEW NODE J */
/* 	ILRSAV(J)- IDENTITY OF DATA RATE AT NEW NODE J */
/* 	LKHD(J)- LIKELIHOOD VALUE FOR NODE J */
/* 	P(J)- COMPUTED POSTERIOR PROB OF PATH ENDING AT NEW NODE J */
/* 	PSELEM(K)-COMPUTED POSTERIOR PROB OF ELEM K */
/* 	SPDHAT -COND MEAN ESTIMATE OF INSTANT DATA RATE */
/* 	PX- POSTERIOR PROB THAT KEYSTATE EQUALS 1 */

/* 	THE FOLLOWING SUBROUTINES ARE UTILIZED: */
/* 	TRPROB- COMPUTES TRANSITION PROBABILITIES */
/* 	PATH- COMPUTES IDENTITY OF NEW PATHS */
/* 	LIKHD- COMPUTES THE LIKELIHOOD OF EACH PATH EXTENSION */
/* 	PROBP- COMPUTES POSTERIOR PROBS OF EACH NEW PATH */
/* 	SPROB- COMPUTES POSTERIOR PROBS OF EACH STATE */
/* 	SAVE- SAVES THE HIGHEST PROB PATHS */
/* 	TRELIS- FORMS A TRELLIS OF SAVED PATHS */
/* 	TRANSL- TRANSLATES THE LETTER ESTIMATE */

/* 	ALL TABLES OF CONSTANTS ARE STORED IN COMMON. */
/* 	FOR EACH SAVED PATH, COMPUTE: */
/* 	TRANSITION PROBABILITY TO NEW STATE (TRPROB); */
/* 	IDENTITY OF EACH NEW PATH EXTENDED (PATH); */
/* 	LIKELIHOOD OF EACH STATE EXTENSION (LIKHD); */


	if (init) {
		for(int i=0;i<25;i++) {
			lambda[i] = 5;
			ilrate[i]= ((i/5+1)*10);
			dur[i]=9e3f;
			pathsv[i]=5;
		}
		for(int i=0;i<750;i++) {
			p[i]=1.f;
			lamsav[i]=5;
			dursav[i]=0.f;
			ilrsav[i]=20;
		}
		init = 0;
	}

    integer i1 = isave;
    for (int i = 1; i <= i1; ++i) {
		integer ipath = i;
		trprob_(&ipath, &lambda[i - 1], &dur[i - 1], &ilrate[i - 1], pin);
		path_(&ipath, &lambda[i - 1], &dur[i - 1], &ilrate[i - 1],lamsav, dursav, ilrsav);
		likhd_(z, rn, &ipath, &lambda[i - 1], &dur[i - 1], &ilrate[i- 1], pin, lkhd);
    }
/* 	HAVING OBTAINED ALL NEW PATHS, COMPUTE: */
/* 	POSTERIOR PROBABILITY OF EACH NEW PATH(PROBP); */
/* 	POSTERIOR PROBABILITY OF KEYSTATE, ELEM STATE, */
/* 	CONDITIONAL MEAN ESTIMATE OF SPEED(SPROB); */
    probp_(p, pin, &isave, lkhd);
    sprob_(p, &isave, ilrsav, &pelm, elmhat, spdhat, px);

    *xhat = 0;
    if (*px > .5f) {
		*xhat = 1;
    }
/* 	SAVE THE PATHS WITH HIGHEST PROBABILITY, AND */
/* 	STORE THE VALUES CORRESPONDING TO THESE PATHS: */

    savep_(p, pathsv, &isave, imax, lamsav, dursav, ilrsav, lambda, dur, ilrate, sort, pmax);

/* 	UPDATE TRELLIS WITH NEW SAVED NODES, AND */
/* 	OBTAIN LETTER STATE ESTIMATE: */

    integer retstat=trelis_(isave, pathsv, lambda, *imax, &ipmax);

    return retstat;
} /* proces_ */

