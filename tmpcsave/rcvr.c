/* rcvr.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include "common.h"

/* Common Block Declarations */


/* Table of constant values */

static real c_b2 = 6.28319f;

/* Subroutine */ int rcvr_(real *zin, real *zout)
{
    /* Initialized data */

    static real theta = 0.f;
    static real thetlo = 0.f;

    /* Builtin functions */
    double r_mod(real *, real *), cos(doublereal), sin(doublereal);

    /* Local variables */
    static real z1, zq, zglp, zilp, zqlp;


/*       THIS SUBROUTINE CONVERTS THE INPUT SIGNAL AT */
/*       RADIAN FREQ  WC TO 1000 Hz. */


    theta += blk2.wc * blk1.tau;
    theta = r_mod(&theta, &c_b2);
    z1 = *zin * cos(theta);
    zq = *zin * sin(theta);
    zilp += (z1 - zilp) * .07f;
    zqlp += (zq - zqlp) * .07f;
    thetlo += blk1.tau * 6283.2f;
    thetlo = r_mod(&thetlo, &c_b2);
    *zout = zilp * cos(thetlo) + zglp * sin(thetlo);
    return 0;
} /* rcvr_ */

