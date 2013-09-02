/* morse.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Main program */ int MAIN__()
{
    /* Initialized data */

    static real rn = (float).1;
    static integer np = 0;

    /* Format strings */
    static char fmt_4[] = "(\002MORSE:, X, PMAX, IPMAX, ELMHAT, ZSIG,RN,SPDH\
AT,ZDET\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer n;
    static real x, z__;
    static integer n1, n2;
    static real s1[512], s2[512], s3[512], s4[512], px;
    static integer imax, xhat;
    static real pmax, zdet, zsig;
    extern /* Subroutine */ int rcvr_(real *, real *);
    static real zrcv, zout;
    extern /* Subroutine */ int noise_(real *, real *, real *), initl_(), 
	    stats_(real *, real *, real *, integer *, real *, real *, real *, 
	    real *, integer *), bpfdet_(real *, real *);
    static integer elmhat;
    static real spdhat;
    extern /* Subroutine */ int simsgi_(real *, real *), proces_(real *, real 
	    *, integer *, real *, integer *, integer *, real *, integer *, 
	    real *);
    static integer ltrhat;
    extern /* Subroutine */ int inputl_();

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, fmt_4, 0 };


/* 	CALL SRAND(86456) */
    initl_();
    inputl_();
    s_wsfe(&io___3);
    e_wsfe();
L1:
    for (n1 = 1; n1 <= 512; ++n1) {
	for (n2 = 1; n2 <= 18; ++n2) {
	    simsgi_(&x, &zsig);
	    rcvr_(&zsig, &zrcv);
	    bpfdet_(&zrcv, &zdet);
	    ++np;
/* 	DECIMATE 4 kHz by 40  down to 100Hz - 5 ms sample time for PROCES */
	    if (np < 40) {
		goto L3;
	    }
	    np = 0;
	    noise_(&zdet, &rn, &zout);
/* 	RN = RAND() */
	    rn = (float).01;
	    proces_(&x, &rn, &xhat, &px, &elmhat, &ltrhat, &spdhat, &imax, &
		    pmax);
L3:
	    ;
	}
	n = n1;
/* 	PRINT 5, X,  PMAX, IMAX, ELMHAT,ZSIG, RN, SPDHAT, ZDET */
/* 5	FORMAT('MORSE:',2(F10.3,2X),2(I7,1X),4(F10.3,2X) ) */
	stats_(&zdet, &z__, &px, &xhat, s1, s2, s3, s4, &n);
/* L2: */
    }
/* 	CALL DISPLA(S1,S2,S3,S4) */
    goto L1;
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int morse_ () { MAIN__ (); return 0; }
#ifdef __cplusplus
	}
#endif
