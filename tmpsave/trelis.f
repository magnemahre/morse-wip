	SUBROUTINE TRELIS (ISAVE, PATHSV, LAMBDA, IMAX)

C
C    THIS SUBROUTINE STORES THE SAVED NODES AT EACH
C    STAGE AND FORMS THE TREE OF SAVED PATHS LINKING
C    THE NODES. DECODING IS ACCOMPLISHED BY FINDING
C    THE CONVERGENT PATH IF IT OCCURS WITHIN A MAXIMUM
C    DELAY SET BY THE PARAMETER NDELAY. IF CONVERGENCE 
C    TO A SINGLE PATH DOES NOT OCCURS, THEN DECODING IS 
C    DONE BY READING THE LETTER ON THE PATH WITH HIGHEST 
C    PROBABILITY


	INTEGER PATHSV
	DIMENSION PATHSV (25), LAMBDA (25), PTHTRL (200,25)
	DIMENSION LMDSAV (200,25), IPNOD (25), LTRSV (200)
	COMMON / BLKEND / IEND

	DATA PTHTRL/5000*5/, LMDSAV/ 5000*5/
	DATA N/0/, NDELAY /200/
	DATA IPNOD/25*1/, NCALL/0/, NMAX/0/, MMAX/0/

C	KEEP AVERAGE OF ISAVE, NDEL FOR DATA ANALYSIS:

	NCALL = NCALL + 1
	IF (IEND. NE. 1) GO TO 10
	ISAVG = XSAVG
	NDLAVG = XDLAVG
	IEND = 0

	PRINT 2000, ISAVG, NDLAVG
2000	FORMAT (1X, 'AVG NO OF PATHS SAVED: ',I2, 2X,'AVG DECODE DELAY: ',
     &I3)
C	2 I3)

	PRINT 330, XMMAX, XNMAX
330	FORMAT (1X, 'PERCENT OF TIME PATHS = 25: ', F3.2, 2X,
     &'PERCENT OF TIME DELAY = 200: ', F3.2)
C	2 'PERCENT OF TIME DELAY = 200: ', F3.2)
	READ 2000, WAIT

10	XSAVG = (XSAVG * (NCALL - 1) + ISAVE)/ NCALL
	XDLAVG = (XDLAVG * (NCALL - 1) + NDEL)/ NCALL
	IF (NDEL. NE. NDELAY) GO TO 20
	NMAX = NMAX + 1
	XNMAX = NMAX
	XNMAX = XNMAX/ NCALL

20	IF (ISAVE. NE. 25) GO TO 30
	MMAX = MMAX + 1
	XMMAX = MMAX
	XMMAX = XMMAX/ NCALL

30	CONTINUE

C	STORE PATHSV AND CORRESPONDING LAMBDA IN THE
C	TRELLIS USING A CIRCULAR BUFFER OF LENGTH NDELAY :

	N = N + 1
	IF (N. EQ. NDELAY +1) N = 1

	DO 100 I = 1, ISAVE
	PTHTRL(N, I)=PATHSV(I)
	LMDSAV(N,I)=LAMBDA(I)
100	CONTINUE


C	PERFORM DYNAMIC PROGRAM ROUTINE TO FIND CONVERGENT PATH:

	K=0
	DO 180 I=1,ISAVE
	IPNOD(I)=I
180	CONTINUE		

190	K=K+1
	IF(K.EQ.NDELAY) GO TO 700

	DO 200 IP=1,ISAVE
	I=N-K+1
	IF (I .LE. 0) I=NDELAY+I
	IPNOD(IP)=PTHTRL(I,IPNOD(IP))
	IF (IP .EQ. IMAX) IPMAX=IPNOD(IP)
200	CONTINUE


C	IF ALL NODES ARE EQUAL,THEN PATHS CONVERGE:

	DO 300 IEQ=2,ISAVE
	IF (IPNOD(1).NE. IPNOD(IEQ)) GO TO 190
300	CONTINUE

C	PATHS CONVERGE; SET NDEL:

	NDEL=K+1

C	IF POINT OF CONVERGENCE IS SAME AS IT WAS ON
C	LAST CALL, THEN NO NEED TO RE-DECODE SAME NODE:

	IF (NDEL.EQ. NDELST+1) GO TO 800

C	IF POINT OF CONVERGENCE OCCURS AT SAME DELAY AS LAST CALL, THEN TRANSLATE:

	IF(NDEL.NE. NDELST) GO TO 350
	I=N-NDEL+1
	IF(I.LE. 0) I=NDELAY+1
	LTR=LMDSAV(I,IPNOD(1))
	CALL TRANSL(LTR)
	GO TO 800

C	OTHERWISE,POINT OF CONVERGENCE HAS OCCURED
C	EARLIER ON THIS CALL, SO NEED TO TRANSLATE
C	EVERYTHING ON THE CONVERGENT PATH FROM
C	PREVIOUS POINT OF CONVERGENCE TO THIS POINT:


350	K0=0
	IP=IPNOD(1)	

	DO 400 K=NDEL,NDELST
	KD=KD+1
	I=N-K+1
	IF (I.LE. 0) I=NDELAY+I
	LTRSV(KD) = LMDSAV(I,IP)
	IP=PTHTRL(I,IP)
400	CONTINUE

C	REVERSE ORDER OF DECODED LETTERS, SINCE THEY
C	WERE OBTAINED FROM THE TRELLIS IN REVERSE;
C	TRANSLATE EACH:

	DO 500 I=1,KD
	LTR=LTRSV(KD-I+1)
	CALL TRANSL(LTR)
500	CONTINUE
	GO TO 800

700	CONTINUE
C	PATHS HAVE NOT CONVERGED AT MAXIMUM ALLOWABLE
C	DELAY, SO TRANSLATE WHAT IS ON HIGHEST
C	PROBABILITY PATH:

	NDEL=NDELAY
	I=N-NDELAY+1
	IF(I .LE. 0) I=NDELAY+I
	LTR=LMDSAV(I,IPMAX)
	CALL TRANSL(LTR)

C	PRUNE AWAY NODES WHICH ARE NOT ON THIS PATH:

	DO 750 K=1,ISAVE
	IF (IPNOD(K) .EQ. IPMAX) GO TO 750
	LAMBDA(K)=0
750	CONTINUE


800	NDELST=NDEL
	RETURN
	END


