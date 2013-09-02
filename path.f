	SUBROUTINE PATH(IP,LAMBDA,DUR,ILRATE,LAMSAV,DURSAV,ILRSAV)

C
C  PATH COMPUTES THE LTR STATE, DURATION, AND DATA RATE OF
C  EACH NEW PATH EXTENDED TO STATE N
C  
C  VARIABLES:
C  IP-		SAVED PATH IDENTITY
C  LAMBDA-	LTR STATE OF SAVED PATH
C  DUR-		DURATION OF ELEMENT ON SAVED PATH
C  ILRATE-	DATA RATE OF ELEMENT ON SAVED PATH
C  LAMSAV-	NEW LTR STATES FOR EACH PATH EXTENSION
C  DURSAV-	NEW ELEM DURATIONS FOR EACH PATH EXTENSION
C  ILRSAV-	NEW DATA RATES FOR EACH PATH EXTENSION
C  J-		NEW PATH IDENTITY

C  THE LETTER TRANSITION TABLE, MEMFCN, IS STORED IN COMMON.
C 

	DIMENSION LAMSAV(750),DURSAV(750),ILRSAV( 750)
	DIMENSION MEMFCN(400,6),IELMST(400),ILAMI(16)
	DIMENSION ILAMX(6),ISX(6),MEMDEL(6,6)

	COMMON/BLKLAM/IELMST,ILAMI,ILAMX
	COMMON/BLKMEM/MEMFCN
	COMMON/BLKS/ISX
	COMMON/BLKRAT/MEMDEL



C  FOR EACH ELEM STATE K, AND EACH SPEED I, COMPUTE:
C 
	DO 100 K=1,6
	DO 100 I=1,5


C  STATE IDENTITY N:
C 

	N=(I-1)*6+K


C  NEW PATH IDENTITY:
C 

	J=(IP-1)*30+N

C  NEW LTR STATE:
C 


	IF(LAMBDA .NE. 0) GO TO 50
	LAMSAV(J)=0    
	GO TO 100

50	LAMSAV(J)=MEMFCN(LAMBDA,K)
	IF(LAMSAV(J).EQ. 0) GO TO 100

C  NEW DURATION:
C  OBTAIN KEYSTATE OF SAVED PATH AND NEW STATE:


	ILELM=ILAMI(IELMST(LAMBDA))
	IXL=ILAMX(ILELM)
	IXS=ISX(K)

C CALCULATE DURATION - ADD SAMPLE DURATION 5 ms FOR EACH VALID PATH 

	DURSAV(J)=DUR*(1-IXS-IXL+2*IXS*IXL)+5.

C	NEW DATA RATE:

	ILRSAV(J)=ILRATE+(I-3)*MEMDEL(ILELM,K)
	GOTO 100

	PRINT 75, J, ILRSAV(J),LAMBDA, DUR, DURSAV(J) 
75	FORMAT('PATH:',3(I3,2X), 2(F8.3,2X))

100	CONTINUE
200	RETURN
	END

