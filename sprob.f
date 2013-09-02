	SUBROUTINE SPROB (P,ISAVE,ILRSAV,PELM,KHAT,SPDHAT,PX)


C
C    SPROB COMPUTES THE POSTERIOR PROBS OF THE ELEMENT
C    STATES, DATA RATE STATES, AND KEYSTATES BY SUMMING
C    OVER THE APPROPRIETE PATHS.

C
C    VARIABLE:

C    P-		INPUT PATH PROBABILITIES
C    ISAVE- 	NUMBER OF PATHS SAVED
C    PSELEM-	OUTPUT ELEMENT PROB
C    SPDHAT-	OUTPUT SPEED ESTIMATE (DATA RATE WPM)
C    PX- 	OUTPUT KEYSTATE PROBABILITY
C



	DIMENSION P(750),PSELEM(6),ILRSAV(750)



C	INITIALIZE:

	SPDHAT=0.
	PX=0.

C	FOR EACH STATE EXTENSION OF PATH M:
C	OBTAIN ELEMENT STATE PROBS,KEYSTATE PROBS,SPEED EST:

	DO 100 K=1,6
	PSELEM(K)=0.
	
	DO 100 I=1,5
	N=(I-1)*6+K

	DO 100 M=1,ISAVE
	J=(M-1)*30+N
	PSELEM(K)=PSELEM(K)+P(J)
	SPDHAT=SPDHAT+ILRSAV(J)*P(J)
	IF(K.GT.2) GO TO 100
	PX=PX+P(J)
100	CONTINUE

	PELM=0.

	DO 200 K=1,6
	IF(PSELEM(K).LT.PELM) GO TO 200
	PELM=PSELEM(K)
	KHAT=K

200	CONTINUE

	RETURN
	END
