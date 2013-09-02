	SUBROUTINE RCVR (ZIN, ZOUT)

C
C       THIS SUBROUTINE CONVERTS THE INPUT SIGNAL AT
C       RADIAN FREQ  WC TO 1000 Hz.
C
C

	COMMON/BLK1/TAU/BLK2/WC
	DATA THETA/0./,THETLO/0./

	THETA = THETA + WC*TAU
	THETA = AMOD(THETA,6.28319)

	ZI = ZIN*COS(THETA)
	ZQ = ZIN*SIN(THETA)
	ZILP = ZILP+ .070* (ZI-ZILP)
	ZQLP = ZQLP+ .070* (ZQ-ZQLP)

	THETLO = THETLO+6283.2* TAU
	THETLO = AMOD(THETLO, 6.28319)

	ZOUT = ZILP*COS(THETLO) + ZQLP* SIN(THETLO)

	RETURN
	END
