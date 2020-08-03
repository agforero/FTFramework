C
C     Include file oscltx.inc
C

C$ Abstract
C
C     This file declares public parameters related to 
C     the SPICELIB routine OSCLTX.
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 25-JAN-2017 (NJB)
C
C-&


C
C     Size of the OSCLTX output argument array ELTS.
C
C     The output array ELTS is intended to contain unused space to hold
C     additional elements that may be added in a later version of this
C     routine. In order to maintain forward compatibility, user
C     applications should declare ELTS as follows:
C
C     DOUBLE PRECISION   ELTS( OSCXSZ )
C

      INTEGER               OSCXSZ
      PARAMETER           ( OSCXSZ = 20 )


C
C     End of include file oscltx.inc
C
