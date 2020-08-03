C$Procedure ZZDSCM ( SGP4 deep space common calculations )

      SUBROUTINE ZZDSCM ( EPOCH,  ECCP,  ARGPP,   TC,     INCLP,
     .                    NODEP,  NP,    SNODM,   CNODM,
     .                    SINIM,  COSIM, SINOMM,  COSOMM,
     .                    DAY,    E3,    EE2,     ECCM,   EMSQ, GAM,
     .                    PEO,    PGHO,  PHO,     PINCO,  PLO,
     .                    RTEMSQ, SE2,   SE3,     SGH2,   SGH3, SGH4,
     .                    SH2,    SH3,   SI2,     SI3,    SL2,  SL3,
     .                    SL4,    S1,    S2,      S3,     S4,   S5,
     .                    S6,     S7,    SS1,     SS2,    SS3,  SS4,
     .                    SS5,    SS6,   SS7,     SZ1,    SZ2,  SZ3,
     .                    SZ11,   SZ12,  SZ13,    SZ21,   SZ22, SZ23,
     .                    SZ31,   SZ32,  SZ33,    XGH2,   XGH3, XGH4,
     .                    XH2,    XH3,   XI2,     XI3,    XL2,  XL3,
     .                    XL4,    XN,    Z1,      Z2,     Z3,   Z11,
     .                    Z12,    Z13,   Z21,     Z22,    Z23,  Z31,
     .                    Z32,    Z33,   ZMOL,    ZMOS )

C$ Abstract
C
C     This subroutine provides deep space common items used by both the
C     secular and periodics subroutines.
C
C     This routine previously had the name DPPER, but the functions
C     inside weren't well organized.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     None.
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION         EPOCH
      DOUBLE PRECISION         ECCP
      DOUBLE PRECISION         ARGPP
      DOUBLE PRECISION         TC
      DOUBLE PRECISION         INCLP
      DOUBLE PRECISION         NODEP
      DOUBLE PRECISION         NP
      DOUBLE PRECISION         SNODM
      DOUBLE PRECISION         CNODM
      DOUBLE PRECISION         SINIM
      DOUBLE PRECISION         COSIM
      DOUBLE PRECISION         SINOMM
      DOUBLE PRECISION         COSOMM
      DOUBLE PRECISION         DAY
      DOUBLE PRECISION         E3
      DOUBLE PRECISION         EE2
      DOUBLE PRECISION         ECCM
      DOUBLE PRECISION         EMSQ
      DOUBLE PRECISION         GAM
      DOUBLE PRECISION         RTEMSQ
      DOUBLE PRECISION         SE2
      DOUBLE PRECISION         PEO
      DOUBLE PRECISION         PGHO
      DOUBLE PRECISION         PHO
      DOUBLE PRECISION         PINCO
      DOUBLE PRECISION         PLO
      DOUBLE PRECISION         SE3
      DOUBLE PRECISION         SGH2
      DOUBLE PRECISION         SGH3
      DOUBLE PRECISION         SGH4
      DOUBLE PRECISION         SH2
      DOUBLE PRECISION         SH3
      DOUBLE PRECISION         SI2
      DOUBLE PRECISION         SI3
      DOUBLE PRECISION         SL2
      DOUBLE PRECISION         SL3
      DOUBLE PRECISION         SL4
      DOUBLE PRECISION         S1
      DOUBLE PRECISION         S2
      DOUBLE PRECISION         S3
      DOUBLE PRECISION         S4
      DOUBLE PRECISION         S5
      DOUBLE PRECISION         S6
      DOUBLE PRECISION         S7
      DOUBLE PRECISION         SS1
      DOUBLE PRECISION         SS2
      DOUBLE PRECISION         SS3
      DOUBLE PRECISION         SS4
      DOUBLE PRECISION         SS5
      DOUBLE PRECISION         SS6
      DOUBLE PRECISION         SS7
      DOUBLE PRECISION         SZ1
      DOUBLE PRECISION         SZ2
      DOUBLE PRECISION         SZ3
      DOUBLE PRECISION         SZ11
      DOUBLE PRECISION         SZ12
      DOUBLE PRECISION         SZ13
      DOUBLE PRECISION         SZ21
      DOUBLE PRECISION         SZ22
      DOUBLE PRECISION         SZ23
      DOUBLE PRECISION         SZ31
      DOUBLE PRECISION         SZ32
      DOUBLE PRECISION         SZ33
      DOUBLE PRECISION         XGH2
      DOUBLE PRECISION         XGH3
      DOUBLE PRECISION         XGH4
      DOUBLE PRECISION         XH2
      DOUBLE PRECISION         XH3
      DOUBLE PRECISION         XI2
      DOUBLE PRECISION         XI3
      DOUBLE PRECISION         XL2
      DOUBLE PRECISION         XL3
      DOUBLE PRECISION         XL4
      DOUBLE PRECISION         XN
      DOUBLE PRECISION         Z1
      DOUBLE PRECISION         Z2
      DOUBLE PRECISION         Z3
      DOUBLE PRECISION         Z11
      DOUBLE PRECISION         Z12
      DOUBLE PRECISION         Z13
      DOUBLE PRECISION         Z21
      DOUBLE PRECISION         Z22
      DOUBLE PRECISION         Z23
      DOUBLE PRECISION         Z31
      DOUBLE PRECISION         Z32
      DOUBLE PRECISION         Z33
      DOUBLE PRECISION         ZMOL
      DOUBLE PRECISION         ZMOS

C$ Brief_I/O
C
C    See Detailed_input and Detailed_Output.
C     
C$ Detailed_Input
C
C    EPOCH       Epcoh of TLE set as Julian day value.
C
C    EP          Eccentricity
C
C    ARGPP       Argument of perigee
C
C    TC          Minutes past EPOCH, nominally zero.
C
C    INCLP       Inclination
C
C    NODEP       Right ascension of ascending node
C
C    NP          Mean motion
C
C$ Detailed_Output
C
C    SINIM       SIN of mean inclination.   
C
C    COSIM       COS of mean inclination.
C
C    SINOMM      Internal SGD4 parameter.
C
C    COSOMM      Internal SGD4 parameter.
C
C    SNODM       Internal SGD4 parameter.
C
C    CNODM       Internal SGD4 parameter.
C
C    DAY         Internal SGD4 parameter.
C
C    E3          Internal SGD4 parameter.
C
C    EE2         Internal SGD4 parameter.
C
C    EM          Eccentricity.
C
C    EMSQ        Eccentricity squared.
C
C    GAM         Internal SGD4 parameter.
C
C    PEO         Internal SGD4 parameter.
C
C    PGHO        Internal SGD4 parameter.
C
C    PHO         Internal SGD4 parameter.
C
C    PINCO       Internal SGD4 parameter.
C
C    PLO         Internal SGD4 parameter.
C
C    RTEMSQ      Internal SGD4 parameter.
C
C    SE2         Internal SGD4 parameter.
C
C    SE3         Internal SGD4 parameter.
C
C    SGH2        Internal SGD4 parameter.
C
C    SGH3        Internal SGD4 parameter.
C
C    SGH4        Internal SGD4 parameter.
C
C    SH2         Internal SGD4 parameter.
C
C    SH3         Internal SGD4 parameter.
C
C    SI2         Internal SGD4 parameter.
C
C    SI3         Internal SGD4 parameter.
C
C    SL2         Internal SGD4 parameter.
C
C    SL3         Internal SGD4 parameter.
C
C    SL4         Internal SGD4 parameter.
C
C    S1          S coeffcients
C
C    S2             ...
C
C    S3             ...
C
C    S4             ...
C
C    S5             ...
C
C    S6             ...
C
C    S7             ...
C
C    SS1            ...
C
C    SS2            ...
C
C    SS3            ...
C
C    SS4            ...
C
C    SS5            ...
C
C    SS6            ...
C
C    SS7            ...
C
C    SZ1            ...
C
C    SZ2            ...
C
C    SZ3            ...
C
C    SZ11           ...
C
C    SZ12           ...
C
C    SZ13           ...
C
C    SZ21           ...
C
C    SZ22           ...
C
C    SZ23           ...
C
C    SZ31           ...
C
C    SZ32           ...
C
C    SZ33           ...
C
C    XGH2        Internal SGD4 parameter.
C
C    XGH3        Internal SGD4 parameter.
C
C    XGH4        Internal SGD4 parameter.
C
C    XH2         Internal SGD4 parameter.
C
C    XH3         Internal SGD4 parameter.
C
C    XI2         Internal SGD4 parameter.
C
C    XI3         Internal SGD4 parameter.
C
C    XL2         Internal SGD4 parameter.
C
C    XL3         Internal SGD4 parameter.
C
C    XL4         Internal SGD4 parameter.
C
C    NM          Mean motion
C
C    Z1          Z coeffcients
C
C    Z2             ...
C
C    Z3             ...
C
C    Z11            ...
C
C    Z12            ...
C
C    Z13            ...
C
C    Z21            ...
C
C    Z22            ...
C
C    Z23            ...
C
C    Z31            ...
C
C    Z32            ...
C
C    Z33            ...
C
C    ZMOL        Internal SGD4 parameter.
C
C    ZMOS        Internal SGD4 parameter.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     None.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is based on the DSCOM code by David Vallado
C     corresponding to "Revisiting Spacetrack Report #3" [4].
C     The intent is to maintain the original Vallado algorithm,
C     changing code only to meet NAIF format standards and to
C     integrate with SPICELIB.
C
C        Capitalize all variables.
C
C        ENDIF replaced with END IF.
C
C        ENDDO replaced with END DO.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C   [1] Hoots, F. R., and Roehrich, R. L. 1980. "Models for
C       Propagation of the NORAD Element Sets." Spacetrack Report #3.
C       U.S. Air Force: Aerospace Defense Command.
C
C   [2] Hoots, Felix R. "Spacetrack Report #6: Models for Propagation
C       of Space Command Element Sets." Space Command,
C       U. S. Air Force, CO.
C
C   [3] Hoots, Felix R., P. W. Schumacher, and R. A. Glover. 2004.
C       History of Analytical Orbit Modeling in the U. S. Space
C       Surveillance System. Journal of Guidance, Control, and
C       Dynamics. 27(2):174-185.
C
C   [4] Vallado, David, Crawford, Paul, Hujsak, Richard,
C       and Kelso, T.S. 2006. Revisiting Spacetrack Report #3. Paper
C       AIAA 2006-6753 presented at the AIAA/AAS Astrodynamics
C       Specialist Conference, August 21-24, 2006. Keystone, CO.
C
C$ Author_and_Institution
C
C     David Vallado   (AGI)
C     E. D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 11-NOV-2014 (EDW)
C
C        Based on routine DSCOM, 28-JUN-2005, Vallado 2006 [4].
C
C-&

C$ Index_Entries
C
C   SGP4
C
C-&

C
C     Local Variables
C
      DOUBLE PRECISION         ZES
      DOUBLE PRECISION         A1
      DOUBLE PRECISION         A10
      DOUBLE PRECISION         A2
      DOUBLE PRECISION         A3
      DOUBLE PRECISION         A4
      DOUBLE PRECISION         A5
      DOUBLE PRECISION         A6
      DOUBLE PRECISION         A7
      DOUBLE PRECISION         A8
      DOUBLE PRECISION         A9
      DOUBLE PRECISION         BETASQ
      DOUBLE PRECISION         C1L
      DOUBLE PRECISION         C1SS
      DOUBLE PRECISION         CC
      DOUBLE PRECISION         CTEM
      DOUBLE PRECISION         STEM
      DOUBLE PRECISION         X1
      DOUBLE PRECISION         X2
      DOUBLE PRECISION         X3
      DOUBLE PRECISION         X4
      DOUBLE PRECISION         X5
      DOUBLE PRECISION         X6
      DOUBLE PRECISION         X7
      DOUBLE PRECISION         X8
      DOUBLE PRECISION         XNODCE
      DOUBLE PRECISION         XNOI
      DOUBLE PRECISION         ZCOSG
      DOUBLE PRECISION         ZCOSGL
      DOUBLE PRECISION         ZCOSGS
      DOUBLE PRECISION         ZCOSH
      DOUBLE PRECISION         ZCOSHL
      DOUBLE PRECISION         ZCOSI
      DOUBLE PRECISION         ZCOSIL
      DOUBLE PRECISION         ZCOSIS
      DOUBLE PRECISION         ZEL
      DOUBLE PRECISION         ZSING
      DOUBLE PRECISION         ZSINGL
      DOUBLE PRECISION         ZSINGS
      DOUBLE PRECISION         ZSINH
      DOUBLE PRECISION         ZSINHL
      DOUBLE PRECISION         ZSINI
      DOUBLE PRECISION         ZSINIL
      DOUBLE PRECISION         ZSINIS
      DOUBLE PRECISION         ZX
      DOUBLE PRECISION         ZY

      INTEGER                  LSFLG

C
C     SPICELIB routines.
C
      DOUBLE PRECISION         TWOPI
      LOGICAL                  RETURN

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSCM' )

C
C     Constants
C
      ZES    =  0.01675D0
      ZEL    =  0.05490D0
      C1SS   =  2.9864797D-6
      C1L    =  4.7968065D-7
      ZSINIS =  0.39785416D0
      ZCOSIS =  0.91744867D0
      ZCOSGS =  0.1945905D0
      ZSINGS = -0.98088458D0

C
C     DEEP SPACE PERIODICS INITIALIZATION
C
      XN     = NP
      ECCM   = ECCP
      SNODM  = DSIN(NODEP)
      CNODM  = DCOS(NODEP)
      SINOMM = DSIN(ARGPP)
      COSOMM = DCOS(ARGPP)
      SINIM  = DSIN(INCLP)
      COSIM  = DCOS(INCLP)
      EMSQ   = ECCM*ECCM
      BETASQ = 1.0D0-EMSQ
      RTEMSQ = DSQRT(BETASQ)

C
C     INITIALIZE LUNAR SOLAR TERMS
C
C     Note, EPOCH + 18261.5D0 corresponds to JD days since 
C     1899-12-31 12:00:00.
C
      PEO    = 0.0D0
      PINCO  = 0.0D0
      PLO    = 0.0D0
      PGHO   = 0.0D0
      PHO    = 0.0D0
      DAY    = EPOCH + 18261.5D0 + TC/1440.0D0
      XNODCE = DMOD(4.5236020D0 - 9.2422029D-4*DAY,TWOPI())
      STEM   = DSIN(XNODCE)
      CTEM   = DCOS(XNODCE)
      ZCOSIL = 0.91375164D0 - 0.03568096D0*CTEM
      ZSINIL = DSQRT(1.0D0 - ZCOSIL*ZCOSIL)
      ZSINHL = 0.089683511D0*STEM / ZSINIL
      ZCOSHL = DSQRT(1.0D0 - ZSINHL*ZSINHL)
      GAM    = 5.8351514D0 + 0.0019443680D0*DAY
      ZX     = 0.39785416D0*STEM/ZSINIL
      ZY     = ZCOSHL*CTEM + 0.91744867D0*ZSINHL*STEM
      ZX     = DATAN2(ZX,ZY)
      ZX     = GAM + ZX - XNODCE
      ZCOSGL = DCOS(ZX)
      ZSINGL = DSIN(ZX)

C
C     DO SOLAR TERMS
C
      ZCOSG = ZCOSGS
      ZSING = ZSINGS
      ZCOSI = ZCOSIS
      ZSINI = ZSINIS
      ZCOSH = CNODM
      ZSINH = SNODM
      CC    = C1SS
      XNOI  = 1.0D0 / XN

C
C     Loop over the lunar and solar term flags.
C
      DO LSFLG = 1,2

         A1 =   ZCOSG*ZCOSH + ZSING*ZCOSI*ZSINH
         A3 =  -ZSING*ZCOSH + ZCOSG*ZCOSI*ZSINH
         A7 =  -ZCOSG*ZSINH + ZSING*ZCOSI*ZCOSH
         A8 =   ZSING*ZSINI
         A9 =   ZSING*ZSINH + ZCOSG*ZCOSI*ZCOSH
         A10=   ZCOSG*ZSINI
         A2 =   COSIM*A7 + SINIM*A8
         A4 =   COSIM*A9 + SINIM*A10
         A5 =  -SINIM*A7 + COSIM*A8
         A6 =  -SINIM*A9 + COSIM*A10

         X1 =  A1*COSOMM + A2*SINOMM
         X2 =  A3*COSOMM + A4*SINOMM
         X3 = -A1*SINOMM + A2*COSOMM
         X4 = -A3*SINOMM + A4*COSOMM
         X5 =  A5*SINOMM
         X6 =  A6*SINOMM
         X7 =  A5*COSOMM
         X8 =  A6*COSOMM

         Z31= 12.0D0*X1*X1 - 3.0D0*X3*X3
         Z32= 24.0D0*X1*X2 - 6.0D0*X3*X4
         Z33= 12.0D0*X2*X2 - 3.0D0*X4*X4
         Z1 =  3.0D0* (A1*A1 + A2*A2) + Z31*EMSQ
         Z2 =  6.0D0* (A1*A3 + A2*A4) + Z32*EMSQ
         Z3 =  3.0D0* (A3*A3 + A4*A4) + Z33*EMSQ
         Z11= -6.0D0*A1*A5 + EMSQ* (-24.0D0*X1*X7-6.0D0*X3*X5)
         Z12= -6.0D0* (A1*A6 + A3*A5) + EMSQ*
     .           ( -24.0D0*(X2*X7+X1*X8) - 6.0D0*(X3*X6+X4*X5) )
         Z13= -6.0D0*A3*A6 + EMSQ*(-24.0D0*X2*X8 - 6.0D0*X4*X6)
         Z21=  6.0D0*A2*A5 + EMSQ*(24.0D0*X1*X5-6.0D0*X3*X7)
         Z22=  6.0D0* (A4*A5 + A2*A6) + EMSQ*
     .           (  24.0D0*(X2*X5+X1*X6) - 6.0D0*(X4*X7+X3*X8) )
         Z23=  6.0D0*A4*A6 + EMSQ*(24.0D0*X2*X6 - 6.0D0*X4*X8)
         Z1 = Z1 + Z1 + BETASQ*Z31
         Z2 = Z2 + Z2 + BETASQ*Z32
         Z3 = Z3 + Z3 + BETASQ*Z33
         S3 = CC*XNOI
         S2 = -0.5D0*S3 / RTEMSQ
         S4 = S3*RTEMSQ
         S1 = -15.0D0*ECCM*S4
         S5 = X1*X3 + X2*X4
         S6 = X2*X3 + X1*X4
         S7 = X2*X4 - X1*X3

C
C        DO LUNAR TERMS
C
         IF (LSFLG.EQ.1) THEN
            SS1   = S1
            SS2   = S2
            SS3   = S3
            SS4   = S4
            SS5   = S5
            SS6   = S6
            SS7   = S7
            SZ1   = Z1
            SZ2   = Z2
            SZ3   = Z3
            SZ11  = Z11
            SZ12  = Z12
            SZ13  = Z13
            SZ21  = Z21
            SZ22  = Z22
            SZ23  = Z23
            SZ31  = Z31
            SZ32  = Z32
            SZ33  = Z33
            ZCOSG = ZCOSGL
            ZSING = ZSINGL
            ZCOSI = ZCOSIL
            ZSINI = ZSINIL
            ZCOSH = ZCOSHL*CNODM+ZSINHL*SNODM
            ZSINH = SNODM*ZCOSHL-CNODM*ZSINHL
            CC    = C1L
         END IF

      END DO

      ZMOL  = DMOD( 4.7199672D0 + 0.22997150D0*DAY-GAM, TWOPI() )
      ZMOS  = DMOD( 6.2565837D0 + 0.017201977D0*DAY,    TWOPI() )

C
C     DO SOLAR TERMS
C
      SE2 =   2.0D0*SS1*SS6
      SE3 =   2.0D0*SS1*SS7
      SI2 =   2.0D0*SS2*SZ12
      SI3 =   2.0D0*SS2*(SZ13-SZ11)
      SL2 =  -2.0D0*SS3*SZ2
      SL3 =  -2.0D0*SS3*(SZ3-SZ1)
      SL4 =  -2.0D0*SS3*(-21.0D0-9.0D0*EMSQ)*ZES
      SGH2=   2.0D0*SS4*SZ32
      SGH3=   2.0D0*SS4*(SZ33-SZ31)
      SGH4= -18.0D0*SS4*ZES
      SH2 =  -2.0D0*SS2*SZ22
      SH3 =  -2.0D0*SS2*(SZ23-SZ21)

C
C     DO LUNAR TERMS
C
      EE2 =   2.0D0*S1*S6
      E3  =   2.0D0*S1*S7
      XI2 =   2.0D0*S2*Z12
      XI3 =   2.0D0*S2*(Z13-Z11)
      XL2 =  -2.0D0*S3*Z2
      XL3 =  -2.0D0*S3*(Z3-Z1)
      XL4 =  -2.0D0*S3*(-21.0D0-9.0D0*EMSQ)*ZEL
      XGH2=   2.0D0*S4*Z32
      XGH3=   2.0D0*S4*(Z33-Z31)
      XGH4= -18.0D0*S4*ZEL
      XH2 =  -2.0D0*S2*Z22
      XH3 =  -2.0D0*S2*(Z23-Z21)

      CALL CHKOUT ( 'ZZDSCM' )

      RETURN
      END

