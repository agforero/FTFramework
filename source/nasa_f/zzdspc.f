C$Procedure ZZDSPC ( SGP4 deep space routine )

      SUBROUTINE ZZDSPC( IREZ,  D2201, D2211, D3210,   D3222,
     .                   D4410, D4422, D5220, D5232,   D5421, D5433,
     .                   DEDT,  DEL1,  DEL2,  DEL3,    DIDT,  DMDT,
     .                   DNODT, DOMDT, ARGPO, ARGPDOT, T,     TC,
     .                   GSTO,  XFACT, XLAMO, NO,      ATIME, ECCM,
     .                   ARGPM, INCLM, XLI,   MM,      XNI,   NODEM,
     .                   DNDT,  XN  )

C$ Abstract
C
C     This subroutine provides deep space contributions to mean
C     elements for perturbing third body. These effects have been
C     averaged over one revolution of the sun and moon. For earth
C     resonance effects, the effects have been averaged over NO
C     revolutions of the satellite.
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

      INTEGER                  IREZ

      DOUBLE PRECISION         D2201
      DOUBLE PRECISION         D2211
      DOUBLE PRECISION         D3210
      DOUBLE PRECISION         D3222
      DOUBLE PRECISION         D4410
      DOUBLE PRECISION         D4422
      DOUBLE PRECISION         D5220
      DOUBLE PRECISION         D5232
      DOUBLE PRECISION         D5421
      DOUBLE PRECISION         D5433
      DOUBLE PRECISION         DEDT
      DOUBLE PRECISION         DEL1
      DOUBLE PRECISION         DEL2
      DOUBLE PRECISION         DEL3
      DOUBLE PRECISION         DIDT
      DOUBLE PRECISION         DMDT
      DOUBLE PRECISION         DNODT
      DOUBLE PRECISION         DOMDT
      DOUBLE PRECISION         ARGPO
      DOUBLE PRECISION         ARGPDOT
      DOUBLE PRECISION         T
      DOUBLE PRECISION         TC
      DOUBLE PRECISION         GSTO
      DOUBLE PRECISION         XFACT
      DOUBLE PRECISION         XLAMO
      DOUBLE PRECISION         NO
      DOUBLE PRECISION         ATIME
      DOUBLE PRECISION         ECCM
      DOUBLE PRECISION         ARGPM
      DOUBLE PRECISION         INCLM
      DOUBLE PRECISION         XLI
      DOUBLE PRECISION         MM
      DOUBLE PRECISION         XNI
      DOUBLE PRECISION         NODEM
      DOUBLE PRECISION         DNDT
      DOUBLE PRECISION         XN

C$ Brief_I/O
C
C    See Detailed_input and Detailed_Output.
C
C$ Detailed_Input
C
C     D2201      D coeffcients
C
C     D2211         ...
C
C     D3210         ...
C
C     D3222         ...
C
C     D4410         ...
C
C     D4422         ...
C
C     D5220         ...
C
C     D5232         ...
C
C     D5421         ...
C
C     D5433         ...
C
C     DEDT       Internal SGD4 parameter.
C
C     DEL1       Internal SGD4 parameter.
C
C     DEL2       Internal SGD4 parameter.
C
C     DEL3       Internal SGD4 parameter.
C
C     DIDT       Internal SGD4 parameter.
C
C     DMDT       Internal SGD4 parameter.
C
C     DNODT      Internal SGD4 parameter.
C
C     DOMDT      Internal SGD4 parameter.
C
C     IREZ       Flag for resonance: 0-none, 1-one day, 2-half day
C.
C     ARGPO      Argument of perigee
C
C     ARGPDOT    Argument of perigee dot (rate)
C
C     T          Time of evaluation.
C
C     TC         Internal SGD4 parameter.
C
C     GSTO       Grenwich Sidereal Time
C
C     XFACT      Internal SGD4 parameter.
C
C     XLAMO      Internal SGD4 parameter.
C
C     NO         Mean motion.
C
C     ATIME      Internal SGD4 parameter.
C
C     EM         Eccentricity.
C
C     FT         Internal SGD4 parameter.
C
C     ARGPM      Argument of perigee
C
C     INCLM      Inclination.
C
C     XLI        Internal SGD4 parameter.
C
C     MM         Mean anomaly.
C
C     XNI        Mean motion.
C
C     NODEM      Right ascension of ascending node
C
C$ Detailed_Output
C
C     ATIME      Internal SGD4 parameter.
C
C     EM         Calculated mean eccentricity.
C
C     ARGPM      Calculated mean argument of perigee.
C
C     INCLM      Calculated mean inclination.
C
C     XLI        Internal SGD4 parameter.
C
C     MM         Calculated mean anomaly.
C
C     XNI        Calculated  mean motion.
C
C     NODEM      Calculated mean right ascension of 
C                ascending node.
C
C     DNDT       Value XN-NO.
C
C     NM         Calculated mean motion.
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
C     This routine is based on the DSPACE code by David Vallado
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
C-    SPICELIB Version 1.0.0, OCT-14-2014 (EDW)
C
C        Based on routine DPSPACE, 28-JUN-2005, Vallado 2006 [4].
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

      DOUBLE PRECISION         DELT
      DOUBLE PRECISION         FASX2
      DOUBLE PRECISION         FASX4
      DOUBLE PRECISION         FASX6
      DOUBLE PRECISION         FT
      DOUBLE PRECISION         G22
      DOUBLE PRECISION         G32
      DOUBLE PRECISION         G44
      DOUBLE PRECISION         G52
      DOUBLE PRECISION         G54
      DOUBLE PRECISION         RPTIM
      DOUBLE PRECISION         STEP2
      DOUBLE PRECISION         STEPN
      DOUBLE PRECISION         STEPP
      DOUBLE PRECISION         THETA
      DOUBLE PRECISION         X2LI
      DOUBLE PRECISION         X2OMI
      DOUBLE PRECISION         XL
      DOUBLE PRECISION         XLDOT
      DOUBLE PRECISION         XNDDT
      DOUBLE PRECISION         XNDT
      DOUBLE PRECISION         XOMI

      INTEGER                  IRET
      INTEGER                  IRETN


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

      CALL CHKIN ( 'ZZDSPC' )

C
C     Constants
C
      FASX2 = 0.13130908D0
      FASX4 = 2.8843198D0
      FASX6 = 0.37448087D0
      G22   = 5.7686396D0
      G32   = 0.95240898D0
      G44   = 1.8014998D0
      G52   = 1.0508330D0
      G54   = 4.4108898D0
      RPTIM = 4.37526908801129966D-3
      STEPP =    720.0D0
      STEPN =   -720.0D0
      STEP2 = 259200.0D0

C
C     Calculate deep space resonance effects.
C
      DNDT  = 0.0D0
      THETA = DMOD(GSTO + TC*RPTIM,TWOPI())
      ECCM  = ECCM + DEDT*T

      INCLM = INCLM + DIDT*T
      ARGPM = ARGPM + DOMDT*T
      NODEM = NODEM + DNODT*T
      MM    = MM + DMDT*T

C
C   sgp4fix for negative inclinations
C   the following if statement should be commented out
C
C        IF( INCLM .LT. 0.0D0) THEN
C            INCLM  = -INCLM
C            ARGPM  = ARGPM-PI
C            NODEM  = NODEM+PI
C        END IF
C

C
C     sgp4fix for propagator problems
C
C     The following integration works for negative time steps and
C     periods. The specific changes are unknown because the original
C     code was so convoluted
C
C     sgp4fix Take out atime = 0.0 and fix for faster operation
C

C
C     Just in case - should be set in loops if used.
C
      FT    = 0.0D0

      IF (IREZ .ne. 0) THEN

C
C     UPDATE RESONANCES : NUMERICAL (EULER-MACLAURIN) INTEGRATION
C
C     EPOCH RESTART
C

C
C        sgp4fix streamline check
C
         IF ( (    ATIME .EQ. 0.0D0) .OR.
     .        (T * ATIME .LE. 0.0D0) .OR.
     .        (  DABS(T) .LT. DABS(ATIME))  ) THEN
            ATIME  = 0.0D0
            XNI    = NO
            XLI    = XLAMO
         END IF

C
C        sgp4fix move check outside loop
C
         IF (T .GT. 0.0D0) THEN
            DELT = STEPP
         ELSE
            DELT = STEPN
         END IF

C
C        ADDED FOR DO LOOP
C
         IRETN = 381

C
C        ADDED FOR LOOP
C
         IRET  =   0

         DO WHILE (IRETN.EQ.381)

C
C           DOT TERMS CALCULATED
C
C           NEAR - SYNCHRONOUS RESONANCE TERMS
C
            IF (IREZ .NE. 2) THEN

               XNDT  = DEL1*DSIN(XLI-FASX2) +
     .                  DEL2*DSIN(2.0D0*(XLI-FASX4)) +
     .                  DEL3*DSIN(3.0D0*(XLI-FASX6))
               XLDOT = XNI + XFACT
               XNDDT = DEL1*DCOS(XLI-FASX2) +
     .            2.0D0*DEL2*DCOS(2.0D0*(XLI-FASX4)) +
     .            3.0D0*DEL3*DCOS(3.0D0*(XLI-FASX6))
               XNDDT = XNDDT*XLDOT

            ELSE

C
C              NEAR - HALF-DAY RESONANCE TERMS
C

               XOMI = ARGPO + ARGPDOT*ATIME
               X2OMI= XOMI + XOMI
               X2LI = XLI + XLI
               XNDT = D2201*DSIN(X2OMI+XLI-G22) +
     .                 D2211*DSIN(XLI-G22) +
     .                 D3210*DSIN( XOMI+XLI-G32) +
     .                 D3222*DSIN(-XOMI+XLI-G32) +
     .                 D4410*DSIN(X2OMI+X2LI-G44)+
     .                 D4422*DSIN(X2LI-G44)+
     .                 D5220*DSIN( XOMI+XLI-G52) +
     .                 D5232*DSIN(-XOMI+XLI-G52) +
     .                 D5421*DSIN( XOMI+X2LI-G54)+
     .                 D5433*DSIN(-XOMI+X2LI-G54)
               XLDOT = XNI+XFACT
               XNDDT = D2201*DCOS(X2OMI+XLI-G22) +
     .                  D2211*DCOS(XLI-G22)+
     .                  D3210*DCOS( XOMI+XLI-G32) +
     .                  D3222*DCOS(-XOMI+XLI-G32) +
     .                  D5220*DCOS( XOMI+XLI-G52) +
     .                  D5232*DCOS(-XOMI+XLI-G52) +
     .                  2.0D0*(D4410*DCOS(X2OMI+X2LI-G44) +
     .                  D4422*DCOS(X2LI-G44) +
     .                  D5421*DCOS( XOMI+X2LI-G54) +
     .                  D5433*DCOS(-XOMI+X2LI-G54))
               XNDDT = XNDDT*XLDOT

            END IF

C
C           INTEGRATOR
C
C           sgp4fix move end checks to end of routine
C

            IF (DABS(T-ATIME).GE.STEPP) THEN
               IRET  = 0
               IRETN = 381
            ELSE
               FT    = T-ATIME
               IRETN = 0
            END IF

            IF (IRETN.EQ.381) THEN
               XLI   = XLI + XLDOT*DELT + XNDT*STEP2
               XNI   = XNI + XNDT*DELT + XNDDT*STEP2
               ATIME = ATIME + DELT
            END IF

         END DO

         XN = XNI + XNDT*FT  + XNDDT*FT*FT*0.5D0
         XL = XLI + XLDOT*FT + XNDT*FT*FT*0.5D0

         IF(IREZ .NE. 1) THEN
            MM   = XL-2.0D0*NODEM+2.0D0*THETA
            DNDT = XN-NO
         ELSE
            MM   = XL-NODEM-ARGPM+THETA
            DNDT = XN-NO
         END IF

         XN = NO + DNDT

      END IF

      CALL CHKOUT ( 'ZZDSPC' )

      RETURN
      END

