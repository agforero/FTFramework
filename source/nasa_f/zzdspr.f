C$Procedure ZZDSPR ( SGP4 deep space long period )

      SUBROUTINE ZZDSPR (  OPMODE, 
     .                     E3,   EE2,   PEO,   PGHO,  PHO,  PINCO,
     .                     PLO,  SE2,   SE3,   SGH2,  SGH3, SGH4,
     .                     SH2,  SH3,   SI2,   SI3,   SL2,  SL3,
     .                     SL4,  T,     XGH2,  XGH3,  XGH4, XH2,
     .                     XH3,  XI2,   XI3,   XL2,   XL3,  XL4,
     .                     ZMOL, ZMOS,  INCLO, DOINIT,
     .                     ECCP, INCLP, NODEP, ARGPP, MP)

C$ Abstract
C
C     This subroutine provides deep space long period periodic
C     contributions to the mean elements. By design, these periodics
C     are zero at epoch. This used to be dscom which included
C     initialization, but it's really a recurring function.
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

      INCLUDE 'zzsgp4.inc'

      INTEGER               OPMODE
      DOUBLE PRECISION      E3
      DOUBLE PRECISION      EE2
      DOUBLE PRECISION      PEO
      DOUBLE PRECISION      PGHO
      DOUBLE PRECISION      PHO
      DOUBLE PRECISION      PINCO
      DOUBLE PRECISION      PLO
      DOUBLE PRECISION      SE2
      DOUBLE PRECISION      SE3
      DOUBLE PRECISION      SGH2
      DOUBLE PRECISION      SGH3
      DOUBLE PRECISION      SGH4
      DOUBLE PRECISION      SH2
      DOUBLE PRECISION      SH3
      DOUBLE PRECISION      SI2
      DOUBLE PRECISION      SI3
      DOUBLE PRECISION      SL2
      DOUBLE PRECISION      SL3
      DOUBLE PRECISION      SL4
      DOUBLE PRECISION      T
      DOUBLE PRECISION      XGH2
      DOUBLE PRECISION      XGH3
      DOUBLE PRECISION      XGH4
      DOUBLE PRECISION      XH2
      DOUBLE PRECISION      XH3
      DOUBLE PRECISION      XI2
      DOUBLE PRECISION      XI3
      DOUBLE PRECISION      XL2
      DOUBLE PRECISION      XL3
      DOUBLE PRECISION      XL4
      DOUBLE PRECISION      ZMOL
      DOUBLE PRECISION      ZMOS
      DOUBLE PRECISION      INCLO
      LOGICAL               DOINIT
      DOUBLE PRECISION      ECCP
      DOUBLE PRECISION      INCLP
      DOUBLE PRECISION      NODEP
      DOUBLE PRECISION      ARGPP
      DOUBLE PRECISION      MP

C$ Brief_I/O
C
C
C    See Detailed_input and Detailed_Output.
C
C$ Detailed_Input
C
C     OPMODE     Flag indicating which technique
C                to use to calculate sidereal time.
C 
C     E3         Internal SGD4 parameter.
C
C     EE2        Internal SGD4 parameter.
C
C     PEO        Internal SGD4 parameter.
C
C     PGHO       Internal SGD4 parameter.
C
C     PHO        Internal SGD4 parameter.
C
C     PINCO      Internal SGD4 parameter.
C
C     PLO        Internal SGD4 parameter.
C
C     SE2        Internal SGD4 parameter.
C
C     SE3        Internal SGD4 parameter.
C
C     SGH2       Internal SGD4 parameter.
C
C     SGH3       Internal SGD4 parameter.
C
C     SGH4       Internal SGD4 parameter.
C
C     SH2        Internal SGD4 parameter.
C
C     SH3        Internal SGD4 parameter.
C
C     SI2        Internal SGD4 parameter.
C
C     SI3        Internal SGD4 parameter.
C
C     SL2        Internal SGD4 parameter.
C
C     SL3        Internal SGD4 parameter.
C
C     SL4        Internal SGD4 parameter.
C
C     T          Time for state evaluation.
C
C     XGH2       Internal SGD4 parameter.
C
C     XGH3       Internal SGD4 parameter.
C
C     XGH4       Internal SGD4 parameter.
C
C     XH2        Internal SGD4 parameter.
C
C     XH3        Internal SGD4 parameter.
C
C     XI2        Internal SGD4 parameter.
C
C     XI3        Internal SGD4 parameter.
C
C     XL2        Internal SGD4 parameter.
C
C     XL3        Internal SGD4 parameter.
C
C     XL4        Internal SGD4 parameter.
C
C     ZMOL       Internal SGD4 parameter.
C
C     ZMOS       Internal SGD4 parameter.
C
C     INCLO      Unused argument. Maintained in call for historical
C                referene.
C
C     DOINIT     Flag indicating initialization state. True to
C                initialize, false otherwise.
C
C     ECCP       Eccentricity.
C
C     INCLP      Inclination.
C
C     NODEP      Right ascension of  ascending node.
C
C     ARGPP      Argument of periapsis.
C
C     MP         Mean anomoly.
C
C$ Detailed_Output
C
C     ECCP       Calculated perturbed eccentricity.
C
C     INCLP      Calculated perturbed inclination.
C
C     NODEP      Calculated perturbed right ascension of 
C                ascending node.
C
C     ARGPP      Calculated perturbed argument of periapsis.
C
C     MP         Calculated perturbed mean anomoly.
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
C     This routine is based on the DPPER code by David Vallado
C     corresponding to "Revisiting Spacetrack Report #3" [4].
C     The intent is to maintain the original Vallado algorithm,
C     changing code only to meet NAIF format standards and to
C     integrate with SPICELIB.
C
C        Capitalize all variables.
C
C        ENDIF replaced with END IF.
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
C-    SPICELIB Version 1.0.0, SEP-15-2014 (EDW)
C
C        Based on routine DDPER, 28-JUN-2005, Vallado 2006 [4].
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
      DOUBLE PRECISION      ALFDP
      DOUBLE PRECISION      BETDP
      DOUBLE PRECISION      COSIP
      DOUBLE PRECISION      COSOP
      DOUBLE PRECISION      DALF
      DOUBLE PRECISION      DBET
      DOUBLE PRECISION      DLS
      DOUBLE PRECISION      F2
      DOUBLE PRECISION      F3
      DOUBLE PRECISION      PE
      DOUBLE PRECISION      PGH
      DOUBLE PRECISION      PH
      DOUBLE PRECISION      PINC
      DOUBLE PRECISION      PL
      DOUBLE PRECISION      SEL
      DOUBLE PRECISION      SES
      DOUBLE PRECISION      SGHL
      DOUBLE PRECISION      SGHS
      DOUBLE PRECISION      SHL
      DOUBLE PRECISION      SHS
      DOUBLE PRECISION      SIL
      DOUBLE PRECISION      SINIP
      DOUBLE PRECISION      SINOP
      DOUBLE PRECISION      SINZF
      DOUBLE PRECISION      SIS
      DOUBLE PRECISION      SLL
      DOUBLE PRECISION      SLS
      DOUBLE PRECISION      XLS
      DOUBLE PRECISION      XNOH
      DOUBLE PRECISION      ZEL
      DOUBLE PRECISION      ZES
      DOUBLE PRECISION      ZF
      DOUBLE PRECISION      ZM
      DOUBLE PRECISION      ZNL
      DOUBLE PRECISION      ZNS

C
C     SPICELIB routines.
C
      DOUBLE PRECISION         TWOPI
      DOUBLE PRECISION         PI
      LOGICAL                  RETURN

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSPR' )

C
C     Local constants.
C
      ZES  = 0.01675D0
      ZEL  = 0.05490D0
      ZNS  = 1.19459D-5
      ZNL  = 1.5835218D-4

C
C     Calculate time varying periodics.
C
      ZM   = ZMOS + ZNS*T

      IF (DOINIT) ZM = ZMOS

      ZF   = ZM + 2.0D0*ZES*DSIN(ZM)
      SINZF= DSIN(ZF)
      F2   =  0.5D0*SINZF*SINZF - 0.25D0
      F3   = -0.5D0*SINZF*DCOS(ZF)
      SES  = SE2*F2 + SE3*F3
      SIS  = SI2*F2 + SI3*F3
      SLS  = SL2*F2 + SL3*F3 + SL4*SINZF
      SGHS = SGH2*F2 + SGH3*F3 + SGH4*SINZF
      SHS  = SH2*F2 + SH3*F3
      ZM   = ZMOL + ZNL*T

      IF (DOINIT) ZM = ZMOL

      ZF   = ZM + 2.0D0*ZEL*DSIN(ZM)
      SINZF= DSIN(ZF)
      F2   =  0.5D0*SINZF*SINZF - 0.25D0
      F3   = -0.5D0*SINZF*DCOS(ZF)
      SEL  = EE2*F2 + E3*F3
      SIL  = XI2*F2 + XI3*F3
      SLL  = XL2*F2 + XL3*F3 + XL4*SINZF
      SGHL = XGH2*F2 + XGH3*F3 + XGH4*SINZF
      SHL  = XH2*F2 + XH3*F3
      PE   = SES + SEL
      PINC = SIS + SIL
      PL   = SLS + SLL
      PGH  = SGHS + SGHL
      PH   = SHS + SHL

      IF ( .NOT. DOINIT) THEN

         PE    = PE   - PEO
         PINC  = PINC - PINCO
         PL    = PL   - PLO
         PGH   = PGH  - PGHO
         PH    = PH   - PHO
         INCLP = INCLP  + PINC
         ECCP  = ECCP   + PE
         SINIP = DSIN(INCLP)
         COSIP = DCOS(INCLP)

C
C        Apply periodics directly.
C
C        sgp4fix for lyddane choice
C
C        strn3 used original inclination - this is technically
C        feasible
C
C        gsfc used perturbed inclination - also technically feasible
C        probably best to readjust the 0.2 limit value and limit
C        discontinuity
C
C        0.2 rad = 11.45916 deg
C
C        use next line for original strn3 approach and original
C        inclination
C
C            IF (INCLO.GE.0.2D0) THEN
C
C        use next line for gsfc version and perturbed inclination
C
         IF (INCLP.GE.0.2D0) THEN

            PH     = PH/SINIP
            PGH    = PGH - COSIP*PH
            ARGPP  = ARGPP + PGH
            NODEP  = NODEP + PH
            MP     = MP + PL

         ELSE

C
C           Apply periodics with Lyddane modification.
C
            SINOP  = DSIN(NODEP)
            COSOP  = DCOS(NODEP)
            ALFDP  = SINIP*SINOP
            BETDP  = SINIP*COSOP
            DALF   =  PH*COSOP + PINC*COSIP*SINOP
            DBET   = -PH*SINOP + PINC*COSIP*COSOP
            ALFDP  = ALFDP + DALF
            BETDP  = BETDP + DBET
            NODEP  = DMOD(NODEP,TWOPI())

C
C           sgp4fix for afspc written intrinsic functions
C           NODEP used without a trigonometric function ahead
C

            IF ( (NODEP  .LT. 0.0D0)      .AND.
     .           (OPMODE .EQ. AFSPC)  )   THEN
               NODEP = NODEP + TWOPI()
            END IF

            XLS    = MP + ARGPP + COSIP*NODEP
            DLS    = PL + PGH - PINC*NODEP*SINIP
            XLS    = XLS + DLS
            XNOH   = NODEP
            NODEP  = DATAN2(ALFDP,BETDP)

C
C           sgp4fix for afspc written intrinsic functions
C           NODEP used without a trigonometric function ahead
C

            IF ( (NODEP .LT. 0.0D0)  .AND.
     .           (OPMODE .EQ. AFSPC)  )  THEN
               NODEP = NODEP + TWOPI()
            END IF

            IF (DABS(XNOH-NODEP) .GT. PI()) THEN

               IF (NODEP .LT. XNOH) THEN
                  NODEP = NODEP+TWOPI()
               ELSE
                  NODEP = NODEP-TWOPI()
               END IF

            END IF

            MP   = MP + PL
            ARGPP=  XLS - MP - COSIP*NODEP

         END IF

      END IF

      CALL CHKOUT ( 'ZZDSPR' )

      RETURN
      END


