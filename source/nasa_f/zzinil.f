C$Procedure ZZINIL ( SGP4 initializer )

      SUBROUTINE ZZINIL ( GEOPHS, OPMODE, ECCO,   
     .                    EPOCH,  INCLO,  NO,     AINV,
     .                    AO,     CON41,  CON42,  COSIO,
     .                    COSIO2, ECCSQ,  OMEOSQ, POSQ,
     .                    RP,     RTEOSQ, SINIO,  GSTO )

C$ Abstract
C
C     This subroutine initializes the SGP4 propagator. All the
C     initialization is consolidated here instead of having multiple
C     loops inside other routines.
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

      DOUBLE PRECISION         GEOPHS    ( 8 )
      INTEGER                  OPMODE

      DOUBLE PRECISION         ECCO
      DOUBLE PRECISION         EPOCH
      DOUBLE PRECISION         INCLO
      DOUBLE PRECISION         NO
      DOUBLE PRECISION         AINV
      DOUBLE PRECISION         AO
      DOUBLE PRECISION         CON41
      DOUBLE PRECISION         CON42
      DOUBLE PRECISION         COSIO
      DOUBLE PRECISION         COSIO2
      DOUBLE PRECISION         ECCSQ
      DOUBLE PRECISION         OMEOSQ
      DOUBLE PRECISION         POSQ
      DOUBLE PRECISION         RP
      DOUBLE PRECISION         RTEOSQ
      DOUBLE PRECISION         SINIO
      DOUBLE PRECISION         GSTO

C$ Brief_I/O
C
C    See Detailed_input and Detailed_Output.
C
C$ Detailed_Input
C
C     GEOPYS         Gephysical constants array.
C
C     ECCO           Eccentricity.
C
C     EPOCH          TLE epoch time in days from
C                    1950-00-00 00:00:00.000000 TDB, i.e.
C                    1949-12-31 00:00:00.000000 TDB.
C
C     INCLO          Inclination
C
C     OPMODE         Flag indicating which technique
C                    to use to calculate sidereal time.
C
C     NO             Mean motion.
C
C$ Detailed_Output
C
C     NO             Mean motion.
C
C     AINV           1.0/A0
C
C     AO             Semi major axis.
C
C     CON41          Value -CON42-COSIO2-COSIO2.
C
C     CON42          1.0 - 5.0*cos(inclination).
C
C     COSIO          Cosine of inclination.
C
C     COSIO2         COSIO squared.
C
C     ECCSQ          Eccentricity squared.
C
C     OMEOSQ         1.0 - ECCO * ECCO.
C
C     POSQ           Semi-parameter squared.
C
C     RP             Radius of perigee.
C
C     RTEOSQ         Square root of (1.0 - ECCO*ECCO).
C
C     SINIO          Sine of inclination.
C
C     GSTO           GST at time of observation.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) SPICE(UNKNOWNMODE) signals when the value of OPMODE
C        does not equal an assigned OPMODE value listed in
C        zzsgp4.inc.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is based on the INITL code by David Vallado
C     corresponding to "Revisiting Spacetrack Report #3" [4].
C     The intent is to maintain the original Vallado algorithm,
C     changing code only to meet NAIF format standards and to
C     integrate with SPICELIB.
C
C        Removed getgravconst call, replaced with GEOPHS array.
C
C        Capitalize all variables.
C
C        ENDIF replaced with END IF.
C
C        RadPerDay     replaced with RADDAY
C
C        Opsmode       replaced with OPMODE
C
C        whichconst    eliminated, function provided by GEOPHS
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
C-    SPICELIB Version 1.0.0, OCT-09-2014 (EDW)
C
C        Based on routine INITL, 28-JUN-2005, Vallado 2006 [4].
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

      DOUBLE PRECISION         ADEL
      DOUBLE PRECISION         AK
      DOUBLE PRECISION         C1
      DOUBLE PRECISION         C1P2P
      DOUBLE PRECISION         D1
      DOUBLE PRECISION         DEL
      DOUBLE PRECISION         FK5R
      DOUBLE PRECISION         J2
      DOUBLE PRECISION         PO
      DOUBLE PRECISION         RADDAY
      DOUBLE PRECISION         TEMP
      DOUBLE PRECISION         TFRAC
      DOUBLE PRECISION         THGR70
      DOUBLE PRECISION         TS70
      DOUBLE PRECISION         TUT1
      DOUBLE PRECISION         X2O3
      DOUBLE PRECISION         XKE

      INTEGER                  IDS70


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

      CALL CHKIN ( 'ZZINIL' )


C
C     This code block replaces the call:
C
C     sgp4fix identify constants and allow alternate values.
C
C     CALL getgravconst( whichconst, tumin,
C     .                  mu, radiusearthkm, xke,
C     .                  j2, j3, j4, j3oj2 )
C
      J2    = GEOPHS(K_J2)
      XKE   = GEOPHS(K_KE)
      X2O3  = 2.0D0/3.0D0

C
C     Calculate auxillary epoch quantities
C
      ECCSQ  = ECCO*ECCO
      OMEOSQ = 1.0D0 - ECCSQ
      RTEOSQ = DSQRT(OMEOSQ)
      COSIO  = DCOS(INCLO)
      COSIO2 = COSIO*COSIO

C
C     Un-KOZAI the mean motion
C
      AK   =  (XKE/NO)**X2O3
      D1   =  0.75D0*J2* (3.0D0*COSIO2-1.0D0) / (RTEOSQ*OMEOSQ)
      DEL  =  D1/(AK*AK)
      ADEL =  AK * ( 1.0D0 - DEL*DEL - DEL*
     .                 (1.0D0/3.0D0 + 134.0D0*DEL*DEL / 81.0D0) )
      DEL  =  D1/(ADEL*ADEL)
      NO   =  NO/(1.0D0 + DEL)

      AO   =  (XKE/NO)**X2O3
      SINIO=  DSIN(INCLO)
      PO   =  AO*OMEOSQ
      CON42=  1.0D0-5.0D0*COSIO2
      CON41=  -CON42-COSIO2-COSIO2
      AINV =  1.0D0/AO
      POSQ =  PO*PO
      RP   =  AO*(1.0D0-ECCO)

C
C     Calculate greenwich location at epoch
C

C
C     sgp4fix Modern approach to finding sidereal time
C

      IF (OPMODE .EQ. IMPRVD ) THEN

C
C        Radians per day, earth rotation, 6.30038809866574D0.
C
         RADDAY = TWOPI() * 1.002737909350795D0
         TEMP   = EPOCH + 2433281.5D0
         TUT1= ( DINT(TEMP-0.5D0) + 0.5D0 - 2451545.0D0 ) / 36525.0D0
         GSTO= 1.75336855923327D0 + 628.331970688841D0*TUT1
     .             + 6.77071394490334D-06*TUT1*TUT1
     .             - 4.50876723431868D-10*TUT1*TUT1*TUT1
     .             + RADDAY*( TEMP-0.5D0-DINT(TEMP-0.5D0) )

      ELSE IF (OPMODE .EQ. AFSPC ) THEN

C
C        sgp4fix Use old way of finding GST
C
C        Count integer number of days from 0 jan 1970
C
         TS70  = EPOCH-7305.0D0
         IDS70 = INT(TS70 + 1.0D-8)
         TFRAC = TS70-IDS70

C
C        Find greenwich location at epoch
C
         C1     = 1.72027916940703639D-2
         THGR70 = 1.7321343856509374D0
         FK5R   = 5.07551419432269442D-15
         C1P2P  = C1+TWOPI()
         GSTO   = THGR70+C1*IDS70+C1P2P*TFRAC+TS70*TS70*FK5R

      ELSE

         CALL SETMSG ( 'Unknown value for OPMODE. Value # not '
     .    //           'coded in zzsgp4.inc.'                         )
         CALL ERRINT ( '#', OPMODE                                    )
         CALL SIGERR ( 'SPICE(UNKNOWNMODE)'                           )
         CALL CHKOUT ( 'ZZINIL'                                       )
         RETURN

      END IF

C
C     Check quadrants
C
      GSTO = DMOD( GSTO,TWOPI() )

      IF ( GSTO .LT. 0.0D0 ) THEN
         GSTO= GSTO + TWOPI()
      END IF

      CALL CHKOUT ( 'ZZINIL' )

      RETURN
      END
