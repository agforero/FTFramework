C$Procedure ZZGRAV ( SGP4 gravitational constants )

      SUBROUTINE ZZGRAV ( GRAV )

C$ Abstract
C
C      Constants assignments equivalent to 21 July 2006 "getgravconst"
C      by David Vallado, AIAA-2006-6753.
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

      INCLUDE               'zzsgp4.inc'

      DOUBLE PRECISION      GRAV( NGRAVS, NGRAVC )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     GRAV       O   Array of geophysical constants.
C
C$ Detailed_Input
C
C     None.
C
C$ Detailed_Output
C
C     GRAV       An array of three geophysical constants vectors.
C
C$ Parameters
C
C     Refer to include file zzsgp4.inc.
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
C     None.
C
C$ Examples
C
C           INCLUDE               'zzsgp4.inc'
C
C           CALL ZZGRAV ( GRAV )
C     
C     C
C     C     Retrieve the constants set based on the WGS ID parameter.
C     C
C
C           MU   = GRAV( P_WGS, P_MU)
C           RAD  = GRAV( P_WGS, P_RAD)
C           XKE  = GRAV( P_WGS, P_XKE)
C           TUMN = GRAV( P_WGS, P_TUMN)
C           J2   = GRAV( P_WGS, P_J2)
C           J3   = GRAV( P_WGS, P_J3)
C           J4   = GRAV( P_WGS, P_J4)
C           J3OJ2= GRAV( P_WGS, P_J3J2)
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
C-    SPICELIB Version 1.0.0 22-JUL-2014 (EDW)
C
C-&

C$ Index_Entries
C
C   SGP4
C
C-&

C
C     The definition of the MU constant is to allow calculation of
C     those constants which are functions of MU in the WGS72 and WGS84
C     assignments.
C


C
C     WGS72 low precision Spacetrack #3 constants. These values should
C     correspond to those in geophysical.ker. If not, something is
C     very wrong.
C

      GRAV( WGS721, P_RAD)  =  6378.135D0

      GRAV( WGS721, P_XKE)  =  0.0743669161D0

      GRAV( WGS721, P_MU)   =  398600.79964D0

      GRAV( WGS721, P_TUMN) =  1.0D0 / GRAV(WGS721, P_XKE)

      GRAV( WGS721, P_J2)   =  0.001082616D0

      GRAV( WGS721, P_J3)   =  -0.00000253881D0

      GRAV( WGS721, P_J4)   =  -0.00000165597D0

      GRAV( WGS721, P_J3J2) =  GRAV(WGS721, P_J3) / GRAV(WGS721, P_J2)


C
C     WGS 72 constants.
C

      GRAV( WGS72, P_MU)   = 398600.8D0

      GRAV( WGS72, P_RAD)  =  6378.135D0

      GRAV( WGS72, P_XKE)  =  60.0D0 /
     .                 DSQRT(GRAV(WGS72, P_RAD)**3/GRAV(WGS72, P_MU))

      GRAV( WGS72, P_TUMN) =  1.0D0 / GRAV(WGS72, P_XKE)

      GRAV( WGS72, P_J2)   =  0.001082616D0

      GRAV( WGS72, P_J3)   =  -0.00000253881D0

      GRAV( WGS72, P_J4)   =  -0.00000165597D0

      GRAV( WGS72, P_J3J2) =  GRAV(WGS72, P_J3) / GRAV(WGS72, P_J2)


C
C     WGS 84 constants.
C

      GRAV( WGS84, P_MU)   =  398600.5D0

      GRAV( WGS84, P_RAD)  =  6378.137D0

      GRAV( WGS84, P_XKE)  =  60.0D0 /
     .                 DSQRT(GRAV(WGS84, P_RAD)**3/GRAV(WGS84, P_MU))

      GRAV( WGS84, P_TUMN) =  1.0D0 / GRAV(WGS84, P_XKE)

      GRAV( WGS84, P_J2)   =  0.00108262998905D0

      GRAV( WGS84, P_J3)   =  -0.00000253215306D0

      GRAV( WGS84, P_J4)   =  -0.00000161098761D0

      GRAV( WGS84, P_J3J2) =  GRAV(WGS84, P_J3) / GRAV(WGS84, P_J2)

      END
      
