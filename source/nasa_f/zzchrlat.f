C$Procedure ZZCHRLAT ( Chord latitude  )
 
      SUBROUTINE ZZCHRLAT ( MIDLAT, DLON, EPTLAT )

C$ Abstract
C
C     Given the latitude of a midpoint of chord on a circle of constant
C     latitude, and given the longitude extent of the chord, compute
C     the latitude of the chord's endpoints. The coordinate system is
C     "latitudinal," aka planetocentric. 
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
C     COORDINATES
C     LATITUDE
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsktol.inc'

      DOUBLE PRECISION      MIDLAT
      DOUBLE PRECISION      DLON
      DOUBLE PRECISION      EPTLAT

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     MIDLAT     I   is the latitude of the midpoint of a chord of a 
C                    latitude circle.
C     DLON       I   is the longitude extent of the chord.
C     EPTLAT     O   is the latitude of the endpoints of the chord.
C
C$ Detailed_Input
C
C     MIDLAT     is the latitude of the midpoint of a chord of a circle
C                of constant latitude. Units are radians. The range of
C                MIDLAT is -pi : pi.
C
C     DLON       is the extent in longitude of the chord. DLON is the
C                difference of the longitudes of the endpoints of the
C                chords, expressed as a non-negative value. Units are
C                radians. DLON must be strictly less than pi.
C
C$ Detailed_Output
C
C     EPTLAT     is the latitude of the circle to which the chord
C                is tangent. Units are radians.
C                
C$ Parameters
C
C     ANGMRG is a tolerance used to determine whether the input latitude
C     is within range. See 'dsktol.inc' for further information.
C
C$ Exceptions
C 
C     1)  If the input longitude extent DLON is negative, or if the
C         extent is greater than or equal to pi radians, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C     2)  If the input latitude is outside the range 
C
C             -pi/2 - ANGMRG  :  pi/2 + ANGMRG
C
C         the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     In the remarks below, the coordinate system is presumed to be
C     latitudinal.
C
C     A "chord" of a circle is a line segment having endpoints on the
C     circle.
C
C     This routine supports partitioning of a tessellated spheroid
C     surface into bands of plates covering specified latitude ranges.
C     Note that, for a plate edge that is a chord of a positive
C     latitude circle, the maximum latitude on the chord is attained at
C     the chord's midpoint. So a set of plates lying on or above the
C     plate containing the latitude circle cannot actually cover the
C     latitude range above the circle: rays emanating from the origin
C     can pass above the latitude circle and miss the set of plates. An
C     analogous situation exists for circles of negative latitude.
C
C     In order to generate a set of plates that cover a positive
C     latitude band, the plates can be arranged so that the maximum
C     latitude of the lowest plate edges is less than or equal to the
C     lower latitude bound of the band. An analogous solution can be
C     derived for bands of negative latitude. This routine can be used
C     to solve these problems: it can compute the latitude of a circle
C     ---where the circle is horizontal and centered on the z axis---
C     in which a regular polygon must be inscribed so that the maximum
C     latitude on the polygon is less than a specified value, or such
C     that the minimum latitude on the polygon is greater than a
C     specified value.
C     
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C
C
C     1)  Find the latitude circle on which a line segment's endpoints
C         must lie, if the longitude extent of the segment is 60
C         degrees, and if the latitude of the segment's midpoint is 30
C         degrees.
C
C
C     Example code begins here. 
C
C
C           PROGRAM EX1
C           IMPLICIT NONE
C
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      DPR
C           DOUBLE PRECISION      RPD
C
C     C
C     C     Local variables
C     C
C           DOUBLE PRECISION      DLON
C           DOUBLE PRECISION      EP     ( 3, 2 )
C           DOUBLE PRECISION      EPTLAT
C           DOUBLE PRECISION      MIDLAT
C           DOUBLE PRECISION      MIDLON
C           DOUBLE PRECISION      MIDRAD
C           DOUBLE PRECISION      MIDPT  ( 3 )
C           DOUBLE PRECISION      OUTLAT
C
C     C
C     C     Let the segment have a longitude extent
C     C     of 60 degrees.
C     C
C           DLON   = 60.D0 * RPD()
C
C     C
C     C     Let the midpoint of the segment have
C     C     a latitude of 30 degrees.
C     C
C           MIDLAT = 30.D0 * RPD()
C
C     C
C     C     Find the latitude of the segment's endpoints.
C     C
C           CALL ZZCHRLAT ( MIDLAT, DLON, EPTLAT )
C
C           WRITE (*,*) 'Endpoint latitude (deg) = ', EPTLAT*DPR()
C
C     C
C     C     Generate the endpoints and the latitude
C     C     of the segment's midpoint.
C     C
C     C     Note that the scale is arbitrary: we can set it
C     C     to 1.0. The longitudes of the endpoints are
C     C     arbitrary as well; only the difference is known.
C     C
C           CALL LATREC ( 1.D0, 0.D0, EPTLAT, EP(1,1) )
C           CALL LATREC ( 1.D0, DLON, EPTLAT, EP(1,2) )
C
C           CALL VLCOM  ( 0.5D0, EP(1,1),
C          .              0.5D0, EP(1,2), MIDPT )
C
C           CALL RECLAT ( MIDPT, MIDRAD, MIDLON, OUTLAT )
C
C           WRITE (*,*) 'Check: difference (deg) = ',
C          .            (MIDLAT - OUTLAT) * DPR()
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit
C     platform, the output was:
C
C
C        Endpoint latitude (deg) =    26.565051177077986
C        Check: difference (deg) =   6.36110936292703354E-015
C
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 20-APR-2016 (NJB) (EDW)
C
C-&
 
C$ Index_Entries
C
C     find latitude of chord's endpoints on latitude circle
C
C-&
 

C
C     SPICELIB functions
C
      DOUBLE PRECISION      BRCKTD
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      PI

      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      MLAT

C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      IF ( ( DLON .LT. 0.D0 ) .OR. ( DLON .GE. PI() ) )  THEN
         
         CALL CHKIN  ( 'ZZCHRLAT'                          )
         CALL SETMSG ( 'The input longitude extent was #; '
     .   //            'this value must be in the range '
     .   //            '[0 : pi ) radians.'                )
         CALL ERRDP  ( '#', DLON                           )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'            )
         CALL CHKOUT ( 'ZZCHRLAT'                          )
         RETURN

      END IF

      IF (  ABS(MIDLAT)  .GT.  ( HALFPI()+ANGMRG )  ) THEN

         CALL CHKIN  ( 'ZZCHRLAT'                          )
         CALL SETMSG ( 'The input latitude was #; this '
     .   //            'value must be in the interval '
     .   //            '-pi/2 : pi/2 (radians).'           )
         CALL ERRDP  ( '#', MIDLAT                         )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'            )
         CALL CHKOUT ( 'ZZCHRLAT'                          )
         RETURN

      END IF

C
C     The input latitude is, at worst, slightly out of range.
C     Bracket it.
C
      MLAT = BRCKTD( MIDLAT, -HALFPI(), HALFPI() )

C
C     The endpoint latitude EPTLAT is defined by 
C                   
C        EPTLAT = atan ( tan(MLAT) * cos( DLON/2 ) )
C
C     For numerical robustness, we'll re-write this using
C     the two-argument arctangent function and well-behaved 
C     trig functions as input arguments:
C     
      EPTLAT = ATAN2 ( SIN(MLAT) * COS(DLON/2),  COS(MLAT) )

      RETURN 
      END
