C$Procedure ZZINVELT ( DSK, in volume element? )

      SUBROUTINE ZZINVELT ( P,      CORSYS, CORPAR, BOUNDS, 
     .                      MARGIN, EXCLUD, INSIDE         )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Test a point represented by a set of Cartesian coordinates for
C     inclusion in a volume element in a specified coordinate system.
C     The volume element is bounded by surfaces on which one coordinate
C     is constant. The test is performed using margins for the element.
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
C     DSK
C
C$ Keywords
C     
C     DSK
C     GEOMETRY
C     INTERSECTION
C     SURFACE
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsktol.inc'
      INCLUDE 'dskdsc.inc'


      DOUBLE PRECISION      P      ( 3 )
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      MARGIN
      INTEGER               EXCLUD
      LOGICAL               INSIDE
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     P          I   Input point.
C     CORSYS     I   Coordinate system code.
C     CORPAR     I   Coordinate system parameters.
C     BOUNDS     I   Coordinate bounds of element.
C     MARGIN     I   Margin used for inclusion testing.
C     EXCLUD     I   Index of coordinate to exclude from test.
C     INSIDE     O   Flag indicating whether point is in element.
C
C$ Detailed_Input
C
C     P          is a point expressed in Cartesian coordinates. The
C                point is to be checked to determine whether it is
C                inside the volume element specified by BOUNDS.
C
C     CORSYS     is an integer parameter identifying the coordinate
C                system in which the bounds are to be computed. See the
C                include file dskdsc.inc for allowed values of CORSYS.
C
C     CORPAR     is an array of parameters associated with the
C                coordinate system. Currently the only supported system
C                that has associated parameters is the planetodetic
C                system. For planetodetic coordinates,
C
C                   CORPAR(1) is the equatorial radius
C
C                   CORPAR(2) is the flattening coefficient. Let RE and
C                   RP represent, respectively, the equatorial and
C                   polar radii of the reference ellipsoid of the
C                   system. Then
C
C                       CORPAR(2) = ( RE - RP ) / RE
C
C
C     BOUNDS     is an 2x3 array containing the bounds of a volume
C                element expressed in a supported coordinate system,
C                specifically one of:
C
C                   LATITUDINAL
C                   PLANETODETIC
C                   RECTANGULAR
C
C
C                BOUNDS defines the volume element used in the
C                comparison. In the element
C
C                   BOUNDS(I,J) 
C
C                J is the coordinate index. I is the bound index.
C
C                   I = 1   ->   lower bound
C                   I = 2   ->   upper bound
C
C                See the routines 
C
C                   ZZINLAT
C                   ZZINPDT
C                   ZZINREC
C
C                for details on the contents of BOUNDS for the
C                respective coordinate systems supported by those
C                routines.
C                
C
C     MARGIN     is a fraction used to expand the volume element for
C                inclusion testing. 
C
C                See the routines 
C
C                   ZZINLAT
C                   ZZINPDT
C                   ZZINREC
C
C                for details regarding the application of MARGIN.
C
C
C     EXCLUD     is either a coordinate index or the parameter NONE.
C
C                If EXCLUD is set to one of
C
C                   { 1, 2, 3 }
C
C                then the indicated coordinate is excluded from
C                comparison with the corresponding volume element
C                boundaries.
C
C                If EXCLUD is set to NONE, all coordinates are
C                compared.
C
C                Exclusion of coordinates is used in cases where a
C                point is known to be on a level surface of a given
C                coordinate. For example, if a point is on the sphere
C                of radius equal to the upper radius bound, radius need
C                not be used in the comparison and in fact can't be
C                meaningfully compared, due to round-off errors.
C
C                See the routines 
C
C                   ZZINLAT
C                   ZZINPDT
C                   ZZINREC
C
C                for details regarding the application of EXCLUD.
C
C$ Detailed_Output
C
C     INSIDE     is a logical flag that is set to .TRUE. if and only if
C                the input coordinates represent a point inside or on
C                the surface of the volume element, according to the
C                comparisons that are performed.
C
C                The value of INSIDE is not affected by the value of
C                any excluded coordinate.
C
C$ Parameters
C
C     See the include files
C
C        dskdsc.inc
C        dsktol.inc
C
C     and the Parameters sections of the routines
C
C        ZZINLAT
C        ZZINPDT
C        ZZINREC
C
C
C$ Exceptions
C
C     1)  If CORSYS is not recognized, the error SPICE(NOTSUPPORTED) 
C         is signaled.
C
C     2)  If MARGIN is negative, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     3)  If EXCLUD is less than 0 or greater than 3, the error will
C         be signaled by a routine in the call tree of this routine.
C
C     4)  If an error occurs while determining the planetodetic
C         coordinates of the input point, the error will be signaled by
C         a routine in the call tree of this routine.
C
C     5)  If any rectangular coordinate upper bound is less than the
C         corresponding lower bound, the error will be signaled by
C         a routine in the call tree of this routine.
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
C     See usage in ZZPTPL02.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by the DSK
C     subsystem.
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
C-    SPICELIB Version 1.0.0, 04-APR-2017 (NJB) 
C
C     03-JUN-2016 (NJB) 
C
C        Original version.
C
C-&
 
C$ Index_Entries
C
C     test point against volume element using margin
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               RETURN
      
C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZINVELT' ) 


      IF ( MARGIN .LT. 0.D0 ) THEN

         CALL SETMSG ( 'Margin must be non-negative but was #.' )
         CALL ERRDP  ( '#', MARGIN                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZINVELT'                               )
         RETURN

      END IF

C
C     Delegate the job to one of the coordinate system-specific
C     routines.
C
      IF ( CORSYS .EQ. LATSYS ) THEN

         CALL ZZINLAT ( P, BOUNDS, MARGIN, EXCLUD, INSIDE )


      ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
 
         CALL ZZINPDT ( P, BOUNDS, CORPAR, MARGIN, EXCLUD, INSIDE )


      ELSE IF ( CORSYS .EQ. RECSYS ) THEN

         CALL ZZINREC ( P, BOUNDS, MARGIN, EXCLUD, INSIDE )


      ELSE

         CALL SETMSG ( 'Coordinate system code # was not recognized.' )
         CALL ERRINT ( '#', CORSYS                                    )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                          )
         CALL CHKOUT ( 'ZZINVELT'                                     )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZINVELT' ) 
      RETURN
      END 

