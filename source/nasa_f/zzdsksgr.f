C$Procedure ZZDSKSGR ( DSK, return segment bounding radius )
 
      DOUBLE PRECISION FUNCTION ZZDSKSGR ( DSKDSC )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Return a radius of a bounding sphere for the spatial region
C     described by a DSK segment descriptor. The radius is measured
C     from the origin of the coordinate system associated from the 
C     segment.
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
C     DAS
C     DSK
C
C$ Keywords
C
C     DSK
C     GEOMETRY
C     UTILITY
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'

      DOUBLE PRECISION      DSKDSC ( * )

      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     DSKDSC     I   is a DSK descriptor.
C 
C     The function returns the maximum radius of the segment.
C
C$ Detailed_Input
C
C     DSKDSC         is a DSK descriptor.
C                    
C$ Detailed_Output
C
C     The function returns the maximum radius attained on the
C     boundary of the coverage region of the segment associated
C     with the input descriptor. The radius is measured from
C     the origin of the coordinate system in which the coverage
C     bounds are given.
C                   
C     Units are km.
C                   
C$ Parameters
C
C     See the INCLUDE file dskdsc.inc for parameter declarations and
C     documentation.
C
C$ Exceptions
C
C     1)  If the coordinate system code in the descriptor is not
C         recognized, the error SPICE(NOTSUPPORTED) is signaled.
C
C     2)  If an invalid radius is found, the error
C         SPICE(INVALIDVALUE) is signaled.
C
C     3)  If an invalid flattening coefficient is found, the error
C         SPICE(INVALIDVALUE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine computes an upper bound on the radius of 
C     the coverage region of a DSK segment. The radius at 
C     a point on the boundary of the coverage region is
C     measured from the origin of the coordinate system.
C
C     This routine can be used to generate a bounding sphere for
C     a surface described by one or more DSK segments. Note that
C     the DSK segments describing the surface of an object don't
C     necessarily have a common central body, so offsets between
C     segments' central bodies must be taken into account when
C     computing a bounding sphere.
C
C$ Examples
C
C
C     1) Get a bounding radius for a DSK segment using the
C        latitudinal coordinate system.
C
C        Example code begins here.
C
C
C                 PROGRAM DSKGSR_EX1
C                 IMPLICIT NONE
C
C                 INCLUDE 'dskdsc.inc'
C
C                 DOUBLE PRECISION      ZZDSKSGR
C                 DOUBLE PRECISION      PI
C
C                 DOUBLE PRECISION      DSKDSC( DSKDSZ )
C
C           C
C           C     Latitudinal coordinates:
C           C
C                 CALL CLEARD( DSKDSZ, DSKDSC )
C
C                 DSKDSC( SYSIDX ) = LATSYS
C
C           C
C           C     Fill in the coordinate bounds in the descriptor.
C           C
C           C     Longitude:
C           C
C                 DSKDSC(MN1IDX) =  - PI()/2
C                 DSKDSC(MX1IDX) =    PI()/2
C           C
C           C     Latitude:
C           C
C                 DSKDSC(MN2IDX) =    0
C                 DSKDSC(MX2IDX) =    PI()/4
C           C
C           C     Radius:
C           C
C                 DSKDSC(MN3IDX) =    1.D3
C                 DSKDSC(MX3IDX) =    1.1D3
C
C                 WRITE (*,*) 'Radius bound for latitudinal section: ',
C                .            ZZDSKSGR( DSKDSC )
C 
C
C                 END
C
C
C        When run on a PC/Linux/gfortran platform, the output
C        from this program was:
C 
C           Radius bound for latitudinal section:    1100.0000000000000
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
C-    SPICELIB Version 1.0.0, 13-FEB-2017 (NJB)
C
C        Added case for planetodetic coordinates; deleted case for
C        cylindrical coordinates. Recoded algorithm for rectangular
C        coordinates.
C
C        22-JUL-2016 (NJB)
C
C           Renamed to ZZDSKSGR. 
C
C        14-MAY-2010 (NJB)
C
C-&
 
C$ Index_Entries
C
C     Return upper bound on dsk segment coverage radius
C     
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      VNORM

C
C     Local variables
C
      DOUBLE PRECISION      BDS    ( 2, 3 )
      DOUBLE PRECISION      F
      DOUBLE PRECISION      MAG    ( 3 )
      DOUBLE PRECISION      MINR
      DOUBLE PRECISION      RE
      DOUBLE PRECISION      RP

      INTEGER               B
      INTEGER               CORSYS
      INTEGER               I
      INTEGER               J


C
C     Use discovery check-in.
C

C
C     Set an initial return value.
C
      ZZDSKSGR = -1.D0

      
C
C     The radius calculation depends on the coordinate system.
C
      CORSYS = NINT( DSKDSC(SYSIDX) )

      IF ( CORSYS .EQ. LATSYS ) THEN
C
C        Fetch the minimum radius from the descriptor.
C
         MINR = DSKDSC(MN3IDX)

         IF ( MINR .LE. 0.D0 ) THEN

            CALL CHKIN  ( 'ZZDSKSGR'                )
            CALL SETMSG ( 'Minimum radius was *.'   )
            CALL ERRDP  ( '*', MINR                 )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'  )
            CALL CHKOUT ( 'ZZDSKSGR'                )
            RETURN

         END IF

C
C        This is as simple as it gets. The radius bounds
C        correspond to the third coordinate in the descriptor.
C
         ZZDSKSGR = DSKDSC( MX3IDX )


      ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C        Fetch the equatorial radius from the descriptor.
C
         RE = DSKDSC(PARIDX)

         IF ( RE .LE. 0.D0 ) THEN

            CALL CHKIN  ( 'ZZDSKSGR'                  )
            CALL SETMSG ( 'Equatorial radius was *.'  )
            CALL ERRDP  ( '*', RE                     )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'    )
            CALL CHKOUT ( 'ZZDSKSGR'                  )
            RETURN

         END IF

C
C        Fetch the flattening coefficient from the descriptor.
C
         F = DSKDSC(PARIDX+1)

         IF ( ( F .GE. 0.D0 ) .AND. ( F .LT. 1.D0 ) ) THEN
C
C           This is the oblate case.
C
C           The maximum radius of an oblate planetodetic boundary
C           occurs on the X-Y plane at the maximum height.
C
            ZZDSKSGR = DSKDSC(MX3IDX) +  RE


         ELSE IF ( F .LT. 0.D0 ) THEN
C
C           This is the prolate case.
C
C           The maximum radius of an prolate planetodetic boundary
C           occurs on the poles at the maximum height.
C
            RE       = DSKDSC(PARIDX)
            RP       = RE  * ( 1.D0 - F )
            ZZDSKSGR = DSKDSC(MX3IDX) +  RP
         
         ELSE
C
C           We have an invalid flattening coefficient.
C
C           If the flattening coefficient is greater than one, the
C           polar radius computed below is negative. If it's equal to
C           one, the polar radius is zero. Either case is a problem, so
C           signal an error and check out.
C
            CALL CHKIN  ( 'ZZDSKSGR'                       )
            CALL SETMSG ( 'Flattening coefficient was *.'  )
            CALL ERRDP  ( '*', F                           )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
            CALL CHKOUT ( 'ZZDSKSGR'                       )
            RETURN

         END IF


      ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C        The bounding cell of the segment has its extreme value at
C        a corner. Just take the maximum of these values.
C
C        First copy the bounds into an appropriately dimensioned
C        array.
C
         CALL MOVED ( DSKDSC(MN1IDX), 6, BDS )


          B = MN1IDX - 1

          DO I = 1, 3

             J      = B + (2*I) - 1

             MAG(I) = MAX(  ABS(DSKDSC(J)), ABS(DSKDSC(J+1))  )

          END DO

          ZZDSKSGR  = VNORM( MAG )


      ELSE
C
C        Never heard of this coordinate system.
C
         CALL CHKIN ( 'ZZDSKSGR'                                       )
         CALL SETMSG( 'The coordinate system code # is not recognized.')
         CALL ERRINT( '#',  CORSYS                                     )
         CALL SIGERR( 'SPICE(NOTSUPPORTED)'                            )
         CALL CHKOUT( 'ZZDSKSGR'                                       )
         RETURN

      END IF

      END
