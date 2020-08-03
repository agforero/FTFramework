C$Procedure ZZELLPLT ( Tessellate an ellipsoid with triangular plates )

      SUBROUTINE ZZELLPLT ( A,    B,  C,     NLON, NLAT,  MAXV, 
     .                      MAXP, NV, VERTS, NP,   PLATES      )
 
C$ Abstract
C
C     Create a set of triangular plates covering a specified triaxial
C     ellipsoid.
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
C     ELLIPSOID
C     PLATE
C     TILE
C     TESSELLATE
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      C
      INTEGER               NLON
      INTEGER               NLAT
      INTEGER               MAXV
      INTEGER               MAXP
      INTEGER               NV
      DOUBLE PRECISION      VERTS  ( 3, * )
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     A          I   Length of ellipsoid semi-axis lying on the x-axis.
C     B          I   Length of ellipsoid semi-axis lying on the y-axis.
C     C          I   Length of ellipsoid semi-axis lying on the z-axis.
C     NLON       I   Number of longitude bands in plate set.
C     NLAT       I   Number of latitude bands in plate set.
C     MAXV       I   Maximum number of vertices to return.
C     MAXP       I   Maximum number of plates to return.
C     NV         O   Number of vertices in output array.
C     VERTS      O   Vertices.
C     NP         O   Number of plates in output array.
C     PLATES     O   Plates.
C
C$ Detailed_Input
C
C     A,
C     B,
C     C          are the lengths of the semi-axes of a triaxial
C                ellipsoid. The ellipsoid is centered at the origin and
C                oriented so that its axes lie on the x, y and z axes.
C                A, B, and C are the lengths of the semi-axes that
C                point in the x, y, and z directions respectively.
C
C
C     NLON       is the number of longitude bands in the output plate
C                set. Each longitude band is bounded by two meridians.
C                All longitude bands have equal angular extent in
C                longitude. The vertices of any plate lie on adjacent
C                longitude band boundaries.
C
C
C     NLAT       is the number of latitude bands in the output plate
C                set. The vertices of each band are bounded by two
C                cones of constant planetocentric latitude. All
C                latitude bands have equal angular extent in
C                planetocentric latitude. The vertices of any plate lie
C                on adjacent latitude band boundaries.
C
C                Each polar "cap" consists of one longitude band.
C
C
C     MAXV       is the maximum number of vertices to return. MAXV must
C                be at least
C
C                   ( NLON * ( NLAT - 1 ) )  +  2
C
C                The array VERTS must have size at least 3*MAXV.
C
C 
C     MAXP       is the maximum number of plates to return. MAXP must
C                be at least
C
C                   2 * NLON * ( NLAT - 1 )
C
C                The array PLATES must have size at least 3*MAXP.
C
C$ Detailed_Output
C
C     NV         is the number of vertices in the output array VERTS.
C
C     VERTS      is an array containing the vertices of the output
C                plate set. There is a vertex at each intersection of a
C                latitude band boundary and a longitude band boundary.
C                The vertices at the north and south poles are at
C                indices NV and NV-1, respectively. The non-polar
C                vertex indices start at 1. Non-polar vertices are
C                indexed in top-down, left-to-right order, with
C                vertices of each latitude band stored contiguously.
C     
C     NP         is the number of plates in the output array PLATES.
C                    
C
C     PLATES     is an array containing the tessellating plate set.
C                Each plate is an array of three vertex indices, where
C                the indices range from 1 to NV.
C
C                The vertices of any plate are ordered in the
C                right-handed sense about the outward normal direction
C                for that plate.
C
C                The non-polar plates---those not having a vertex at
C                the north or south pole---are indexed in top-down,
C                left-to-right order, with plates belonging to each
C                latitude band stored contiguously. Plates constituting
C                the north polar cap follow the non-polar plate set,
C                and plates constituting the south polar cap follow
C                those.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the length of any semi-axis of the ellipsoid is
C         non-positive, the error SPICE(INVALIDAXISLENGTH) is signaled.
C
C     2)  If NLAT is less than 2, the error SPICE(INVALIDCOUNT) is
C         signaled.
C
C     3)  If NLON is less than 3, the error SPICE(INVALIDCOUNT) is
C         signaled.
C
C     4)  If the number of vertices implied by the input values NLON
C         and NLAT exceeds MAXV, the error SPICE(ARRAYTOOSMALL) is
C         signaled.
C
C     5)  If the number of plates implied by the input values NLON
C         and NLAT exceeds MAXP, the error SPICE(ARRAYTOOSMALL) is
C         signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The vertex and plate sets created by this routine are suitable
C     for use in a type 2 DSK segment.
C
C     While the primary purpose of this routine is to support testing,
C     there may be some user applications for which a tessellated plate
C     model is valuable. For example, computing an estimate of the area
C     of a specified surface region may be simplified by using a plate
C     model.
C
C     Note that, for ellipsoids having three distinct radii, the Z
C     components of the vertices on any latitude band boundary (except
C     the poles themselves) will vary with longitude.
C
C     Also note that the horizontal edge of a plate may extend beyond
C     the boundaries of the latitude band containing the plate.
C     
C$ Examples
C
C     See use of C language version in tspice_c test families.
C
C$ Restrictions
C
C     1) For use only by TSPICE.
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
C-    TESTUTIL Version 1.0.0, 30-APR-2014 (NJB)
C
C-&

C$ Index_Entries
C
C     tessellate ellipsoid
C
C-&

C     
C     SPICELIB functions
C
      DOUBLE PRECISION      PI

      LOGICAL               FAILED      
      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      DIR    ( 3 )
      DOUBLE PRECISION      DLAT
      DOUBLE PRECISION      DLON
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      LEVEL
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      S

      INTEGER               BIX
      INTEGER               I
      INTEGER               J
      INTEGER               N
      INTEGER               NNP
      INTEGER               NSP
      INTEGER               PIX
      INTEGER               VIX
      

      IF ( RETURN() ) THEN
         RETURN 
      END IF

      CALL CHKIN ( 'ZZELLPLT' )

C
C     The semi-axes must have positive length.
C
      IF (       ( A .LE. 0.D0 )
     .     .OR.  ( B .LE. 0.D0 )
     .     .OR.  ( C .LE. 0.D0 )   )   THEN
 
         CALL SETMSG ( 'Semi-axis lengths:  A = #, B = #, C = #. ' )
         CALL ERRDP  ( '#', A                                      )
         CALL ERRDP  ( '#', B                                      )
         CALL ERRDP  ( '#', C                                      )
         CALL SIGERR ( 'SPICE(INVALIDAXISLENGTH)'                  )
         CALL CHKOUT ( 'ZZELLPLT'                                    )
         RETURN
 
      END IF

C
C     The longitude and latitude band counts must be realizable.
C
      IF ( NLAT .LT. 2 ) THEN

         CALL SETMSG ( 'The latitude band count must be at least '
     .   //            '2 but was #.'                              )
         CALL ERRINT ( '#', NLAT                                   )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                       )
         CALL CHKOUT ( 'ZZELLPLT'                                    )
         RETURN

      END IF

      IF ( NLON .LT. 3 ) THEN

         CALL SETMSG ( 'The longitude band count must be at least '
     .   //            '3 but was #.'                              )
         CALL ERRINT ( '#', NLON                                   )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                       )
         CALL CHKOUT ( 'ZZELLPLT'                                    )
         RETURN

      END IF

C
C     Compute the vertex and plate counts. Check against available
C     room.
C
C        Vertex count: there are NLAT-2 latitude bands, excluding
C                      the polar caps. These are bounded by NLAT-1 rows
C                      of vertices. Each row of vertices has NLON
C                      members. The caps add two vertices.
C
C        Plate count:  each latitude band, excluding the polar caps,
C                      contains 2*NLON plates. Each cap contains NLON
C                      plates.
C
C
      NV = ( NLON * ( NLAT - 1 ) )  +  2

      NP = 2 * NLON * ( NLAT - 1 )


      IF ( NV .GT. MAXV ) THEN

         CALL SETMSG ( 'The requested plate model requires # '
     .   //            'vertices but the maximum vertex count is #.' )
         CALL ERRINT ( '#', NV                                       )
         CALL ERRINT ( '#', MAXV                                     )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                        )
         CALL CHKOUT ( 'ZZELLPLT'                                      )
         RETURN

      END IF

      IF ( NP .GT. MAXP ) THEN

         CALL SETMSG ( 'The requested plate model requires # '
     .   //            'plates but the maximum plate count is #.' )
         CALL ERRINT ( '#', NP                                    )
         CALL ERRINT ( '#', MAXP                                  )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                     )
         CALL CHKOUT ( 'ZZELLPLT'                                   )
         RETURN

      END IF

C
C     Create the vertex set. The north polar vertex is
C     at index 1; the south vertex is at index NV. It will
C     be convenient to make these the last two vertices.
C
      CALL VPACK ( 0.D0, 0.D0,  C, VERTS(1, NV-1) )
      CALL VPACK ( 0.D0, 0.D0, -C, VERTS(1, NV  ) )

C
C     The latitude bands are equally spaced in planetocentric
C     latitude.
C
      DLAT =     PI() / NLAT 
      DLON = 2 * PI() / NLON

      VIX  = 1

      DO I = 1, NLAT-1 

         LAT = ( PI()/2 ) - ( I * DLAT )

         DO J = 1, NLON

            LON = (J-1) * DLON
C
C           Create a unit direction vector for the current
C           vertex. Scale this vector to make it lie on the
C           ellipsoid's surface; the scaled vector is the
C           current vertex.
C
            CALL LATREC ( 1.D0, LON, LAT, DIR )
            
            LEVEL = (DIR(1)/A)**2 + (DIR(2)/B)**2 + (DIR(3)/C)**2

            S     = 1.D0 / SQRT( LEVEL )

            CALL VSCL ( S, DIR, VERTS(1,VIX) )
C
C           Next vertex.
C
            VIX   = VIX + 1

         END DO

      END DO

C
C     Create the plates for the latitude bounds other than
C     those belonging to the caps.
C
C     The first two inputs are the vertex row and column counts.
C     Next is a logical flag indicating whether longitude wrapping
C     should be used.
C     
      IF ( NLAT .GT. 2 ) THEN

         CALL ZZGRDPLT ( NLAT-1, NLON, .TRUE., N, PLATES )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZELLPLT' )
            RETURN
         END IF

      END IF

C
C     Add the north cap. This is a set of plates; the vertices
C     already have been computed.
C
      PIX = ( NP - 2*NLON ) + 1
      BIX = 0

      CALL ZZCAPPLT ( NLON, .TRUE., .TRUE., 
     .                BIX,  NV-1,   NNP,   PLATES(1,PIX) )

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZELLPLT' )
         RETURN
      END IF

C
C     Add the south cap.  
C
      PIX = PIX +   NLON
      BIX = NV  - ( NLON + 2 )

      CALL ZZCAPPLT ( NLON, .FALSE., .TRUE.,
     .                BIX,  NV,      NSP,    PLATES(1,PIX) )


      CALL CHKOUT ( 'ZZELLPLT' )
      RETURN
      END 


