C$Procedure ZZELLSEC ( Tessellate an ellipsoid section with plates )

      SUBROUTINE ZZELLSEC ( A,      B,      C,     MINLON, MAXLON,
     .                      MINLAT, MAXLAT, NLON,  NLAT,   MAXV, 
     .                      MAXP,   NV,     VERTS, NP,     PLATES )
 
C$ Abstract
C
C     Create a set of triangular plates covering a specified section
C     of the surface of a triaxial ellipsoid. The boundaries of the
C     section are curves of constant planetocentric longitude and
C     latitude.
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
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MAXLAT
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
C     MINLON     I   Minimum longitude of section.
C     MAXLON     I   Minimum longitude of section.
C     MINLAT     I   Minimum latitude of section.
C     MAXLAT     I   Minimum latitude of section.
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
C     MINLON,
C     MAXLON,
C     MINLAT,
C     MAXLAT     are, respectively, the longitude and latitude bounds
C                of a section of the surface the triaxial ellipsoid.
C                The coordinate system is latitudinal. Units are
C                radians.
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
C     MAXV       is the maximum number of vertices to return. The
C                number of vertices created depends on which polar
C                caps are created. In all cases the number does not
C                exceed
C
C                   ( NLON + 1 ) * ( NLAT + 1 )
C
C                The array VERTS must have size at least 3*MAXV.
C
C 
C     MAXP       is the maximum number of plates to return. The
C                number of plates created depends on which polar
C                caps are created and whether longitude wrapping
C                is selected. In all cases the number does not
C                exceed
C
C                   2 * NLON * NLAT
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
C                indices NV and NV-1, respectively, if both polar caps
C                are created. The non-polar vertex indices start at 1.
C                Non-polar vertices are indexed in top-down,
C                left-to-right order, with vertices of each latitude
C                band stored contiguously.
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
C-    TESTUTIL Version 1.0.0, 29-SEP-2014 (NJB)
C
C        Based on ZZELLPLT Version 1.0.0, 30-APR-2014 (NJB)
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
C     Local parameters
C
      DOUBLE PRECISION      TOL
      PARAMETER           ( TOL  = 1.D-12 )

C
C     Local variables
C
      DOUBLE PRECISION      DIR    ( 3 )
      DOUBLE PRECISION      DLAT
      DOUBLE PRECISION      DLON
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      LEVEL
      DOUBLE PRECISION      LMXLON
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      S

      INTEGER               BIX
      INTEGER               I
      INTEGER               J
      INTEGER               LATLB
      INTEGER               LATUB
      INTEGER               N
      INTEGER               NCOLS
      INTEGER               NNP
      INTEGER               NROWS
      INTEGER               NSP
      INTEGER               PIX
      INTEGER               POLIDX
      INTEGER               VIX
      
      LOGICAL               NCAP
      LOGICAL               SCAP
      LOGICAL               WRAP

      IF ( RETURN() ) THEN
         RETURN 
      END IF

      CALL CHKIN ( 'ZZELLSEC' )

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
         CALL CHKOUT ( 'ZZELLSEC'                                  )
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
         CALL CHKOUT ( 'ZZELLSEC'                                  )
         RETURN

      END IF

      IF ( NLON .LT. 3 ) THEN

         CALL SETMSG ( 'The longitude band count must be at least '
     .   //            '3 but was #.'                              )
         CALL ERRINT ( '#', NLON                                   )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                       )
         CALL CHKOUT ( 'ZZELLSEC'                                  )
         RETURN

      END IF

C
C     Decide whether we have two distinct longitude boundaries. First
C     create a local maximum longitude that's greater than the minimum
C     longitude. The logical variable WRAP is .TRUE. if and only if we
C     have 2*pi - TOL radians of longitude coverage, where TOL is a
C     small value.
C
      IF ( MAXLON .GT. MINLON ) THEN
         
         LMXLON = MAXLON
      ELSE
         LMXLON = MAXLON + (2*PI())
      END IF

      WRAP = ( LMXLON - MINLON ) .GT. ( 2*PI() - TOL )  

C
C     Decide whether we have north or south polar caps.
C     
      NCAP = MAXLAT .GT. ( ( PI()/2) - TOL )
      SCAP = MINLAT .LT. ( (-PI()/2) + TOL )

C
C     Compute the vertex counts. 
C
      IF ( WRAP ) THEN
C
C        Vertex count:  When both caps are present, there are NLAT-2
C                       latitude bands, excluding the polar caps. These
C                       are bounded by NLAT-1 rows of vertices. Each
C                       row of vertices has NLON members. The caps add
C                       two vertices.      
C
         IF ( NCAP .AND. SCAP ) THEN
C
C           There are two polar caps.
C
            NV = ( NLON * ( NLAT - 1 ) )  +  2

         ELSE IF ( NCAP .OR. SCAP ) THEN
C
C           There's just one polar cap. Excluding the polar
C           vertex, there are NLAT rows of vertices.
C
            NV = ( NLON * NLAT )  +  1
            
         ELSE
C
C           No polar caps. There are NLAT+1 rows of vertices.
C
            NV = NLON * ( NLAT + 1 )

         END IF

      ELSE

         IF ( NCAP .AND. SCAP ) THEN
C
C           There are two polar caps.
C
            NV = ( (NLON+1) * ( NLAT - 1 ) )  +  2

         ELSE IF ( NCAP .OR. SCAP ) THEN
C
C           There's just one polar cap. Excluding the polar
C           vertex, there are NLAT rows of vertices.
C
            NV = ( (NLON+1) * NLAT )  +  1
            
         ELSE
C
C           No polar caps. There are NLAT+1 rows of vertices.
C
            NV = ( NLON + 1 ) * ( NLAT + 1 )

         END IF

      END IF


C
C     Compute the plate counts. These depend on the set of
C     polar caps.
C
C
C        Plate count:   each latitude band, excluding the polar caps,
C                       contains 2*NLON plates. Each cap contains NLON
C                       plates.
C
      IF ( NCAP .AND. SCAP ) THEN
C
C        There are two polar caps.
C
         NP = 2 * NLON * ( NLAT - 1 )

      ELSE IF ( NCAP .OR. SCAP ) THEN
C
C        There's just one polar cap. Excluding the polar
C        vertex, there are NLAT rows of vertices.
C            
         NP = NLON * ( ( 2 * NLAT ) - 1 )

      ELSE
C
C        No polar caps. There are NLAT+1 rows of vertices.
C
         NP = 2 * NLON * NLAT 

      END IF
 


      IF ( NV .GT. MAXV ) THEN

         CALL SETMSG ( 'The requested plate model requires # '
     .   //            'vertices but the maximum vertex count is #.' )
         CALL ERRINT ( '#', NV                                       )
         CALL ERRINT ( '#', MAXV                                     )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                        )
         CALL CHKOUT ( 'ZZELLSEC'                                      )
         RETURN

      END IF

      IF ( NP .GT. MAXP ) THEN

         CALL SETMSG ( 'The requested plate model requires # '
     .   //            'plates but the maximum plate count is #.' )
         CALL ERRINT ( '#', NP                                    )
         CALL ERRINT ( '#', MAXP                                  )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                     )
         CALL CHKOUT ( 'ZZELLSEC'                                   )
         RETURN

      END IF

C
C     Create the vertex set, excluding the polar caps.
C
C     LATLB will be the index of the first vertex row, excluding
C     the caps. LATUB will be the index of the last vertex row,
C     excluding the caps.
C
      IF ( NCAP .AND. SCAP ) THEN

         LATLB = 2
         LATUB = NLAT 

      ELSE IF ( NCAP ) THEN

         LATLB = 2
         LATUB = NLAT + 1

      ELSE IF ( SCAP ) THEN

         LATLB = 1
         LATUB = NLAT

      ELSE

         LATLB = 1
         LATUB = NLAT + 1

      END IF


C
C     NCOLS is the number of columns of vertices.
C
      IF ( WRAP ) THEN
         NCOLS = NLON
      ELSE
         NCOLS = NLON + 1
      END IF

C
C     The latitude bands are equally spaced in planetocentric latitude.
C
      DLAT = ( MAXLAT - MINLAT ) / NLAT 
      DLON = ( LMXLON - MINLON ) / NLON

      VIX  = 1

      DO I = LATLB, LATUB

         LAT = MAXLAT - ( (I-1) * DLAT )

         DO J = 1, NCOLS

            LON = MINLON  +  ( (J-1) * DLON )
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
C     Create the polar vertices if necessary.
C
      IF ( NCAP .AND. SCAP ) THEN
 
         CALL VPACK ( 0.D0, 0.D0,  C, VERTS(1, NV-1) )
         CALL VPACK ( 0.D0, 0.D0, -C, VERTS(1, NV  ) )

      ELSE IF ( NCAP ) THEN

         CALL VPACK ( 0.D0, 0.D0,  C, VERTS(1, NV) )

      ELSE IF ( SCAP ) THEN

         CALL VPACK ( 0.D0, 0.D0, -C, VERTS(1, NV) )

      END IF



C
C     Create the plates for the latitude bounds other than
C     those belonging to the caps.
C
C     The first two inputs are the vertex row and column counts.
C     Next is a logical flag indicating whether longitude wrapping
C     should be used.
C     
      IF ( LATUB .GT. LATLB ) THEN

         NROWS = LATUB - LATLB + 1

         CALL ZZGRDPLT ( NROWS, NCOLS, WRAP, N, PLATES )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZELLSEC' )
            RETURN
         END IF

      ELSE

         N = 0

      END IF


      IF ( NCAP ) THEN
C
C        Add the north cap. This is a set of plates; the vertices
C        already have been computed.
C
C        PIX is the index of the first cap plate. BIX is the
C        base (predecessor) index of the first vertex in the
C        first vertex row.
C
         PIX = N + 1
         BIX = 0

C
C        POLIDX is the vertex index of the north polar vertex.
C        
         IF ( SCAP ) THEN
            POLIDX = NV - 1
         ELSE
            POLIDX = NV
         END IF

         CALL ZZCAPPLT ( NCOLS, .TRUE., WRAP,  
     .                   BIX,   POLIDX, NNP,  PLATES(1,PIX) )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZELLSEC' )
            RETURN
         END IF

      END IF


      IF ( SCAP ) THEN
C
C        Add the south cap.  
C
         POLIDX = NV

         IF ( NCAP ) THEN

            PIX = PIX +  NNP
            BIX = NV  - ( NCOLS + 2 )

         ELSE

            PIX = N  + 1
            BIX = NV - ( NCOLS + 1 )

         END IF

         CALL ZZCAPPLT ( NCOLS, .FALSE., WRAP,
     .                   BIX,   POLIDX,  NSP,  PLATES(1,PIX) )
      END IF

      CALL CHKOUT ( 'ZZELLSEC' )
      RETURN
      END 


