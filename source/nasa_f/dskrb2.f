C$Procedure DSKRB2 ( DSK, determine range bounds for plate set )

      SUBROUTINE DSKRB2 ( NV,     VRTCES, NP,     PLATES,
     .                    CORSYS, CORPAR, MNCOR3, MXCOR3 )

C$ Abstract
C
C     Determine range bounds for a set of triangular plates to
C     be stored in a type 2 DSK segment.
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
C     DAS
C     DSK
C     FILES
C     PLATE
C     TOPOGRAPHY
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'

      INTEGER               NV
      DOUBLE PRECISION      VRTCES ( 3, * )
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      MNCOR3
      DOUBLE PRECISION      MXCOR3
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     NV         I   Number of vertices.
C     VRTCES     I   Vertices.
C     NP         I   Number of plates.
C     PLATES     I   Plates.
C     CORSYS     I   DSK coordinate system code.
C     CORPAR     I   DSK coordinate system parameters.
C     MNCOR3     O   Lower bound on range of third coordinate.
C     MXCOR3     O   Upper bound on range of third coordinate.
C
C$ Detailed_Input
C
C     NV          is the number of vertices belonging to the input
C                 set of plates.
C
C 
C     VRTCES      is an array of coordinates of the vertices. The Ith
C                 vertex occupies elements (1:3,I) of this array.
C
C
C     NP          is the number of plates in the input plate set.
C
C
C     PLATES      is an array representing the triangular plates of a
C                 shape model. The elements of PLATES are vertex
C                 indices; vertex indices are 1-based. The vertex
C                 indices of the Ith plate occupy elements (1:3,I) of
C                 this array.
C
C     CORSYS      is an integer parameter identifying the coordinate
C                 system in which the bounds are to be computed. The
C                 bounds apply to the third coordinate in each system:
C
C                    Latitudinal:           radius
C                    Planetodetic:          altitude
C                    Rectangular:           Z
C
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
C$ Detailed_Output
C
C     MNCOR3    is a lower bound on the range of the third coordinate
C               of the system identified by CORSYS and CORPAR, taken
C               over all plates.
C
C               For latitudinal and rectangular coordinates, MNCOR3
C               is the greatest lower bound of the third coordinate.
C
C               For planetodetic coordinates, MNCOR3 is an
C               approximation: it is less than or equal to the greatest
C               lower bound.
C
C     MXCOR3    is the least upper bound on the range of the third
C               coordinate of the system identified by CORSYS and
C               CORPAR, taken over all plates.
C
C
C$ Parameters
C
C     See the include file dskdsc.inc for declarations of the public DSK
C     type 2 parameters used by this routine.
C
C$ Exceptions
C
C     1)  If the input coordinate system is not recognized, the
C         error SPICE(NOTSUPPORTED) is signaled.
C
C     2)  If a conversion from rectangular to planetodetic coordinates
C         fails, the error will be signaled by a routine in the call
C         tree of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     Users planning to create DSK files should consider whether the 
C     SPICE DSK creation utility MKDSK may be suitable for their needs.
C
C     This routine supports use of the DSK type 2 segment writer DSKW02
C     by computing bounds on the range of the third coordinates of
C     the input plate set.
C
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C
C     1) Create a three-segment DSK file using plate model data for
C        Phobos. Use latitudinal, rectangular, and planetodetic
C        coordinates in the respective segments. This is not a 
C        realistic example, but it serves to demonstrate use of 
C        the supported coordinate systems.
C
C        For simplicity, use an existing DSK file to provide the 
C        input plate and vertex data. The selected input file has one
C        segment.
C
C
C     C
C     C     Example program for DSKW02, DSKMI2, and DSKRB2
C     C
C     C        Create a three-segment DSK file using plate model data
C     C        for Phobos. Use latitudinal, rectangular, and
C     C        planetodetic coordinates in the respective segments.
C     C
C     C        For simplicity, use an existing DSK file to provide the
C     C        input plate and vertex data. The selected input file has
C     C        one segment.
C     C
C     C           Version 1.0.0 22-JAN-2016 (NJB)
C     C
C           PROGRAM EX1
C           IMPLICIT NONE
C
C           INCLUDE 'dla.inc'
C           INCLUDE 'dskdsc.inc'
C           INCLUDE 'dsk02.inc'
C
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      JYEAR
C           DOUBLE PRECISION      PI
C     C
C     C     Local parameters
C     C
C           INTEGER               FRNMLN
C           PARAMETER           ( FRNMLN = 32 )
C
C           INTEGER               NSEG
C           PARAMETER           ( NSEG   = 3 )
C
C           INTEGER               NAMLEN
C           PARAMETER           ( NAMLEN = 20 )
C
C           INTEGER               FILSIZ
C           PARAMETER           ( FILSIZ = 255 )
C
C           INTEGER               LNSIZE
C           PARAMETER           ( LNSIZE = 80 )
C
C           INTEGER               NCOR
C           PARAMETER           ( NCOR   = 4 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(NAMLEN)    CORNAM ( NCOR )
C           CHARACTER*(FILSIZ)    DSK
C           CHARACTER*(FRNMLN)    FRAME
C           CHARACTER*(FILSIZ)    INDSK
C           CHARACTER*(LNSIZE)    LINE
C     C
C     C     Note: the values of MAXVRT and MAXPLT declared
C     C     in dsk02.inc, and the integer spatial index
C     C     dimension SPAISZ are very large. Smaller buffers
C     C     can be used for most applications.
C     C
C           DOUBLE PRECISION      CORPAR ( NSYPAR )
C           DOUBLE PRECISION      F
C           DOUBLE PRECISION      FINSCL
C           DOUBLE PRECISION      FIRST
C           DOUBLE PRECISION      LAST
C           DOUBLE PRECISION      MNCOR1
C           DOUBLE PRECISION      MNCOR2
C           DOUBLE PRECISION      MNCOR3
C           DOUBLE PRECISION      MXCOR1
C           DOUBLE PRECISION      MXCOR2
C           DOUBLE PRECISION      MXCOR3
C           DOUBLE PRECISION      RE
C           DOUBLE PRECISION      RP
C           DOUBLE PRECISION      SPAIXD ( IXDFIX )
C           DOUBLE PRECISION      VRTCES ( 3, MAXVRT )
C
C           INTEGER               CENTER
C           INTEGER               CORSCL
C           INTEGER               CORSYS
C           INTEGER               DCLASS
C           INTEGER               DLADSC ( DLADSZ )
C           INTEGER               HANDLE
C           INTEGER               INHAN
C           INTEGER               NP
C           INTEGER               NV
C           INTEGER               PLATES ( 3, MAXPLT )
C           INTEGER               SEGNO
C           INTEGER               SPAIXI ( SPAISZ )
C           INTEGER               SURFID
C           INTEGER               VOXPSZ
C           INTEGER               VOXLSZ
C           INTEGER               WORK   ( 2, MAXCEL )
C           INTEGER               WORKSZ
C
C           LOGICAL               FOUND
C     C
C     C     Saved variables
C     C
C     C     Save all large arrays to avoid stack problems.
C     C
C           SAVE
C     C
C     C     Initial values
C     C
C           DATA                  CORNAM / 'radius',
C          .                               'Z-coordinate',
C          .                               'Z-coordinate',
C          .                               'altitude'     /
C
C     C
C     C     Assign names of input and output DSK files.
C     C
C           INDSK = 'phobos_3_3.bds'
C           DSK   = 'phobos_3_3_3seg.bds'
C     C
C     C     Open input DSK for read access; find first segment.
C     C
C           CALL DASOPR ( INDSK, INHAN )
C           CALL DLABFS ( INHAN, DLADSC, FOUND )
C     C
C     C     Fetch vertices and plates from input DSK file.
C     C
C           WRITE (*,*) 'Reading input data...'
C
C           CALL DSKV02 ( INHAN, DLADSC, 1, MAXVRT, NV, VRTCES )
C           CALL DSKP02 ( INHAN, DLADSC, 1, MAXPLT, NP, PLATES )
C
C           WRITE (*,*) 'Done.'
C     C
C     C     Set input array sizes required by DSKMI2.
C     C
C           VOXPSZ = MAXVXP
C           VOXLSZ = MXNVLS
C           WORKSZ = MAXCEL
C     C
C     C     Set fine and coarse voxel scales. (These usually
C     C     need to determined by experimentation.)
C     C
C           FINSCL = 5.D0
C           CORSCL = 4
C     C
C     C     Open a new DSK file.
C     C
C           CALL DSKOPN ( DSK, DSK, 0, HANDLE )
C     C
C     C     Create three segments and add them to the file.
C     C
C           DO SEGNO = 1, NSEG
C     C
C     C        Create spatial index.
C     C
C              WRITE (*,*) 'Creating segment ', SEGNO
C              WRITE (*,*) 'Creating spatial index...'
C
C              CALL DSKMI2 ( NV,     VRTCES, NP,     PLATES, FINSCL,
C          .                 CORSCL, WORKSZ, VOXPSZ, VOXLSZ, .TRUE.,
C          .                 SPAISZ, WORK,   SPAIXD, SPAIXI          )
C
C              WRITE (*,*) 'Done.'
C     C
C     C        Set up inputs describing segment attributes:
C     C
C     C        - Central body: Phobos
C     C        - Surface ID code: user's choice.
C     C          We use the segment number here.
C     C        - Data class: general (arbitrary) shape
C     C        - Body-fixed reference frame
C     C        - Time coverage bounds (TBD)
C     C
C              CENTER = 401
C              SURFID = SEGNO
C              DCLASS = GENCLS
C              FRAME  = 'IAU_PHOBOS'
C
C              FIRST = -50 * JYEAR()
C              LAST  =  50 * JYEAR()
C     C
C     C        Set the coordinate system and coordinate system
C     C        bounds based on the segment index.
C     C
C     C        Zero out the coordinate parameters to start.
C     C
C              CALL CLEARD ( NSYPAR, CORPAR )
C
C              IF ( SEGNO .EQ. 1 ) THEN
C     C
C     C           Use planetocentric latitudinal coordinates. Set
C     C           the longitude and latitude bounds.
C     C
C                 CORSYS = LATSYS
C
C                 MNCOR1 = -PI()
C                 MXCOR1 =  PI()
C                 MNCOR2 = -PI()/2
C                 MXCOR2 =  PI()/2
C
C              ELSE IF ( SEGNO .EQ. 2 ) THEN
C     C
C     C           Use rectangular coordinates. Set the
C     C           X and Y bounds.
C     C
C     C           The bounds shown here were derived from
C     C           the plate data. They lie slightly outside
C     C           of the range spanned by the plates.
C     C
C                 CORSYS = RECSYS
C
C                 MNCOR1 = -1.3D0
C                 MXCOR1 =  1.31D0
C                 MNCOR2 = -1.21D0
C                 MXCOR2 =  1.2D0
C
C              ELSE
C     C
C     C           Set the coordinate system to planetodetic.
C     C
C                 CORSYS    = PDTSYS
C
C                 MNCOR1    = -PI()
C                 MXCOR1    =  PI()
C                 MNCOR2    = -PI()/2
C                 MXCOR2    =  PI()/2
C     C
C     C           We'll use equatorial and polar radii from
C     C           pck00010.tpc. These normally would be fetched
C     C           at run time, but for simplicity, we'll use
C     C           hard-coded values.
C
C                 RE        = 13.0D0
C                 RP        =  9.1D0
C                 F         = ( RE - RP ) / RE
C
C                 CORPAR(1) = RE
C                 CORPAR(2) = F
C
C              END IF
C     C
C     C        Compute plate model radius bounds.
C     C
C              LINE = 'Computing # bounds of plate set...'
C
C              CALL REPMC ( LINE, '#', CORNAM(CORSYS), LINE )
C              WRITE (*,*) LINE
C
C              CALL DSKRB2 ( NV,     VRTCES, NP,     PLATES,
C          .                 CORSYS, CORPAR, MNCOR3, MXCOR3 )
C
C              WRITE (*,*) 'Done.'
C     C
C     C        Write the segment to the file.
C     C
C              WRITE (*,*) 'Writing segment...'
C
C              CALL DSKW02 ( HANDLE,
C          .                 CENTER, SURFID, DCLASS, FRAME,  CORSYS,
C          .                 CORPAR, MNCOR1, MXCOR1, MNCOR2, MXCOR2,
C          .                 MNCOR3, MXCOR3, FIRST,  LAST,   NV,
C          .                 VRTCES, NP,     PLATES, SPAIXD, SPAIXI )
C
C              WRITE (*,*) 'Done.'
C
C           END DO
C     C
C     C     Segregate the data records in the DSK file and
C     C     close the file.
C     C
C           WRITE (*,*) 'Segregating and closing DSK file...'
C
C           CALL DSKCLS ( HANDLE, .TRUE. )
C
C           WRITE (*,*) 'Done.'
C           END
C
C     
C$ Restrictions
C
C     1) For planetodetic coordinates, the computation of the lower
C        altitude bound requires that the surface at altitude MNCOR3 be
C        convex. This is the case for realistic geometries, but can
C        be false if a plate is very large compared to the overall
C        shape model.
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
C       22-JAN-2016 (NJB) 
C
C         Original version.
C-&
 
C$ Index_Entries
C
C     compute range bounds for type 2 dsk segment
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      DPMIN
      DOUBLE PRECISION      VNORM
      DOUBLE PRECISION      VDIST

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      DOUBLE PRECISION      THIRD
      PARAMETER           ( THIRD  = 1.D0/3.D0 )
C
C
C     Local variables
C
      DOUBLE PRECISION      ALT
      DOUBLE PRECISION      CENTER ( 3 )
      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      F
      DOUBLE PRECISION      LAT
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      MAXD
      DOUBLE PRECISION      ORIGIN ( 3 )
      DOUBLE PRECISION      PNEAR  ( 3 )
      DOUBLE PRECISION      RE

      INTEGER               I

C
C     Saved variables
C
      SAVE                  ORIGIN

C
C     Initial values
C
      DATA                  ORIGIN / 3 * 0.D0 /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKRB2' )


      IF ( CORSYS .EQ. LATSYS ) THEN
C
C        The coordinate system is latitudinal.
C
C        Compute radius bounds. Start with the maximum radius.
C        This is simply the maximum norm of the vertices.
C
         MXCOR3 = 0.D0

         DO I = 1, NV
            MXCOR3 = MAX (  VNORM( VRTCES(1,I) ),  MXCOR3  )
         END DO

C
C        Compute the minimum radius of the plate set.
C
         MNCOR3 = DPMAX()

         DO I = 1, NP

            CALL PLTNP ( ORIGIN,
     .                   VRTCES( 1, PLATES(1,I) ),
     .                   VRTCES( 1, PLATES(2,I) ),
     .                   VRTCES( 1, PLATES(3,I) ),
     .                   PNEAR,
     .                   DIST                     )

            MNCOR3 = MIN ( DIST, MNCOR3 )

         END DO
 
         
      ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C        The coordinate system is rectangular. Compute the range
C        of Z-coordinates of the plates.
C
         MNCOR3 = DPMAX()
         MXCOR3 = DPMIN()

         DO I = 1, NV

            MNCOR3 = MIN( MNCOR3, VRTCES(3,I) )
            MXCOR3 = MAX( MXCOR3, VRTCES(3,I) )

         END DO


      ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C        The coordinate system is planetodetic. Compute the range
C        of altitudes of the plates.
C
         RE     = CORPAR(1)
         F      = CORPAR(2)

         MXCOR3 = DPMIN()
         MNCOR3 = DPMAX()

C
C        The maximum altitude is attained at a plate vertex.
C         
         DO I = 1, NV

            CALL RECGEO ( VRTCES(1,I), RE, F, LON, LAT, ALT )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'DSKRB2' )
               RETURN
            END IF

            MXCOR3 = MAX( MXCOR3, ALT )

         END DO

C
C        For the Ith plate, let DMAX(I) be the maximum distance between
C        the plate's center and any of the plate's vertices.
C
C        The minimum altitude is greater than or equal to 
C        the minimum of 
C
C           {altitude of the Ith plate's center - DMAX(I)}
C
C        taken over all plates.
C
         DO I = 1, NP

            CALL VLCOM3 ( THIRD, VRTCES( 1, PLATES(1,I) ),
     .                    THIRD, VRTCES( 1, PLATES(2,I) ),
     .                    THIRD, VRTCES( 1, PLATES(3,I) ),  CENTER )

            MAXD = MAX (  VDIST(  VRTCES( 1, PLATES(1,I) ), CENTER  ),
     .                    VDIST(  VRTCES( 1, PLATES(2,I) ), CENTER  ),
     .                    VDIST(  VRTCES( 1, PLATES(3,I) ), CENTER  )  )

            CALL RECGEO ( CENTER, RE, F, LON, LAT, ALT )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'DSKRB2' )
               RETURN
            END IF

            MNCOR3 = MIN (  MNCOR3,  ( ALT - MAXD )  )

         END DO
         

      ELSE

         CALL SETMSG ( 'Coordinate system # is not supported.' )
         CALL ERRINT ( '#', CORSYS                             )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                   )
         CALL CHKOUT ( 'DSKRB2'                                )
         RETURN

      END IF

      CALL CHKOUT ( 'DSKRB2' )
      RETURN
      END
