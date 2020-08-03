C$Procedure ZZMKSPIN ( Make spatial index of plates )
 
      SUBROUTINE ZZMKSPIN ( NP,     PLATES, VRTCES, VOXSCL, CGSCAL,
     .                      MAXPTR, MXCELL, MAXVXL, CELLS,  NVOX,
     .                      VOXSIZ, VOXORI, NVXTOT, NVXPTR, VXPTR,
     .                      NVXLST, VXLIST, EXTENT, CGRPTR         )
      IMPLICIT NONE

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Create voxel grid data structure and voxel-plate mapping.
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
C     plate voxel index
C
C$ Declarations
  
      INCLUDE               'dsk02.inc'
      
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
      DOUBLE PRECISION      VRTCES ( 3, * )
      DOUBLE PRECISION      VOXSCL
      INTEGER               CGSCAL
      INTEGER               MAXPTR
      INTEGER               MXCELL
      INTEGER               MAXVXL
      INTEGER               CELLS   ( 2, MXCELL )
      INTEGER               NVOX   ( 3 )
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VOXORI ( 3 )
      INTEGER               NVXTOT
      INTEGER               NVXPTR
      INTEGER               VXPTR  ( * )
      INTEGER               NVXLST
      INTEGER               VXLIST ( * )
      DOUBLE PRECISION      EXTENT ( 6 )
      INTEGER               CGRPTR ( * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NP         I   Number of plates.
C     PLATES     I   Array of plates.
C     VRTCES     I   Array of vertices.
C     VOXSCL     I   Voxel scaling factor.
C     CGSCAL     I   Coarse voxel grid scaling factor.
C     MAXPTR     I   Maximum voxel pointer list size.
C     MXCELL     I   Cell array column dimension.
C     MAXVXL     I   Maximum voxel-plate list size.
C     CELLS     I-O  Workspace array for vertex list construction.
C     NVOX       O   Dimensions of voxel grid.
C     VOXSIZ     O   Size of a voxel in model units.
C     VOXORI     O   Origin of voxel grid in model units.
C     NVXTOT     O   Total number of voxels in grid.
C     NVXPTR     O   Number of pointers in VXPTR array.
C     VXPTR      O   Pointer index to the voxel to plates list.
C     NVXLST     O   Length of voxel to plates list.
C     VXLIST     O   Voxel to plates list.
C     EXTENT     O   Vertex coordinate extents in km.
C     CGRPTR     O   Coarse voxel grid pointer array.
C
C$ Detailed_Input
C
C     NP         Total number of plates in model.
C
C     PLATES     Array containing indices in the array VRTCES of
C                the vertices corresponding to each plate. The 
C                elements (1:3,I) of PLATES are the vertex indices
C                of the Ith plate.
C
C     VRTCES     Array containing cartesian coordinates for all
C                vertices in the model reference frame. Elements
C                (1:3,I) of VRTCES are the coordinates of the Ith
C                vertex.
C
C     VOXSCL     The voxel size is determined by the average extent of
C                each plate in XYZ. VOXSCL is a scaling factor that
C                may be used to adjust final voxel size.
C
C     CGSCAL     The coarse voxel grid has voxel edge length
C                larger by factor CGSCAL than the fine voxel
C                edge length.
C
C     MAXPTR     is the maximum number of elements in the output
C                list of voxel pointers.
C
C     MXCELL     is the number of cells in the input cell array. 
C                This is the second dimension of the array.
C
C     MAXVXL     is the maximum number of elements in the output
C                voxel-plate list.
C
C     CELLS      workspace array used to construct the voxel-plate
C                mapping.
C
C$ Detailed_Output
C
C     CELLS      workspace array used to construct the voxel-plate
C                mapping.
C
C     NVOX       Dimensions of the voxel grid in voxel units.
C
C     VOXSIZ     Size of each voxel in model units (km).
C
C     VOXORI     Origin of voxel grid in model units (km).
C
C     NVXTOT     Total number of voxels in grid.
C
C     NVXPTR     Number of pointers in VXPTR array.
C
C     VXPTR      Array of pointers that map voxels to their associated
C                plated IDs in the voxel to plates list. The Nth
C                element of this array contains a pointer to the Mth
C                element of VXLIST, or -1. If VXPTR(N) is -1, the Nth
C                voxel is empty.
C
C     NVXLST     The total number of elements in the plate list.
C
C     VXLIST     The plates list. For the Nth voxel, VXPTR(N) points to
C                VXLIST(M). VXLIST(M) contains the number plates that
C                intersect voxel N. The IDs of those plates are stored
C                in the following elements.
C
C     EXTENT     is an array of extents of vertex coordinates. Elements
C                indexed (2*I)-1 and 2*I are the minimum and maximum 
C                values of the Ith coordinate, taken over all vertices.
C                Units are km.
C
C     CGRPTR     is an of array coarse voxel grid pointers. Null
C                pointers have the value -1. Non-null pointers indicate
C                locations in the voxel pointer array: each non-empty
C                coarse voxel points to a list of pointers for each of
C                the fine voxels it contains.
C                
C$ Parameters
C
C     See the include file dsk02.inc.
C
C$ Exceptions
C
C     1)  If NP is less than 1 or greater than MAXPLT, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C     2)  If the number of coarse voxels exceeds the grid size MAXCGR,
C         the error SPICE(COARSEGRIDOVERFLOW) is signaled.
C
C     3)  If the voxel count is greater than MAXVOX, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C     4)  If the coarse voxel count is less than 1 or greater than
C         MAXCGR, the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     5)  If the coarse voxel scale is less than 1 or more than
C         the cube root of the fine voxel count, the error
C         SPICE(VALUEOUTOFRANGE) will be signaled.
C
C     6)  If the cube of the coarse voxel scale does not divide the
C         fine voxel count evenly, the error SPICE(INCOMPATIBLESCALE)
C         will be signaled.
C
C     7)  If the workspace overflows while this routine accumulates
C         voxel-plate associations, the error will be signaled by 
C         a routine in the call tree of this routine.
C
C     8)  If the voxel-plate association list array overflows while
C         this routine accumulates voxel-plate associations, the error
C         will be signaled by a routine in the call tree of this
C         routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports creation of spatial indexes for DSK type 2
C     segments. This routine determines the fine and coarse voxel 
C     grid dimensions and parameters. It also builds the voxel-plate
C     association data structures.
C
C$ Examples
C
C     See usage in DSKMI2.
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
C     N.J. Bachman  (JPL)
C     J.A. Bytof    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 17-FEB-2017 (NJB)
C
C        Added new error checks.
C
C        29-JAN-2016 (NJB)
C
C           Removed reference to ZZVOXPAD. Renamed argument MAXCEL to
C           MXCELL to accommodate declaration of MAXCEL in dsk02.inc.
C
C           Updated long error messages and header comments.
C
C        26-MAR-2015 (NJB)
C
C           Functional change: now fills in the coarse voxel grid
C           pointer array with valid pointers; this is done in one
C           pass. The old version of the routine only marked the
C           elements of this array as null or non-null.
C
C           The algorithm for computing the set of voxels containing
C           the bounding box of each plate was re-written to improve
C           efficiency.
C
C        30-JUN-2014 (NJB)
C
C           Argument list change: removed work space array PNTRS. The
C           pointer array for the voxel-to-plate list mapping is now
C           constructed in place.
C
C           Changed argument list to include sizes of arrays. Changed
C           error handling to make use of array sizes. Changed call to
C           UNTNGL to accommodate argument list change in that routine.
C           Updated header I/O descriptions.
C
C        13-MAY-2010 (NJB)
C
C           Now accepts input workspace arrays PNTRS and CELLS.
C           Replaced argument VERT with arguments VRTCES and PLATES;
C           this was done to greatly reduce the memory needed by this
C           program.
C
C           Changed INCLUDE file to dsk02.inc.
C
C           Bug fix: determination of empty voxels is now based on
C           voxel pointers, not their target values.
C
C
C        08-OCT-2009 (NJB)
C
C           Re-ordered header sections.
C
C        12-SEP-2004 (EDW)
C
C           Improve/expand comments and descriptions.
C
C        03-FEB-1999 (JAB)
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     spatial index plates voxels
C
C-&

C$ Revisions
C
C        20-MAR-2015 (NJB)
C
C        Functional change: now fills in the coarse voxel grid pointer
C        array with valid pointers; this is done in one pass. The old
C        version of the routine only marked the elements of this array
C        as null or non-null.
C
C        The maximum size of the voxel-plate pointer array has been
C        reduced to MAXVOX/2.
C     
C        The array mapping voxels to the voxel-plate pointer array is 
C        no longer used. Previously, voxel IDs were used as the
C        "A list" elements with which plate lists were associated. 
C        Now, indices in the voxel-plate pointer array are used as
C        "A" values. 
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      BRCKTD
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      DPMIN

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Other functions
C
      INTEGER               ZZVOX2ID

C
C     Local parameters
C

C
C     Fraction of voxel edge length used a tolerance for plate
C     inclusion in voxels:
C
      DOUBLE PRECISION      TOL
      PARAMETER           ( TOL  = 1.D-3 )

C
C     Local Variables
C
      DOUBLE PRECISION      AVEXT
      DOUBLE PRECISION      BXMAX
      DOUBLE PRECISION      BXMIN
      DOUBLE PRECISION      BYMAX
      DOUBLE PRECISION      BYMIN
      DOUBLE PRECISION      BZMAX
      DOUBLE PRECISION      BZMIN
      DOUBLE PRECISION      CVXSIZ
      DOUBLE PRECISION      MDLTOL
      DOUBLE PRECISION      VMOD   ( 3 )
      DOUBLE PRECISION      XEXTNT ( 6 )
      DOUBLE PRECISION      XMAX
      DOUBLE PRECISION      XMIN
      DOUBLE PRECISION      XP     ( 3 )
      DOUBLE PRECISION      XVMAX
      DOUBLE PRECISION      XVMIN
      DOUBLE PRECISION      YMAX
      DOUBLE PRECISION      YMIN
      DOUBLE PRECISION      YP     ( 3 )
      DOUBLE PRECISION      YVMAX
      DOUBLE PRECISION      YVMIN
      DOUBLE PRECISION      ZMAX
      DOUBLE PRECISION      ZMIN
      DOUBLE PRECISION      ZP     ( 3 )
      DOUBLE PRECISION      ZVMAX
      DOUBLE PRECISION      ZVMIN

      INTEGER               CGOF1D
      INTEGER               CGOFF  ( 3 )
      INTEGER               CGRDIM ( 3 )
      INTEGER               CGXYZ  ( 3 )
      INTEGER               CVID
      INTEGER               GXMAX
      INTEGER               GXMIN
      INTEGER               GYMAX
      INTEGER               GYMIN
      INTEGER               GZMAX
      INTEGER               GZMIN
      INTEGER               I
      INTEGER               IX
      INTEGER               IXPTR
      INTEGER               IY
      INTEGER               IZ
      INTEGER               J
      INTEGER               NCELL
      INTEGER               NCGFLG
      INTEGER               NPCG
      INTEGER               NX
      INTEGER               NY
      INTEGER               NZ
      INTEGER               Q
      INTEGER               R
      INTEGER               TO
      INTEGER               VIXYZ  ( 3 )
      INTEGER               VCOORD ( 3 )
      
      LOGICAL               INBOX

C
C     Saved variables
C
C
C     Required for f2c use on Linux, all local variables
C     to static.
C 
      SAVE

C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZMKSPIN' )

C
C     Check NP. 
C
      IF (  ( NP .LT. 1 ) .OR. ( NP .GT. MAXPLT )  ) THEN

         CALL SETMSG ( 'Plate count NP = #; count must ' 
     .   //            'be in the range 1:#.'            )
         CALL ERRINT ( '#',  NP                          )
         CALL ERRINT ( '#',  MAXPLT                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'          )
         CALL CHKOUT ( 'ZZMKSPIN'                        )
         RETURN
         
      END IF

C
C     Make sure the coarse voxel scale is positive. We'll
C     perform additional checks later on. Those checks 
C     require computations that can't be done if the coarse
C     scale is zero.
C
      IF ( CGSCAL .LT. 1 ) THEN

         CALL SETMSG ( 'Coarse voxel scale = #; scale must ' 
     .   //            'be positive.'                       )
         CALL ERRINT ( '#',  CGSCAL                         )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'             )
         CALL CHKOUT ( 'ZZMKSPIN'                           )
         RETURN

      END IF

C
C     Get the average extents of all plates and the
C     overall model extent. The extents have model units; in
C     other words km.
C 
      AVEXT = 0.D0
 
      XMIN =  DPMAX()
      XMAX =  DPMIN()
      YMIN =  DPMAX()
      YMAX =  DPMIN()
      ZMIN =  DPMAX()
      ZMAX =  DPMIN()
 
      DO I = 1, NP
 
         BXMIN = DPMAX()
         BXMAX = DPMIN()
         BYMIN = DPMAX()
         BYMAX = DPMIN()
         BZMIN = DPMAX()
         BZMAX = DPMIN()
 
         XP(1) = VRTCES( 1, PLATES(1,I) )
         XP(2) = VRTCES( 1, PLATES(2,I) )
         XP(3) = VRTCES( 1, PLATES(3,I) )
 
         YP(1) = VRTCES( 2, PLATES(1,I) )
         YP(2) = VRTCES( 2, PLATES(2,I) )
         YP(3) = VRTCES( 2, PLATES(3,I) )
 
         ZP(1) = VRTCES( 3, PLATES(1,I) )
         ZP(2) = VRTCES( 3, PLATES(2,I) )
         ZP(3) = VRTCES( 3, PLATES(3,I) )
 
         DO J = 1, 3
C
C           Determine plate extents.
C
            BXMIN = DMIN1( BXMIN, XP(J) )
            BXMAX = DMAX1( BXMAX, XP(J) )
 
            BYMIN = DMIN1( BYMIN, YP(J) )
            BYMAX = DMAX1( BYMAX, YP(J) )
 
            BZMIN = DMIN1( BZMIN, ZP(J) )
            BZMAX = DMAX1( BZMAX, ZP(J) )

C
C           Determine model extent.
C
            XMIN = DMIN1 ( XMIN, BXMIN )
            XMAX = DMAX1 ( XMAX, BXMAX )
 
            YMIN = DMIN1 ( YMIN, BYMIN )
            YMAX = DMAX1 ( YMAX, BYMAX )
 
            ZMIN = DMIN1 ( ZMIN, BZMIN )
            ZMAX = DMAX1 ( ZMAX, BZMAX )
 
         END DO
 
         EXTENT(1) = XMIN
         EXTENT(2) = XMAX
         EXTENT(3) = YMIN
         EXTENT(4) = YMAX
         EXTENT(5) = ZMIN
         EXTENT(6) = ZMAX
         
C
C        Calculate the cumulative extent of all plates for
C        each degree of freedom.
C
         AVEXT = AVEXT + DABS(BXMAX-BXMIN) +
     .                   DABS(BYMAX-BYMIN) +
     .                   DABS(BZMAX-BZMIN)
 
      END DO

C
C     Calculate the average extent of all plates for
C     and the voxel size, i.e the length of one side
C     of a voxel cube.
C
      AVEXT  = AVEXT / DBLE(3*NP)
      VOXSIZ = VOXSCL * AVEXT
      MDLTOL = VOXSIZ * TOL
      
C
C     Produce a set of vertex extents, extended by MDLTOL,
C     to be used later.
C
      DO I = 1, 5, 2

         XEXTNT(I)   = EXTENT(I)   - MDLTOL
         XEXTNT(I+1) = EXTENT(I+1) + MDLTOL
   
      END DO

C
C     Determine the size of the coarse voxels.
C
      CVXSIZ = VOXSIZ * CGSCAL

C
C     Determine the minima and maxima of the body centered
C     vertex coordinates expressed in coarse voxel units. Scale the 
C     vertices coord values by CVXSIZ: this scales the 
C     axis in the voxel model space producing cubic voxels
C     with length 1 along each edge in voxel space,
C     CVXSIZ along an edge in model space.
C
      XVMIN = XMIN / CVXSIZ
      YVMIN = YMIN / CVXSIZ 
      ZVMIN = ZMIN / CVXSIZ
      XVMAX = XMAX / CVXSIZ
      YVMAX = YMAX / CVXSIZ
      ZVMAX = ZMAX / CVXSIZ
 
C
C     Extend the coarse voxel grid by at least 1/2
C     coarse voxel length along each degree of freedom.
C
      XVMIN = DNINT ( XVMIN - 1.D0 )
      YVMIN = DNINT ( YVMIN - 1.D0 )
      ZVMIN = DNINT ( ZVMIN - 1.D0 )
      XVMAX = DNINT ( XVMAX + 1.D0 )
      YVMAX = DNINT ( YVMAX + 1.D0 )
      ZVMAX = DNINT ( ZVMAX + 1.D0 )
 
C
C     Calculate the coarse voxel grid origin in model units.
C
      VOXORI(1) = XVMIN * CVXSIZ
      VOXORI(2) = YVMIN * CVXSIZ
      VOXORI(3) = ZVMIN * CVXSIZ
 
C
C     Calculate the dimension of the voxel grid in
C     units of (regular) voxels.
C
      NX = NINT( XVMAX - XVMIN ) * CGSCAL
      NY = NINT( YVMAX - YVMIN ) * CGSCAL
      NZ = NINT( ZVMAX - ZVMIN ) * CGSCAL

      NVOX(1) = NX
      NVOX(2) = NY
      NVOX(3) = NZ
 
      NVXTOT = NX*NY*NZ
C
C     Make sure the number of voxels NVXTOT is within range.
C
      IF ( NVXTOT .GT. MAXVOX ) THEN

         CALL SETMSG ( 'Fine voxel count NVXTOT = #; count must ' 
     .   //            'be in the range 1:#.'                    )
         CALL ERRINT ( '#',  NVXTOT                              )
         CALL ERRINT ( '#',  MAXVOX                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                  )
         CALL CHKOUT ( 'ZZMKSPIN'                                )
         RETURN
         
      END IF


C
C     Check the coarse voxel scale. It must be at least 1, and its
C     cube must not exceed the fine voxel count.
C
      IF (      ( CGSCAL .LT. 1                ) 
     .     .OR. ( CGSCAL .GT. NVXTOT**(1.D0/3) )  ) THEN

         CALL SETMSG ( 'Coarse voxel scale = #; scale must ' 
     .   //            'be in the range 1:NVXTOT**3, where '
     .   //            'NVXTOT is the total fine voxel count. '
     .   //            'In this case, NVXTOT = #.'               )
         CALL ERRINT ( '#',  CGSCAL                              )
         CALL ERRINT ( '#',  NVXTOT                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                  )
         CALL CHKOUT ( 'ZZMKSPIN'                                )
         RETURN

      END IF
      
C
C     The cube of the coarse scale must divide the total voxel count
C     evenly. This is a consistency check: the code that derives the
C     voxel grid dimensions should ensure this condition is met.
C
      Q = NVXTOT /   (CGSCAL**3)
      R = NVXTOT - Q*(CGSCAL**3)

      IF ( R .NE. 0 ) THEN

         CALL SETMSG ( 'Coarse voxel scale = #; the cube of '
     .   //            'the scale must divide NVXTOT evenly, ' 
     .   //            'where NVXTOT is the total  fine voxel '
     .   //            'count. In this case, NVXTOT = #.'        )
         CALL ERRINT ( '#',  CGSCAL                              )
         CALL ERRINT ( '#',  NVXTOT                              )
         CALL SIGERR ( 'SPICE(INCOMPATIBLESCALE)'                )
         CALL CHKOUT ( 'ZZMKSPIN'                                )
         RETURN

      END IF

C
C     Check the number of coarse voxels.
C
      NPCG   = CGSCAL ** 3

      NCGFLG = NVXTOT / NPCG
 
      IF ( NCGFLG .GT. MAXCGR ) THEN

         CALL SETMSG ( 'Number of coarse voxels # exceeds '
     .   //            'limit #. Increase coarse voxel scale, '
     .   //            'fine voxel scale, or both.'            )
         CALL ERRINT ( '#', NCGFLG                             )
         CALL ERRINT ( '#', MAXCGR                             )
         CALL SIGERR ( 'SPICE(COARSEGRIDOVERFLOW)'             )
         CALL CHKOUT ( 'ZZMKSPIN'                              )
         RETURN

      END IF

C
C     Enumerate all voxels that each plate might intersect.
C
      CALL ZZINILNK ( MAXPTR, MXCELL, NCELL, VXPTR, CELLS )

C
C     Set the dimensions of the coarse grid.
C
      CGRDIM(1) = NX / CGSCAL
      CGRDIM(2) = NY / CGSCAL
      CGRDIM(3) = NZ / CGSCAL

      CALL CLEARI ( NCGFLG, CGRPTR )


C
C     TO points to the first free location in the VXPTR array.
C
      TO = 1

      DO I = 1, NP
C
C        Find the extents of the Ith plate, where the extents
C        are expanded by TOL in each direction. We truncate
C        the expanded box at a distance of MDLTOL beyond the
C        extents of the vertex set, if necessary.
C
         XP(1) = VRTCES( 1, PLATES(1,I) )
         XP(2) = VRTCES( 1, PLATES(2,I) )
         XP(3) = VRTCES( 1, PLATES(3,I) )
 
         YP(1) = VRTCES( 2, PLATES(1,I) )
         YP(2) = VRTCES( 2, PLATES(2,I) )
         YP(3) = VRTCES( 2, PLATES(3,I) )
 
         ZP(1) = VRTCES( 3, PLATES(1,I) )
         ZP(2) = VRTCES( 3, PLATES(2,I) )
         ZP(3) = VRTCES( 3, PLATES(3,I) )

        

         BXMIN = BRCKTD( MIN( XP(1), XP(2), XP(3) ) - MDLTOL,
     .                   XEXTNT(1),  XEXTNT(2)                )

         BXMAX = BRCKTD( MAX( XP(1), XP(2), XP(3) ) + MDLTOL,
     .                   XEXTNT(1),  XEXTNT(2)                )

         BYMIN = BRCKTD( MIN( YP(1), YP(2), YP(3) ) - MDLTOL,
     .                   XEXTNT(3),  XEXTNT(4)                )

         BYMAX = BRCKTD( MAX( YP(1), YP(2), YP(3) ) + MDLTOL,
     .                   XEXTNT(3),  XEXTNT(4)                )

         BZMIN = BRCKTD( MIN( ZP(1), ZP(2), ZP(3) ) - MDLTOL,
     .                   XEXTNT(5),  XEXTNT(6)                )

         BZMAX = BRCKTD( MAX( ZP(1), ZP(2), ZP(3) ) + MDLTOL,
     .                   XEXTNT(5),  XEXTNT(6)                )

C
C        Find the range of voxel coordinates that contain the bounding
C        box of the plate. All we need look at are the coordinates
C        of the two corners having minimum and maximum coordinates.
C
C        Start with the corner having minimum coordinates:
C       
         CALL VPACK ( BXMIN, BYMIN, BZMIN, VMOD )

         CALL ZZGETVOX ( VOXSIZ, VOXORI, NVOX, VMOD, INBOX, VCOORD )

         IF ( .NOT. INBOX ) THEN
C
C           A corner of the bounding box lies outside the voxel grid.
C           This should never occur.
C
            CALL SETMSG ( 'BUG: bounding box of plate is outside of '
     .      //            'voxel grid. Input coordinates were '
     .      //            '(#, #, #). Plate ID = #.'                  )
            CALL ERRDP  ( '#', VMOD(1)                                )
            CALL ERRDP  ( '#', VMOD(2)                                )
            CALL ERRDP  ( '#', VMOD(3)                                )
            CALL ERRINT ( '#', I                                      )
            CALL SIGERR ( 'SPICE(BUG)'                                )
            CALL CHKOUT ( 'ZZMKSPIN'                                  )
            RETURN         
         
         END IF

C
C        Unpack minimum voxel coordinates from VCOORD.
C
         GXMIN = VCOORD(1)
         GYMIN = VCOORD(2)
         GZMIN = VCOORD(3)

C
C        Now handle the corner having maximum coordinates:
C       
         CALL VPACK ( BXMAX, BYMAX, BZMAX, VMOD )

         CALL ZZGETVOX ( VOXSIZ, VOXORI, NVOX, VMOD, INBOX, VCOORD )

         IF ( .NOT. INBOX ) THEN
C
C           A corner of the bounding box lies outside the voxel grid.
C           This should never occur.
C
            CALL SETMSG ( 'BUG: bounding box of plate is outside of '
     .      //            'voxel grid. Input coordinates were '
     .      //            '(#, #, #). Plate ID = #.'                  )
            CALL ERRDP  ( '#', VMOD(1)                                )
            CALL ERRDP  ( '#', VMOD(2)                                )
            CALL ERRDP  ( '#', VMOD(3)                                )
            CALL ERRINT ( '#', I                                      )
            CALL SIGERR ( 'SPICE(BUG)'                                )
            CALL CHKOUT ( 'ZZMKSPIN'                                  )
            RETURN         
         
         END IF

C
C        Unpack maximum voxel coordinates from VCOORD.
C
         GXMAX = VCOORD(1)
         GYMAX = VCOORD(2)
         GZMAX = VCOORD(3)
                          
C
C        Determine voxels that the bounding box of the plate
C        intersects.
C        
         DO IZ = GZMIN, GZMAX

             DO IY = GYMIN, GYMAX

                 DO IX = GXMIN, GXMAX

                    VIXYZ(1) = IX
                    VIXYZ(2) = IY
                    VIXYZ(3) = IZ
C
C                   Find the coarse voxel containing this voxel, and
C                   compute the offset of this voxel within the coarse
C                   voxel. The output CGXYZ contains the 3-dimensional
C                   coordinates of the coarse voxel within the coarse
C                   grid. The output CGOF1D is the 1-based,
C                   1-dimensional offset of the current voxel (having
C                   coordinates VIXYZ) from the start of the coarse
C                   voxel.
C
                    CALL ZZVOXCVO ( VIXYZ, NVOX,  CGSCAL, 
     .                              CGXYZ, CGOFF, CGOF1D  )

                    IF ( FAILED() ) THEN
                       CALL CHKOUT ( 'ZZMKSPIN' )
                       RETURN
                    END IF

                    CVID = ZZVOX2ID ( CGXYZ, CGRDIM )

                    IF ( CGRPTR(CVID) .EQ. 0 ) THEN
C
C                      The coarse voxel at index CVID is empty so far.
C                      Allocate CGSCAL pointers for it in the VXPTR
C                      array; make the coarse voxel point to the first
C                      element of this sub-array.
C
                       CGRPTR(CVID) = TO
                       TO           = TO + NPCG

                    END IF
C
C                   Let IXPTR be the index in the VXPTR array of the
C                   pointer for the current voxel.
C
                    IXPTR = ( CGRPTR(CVID) - 1 )  +  CGOF1D

                    CALL ZZADDLNK ( IXPTR,  I,     MAXPTR, 
     .                              MXCELL, VXPTR, NCELL, CELLS )

                    IF ( FAILED() ) THEN
                       CALL CHKOUT ( 'ZZMKSPIN' )
                       RETURN
                    END IF


                 END DO

              END DO

          END DO

      END DO

      NVXPTR = TO - 1
          
C
C     Generate two linked lists mapping voxel ID to the plates enclosed
C     within that voxel (if any).
C
C     VXPTR : An array, indexed by voxel ID. For an array element, 
C             VXPTR(VOX_ID), greater than zero, the value identifies an 
C             index in VXLIST, the value of that VXLIST array element 
C             equaling the number of plates contained in the voxel 
C             specified by the ID. The condition VXPTR(VOX_ID) = -1 
C             indicates the voxel contains no plates.
C
C     VXLIST: An array, indexed by the positive entries in VXPTR. The 
C             element, N, identified by a VXPTR value describes the 
C             number of plates contained in the corresponding voxel. 
C
C                 N = VXLIST( VXPTR(VOX_ID) )
C
C             The N elements following VXLIST( VXPTR(VOX_ID) ),
C             contain the IDs of those plates within the voxel.
C
      CALL ZZUNTNGL ( NVXPTR, MXCELL, CELLS, 
     .                MAXVXL, VXPTR,  NVXLST, VXLIST )
     

      CALL CHKOUT ( 'ZZMKSPIN' )
      RETURN
      END
