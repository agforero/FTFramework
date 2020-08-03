C$Procedure DSKMI2 ( DSK, make spatial index for type 2 segment )

      SUBROUTINE DSKMI2 ( NV,     VRTCES, NP,     PLATES, FINSCL,
     .                    CORSCL, WORKSZ, VOXPSZ, VOXLSZ, MAKVTL,
     .                    SPXISZ, WORK,   SPAIXD, SPAIXI         ) 

C$ Abstract
C
C     Make spatial index for a DSK type 2 segment. The index is
C     returned as a pair of arrays, one of type INTEGER and one of type
C     DOUBLE PRECISION. These arrays are suitable for use with the DSK
C     type 2 writer DSKW02.
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
C     DAS
C     DSK
C     FILES
C     PLATE
C     TOPOGRAPHY
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsk02.inc'
      
      INTEGER               NV
      DOUBLE PRECISION      VRTCES ( 3, * )
      INTEGER               NP
      INTEGER               PLATES ( 3, * )
      DOUBLE PRECISION      FINSCL
      INTEGER               CORSCL
      INTEGER               WORKSZ
      INTEGER               VOXPSZ
      INTEGER               VOXLSZ
      LOGICAL               MAKVTL
      INTEGER               SPXISZ
      INTEGER               WORK   ( 2, WORKSZ )
      DOUBLE PRECISION      SPAIXD ( IXDFIX )
      INTEGER               SPAIXI ( SPXISZ )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     IXDFIX     P   Size of fixed-size portion of d.p. index component.
C     IXIFIX     P   Size of fixed-size portion of integer index 
C                    component.
C     NV         I   Number of vertices.
C     VRTCES     I   Vertices.
C     NP         I   Number of plates.
C     PLATES     I   Plates.
C     FINSCL     I   Fine voxel scale.
C     CORSCL     I   Coarse voxel scale.
C     WORKSZ     I   Workspace size.
C     VOXPSZ     I   Voxel-plate pointer array size.
C     VOXLSZ     I   Voxel-plate list array size.
C     MAKVTL     I   Vertex-plate list flag.
C     SPXISZ     I   Spatial index integer component size.
C     WORK       I   Workspace.
C     SPAIXD     I   Double precision component of spatial index.
C     SPAIXI     I   Integer component of spatial index.
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
C
C     FINSCL      is the fine voxel scale. This scale determines the
C                 edge length of the cubical voxels comprising the fine
C                 voxel grid: the edge length VOXSIZ is approximately
C                 
C                     FINSCL * {average plate extent}
C
C                 where the extents of a plate are the respective
C                 differences between the maximum and minimum
C                 coordinate values of the plate's vertices.
C                 
C                 The relationship between VOXSIZ and the average plate
C                 extent is approximate because the VOXSIZ is adjusted
C                 so that each dimension of the fine voxel grid is an
C                 integer multiple of the coarse voxel scale.
C
C                 See the Particulars section below for further 
C                 information on voxel scales.
C      
C     
C     CORSCL      is the coarse voxel scale. This integer scale is the 
C                 ratio of the edge length of coarse voxels to
C                 that of fine voxels. The coarse scale must be
C                 large enough so that the total number of coarse
C                 voxels does not exceed MAXCGR (see the Parameters
C                 section below).
C
C
C     WORKSZ      is the second dimension of the workspace array WORK.
C                 WORKSZ must be at least as large as the greater of
C                 
C                    - the number of fine voxel-plate associations
C
C                      This number is equal to
C
C                         NP * {average number of fine voxels 
C                               intersected by each plate}
C
C                    - the number of vertex-plate associations, if
C                      the vertex-plate mapping is constructed.
C 
C                      This number is equal to
C
C                         NV + ( 3 * NP )
C                      
C
C     VOXPSZ      is the size of the fine voxel-plate pointer array.
C                 This array maps fine voxels to lists of plates that
C                 intersect those voxels. VOXPSZ must be at least as
C                 large as
C
C                          3
C                    CORSCL  * {number of non-empty coarse voxels}
C
C
C     VOXLSZ      is the size of the fine voxel-plate list array. This
C                 array contains, for each non-empty fine voxel, the
C                 count of plates that intersect that voxel and the
C                 IDs of those plates. VOXLSZ must be at least as large
C                 as
C
C                         NP * {average number of fine voxels 
C                               intersected by each plate}
C                    
C                     +   {number of non-empty fine voxels}
C
C
C     MAKVTL      is a logical flag that, when set to .TRUE., indicates
C                 that a  vertex-plate association list is to be
C                 constructed.
C
C                 The amount of workspace that is needed may depend on
C                 whether a vertex-plate association list is
C                 constructed. When this list is constructed, the size
C                 of the integer component of the spatial index is
C                 increased by the size of the list and the size of a
C                 vertex-plate pointer array; the total of these sizes
C                 is
C
C                    ( 2 * NV ) + ( 3 * NP )
C
C
C     SPXISZ      is the declared size of the output array SPAIXI. This
C                 size must be at least as large as the sum of
C
C                    - the fixed-size part of the integer component of
C                      the index, which includes the coarse voxel grid;
C                      this value is 
C
C                         IXIFIX
C
C                    - the size VOXPSZ of the voxel-plate pointer array
C
C                    - the size VOXLSZ of the voxel-plate association
C                      list
C      
C                 plus, if the vertex-plate association list is
C                 constructed,
C
C                    - the size NV of the vertex-plate pointer array
C       
C                    - the size of the vertex-plate association list; 
C                      this size is
C
C                         NV + ( 3 * NP )
C                 
C                  
C     WORK        is the workspace array. The array should be declared
C                 with dimensions
C
C                    (2, WORKSZ)
C
C                 See the description of WORKSZ above.
C
C
C$ Detailed_Output
C
C     WORK        is the workspace array, modified by the operations
C                 performed by this routine.
C
C     SPAIXD,
C     SPAIXI      are, respectively, the double precision and integer
C                 components of the spatial index of the segment.
C
C                 SPAIXD must be declared with size at least IXDFIX.
C                 SPAIXI must be declared with size at least SPXISZ.
C
C$ Parameters
C
C     IXDFIX      is the size of the double precision component of
C                 the spatial index.
C
C     IXIFIX      is the size of the fixed-size portion of the integer
C                 component of the spatial index.
C
C     See the include file dsk02.inc for declarations of the public DSK
C     type 2 parameters used by this routine.
C
C$ Exceptions
C
C     1)  If the fine voxel scale is non-positive, the error
C         SPICE(BADFINEVOXELSCALE) is signaled.
C
C     2)  If the coarse voxel scale is less than 1, the error
C         SPICE(BADCOARSEVOXSCALE) is signaled.
C
C     3)  If NV is less than 3 or greater than MAXVRT, the error
C         SPICE(BADVERTEXCOUNT) is signaled.
C
C     4)  If NP is less than 1 or greater than MAXPLT, the error
C         SPICE(BADPLATECOUNT) is signaled.
C
C     5)  If the workspace size WORKSZ is less than NP+1, the error
C         SPICE(WORKSPACETOOSMALL) is signaled. This is merely a 
C         sanity check; normally the workspace will need to be 
C         substantially larger than this reference value. See the 
C         description of WORKSZ in the header section Detailed_Input
C         above.
C
C     6)  If the voxel-plate pointer array size VOXPSZ is less than 1,
C         the error SPICE(PTRARRAYTOOSMALL) is signaled. This is merely
C         a sanity check; normally this pointer array will need to be
C         substantially larger than this reference value. See the
C         description of VOXPSZ in the header section Detailed_Input
C         above.
C
C     7)  If the voxel-plate list array size VOXLSZ is less than NP+1,
C         the error SPICE(PLATELISTTOOSMALL) is signaled. This is
C         merely a sanity check; normally this array will need to be
C         substantially larger than this reference value. See the
C         description of VOXLSZ in the header section Detailed_Input
C         above.
C
C     8)  If the size SPXISZ of the integer array SPAIXI is too small
C         to contain its constituent structures, where the sizes
C         of these structures are derived from the inputs
C
C             NV, NP, VOXPSZ, VOXLSZ
C             
C         the error SPICE(INTINDEXTOOSMALL) will be signaled.
C
C     9)  If there is insufficient room to create any of the data
C         structures contained in the spatial index, the error
C         will be diagnosed and signaled by a routine in the call
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
C     by creating the "spatial index" arrays required as inputs to that
C     routine.
C
C     A spatial index is a group of data structures that facilitates
C     rapid high-level computations involving sets of plates. The data
C     structures created by this routine are aggregated into arrays
C     of type INTEGER and type DOUBLE PRECISION. 
C
C
C     Voxel grids
C     ===========
C
C     A key geometric computation---probably the most important, as it
C     serves as a foundation for other high-level computations---is
C     finding the intersection of a ray with the plate set. DSK type 2
C     segments use data structures called "voxel grids" as part of
C     their indexing mechanism. There is a "coarse grid": a box that
C     completely encloses a DSK type 2 segment's plate set, and which
C     is composed of identically-sized cubes called "coarse voxels."
C     Each coarse voxel in composed of smaller cubes called "fine
C     voxels." When the term "voxel" is used without qualification, it
C     refers to fine voxels.
C
C     Type 2 DSK segments contain data structures that associate plates
C     with the fine voxels intersected by those plates. These
C     structures enable the type 2 DSK software to rapidly find plates
C     in a given region of space.
C
C     Voxel scales
C     ============
C     
C     There are two voxel scales:
C
C        - The coarse voxel scale is the integer ratio of the
C          edge length of a coarse voxel to the edge length of
C          a fine voxel
C
C        - The fine voxel scale is the double precision ratio
C          of the edge length of a fine voxel to the average
C          extent of the plates in the input plate set. "Extents"
C          of a plate are the absolute values of the differences 
C          between the respective maximum and minimum X, Y, and Z
C          coordinates of the plate's vertices.
C
C     Voxel scales determine the resolution of the voxel grid. 
C     Voxel scales must be chosen to satisfy size constraints and
C     provide reasonable plate lookup performance.
C
C     The following considerations apply to spatial indexes of
C     type 2 DSK segments:
C
C        1)  The maximum number of coarse voxels is fixed at MAXCGR
C            (declared in dsk02.inc).
C
C        2)  If there are too few fine voxels, the average number of
C            plates per fine voxel will be very large. This largely
C            negates the performance improvement afforded by having an
C            index. Also, the number of plates per voxel may exceed
C            limits imposed by DSK subroutines that use static arrays.
C
C        3)  If there are too many fine voxels, the average number of
C            voxels intersected by a given plate may be too large for
C            all the plate-voxel associations to be stored. In
C            addition, the time needed to examine the plate lists for
C            each voxel (including the empty ones) may become quite
C            large, again negating the value of the index.
C          
C     In many cases, voxel scales yielding optimum performance must be
C     determined by experiment. However, the following heuristics can
C     provide reasonable starting values:
C
C        Let NP be the number of plates. Let FS be the fine voxel
C        scale. Then a reasonable value of FS may be
C
C                   (0.25D0)
C           FS =  NP       / 8.D0
C
C        In general, FS should not smaller than 1.
C        
C
C$ Examples
C
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
C-    SPICELIB Version 1.0.0, 13-DEC-2016 (NJB)
C
C        Updated check on NV.
C
C        16-MAR-2016 (NJB)
C
C        Now zeros out the size of the vertex-plate list
C        when the list is not created.
C
C        23-JAN-2016 (NJB)
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     make spatial index for type 2 dsk segment
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C

C
C     Local variables
C
      INTEGER               I
      INTEGER               J
      INTEGER               NSHIFT
      INTEGER               NVXTOT
      INTEGER               REQSIZ
      INTEGER               VTLIDX
      INTEGER               VTPIDX
      INTEGER               VTXLSZ
      INTEGER               VXLIDX
      INTEGER               VXPIDX


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKMI2' )

C
C     Perform error checks on inputs.
C
      IF (  FINSCL .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Fine voxel scale = #; scale must ' 
     .   //            'be positive. Usually scale should '
     .   //            'be > 1.0.'                         )
         CALL ERRDP  ( '#',  FINSCL                        )
         CALL SIGERR ( 'SPICE(BADFINEVOXELSCALE)'          )
         CALL CHKOUT ( 'DSKMI2'                            )
         RETURN
         
      END IF

      IF (  CORSCL .LT. 1 ) THEN

         CALL SETMSG ( 'Coarse voxel scale = #; scale must ' 
     .   //            'be >= 1.'                           )
         CALL ERRINT ( '#',  CORSCL                         )
         CALL SIGERR ( 'SPICE(BADCOARSEVOXSCALE)'           )
         CALL CHKOUT ( 'DSKMI2'                             )
         RETURN
         
      END IF


      IF (  ( NV .LT. 3 ) .OR. ( NV .GT. MAXVRT )  ) THEN

         CALL SETMSG ( 'Vertex count NV = #; count must ' 
     .   //            'be in the range 3:#.'            )
         CALL ERRINT ( '#',  NV                          )
         CALL ERRINT ( '#',  MAXVRT                      )
         CALL SIGERR ( 'SPICE(BADVERTEXCOUNT)'           )
         CALL CHKOUT ( 'DSKMI2'                          )
         RETURN
         
      END IF

      IF (  ( NP .LT. 1 ) .OR. ( NP .GT. MAXPLT )  ) THEN

         CALL SETMSG ( 'Plate count NP = #; count must ' 
     .   //            'be in the range 1:#.'            )
         CALL ERRINT ( '#',  NP                          )
         CALL ERRINT ( '#',  MAXPLT                      )
         CALL SIGERR ( 'SPICE(BADPLATECOUNT)'            )
         CALL CHKOUT ( 'DSKMI2'                          )
         RETURN
         
      END IF

      IF ( WORKSZ .LT. (NP+1)  ) THEN

         CALL SETMSG ( 'Workspace size = #; size is too ' 
     .   //            'small to hold all voxel-plate '
     .   //            'associations. Size should be at '
     .   //            'least # * (average number of '
     .   //            'voxels intersected by each plate).' )
         CALL ERRINT ( '#',  WORKSZ                         )
         CALL ERRINT ( '#',  NP                             )
         CALL SIGERR ( 'SPICE(WORKSPACETOOSMALL)'           )
         CALL CHKOUT ( 'DSKMI2'                             )
         RETURN
         
      END IF

      IF ( VOXPSZ .LT. 1 ) THEN

         CALL SETMSG ( 'Voxel-pointer array size = #; size is ' 
     .   //            'too small to hold all voxel-plate '
     .   //            'list pointers. Size should be at '
     .   //            'least # * (number of non-empty '
     .   //            'coarse voxels).'                      )
         CALL ERRINT ( '#',  VOXPSZ                           )
         CALL ERRINT ( '#',  CORSCL**3                        )
         CALL SIGERR ( 'SPICE(PTRARRAYTOOSMALL)'              )
         CALL CHKOUT ( 'DSKMI2'                               )
         RETURN
         
      END IF


      IF ( VOXLSZ .LT. (NP+1)  ) THEN

         CALL SETMSG ( 'Voxel-plate list array size = #; size ' 
     .   //            'is too small to hold all voxel-plate '
     .   //            'associations. Size should be at '
     .   //            'least # * (average number of '
     .   //            'voxels intersected by each plate).'   )
         CALL ERRINT ( '#',  VOXLSZ                           )
         CALL ERRINT ( '#',  NP                               )
         CALL SIGERR ( 'SPICE(PLATELISTTOOSMALL)'             )
         CALL CHKOUT ( 'DSKMI2'                               )
         RETURN
         
      END IF

C
C     Check the size of the integer spatial index array. The
C     declared size must be large enough to hold:
C
C        - the fixed-size part of the index, which includes
C          the coarse voxel grid
C
C        - the voxel-plate pointer array
C
C        - the voxel-plate association list
C      
C     plus, if the vertex-plate association list is constructed,
C
C        - the vertex-plate pointer array
C       
C        - the vertex-plate association list
C        
C
      REQSIZ = IXIFIX + VOXPSZ + VOXLSZ 

      IF ( MAKVTL ) THEN
C
C        Add on the sizes of the vertex-plate pointer array (NV)
C        and the vertex-plate list array (NV + 3*NP).
C
         VTXLSZ = NV  +  3*NP

         REQSIZ = REQSIZ + NV + VTXLSZ

      ELSE
         VTXLSZ = 0
      END IF


      IF ( SPXISZ .LT. REQSIZ ) THEN

         CALL SETMSG ( 'Integer spatial index size = #; size '     
     .   //            'must be at least #.'                   )
         CALL ERRINT ( '#',  SPXISZ                            )
         CALL ERRINT ( '#',  REQSIZ                            )
         CALL SIGERR ( 'SPICE(INTINDEXTOOSMALL)'               )
         CALL CHKOUT ( 'DSKMI2'                                )
         RETURN

      END IF


C
C     Set known values in spatial index arrays.
C
      SPAIXI(SICGSC) = CORSCL
     
C
C     Prepare indices in the spatial index arrays.
C
C        VXPIDX is the start index of the voxel pointer array.
C
      VXPIDX = SICGRD + MAXCGR
C
C        VXLIDX is the start index of the voxel-plate list. This
C        list is offset from the start of the pointer array by
C        the input size given for that array. The size is the
C        total room available, not the room actually used.
C
      VXLIDX = VXPIDX + VOXPSZ

C
C     Create spatial index for plates.
C
      CALL ZZMKSPIN ( NP,      PLATES,  VRTCES, 
     .                FINSCL,  CORSCL,  VOXPSZ, 
     .                WORKSZ,  VOXLSZ,  WORK,
     .                SPAIXI(SIVGRX),
     .                SPAIXD(SIVXSZ),
     .                SPAIXD(SIVXOR),
     .                NVXTOT,  
     .                SPAIXI(SIVXNP),  
     .                SPAIXI(VXPIDX),  
     .                SPAIXI(SIVXNL),  
     .                SPAIXI(VXLIDX),  
     .                SPAIXD(SIVTBD), 
     .                SPAIXI(SICGRD)           )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'DSKMI2' )
         RETURN
      END IF

C
C     At this point the voxel plate list is offset from the
C     start of the voxel pointer array by the allocated size
C     of the array. We need to shift the plate list so that
c     it starts right after the end of the pointer array.
C
      NSHIFT = VOXPSZ - SPAIXI(SIVXNP) 

      DO I = 1, SPAIXI(SIVXNL) 

         J                = VXLIDX - 1 + I

         SPAIXI(J-NSHIFT) = SPAIXI(J)

      END DO

C
C     Update the voxel list start index to reflect the shift.
C     
      VXLIDX = VXLIDX - NSHIFT

C
C     Create vertex-plate mapping, if requested, as indicated
C     by the vertex-plate list size.
C
      IF ( MAKVTL ) THEN
C
C           VTPIDX is the start index of the vertex pointer array.
C
         VTPIDX = VXLIDX + SPAIXI(SIVXNL)
C
C           VXLIDX is the start index of the vertex-plate list. The
C           list start is offset from the vertex pointer array by
C           the size of the array, which is always NV.
C
         VTLIDX = VTPIDX + NV

         CALL ZZVRTPLT ( NV,     NP,     PLATES, 
     .                   WORKSZ, VTXLSZ, WORK, 
     .                   SPAIXI(VTPIDX), 
     .                   SPAIXI(SIVTNL), 
     .                   SPAIXI(VTLIDX)         )

      ELSE
C
C        Zero out the size of the vertex-plate list.
C
         SPAIXI(SIVTNL) = 0

      END IF

      CALL CHKOUT ( 'DSKMI2' )
      RETURN
      END
