C$Procedure DSKW02 ( DSK, write type 2 segment )
 
      SUBROUTINE DSKW02 ( HANDLE, 
     .                    CENTER, SURFID, DCLASS, FRAME,  CORSYS,
     .                    CORPAR, MNCOR1, MXCOR1, MNCOR2, MXCOR2,
     .                    MNCOR3, MXCOR3, FIRST,  LAST,   NV,      
     .                    VRTCES, NP,     PLATES, SPAIXD, SPAIXI ) 

      IMPLICIT NONE

C$ Abstract
C
C     Write a type 2 segment to a DSK file.
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

      INCLUDE  'dskdsc.inc'
      INCLUDE  'dsktol.inc'
      INCLUDE  'dsk02.inc'

      INTEGER               HANDLE
      INTEGER               CENTER
      INTEGER               SURFID
      INTEGER               DCLASS
      CHARACTER*(*)         FRAME
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      MNCOR1
      DOUBLE PRECISION      MXCOR1
      DOUBLE PRECISION      MNCOR2
      DOUBLE PRECISION      MXCOR2
      DOUBLE PRECISION      MNCOR3
      DOUBLE PRECISION      MXCOR3
      DOUBLE PRECISION      FIRST
      DOUBLE PRECISION      LAST
      INTEGER               NV
      DOUBLE PRECISION      VRTCES ( 3, NV )
      INTEGER               NP
      INTEGER               PLATES ( 3, NP )
      DOUBLE PRECISION      SPAIXD ( * )
      INTEGER               SPAIXI ( * )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle assigned to the opened DSK file.
C     CENTER     I   Central body ID code.
C     SURFID     I   Surface ID code.
C     DCLASS     I   Data class.
C     FRAME      I   Reference frame.
C     CORSYS     I   Coordinate system code.
C     CORPAR     I   Coordinate system parameters.
C     MNCOR1     I   Minimum value of first coordinate.
C     MXCOR1     I   Maximum value of first coordinate.
C     MNCOR2     I   Minimum value of second coordinate.
C     MXCOR2     I   Maximum value of second coordinate.
C     MNCOR3     I   Minimum value of third coordinate.
C     MXCOR3     I   Maximum value of third coordinate.
C     FIRST      I   Coverage start time.
C     LAST       I   Coverage stop time.
C     NV         I   Number of vertices.
C     VRTCES     I   Vertices.
C     NP         I   Number of plates.
C     PLATES     I   Plates.
C     SPAIXD     I   Double precision component of spatial index.
C     SPAIXI     I   Integer component of spatial index.
C     ANGMRG     P   Angular round-off margin.
C     GENCLS     P   General surface DSK class.
C     SVFCLS     P   Single-valued function DSK class.
C     NSYPAR     P   Maximum number of coordinate system parameters in
C                    a DSK descriptor.
C     MAXCGR     P   Maximum DSK type 2 coarse voxel count.
C     MAXPLT     P   Maximum DSK type 2 plate count.
C     MAXVOX     P   Maximum DSK type 2 voxel count.
C     MAXVRT     P   Maximum DSK type 2 vertex count.
C     
C$ Detailed_Input
C
C     HANDLE      is the DAS file handle associated with a DSK file.
C                 The file must be open for write access.
C
C     CENTER      is the ID code of the body whose surface is described
C                 by the input plate model. CENTER refers to an
C                 ephemeris object.
C
C     SURFID      is the ID code of the surface described by the input
C                 plate model. Multiple surfaces (for example, surfaces
C                 having different resolutions) may be associated with
C                 a given body.
C
C     DCLASS      is the data class of the input data set. See the 
C                 INCLUDE file dskdsc.inc for values and meanings.
C
C     FRAME       is the name of the reference frame with respect to
C                 which the input data are expressed.
C
C     CORSYS      is the coordinate system in which the spatial coverage
C                 of the input data is expressed. CORSYS is an integer
C                 code. The supported values of CORPAR are
C
C                    Parameter name      Coordinate system
C                    --------------      -----------------
C                    LATSYS              Planetocentric latitudinal
C                    RECSYS              Rectangular (Cartesian)
C                    PDTSYS              Planetodetic
C                  
C                 See the INCLUDE file dskdsc.inc for parameter 
C                 declarations.
C              
C          
C     CORPAR      is an array of parameters associated with the input
C                 coordinate system. 
C
C                 For latitudinal and rectangular coordinates, CORPAR
C                 is ignored.
C                
C                 For planetodetic coordinates, the contents of CORPAR
C                 are:
C
C                    Element         Contents
C                    ---------       -----------------------------------
C                    CORPAR(1)       Equatorial radius of reference 
C                                    spheroid.
C
C                    CORPAR(2)       Flattening coefficient. The polar
C                                    radius of the reference spheroid
C                                    is given by
C
C                                       CORPAR(1) * ( 1 - CORPAR(2) )
C
C                    CORPAR(3)...
C                    CORPAR(NSYPAR)  Unused.
C                             
C        
C     MNCOR1,
C     MXCOR1,
C     MNCOR2,
C     MXCOR2,
C     MNCOR3,
C     MXCOR3      are, respectively, lower and upper bounds of
C                 each of the coordinates of the input data, where the
C                 coordinate system is defined by CORSYS and CORPAR.
C                 These bounds define the region for which the segment
C                 provides data.
C                 
C                 Distance units are km. Angular units are radians.
C
C                 The interpretation of these bounds depends on the data
C                 class; see DCLASS above.
C
C                    Single-valued surfaces
C                    ----------------------
C
C                    If the segment has data class SVFCLS (see
C                    dskdsc.inc), the segment defines a surface as a
C                    single-valued function of its domain coordinates:
C                    for example, it may define the radius of the
C                    surface as a function of planetocentric longitude
C                    and latitude, or Z as a function of X and Y.
C
C                    In this case, the input data must cover a
C                    rectangle in dimensions 1 and 2 of the input
C                    coordinate system: the set of points
C
C                       R = { (x,y): MNCOR1 < x < MXCOR1; 
C                                    MNCOR2 < y < MXCOR2  }
C
C                    must be completely covered by the input data. In
C                    other words, for each point (x,y) of R, there must
C                    be some plate containing a point whose first two
C                    coordinates are (x,y).
C
C                    The plate set may extend beyond the coordinate
C                    range defined by the bounds on the domain.
C
C                    Normally for single-valued surfaces, MNCOR3 and
C                    MXCOR3 are the minimum and maximum values of the 
C                    function attained on the domain.
C     
C                    
C                    General surfaces
C                    ----------------
C
C                    If the segment has data class GENCLS (see
C                    dskdsc.inc), the segment simply contains a
C                    collection of plates: no guarantees are made about
C                    the topology of the surface. The coordinate bounds
C                    simply indicate the spatial region for which the
C                    segment provides data.
C
C                    Note that shapes of small bodies such as asteroids
C                    and comet nuclei may fall into the "general
C                    surface" category. Surface features such as cliffs,
C                    caves, and arches can prevent a surface from being
C                    represented as a single-valued function of latitude
C                    and longitude, for example.
C
C
C                 Longitude interpretation and restrictions
C                 -----------------------------------------
C
C                 The following rules apply to longitudes provided in
C                 the arguments
C
C                    MNCOR1
C                    MXCOR1
C
C                 All angles have units of radians. The tolerance
C                 ANGMRG is used for the comparisons shown below.
C                
C                    1) Longitudes must be in the range 
C
C                          -2*pi  :  2*pi 
C
C                       Values that are out of range by no more than
C                       ANGMRG are bracketed to be in range.
C
C
C                    2) It is acceptable for the longitude bounds to be
C                       out of order. If
C
C                          MXCOR1 < MNCOR1
C                                 
C                       then either MXCOR1 is treated by the DSK
C                       subsystem as though it were MXCOR1 + 2*pi, or
C                       MNCOR1 is treated as MNCOR1 - 2*pi: whichever
C                       shift puts the bounds in the allowed range is
C                       made.
C
C                       The input longitude bounds must not be equal.
C                       If the lower bound is greater than the upper
C                       bound, the difference between the bounds must
C                       not be an integer multiple of 2*pi.
C
C                       Aside from any small changes made to move the
C                       input values of MNCOR1 or MXCOR1 into range,
C                       the values are stored in the DSK segment as is.
C
C
C                    3) MXCOR1 must not exceed MNCOR1 by more than 2*pi.
C                       Values that are out of range by no more than
C                       ANGMRG are bracketed to be in range.
C
C
C     FIRST,
C     LAST        are the endpoints of the time interval over which
C                 this data set is applicable. These times are 
C                 expressed as seconds past J2000 TDB.
C
C     NV          is the number of vertices belonging to the plate
C                 model.
C
C     VRTCES      is an array of coordinates of the vertices.
C                 The Ith vertex occupies elements (1:3,I) of
C                 this array.
C
C     NP          is the number of plates in the plate model.
C
C     PLATES      is an array representing the plates of the model.
C                 The elements of PLATES are vertex indices. The vertex
C                 indices of the Ith plate occupy elements (1:3,I) of
C                 this array.
C
C     SPAIXD,
C     SPAIXI      are, respectively, the double precision and integer
C                 components of the spatial index of the segment.
C
C                 It is strongly recommended that the helper routine
C                 DSKMI2 be used to create these arrays. See the 
C                 Examples section below.
C
C     
C$ Detailed_Output
C
C     None. This routine operates by side effects.
C 
C$ Parameters
C
C     See the SPICELIB include files
C
C        dsk02.inc
C        dskdsc.inc
C        dsktol.inc
C
C     for declarations and detailed descriptions of the parameters
C     referenced in this header.
C
C$ Exceptions
C
C     1)  If the reference frame name FRAME could not be mapped to
C         an ID code, the error SPICE(FRAMEIDNOTFOUND) is signaled.
C
C     2)  If the segment stop time precedes the start time, the
C         error SPICE(TIMESOUTOFORDER) is signaled.
C
C     3)  If an input longitude value is outside the range
C         
C            [ -2*pi - ANGMRG,   2*pi + ANGMRG ]
C
C         the error SPICE(VALUEOUTOFRANGE) will be signaled. Longitudes
C         outside of the range by a smaller amount than ANGMRG will be
C         truncated to lie in the interval [-2*pi, 2*pi].
C
C     4)  If the absolute value of the difference between the input
C         maximum longitude and the minimum longitude is more than 2*pi
C         + ANGMRG, the error SPICE(INVALIDLONEXTENT) will be signaled.
C         If either longitude bound exceeds the other by an amount
C         between 2*pi and 2*pi+ANGMRG, the larger value will be
C         truncated to the smaller value plus 2*pi.
C
C     5)  If an input latitude value is outside the range
C         
C            [ -pi/2 - ANGMRG,   pi/2 + ANGMRG ]
C
C         the error SPICE(VALUEOUTOFRANGE) will be signaled. Latitudes
C         outside of the range by a smaller amount than ANGMRG will be
C         truncated to lie in the interval [-pi/2, pi/2].
C    
C     6)  If the coordinate system is latitudinal and the lower radius
C         bound is negative, or if the upper radius bound is
C         non-positive, the error SPICE(VALUEOUTOFRANGE) will be
C         signaled.
C
C     7)  If the coordinate system is latitudinal or planetodetic and
C         the bounds of the latitude, radius or altitude coordinate are
C         out of order, the error SPICE(BOUNDSOUTOFORDER) will be
C         signaled.
C
C     8)  If the coordinate system is latitudinal or planetodetic and
C         the lower and upper bounds of the longitude, latitude, radius
C         or altitude coordinate, respectively, are equal, the error
C         SPICE(ZEROBOUNDSEXTENT) will be signaled. If the lower
C         longitude bound is greater than the upper bound, and if the
C         difference between the bounds is an integer multiple of 2*pi,
C         the same error will be signaled.
C
C     9)  If the coordinate system is planetodetic and the input
C         equatorial radius is non-positive, the error
C         SPICE(VALUEOUTOFRANGE) will be signaled.
C
C    10)  If the coordinate system is planetodetic and the input
C         flattening coefficient is greater than or equal to 1, the
C         error SPICE(VALUEOUTOFRANGE) will be signaled.
C
C    11)  If the coordinate system is planetodetic, and if the minimum
C         altitude is less than the maximum of
C
C                    2           2
C              {  -(B / A),   -(A / B)  }
C
C         where A and B are the semi-major and semi-minor axis lengths
C         of the reference ellipsoid, the error
C         SPICE(DEGENERATESURFACE) will be signaled.
C
C    12)  If the coordinate system is rectangular and any coordinate
C         lower bound is greater than or equal to the corresponding
C         upper bound, the error SPICE(BOUNDSOUTOFORDER) will be
C         signaled.
C
C    13)  If the coordinate system code is not recognized, the error
C         SPICE(NOTSUPPORTED) will be signaled.
C
C    14)  If any vertex index belonging to an input plate is outside
C         of the range 1:NV, the error SPICE(BADVERTEXINDEX) will be
C         signaled.
C
C    15)  If NV is less than 1 or greater than MAXVRT, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C    16)  If NP is less than 1 or greater than MAXPLT, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C    17)  If any voxel grid extent is less than 1 or greater than
C         MAXVOX, the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C    18)  If the voxel count is greater than MAXVOX, the error
C         SPICE(VALUEOUTOFRANGE) is signaled.
C
C    19)  If the coarse voxel count is less than 1 or greater than
C         MAXCGR, the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C    20)  If the coarse voxel scale is less than 1 or more than
C         the cube root of the fine voxel count, the error
C         SPICE(VALUEOUTOFRANGE) will be signaled.
C
C    21)  If the cube of the coarse voxel scale does not divide the
C         fine voxel count evenly, the error SPICE(INCOMPATIBLESCALE)
C         will be signaled.
C
C    22)  If the input data class is not recognized, the error
C         SPICE(NOTSUPPORTED) is signaled.
C
C    23)  If an error occurs while writing the segment to the output
C         DSK file, the error will be diagnosed by routines in the call
C         tree of this routine.
C
C$ Files
C
C     See argument HANDLE.
C
C$ Particulars
C
C     This routine writes a type 2 segment to a DSK file that 
C     has been opened for write access.
C
C     Users planning to create DSK files should consider whether the 
C     SPICE DSK creation utility MKDSK may be suitable for their needs.
C
C     This routine is supported by the routines DSKMI2 and DSKRB2
C     DSKMI2 simplifies use of this routine by creating the "spatial
C     index" arrays required as inputs by this routine. DSKRB2 computes
C     bounds on the third coordinate of the input plate set.
C
C     Spatial Indexes
C     ===============
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
C-    SPICELIB Version 1.0.0, 04-MAR-2017 (NJB)
C
C        Fixed some comment typos.
C
C     10-OCT-2016 (NJB)
C
C        New error checks on inputs were added.
C
C     07-MAR-2016 (NJB)
C
C        New error checks on inputs were added.
C
C        Argument list change: spatial index is now passed in
C        as two arrays: SPAIXD and SPAIXI.
C
C        Argument CORPAR was added.
C
C        Double precision data are now written to the output
C        segment before integer data.
C
C        22-AUG-2012 (NJB)
C
C           Bug fix: corrected upper bound in test for
C           vertex count.
C     
C        13-MAY-2010 (NJB)
C
C           Updated to reflect new type 2 segment design: normal
C           vectors, plate centers, and lengths of longest plate sides
C           are no longer stored in these segments.
C     
C        03-APR-2010 (NJB)
C
C           New interface; general coordinates are supported. Time
C           bounds, surface ID, data class, and bounds of third
C           coordinate have been added. Albedo inputs have been
C           deleted.
C
C        09-OCT-2009 (NJB)
C
C           Header was added.
C
C        31-OCT-2006 (NJB)
C
C           Input arguments CGSCAL and VTXBDS were removed. 
C
C        27-JUN-2006 (NJB)
C
C           Initial version.
C-&
 
C$ Index_Entries
C
C     write a type 2 dsk segment
C
C-&
 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPR
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      TWOPI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C


C
C     Local variables
C
      DOUBLE PRECISION      A
      DOUBLE PRECISION      ALTLIM
      DOUBLE PRECISION      B
      DOUBLE PRECISION      DESCR  ( DSKDSZ )
      DOUBLE PRECISION      R
      DOUBLE PRECISION      SEGBDS ( 2, 2 )
      DOUBLE PRECISION      VOXORI ( 3 )
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VTXBDS ( 2, 3 )

      INTEGER               CGRSCL
      INTEGER               FRMCDE
      INTEGER               I
      INTEGER               J
      INTEGER               K
      INTEGER               NCGR
      INTEGER               NVXTOT
      INTEGER               PVOXPL
      INTEGER               PVOXPT
      INTEGER               PVTXPL
      INTEGER               PVTXPT
      INTEGER               Q
      INTEGER               VGREXT ( 3 )
      INTEGER               VOXNPT
      INTEGER               VOXNPL
      INTEGER               VTXNPL
      

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKW02' )

C
C     Map the input reference frame name to an ID code.
C
      CALL NAMFRM( FRAME, FRMCDE )

      IF ( FRMCDE .EQ. 0  ) THEN

         CALL SETMSG ( 'Input reference frame # could not '
     .   //            'be mapped to an ID code. The frame '
     .   //            'name might be misspelled, or possibly '
     .   //            'a required frame kernel was not loaded. ' )
         CALL ERRCH  ( '#',  FRAME                                )
         CALL SIGERR ( 'SPICE(FRAMEIDNOTFOUND)'                   )
         CALL CHKOUT ( 'DSKW02'                                   )
         RETURN

      END IF

C
C     Make sure the time bounds are in order.
C
      IF ( LAST .LE. FIRST ) THEN

         CALL SETMSG ( 'Segment time bounds must be increasing; '
     .   //            'bounds were #:#.'                         )
         CALL ERRDP  ( '#',  FIRST                                )
         CALL ERRDP  ( '#',  LAST                                 )
         CALL SIGERR ( 'SPICE(TIMESOUTOFORDER)'                   )
         CALL CHKOUT ( 'DSKW02'                                   )
         RETURN

      END IF

C
C     If applicable, check segment boundaries. Check the 
C     coordinate system as well.
C
      IF ( ( CORSYS .EQ. LATSYS ) .OR. ( CORSYS .EQ. PDTSYS )  ) THEN
C
C        Reject invalid latitudes and longitudes. Move
C        values that are slightly out of range into range.
C
C        Longitude bounds must be distinct.
C
         IF ( MNCOR1 .EQ. MXCOR1 ) THEN
 
            CALL SETMSG ( 'Minimum longitude # radians (# degrees) '
     .      //            'was equal to maximum longitude. Longitude '
     .      //            'bounds must be distinct.'                 )
            CALL ERRDP  ( '#', MNCOR1                                )
            CALL ERRDP  ( '#', MNCOR1 * DPR()                        )
            CALL SIGERR ( 'SPICE(ZEROBOUNDSEXTENT)'                  )
            CALL CHKOUT ( 'DSKW02'                                   )
            RETURN

         END IF

C
C        Check minimum longitude.
C
         IF (      ( MNCOR1 .LT. -TWOPI()-ANGMRG ) 
     .        .OR. ( MNCOR1 .GT.  TWOPI()-ANGMRG ) ) THEN

            CALL SETMSG ( 'Minimum longitude # radians (# degrees) '
     .      //            'was outside of valid range '
     .      //            '[-2*pi, 2*pi - ANGMRG]'                   )
            CALL ERRDP  ( '#', MNCOR1                                )
            CALL ERRDP  ( '#', MNCOR1 * DPR()                        )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                   )
            CALL CHKOUT ( 'DSKW02'                                   )
            RETURN

         END IF 

C
C        The minimum longitude is too small by ANGMRG, at worst.
C        Make it greater than or equal to -2*pi.
C
         SEGBDS(1,1) = MAX( -TWOPI(), MNCOR1 )

C
C        Check maximum longitude.
C
         IF (      ( MXCOR1 .LT. -TWOPI()+ANGMRG ) 
     .        .OR. ( MXCOR1 .GT.  TWOPI()+ANGMRG ) ) THEN

            CALL SETMSG ( 'Maximum longitude # radians (# degrees) '
     .      //            'was outside of valid range '
     .      //            '[-2*pi+ANGMRG, 2*pi]'                     )
            CALL ERRDP  ( '#', MXCOR1                                )
            CALL ERRDP  ( '#', MXCOR1 * DPR()                        )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                   )
            CALL CHKOUT ( 'DSKW02'                                   )
            RETURN

         END IF 

C
C        The maximum longitude is too large by ANGMRG, at worst.
C        Make it less than or equal to 2*pi.
C
         SEGBDS(2,1) = MIN( TWOPI(), MXCOR1 )

C
C        The longitude extent cannot exceed 2*pi.
C
         IF ( MXCOR1 .GT. (MNCOR1 + TWOPI() + ANGMRG) ) THEN

            CALL SETMSG ( 'Longitude bounds #:# radians '
     .      //            '(#:# degrees) are too far apart.' )
            CALL ERRDP  ( '#', MXCOR1                        )
            CALL ERRDP  ( '#', MXCOR2                        )
            CALL ERRDP  ( '#', MXCOR1 * DPR()                )
            CALL ERRDP  ( '#', MXCOR2 * DPR()                )
            CALL SIGERR ( 'SPICE(INVALIDLONEXTENT)'          )
            CALL CHKOUT ( 'DSKW02'                           )
            RETURN

         END IF

         IF ( MXCOR1 .LT. (MNCOR1 - TWOPI() - ANGMRG) ) THEN

            CALL SETMSG ( 'Longitude bounds #:# radians '
     .      //            '(#:# degrees) are too far apart.' )
            CALL ERRDP  ( '#', MXCOR1                        )
            CALL ERRDP  ( '#', MXCOR2                        )
            CALL ERRDP  ( '#', MXCOR1 * DPR()                )
            CALL ERRDP  ( '#', MXCOR2 * DPR()                )
            CALL SIGERR ( 'SPICE(INVALIDLONEXTENT)'          )
            CALL CHKOUT ( 'DSKW02'                           )
            RETURN

         END IF

         
         IF ( SEGBDS(2,1) .GT. SEGBDS(1,1) ) THEN
C
C           The upper bound exceeds the lower by at most 2*pi + ANGMRG.
C           Trim the upper bound to make the difference no more than
C           2*pi.
C         
            SEGBDS(2,1) = MIN( SEGBDS(2,1), SEGBDS(1,1)+TWOPI() )

         ELSE IF ( SEGBDS(2,1) .LT. SEGBDS(1,1) ) THEN
C
C           The lower bound exceeds the upper by at most 2*pi + ANGMRG.
C           Advance the upper bound, if necessary, to make the
C           difference no more than 2*pi.
C         
            SEGBDS(2,1) = MAX( SEGBDS(2,1), SEGBDS(1,1)-TWOPI() )

         END IF

C
C        Make sure the adjusted longitude bounds don't describe an
C        interval that could be interpreted as having length zero,
C        if the bounds were placed in order. If the lower bound is
C        greater than the upper bound, then the difference between
C        the bounds must not be an integer multiple of 2*pi.
C
         IF (      ( SEGBDS(2,1) .EQ. SEGBDS(1,1)           ) 
     .        .OR. ( SEGBDS(2,1) .EQ. SEGBDS(1,1) - TWOPI() )  ) THEN

            CALL SETMSG ( 'After adjustment, minimum longitude '
     .      //            '# radians (# degrees) was equal to '
     .      //            'maximum longitude. Longitude bounds must '
     .      //            'be distinct.'                             )
            CALL ERRDP  ( '#', SEGBDS(1,1)                           )
            CALL ERRDP  ( '#', MNCOR1 * DPR()                        )
            CALL SIGERR ( 'SPICE(ZEROBOUNDSEXTENT)'                  )
            CALL CHKOUT ( 'DSKW02'                                   )
            RETURN

         END IF

C
C        Check minimum latitude.
C
         IF (      ( MNCOR2 .LT. -HALFPI()-ANGMRG ) 
     .        .OR. ( MNCOR2 .GT.  HALFPI()-ANGMRG ) ) THEN

            CALL SETMSG ( 'Minimum latitude # radians (# degrees) '
     .      //            'was outside of valid range '
     .      //            '[-pi/2, pi/2 - ANGMRG]'                   )
            CALL ERRDP  ( '#', MNCOR2                                )
            CALL ERRDP  ( '#', MNCOR2 * DPR()                        )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                   )
            CALL CHKOUT ( 'DSKW02'                                   )
            RETURN

         END IF 

C
C        Trim the lower latitude bound to make it at least -pi/2.
C
         SEGBDS(1,2) = MAX( -HALFPI(), MNCOR2 )

C
C        Check maximum latitude.
C
         IF (      ( MXCOR2 .LT. -HALFPI()+ANGMRG ) 
     .        .OR. ( MXCOR2 .GT.  HALFPI()+ANGMRG ) ) THEN

            CALL SETMSG ( 'Maximum latitude # radians (# degrees) '
     .      //            'was outside of valid range '
     .      //            '[-pi/2+ANGMRG, pi/2]'                     )
            CALL ERRDP  ( '#', MXCOR2                                )
            CALL ERRDP  ( '#', MXCOR2 * DPR()                        )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                   )
            CALL CHKOUT ( 'DSKW02'                                   )
            RETURN

         END IF 

C
C        Trim the upper latitude bound to make it no more than -pi/2.
C
         SEGBDS(1,2) = MAX( -HALFPI(), MNCOR2 )
         SEGBDS(2,2) = MIN(  HALFPI(), MXCOR2 )

C
C        The latitude bounds must be in order.
C        
         IF ( MXCOR2 .LT. MNCOR2 ) THEN

            CALL SETMSG ( 'Latitude bounds # and # are out of order.' )
            CALL ERRDP  ( '#',  MNCOR2                                )
            CALL ERRDP  ( '#',  MXCOR2                                )
            CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)'                   )
            CALL CHKOUT ( 'DSKW02'                                    )
            RETURN

         END IF

         
         IF ( CORSYS .EQ. LATSYS ) THEN
C
C           The coordinate system is latitudinal. Check radius
C           bounds.
C
            IF ( MNCOR3 .LT. 0.D0 ) THEN

               CALL SETMSG ( 'Radius lower bound must be '
     .         //            'non-negative but was #.'     )
               CALL ERRDP  ( '#',  MNCOR3                  )
               CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'      )
               CALL CHKOUT ( 'DSKW02'                      )
               RETURN

            END IF

            IF ( MXCOR3 .LE. 0.D0 ) THEN

               CALL SETMSG ( 'Radius upper bound must be '
     .         //            'strictly positive but was #.'  )
               CALL ERRDP  ( '#',  MXCOR3                    )
               CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'        )
               CALL CHKOUT ( 'DSKW02'                        )
               RETURN

            END IF

         END IF


         IF ( CORSYS .EQ. PDTSYS ) THEN
C          
C           The coordinate system is planetodetic. Check the coordinate
C           parameters as well.
C
            IF ( CORPAR(1) .LE. 0.D0 ) THEN

               CALL SETMSG ( 'Equatorial radius was #; this radius '
     .         //            'must be strictly positive.'            )
               CALL ERRDP  ( '#', CORPAR(1)                          )
               CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                )
               CALL CHKOUT ( 'DSKW02'                                )
               RETURN

            END IF

          
            IF ( CORPAR(2) .GE. 1.D0 ) THEN

               CALL SETMSG ( 'Flattening coefficient was #; this '
     .         //            'value must be strictly less than 1.' )
               CALL ERRDP  ( '#', CORPAR(2)                        )
               CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'              )
               CALL CHKOUT ( 'DSKW02'                              )
               RETURN

            END IF

C
C           Make sure the surface of minimum altitude is smooth and
C           non-self-intersecting.
C
            A      = CORPAR(1)
            B      = A * ( 1 - CORPAR(2) )

            ALTLIM = MAX( -(A**2)/B, -(B**2)/A )

            IF (  MNCOR3  .LE.  ALTLIM  ) THEN

               CALL SETMSG ( 'Reference ellipsoid has semi-axis '
     .         //            'lengths # and #. The minimum '
     .         //            'altitude was #. The minimum altitude '
     .         //            'is required to be greater than the '
     .         //            'maximum of {-(A**2)/B, -(B**2)/A}, '
     .         //            'which is #.'                         )
               CALL ERRDP  ( '#', A                                )
               CALL ERRDP  ( '#', B                                )
               CALL ERRDP  ( '#', MNCOR3                           )
               CALL ERRDP  ( '#', ALTLIM                           )
               CALL SIGERR ( 'SPICE(DEGENERATESURFACE)'            )
               CALL CHKOUT ( 'DSKW02'                              )
               RETURN

            END IF

         END IF

C
C        The bounds of the third coordinate, whether radius or altitude,
C        must be in order and must have positive extent.
C
         IF ( MXCOR3 .LT. MNCOR3 ) THEN

            IF ( CORSYS .EQ. LATSYS ) THEN

               CALL SETMSG ( 'Radius bounds # and # are out of order' )
            ELSE
               CALL SETMSG ( 'Altitude bounds # and # are out of '
     .         //            'order.'                                 )
            END IF

            CALL ERRDP  ( '#',  MNCOR3              )
            CALL ERRDP  ( '#',  MXCOR3              )
            CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)' )
            CALL CHKOUT ( 'DSKW02'                  )
            RETURN

         END IF

         IF ( MXCOR3 .EQ. MNCOR3 ) THEN

            IF ( CORSYS .EQ. LATSYS ) THEN

               CALL SETMSG ( 'Radius bounds # and # must have '
     .         //            'positive extent but are equal.'  )
            ELSE
               CALL SETMSG ( 'Altitude bounds # and # must have '
     .         //            'positive extent but are equal.'  )
            END IF

            CALL ERRDP  ( '#',  MNCOR3              )
            CALL ERRDP  ( '#',  MXCOR3              )
            CALL SIGERR ( 'SPICE(ZEROBOUNDSEXTENT)' )
            CALL CHKOUT ( 'DSKW02'                  )
            RETURN

         END IF

 

      ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C        All coordinate bounds must be in strictly increasing order.
C
         IF (      ( MXCOR1 .LE. MNCOR1 )
     .        .OR. ( MXCOR2 .LE. MNCOR2 )
     .        .OR. ( MXCOR3 .LE. MNCOR3 )  ) THEN

            CALL SETMSG ( 'Rectangular coordinate bounds must be '
     .      //            'strictly increasing in each dimension. '
     .      //            'The bounds were:  X = #:#; Y = #:#; '
     .      //            'Z = #:#.'                               )
            CALL ERRDP  ( '#', MNCOR1                              )
            CALL ERRDP  ( '#', MXCOR1                              )
            CALL ERRDP  ( '#', MNCOR2                              )
            CALL ERRDP  ( '#', MXCOR2                              )
            CALL ERRDP  ( '#', MNCOR3                              )
            CALL ERRDP  ( '#', MXCOR3                              )
            CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)'                )
            CALL CHKOUT ( 'DSKW02'                                 )
            RETURN

         END IF

         SEGBDS(1,1) = MNCOR1
         SEGBDS(2,1) = MXCOR1
         SEGBDS(1,2) = MNCOR2
         SEGBDS(2,2) = MXCOR2

      ELSE

         CALL SETMSG ( 'Coordinate system code # is not recognized.' )
         CALL ERRINT ( '#', CORSYS                                   )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                         )
         CALL CHKOUT ( 'DSKW02'                                      )
         RETURN
 
      END IF

C
C     Check the data class.
C
      IF (  ( DCLASS .LT. 1 ) .OR. ( DCLASS .GT. 2 )  ) THEN

         CALL SETMSG ( 'Data class # is not recognized.' )
         CALL ERRINT ( '#', DCLASS                       )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'             )
         CALL CHKOUT ( 'DSKW02'                          )
         RETURN
         
      END IF

C
C     Check NV and NP. 
C
C     Note that we don't apply Euler's law, since the data
C     set need not represent a complete surface.
C
      IF (  ( NV .LT. 1 ) .OR. ( NV .GT. MAXVRT )  ) THEN

         CALL SETMSG ( 'Vertex count NV = #; count must ' 
     .   //            'be in the range 1:#.'            )
         CALL ERRINT ( '#',  NV                          )
         CALL ERRINT ( '#',  MAXVRT                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'          )
         CALL CHKOUT ( 'DSKW02'                          )
         RETURN
         
      END IF

      IF (  ( NP .LT. 1 ) .OR. ( NP .GT. MAXPLT )  ) THEN

         CALL SETMSG ( 'Plate count NP = #; count must ' 
     .   //            'be in the range 1:#.'            )
         CALL ERRINT ( '#',  NP                          )
         CALL ERRINT ( '#',  MAXPLT                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'          )
         CALL CHKOUT ( 'DSKW02'                          )
         RETURN
         
      END IF

C
C     Check the vertex indices in the plates.
C
      DO I = 1, NP

         DO J = 1, 3

            K = PLATES(J,I) 

            IF ( ( K .LT. 1 ) .OR. ( K .GT. NV ) ) THEN

               CALL SETMSG ( 'Vertex index # of plate # was #; '
     .         //            'vertex indices must be in the '
     .         //            'range 1:NV. The input NV = #.'   )
               CALL ERRINT ( '#', J                            )
               CALL ERRINT ( '#', I                            )
               CALL ERRINT ( '#', K                            )
               CALL ERRINT ( '#', NV                           )
               CALL SIGERR ( 'SPICE(BADVERTEXINDEX)'           )
               CALL CHKOUT ( 'DSKW02'                          )
               RETURN

            END IF

         END DO

      END DO

C
C     Locate the spatial index elements. Some of the elements are at
C     fixed addresses; for others the addresses must be calculated.
C
C     The two components of the spatial index together contain the
C     following items:
C
C        VGREXT      is an array containing the extents of the voxel
C                    grid in the X, Y, and Z directions of the
C                    body-fixed frame. The extents are measured as
C                    voxel counts.
C
C        ORIGIN      is the position, in the body-fixed, body-centered
C                    reference frame associated with BODY, of the
C                    origin of the both the fine and coarse voxel grids
C                    associated with this model.
C
C        VOXSIZ      is the voxel edge length in km.
C
C        CGRSCL      is the coarse voxel grid scale factor: the edge
C                    length of each coarse voxel is scaled up from the
C                    length of a fine voxel edge by this factor.
C
C        CGRPTR      is an array of pointers associated with this
C                    model's coarse voxel grid; these pointers map
C                    one-dimensional coarse voxel indices to start
C                    indices in the fine voxel pointer list.
C
C        VOXNPT      is the cardinality of the fine voxel pointer list.
C
C        VOXPTR      is the fine voxel pointer list. For each fine
C                    voxel belonging to a non-empty coarse voxel, there
C                    is a pointer in this list that identifies the
C                    start index in VOXPLT of the list of plate indices
C                    associated with this fine voxel.
C
C                    The start index in VOXPTR of the set of pointers
C                    associated with a coarse voxel is given by the
C                    element of CGRPTR associated with that coarse
C                    voxel.
C
C                    Within a given coarse voxel, each fine voxel has
C                    an associated one-dimensional offset from the
C                    corner of the coarse voxel closest to the origin
C                    of the voxel grids. This offset gives the location
C                    in VOXPTR of the plate list pointer for the fine
C                    voxel.
C
C        VOXNPL      is the cardinality of the plate list of the fine
C                    voxel-plate mapping.
C
C        VOXPLT      is the plate list of the fine voxel-plate mapping.
C
C        VTXPTR      is the vertex pointer list.
C
C        VTXNPL      is the cardinality of the plate list of the
C                    vertex-plate mapping.
C
C
C
C     Extract double precision elements of the spatial index.
C
      CALL MOVED ( SPAIXD(SIVTBD), 6, VTXBDS )
      CALL VEQU  ( SPAIXD(SIVXOR),    VOXORI )

      VOXSIZ = SPAIXD(SIVXSZ)       

C
C     Extract scalars and small fixed-size arrays from the integer
C     component of the spatial index.
C
C     Fetch grid extents (in units of whole voxels):
C
      CALL MOVEI ( SPAIXI(SIVGRX), 3, VGREXT )

C
C     Fetch coarse grid scale, voxel pointer count, and voxel-plate
C     list count.
C
      CGRSCL = SPAIXI(SICGSC)
      VOXNPT = SPAIXI(SIVXNP)
      VOXNPL = SPAIXI(SIVXNL)

C
C     Create a pointer to the voxel-plate pointer array.
C
      PVOXPT = SICGRD + MAXCGR

C
C     Create a pointer to the voxel-plate list array.
C
      PVOXPL = PVOXPT + VOXNPT

C
C     Create a pointer to the vertex pointer array.
C
      PVTXPT = PVOXPL + VOXNPL

C
C     Create a pointer to the vertex-plate list array.
C
      PVTXPL = PVTXPT + NV

C
C     Fetch vertex-plate list size.
C
      VTXNPL = SPAIXI(SIVTNL)

C
C     Check the input parameters.
C
C
C
C     Make sure the voxel grid extents are within range.
C
      DO I = 1, 3

         IF (      ( VGREXT(I) .LT. 1      ) 
     .        .OR. ( VGREXT(I) .GT. MAXVOX )  ) THEN

            CALL SETMSG ( 'Voxel grid extents are = (#, #, #); ' 
     .   //               'all be in the range 1:#.'             )
            CALL ERRINT ( '#',  VGREXT(1)                        )
            CALL ERRINT ( '#',  VGREXT(2)                        )
            CALL ERRINT ( '#',  VGREXT(3)                        )
            CALL ERRINT ( '#',  MAXVOX                           )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'               )
            CALL CHKOUT ( 'DSKW02'                               )
            RETURN

         END IF

      END DO
            
C
C     Make sure the number of voxels NVXTOT is within range.
C
      NVXTOT = VGREXT(1) * VGREXT(2) * VGREXT(3) 

      IF ( NVXTOT .GT. MAXVOX ) THEN

         CALL SETMSG ( 'Fine voxel count NVXTOT = #; count must ' 
     .   //            'be in the range 1:#.'                    )
         CALL ERRINT ( '#',  NVXTOT                              )
         CALL ERRINT ( '#',  MAXVOX                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                  )
         CALL CHKOUT ( 'DSKW02'                                  )
         RETURN
         
      END IF
      

C
C     Check the coarse voxel scale. It must be at least 1, and its
C     cube must not exceed the fine voxel count.
C
      IF (      ( CGRSCL .LT. 1                ) 
     .     .OR. ( CGRSCL .GT. NVXTOT**(1.D0/3) )  ) THEN

         CALL SETMSG ( 'Coarse voxel scale = #; scale must ' 
     .   //            'be in the range 1:NVXTOT**3, where '
     .   //            'NVXTOT is the total fine voxel count. '
     .   //            'In this case, NVXTOT = #.'               )
         CALL ERRINT ( '#',  CGRSCL                              )
         CALL ERRINT ( '#',  NVXTOT                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                  )
         CALL CHKOUT ( 'DSKW02'                                  )
         RETURN

      END IF

C
C     The cube of the coarse scale must divide the total voxel count
C     evenly.
C
      Q = NVXTOT /   (CGRSCL**3)
      R = NVXTOT - Q*(CGRSCL**3)

      IF ( R .NE. 0 ) THEN

         CALL SETMSG ( 'Coarse voxel scale = #; the cube of '
     .   //            'the scale must divide NVXTOT evenly, ' 
     .   //            'where NVXTOT is the total  fine voxel '
     .   //            'count. In this case, NVXTOT = #.'        )
         CALL ERRINT ( '#',  CGRSCL                              )
         CALL ERRINT ( '#',  NVXTOT                              )
         CALL SIGERR ( 'SPICE(INCOMPATIBLESCALE)'                )
         CALL CHKOUT ( 'DSKW02'                                  )
         RETURN

      END IF


C                
C     NCGR        is the number of voxels in the coarse voxel grid 
C                 associated with this model. Since each coarse voxel
C                 is a cube containing an integer number of fine
C                 voxels, this number is determined by NVXTOT and
C                 CGRSCL.
C
      NCGR = NVXTOT / ( CGRSCL**3 )

C
C     Make sure NCGR is within range.
C
      IF (  ( NCGR .LT. 1 ) .OR. ( NCGR .GT. MAXCGR )  ) THEN

         CALL SETMSG ( 'Coarse voxel count = #; count must ' 
     .   //            'be in the range 1:#.'                )
         CALL ERRINT ( '#',  NCGR                            )
         CALL ERRINT ( '#',  MAXCGR                          )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'              )
         CALL CHKOUT ( 'DSKW02'                              )
         RETURN
         
      END IF


C
C     Start a new DLA segment.
C
      CALL DLABNS ( HANDLE )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'DSKW02' )
         RETURN
      END IF
     
C
C     Write the d.p. data to the segment first. In the segregated
C     segment, d.p. data will precede integer data, so less
C     rearrangement will occur this way.
C
C
C     First, fill in the DSK segment descriptor.
C
      CALL CLEARD ( DSKDSZ, DESCR )

      DESCR ( SRFIDX ) = SURFID
      DESCR ( CTRIDX ) = CENTER
      DESCR ( CLSIDX ) = DCLASS
      DESCR ( TYPIDX ) = 2
      DESCR ( FRMIDX ) = FRMCDE
      DESCR ( SYSIDX ) = CORSYS

      CALL MOVED ( CORPAR, NSYPAR, DESCR( PARIDX ) )

      DESCR ( MN1IDX ) = SEGBDS(1,1)
      DESCR ( MX1IDX ) = SEGBDS(2,1)
      DESCR ( MN2IDX ) = SEGBDS(1,2)
      DESCR ( MX2IDX ) = SEGBDS(2,2)
      DESCR ( MN3IDX ) = MNCOR3
      DESCR ( MX3IDX ) = MXCOR3
      DESCR ( BTMIDX ) = FIRST
      DESCR ( ETMIDX ) = LAST

C
C     Now write the descriptor into the segment.
C
      CALL DASADD ( HANDLE, DSKDSZ, DESCR )
 
C
C     Add the voxel grid origin and voxel size.
C     Finish with the vertex data. 
C
      CALL DASADD ( HANDLE, 6,      VTXBDS )
      CALL DASADD ( HANDLE, 3,      VOXORI )
      CALL DASADD ( HANDLE, 1,      VOXSIZ )
      CALL DASADD ( HANDLE, 3*NV,   VRTCES )
  
C
C     Next add the integer data to the segment.
C
C     NV is the number of vertices.
C     NP is the number of plates.
C     NVXTOT is the number of voxels in the spatial index.
C
      CALL DASADI ( HANDLE, 1,      NV     )
      CALL DASADI ( HANDLE, 1,      NP     )
      CALL DASADI ( HANDLE, 1,      NVXTOT )
      CALL DASADI ( HANDLE, 3,      VGREXT )
      CALL DASADI ( HANDLE, 1,      CGRSCL )
      CALL DASADI ( HANDLE, 1,      VOXNPT )
      CALL DASADI ( HANDLE, 1,      VOXNPL )
      CALL DASADI ( HANDLE, 1,      VTXNPL )
      CALL DASADI ( HANDLE, 3*NP,   PLATES )
      CALL DASADI ( HANDLE, VOXNPT, SPAIXI(PVOXPT) )
      CALL DASADI ( HANDLE, VOXNPL, SPAIXI(PVOXPL) )
      CALL DASADI ( HANDLE, NV,     SPAIXI(PVTXPT) )
      CALL DASADI ( HANDLE, VTXNPL, SPAIXI(PVTXPL) )
      CALL DASADI ( HANDLE, NCGR,   SPAIXI(SICGRD) )

C
C     End the segment.
C
      CALL DLAENS ( HANDLE )

      CALL CHKOUT ( 'DSKW02' )
      RETURN
      END 

