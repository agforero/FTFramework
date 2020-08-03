C$Procedure DSKX02 ( DSK, ray-surface intercept, type 2 )

      SUBROUTINE DSKX02 ( HANDLE, DLADSC, VERTEX, 
     .                    RAYDIR, PLID,   XPT,    FOUND )
 
C$ Abstract
C
C     Determine the plate ID and body-fixed coordinates of the
C     intersection of a specified ray with the surface defined by a
C     type 2 DSK plate model.
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
C      None.
C
C$ Keywords
C
C      GEOMETRY
C      SHAPES
C
C$ Declarations
 
      IMPLICIT NONE
 
      INCLUDE              'dla.inc'
      INCLUDE              'dskdsc.inc'
      INCLUDE              'dsk02.inc'
      INCLUDE              'dsktol.inc'

 
      INTEGER               HANDLE
      INTEGER               DLADSC   ( * )
      DOUBLE PRECISION      VERTEX   ( 3 )
      DOUBLE PRECISION      RAYDIR   ( 3 )
      INTEGER               PLID
      DOUBLE PRECISION      XPT      ( 3 )
      LOGICAL               FOUND
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK kernel containing plate model.
C     DLADSC     I   DLA descriptor of plate model segment.
C     VERTEX     I   Ray's vertex in the  body fixed frame.
C     RAYDIR     I   Ray direction in the body fixed frame.
C     PLID       O   ID code of the plate intersected by the ray.
C     XPT        O   Intercept.
C     FOUND      O   Flag indicating whether intercept exists.
C     XFRACT     P   Plate expansion fraction.
C
C$ Detailed_Input
C
C     HANDLE     is the file handle of a DSK file containing a shape
C                model for a target body.  The shape model is stored
C                in a type 2 DSK segment.
C   
C     DLADSC     is the DLA descriptor of a type 2 DSK segment
C                containing plate model data representing the surface of
C                the target body.  The caller should declare DLADSC
C                with size DLADSZ; this size parameter is defined in
C                the INCLUDE file dla.inc.  Normally this descriptor
C                will be obtained by a search through a DSK file
C                using the DLA search routines; see the Examples
C                header section below for a working code example
C                illustrating a simple search.
C
C     VERTEX     is the vertex of a ray.  VERTEX is expressed relative
C                to the body fixed reference frame associated with the
C                target body.  This reference frame is the same frame
C                relative to which the vertices of the plate model
C                are expressed.  Units are km.
C
C                The vertex is required to be outside the target
C                body.
C
C     RAYDIR     is the ray's direction vector.  RAYDIR is expressed
C                relative to the body fixed reference frame associated
C                with the target body.
C
C$ Detailed_Output
C
C     PLID       is the ID of the plate closest to the input ray's
C                vertex at which a ray-surface intercept exists.
C                If no intercept exists, PLID is undefined.
C
C     XPT        is the ray-target intercept closest to the ray's
C                vertex, if an intercept exists.  XPT is expressed
C                relative to the body-fixed reference frame associated
C                with the target body.  Units are km.
C
C                If no intercept exists, XPT is undefined.
C
C     FOUND      is a logical flag that indicates whether or not the
C                ray does indeed intersect the target.  If the ray
C                intersects a plate FOUND is .TRUE.  Otherwise FOUND is
C                FALSE.
C
C$ Parameters
C
C
C     XFRACT     is the default plate expansion fraction. This
C                parameter can be overridden. See the include file
C
C                   dsktol.inc
C
C                for details.
C
C
C     See the include files
C
C        dla.inc
C        dsk02.inc
C        dskdsc.inc
C        dsktol.inc
C
C     for declarations of parameters used by this routine.
C
C$ Exceptions
C
C     1)  If the input handle is invalid, the error will be diagnosed by
C         routines in the call tree of this routine.
C
C     2)  If a file read error occurs, the error will be diagnosed by
C         routines in the call tree of this routine.
C
C     3)  If the input DLA descriptor is invalid, the effect of this
C         routine is undefined. The error *may* be diagnosed by routines
C         in the call  tree of this routine, but there are no
C         guarantees.
C
C     4)  If an error occurs while trying to look up any component
C         of the shape model, the error will be diagnosed by routines
C         in the call tree of this routine.
C
C     5)  If the input ray direction is the zero vector, the error
C         SPICE(ZEROVECTOR) is signaled.
C
C     6)  If the coarse voxel grid scale of the shape model is less than
C         1, the error SPICE(VALUEOUTOFRANGE) is signaled.
C    
C     7)  If the coarse voxel grid of the shape model contains more
C         than MAXCGR (see pltmax.inc) voxels, the error
C         SPICE(GRIDTOOLARGE) is signaled.
C
C     8)  If the plate list for any intersected voxel is too large
C         for this routine to buffer, the error SPICE(ARRAYTOOSMALL)
C         is signaled.
C
C     9)  Due to round-off errors, results from this routine may
C         differ across platforms.  Results also may differ from
C         those expected---and not necessarily by a small amount.
C         For example, a ray may miss a plate it was expected to
C         hit and instead hit another plate considerably farther
C         from the ray's vertex, or miss the target entirely.
C
C     10) In the event that an intercept point lies on multiple 
C         plates (that is, the point is on an edge or vertex), 
C         a plate will be selected.  Due to round-off error, the
C         selection may vary across platforms.
C
C$ Files
C
C     See the description of the input argument HANDLE.
C
C$ Particulars
C    
C     This routine solves the ray-surface intercept problem for 
C     a specified ray and a surface represented by triangular plate 
C     model.  The surface representation is provided by data in a
C     type 2 segment of a DSK file.
C
C     This routine does not assume that the segment from which the
C     surface model data are read represents the entire surface of
C     the target body.  A program could call this routine repeatedly
C     to find the surface intercept of a ray and a shape model
C     partitioned into multiple segments.
C
C     In general, this routine should be expected to run faster
C     when used with smaller shape models.  
C
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as input,
C     the compiler and supporting libraries, and the machine specific
C     arithmetic implementation.
C
C     1) Find the surface intercept points corresponding to a latitude/
C        longitude grid of a specified resolution, for a specified
C        target body.  This simple program assumes the shape model for
C        the target body is stored in a single type 2 DSK segment, and
C        that this segment is the first one in the DSK file to which it
C        belongs.
C
C        In the following example program, the file
C
C           phobos_3_3.bds
C
C        is a DSK file containing a type 2 segment that provides a
C        plate model representation of the surface of Phobos. In order
C        to duplicate the example output, this kernel name should be
C        supplied at the prompt.
C
C
C     C
C     C     PROGRAM DSKX02_EX1
C     C     IMPLICIT NONE
C     C
C           INCLUDE 'dla.inc'
C           INCLUDE 'dskdsc.inc'
C           INCLUDE 'dsk02.inc'
C     C
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      DSKSGR
C           DOUBLE PRECISION      RPD
C     C
C     C
C     C     Local parameters
C     C
C           INTEGER               FILSIZ
C           PARAMETER           ( FILSIZ = 255 )
C
C           INTEGER               NLAT
C           PARAMETER           ( NLAT   = 9 )
C
C           INTEGER               NLON
C           PARAMETER           ( NLON   = 9 )
C
C     C
C     C     Local parameters
C     C
C           DOUBLE PRECISION      TOL
C           PARAMETER           ( TOL   =  1.D-12 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(FILSIZ)    DSK
C
C           DOUBLE PRECISION      DSKDSC ( DSKDSZ )
C           DOUBLE PRECISION      LAT
C           DOUBLE PRECISION      LON
C           DOUBLE PRECISION      MAXR
C           DOUBLE PRECISION      R
C           DOUBLE PRECISION      RAYDIR ( 3 )
C           DOUBLE PRECISION      VERTEX ( 3 )
C           DOUBLE PRECISION      XLAT
C           DOUBLE PRECISION      XLON
C           DOUBLE PRECISION      XPT    ( 3 )
C           DOUBLE PRECISION      XR
C
C           INTEGER               DLADSC ( DLADSZ )
C           INTEGER               HANDLE
C           INTEGER               I
C           INTEGER               J
C           INTEGER               PLID
C
C           LOGICAL               FOUND
C
C     C
C     C     Prompt for the name of the DSK to read.
C     C
C           CALL PROMPT ( 'Enter DSK name > ', DSK )
C     C
C     C     Open the DSK file for read access.
C     C     We use the DAS-level interface for
C     C     this function.
C     C
C           CALL DASOPR ( DSK, HANDLE )
C     C
C     C     Begin a forward search through the
C     C     kernel, treating the file as a DLA.
C     C     In this example, it's a very short
C     C     search.
C     C
C           CALL DLABFS ( HANDLE, DLADSC, FOUND )
C
C           IF ( .NOT. FOUND ) THEN
C     C
C     C        We arrive here only if the kernel
C     C        contains no segments.  This is
C     C        unexpected, but we're prepared for it.
C     C
C              CALL SETMSG ( 'No segments found '
C          .   //            'in DSK file #.'    )
C              CALL ERRCH  ( '#',  DSK           )
C              CALL SIGERR ( 'SPICE(NODATA)'     )
C
C           END IF
C
C     C
C     C     If we made it this far, DLADSC is the
C     C     DLA descriptor of the first segment.
C     C
C     C     We're going to generate the intercept points
C     C     using a set of rays which point toward the
C     C     origin and whose vertices are on a 
C     C     specified lat/lon grid.  To start out we
C     C     must pick a reasonable range from the origin
C     C     for the vertices:  the range must be large
C     C     enough so that the vertices are guaranteed
C     C     to be outside the target body but small
C     C     enough that we don't lose too much precision
C     C     in the surface intercept computation.
C     C
C     C     We'll look up the upper bound for the target
C     C     radius, then use 2 times this value as the
C     C     vertex magnitude.
C     C
C           CALL DSKGD ( HANDLE, DLADSC, DSKDSC )
C
C           MAXR = DSKSGR( DSKDSC )
C           R    = 2.D0 * MAXR
C
C     C
C     C     Now generate the intercept points.  We generate
C     C     intercepts along latitude bounds, working from
C     C     north to south. Latitude ranges from +80 to
C     C     -80 degrees. Longitude ranges from 0 to 320 
C     C     degrees. The increment is 20 degrees for 
C     C     latitude and 40 degrees for longitude.
C     C
C           DO I = 1, NLAT
C
C              LAT = RPD() * ( 100.D0 - 20.D0*I )
C
C              DO J = 1, NLON
C
C                 LON = RPD() * 40.D0 * (J-1)
C     C
C     C           Produce a ray vertex for the current
C     C           lat/lon value.  Negate the vertex to
C     C           produce the ray's direction vector.
C     C
C                 CALL LATREC ( R, LON, LAT, VERTEX )
C                 CALL VMINUS ( VERTEX,      RAYDIR )
C     C
C     C           Find the surface intercept for this
C     C           ray.
C     C
C                 CALL DSKX02 ( HANDLE, DLADSC, VERTEX,
C          .                    RAYDIR, PLID,   XPT,    FOUND  )
C     C
C     C           Since the ray passes through the origin on
C     C           the body-fixed frame associated with the
C     C           target body, we'd rarely expect to find that
C     C           the ray failed to intersect the surface.
C     C           For safety, we check the FOUND flag.  (A
C     C           "not found" condition could be a sign of
C     C           a bug.)
C     C
C                 IF ( .NOT. FOUND ) THEN
C
C                    WRITE(*,*) ' '
C                    WRITE(*,*) 'Intercept not found!'
C                    WRITE(*,*) '   Ray vertex:'
C                    WRITE(*,*) '   Longitude (deg): ', LON/RPD()
C                    WRITE(*,*) '   Latitude  (deg): ', LAT/RPD()
C                    WRITE(*,*) '   Range      (km): ', R
C                    WRITE(*,*) ' '
C
C                 ELSE
C     C
C     C              This is the normal case.  Display the
C     C              intercept plate ID and the intercept
C     C              point in both Cartesian and latitudinal
C     C              coordinates.  Show the corresponding ray
C     C              vertex to facilitate validation of results.
C     C
C     C              Use RECRAD rather than RECLAT to produce
C     C              non-negative longitudes.
C     C
C                    CALL RECRAD ( XPT, XR, XLON, XLAT )
C
C                    WRITE(*,*) ' '
C                    WRITE(*,*) 'Intercept found:'
C                    WRITE(*,*) '   Plate ID:             ', PLID
C                    WRITE(*, '(1X,A,3E14.6)' )
C          .         '   Cartesian coordinates:', XPT
C                    WRITE(*,*) '   Latitudinal coordinates:'
C                    WRITE(*,*) '   Longitude (deg): ', XLON/RPD()
C                    WRITE(*,*) '   Latitude  (deg): ', XLAT/RPD()
C                    WRITE(*,*) '   Range      (km): ', XR
C                    WRITE(*,*)
C                    WRITE(*,*) '   Ray vertex:'
C                    WRITE(*,*) '   Longitude (deg): ', LON/RPD()
C                    WRITE(*,*) '   Latitude  (deg): ', LAT/RPD()
C                    WRITE(*,*) '   Range      (km): ', R
C                    WRITE(*,*) ' '
C
C                 END IF
C
C              END DO
C
C           END DO
C     C
C     C     Close the kernel.  This isn't necessary in a stand-
C     C     alone program, but it's good practice in subroutines
C     C     because it frees program and system resources.
C     C
C           CALL DASCLS ( HANDLE )
C           END
C
C
C
C     When this program was executed on a PC/Linux/gfortran platform,
C     the output for the first 3 points (the rest of the output
C     is not shown due to its large volume) was:
C
C
C  Enter DSK name    > phobos_3_3.bds
C
C    Intercept found:
C       Plate ID:                   306940
C       Cartesian coordinates:  0.152088E+01  0.000000E+00  0.862533E+01
C       Latitudinal coordinates:
C       Longitude (deg):    0.0000000000000000
C       Latitude  (deg):    79.999999999999986
C       Range      (km):    8.7583866856211507
C
C       Ray vertex:
C       Longitude (deg):    0.0000000000000000
C       Latitude  (deg):    80.000000000000000
C       Range      (km):    28.023536291251524
C
C
C    Intercept found:
C       Plate ID:                   317112
C       Cartesian coordinates:  0.118970E+01  0.998280E+00  0.880777E+01
C       Latitudinal coordinates:
C       Longitude (deg):    40.000000000000092
C       Latitude  (deg):    80.000000000000000
C       Range      (km):    8.9436459265318593
C
C       Ray vertex:
C       Longitude (deg):    40.000000000000000
C       Latitude  (deg):    80.000000000000000
C       Range      (km):    28.023536291251524
C
C
C    Intercept found:
C       Plate ID:                   324141
C       Cartesian coordinates:  0.277775E+00  0.157534E+01  0.907203E+01
C       Latitudinal coordinates:
C       Longitude (deg):    80.000000000000043
C       Latitude  (deg):    79.999999999999986
C       Range      (km):    9.2119797003191035
C
C       Ray vertex:
C       Longitude (deg):    80.000000000000000
C       Latitude  (deg):    80.000000000000000
C       Range      (km):    28.023536291251524
C
C
C$ Restrictions
C
C     1) This is prototype code.  When this code is graduated into
C        SPICELIB, the functionality or interface could change.
C
C$ Literature_References
C
C     Woo, A. "Fast Ray-Box Intersection", Graphic Gems I, 395-396.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     J.A. Bytof      (JPL)
C     W.L. Taber      (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0 04-APR-2017 (NJB) 
C
C        Added test for containment of intersection point
C        within segment boundaries. Updated logic for saving
C        segment attributes so that errors won't cause saved
C        values to be improperly re-used on subsequent calls.
C
C     24-FEB-2016 (NJB) 
C
C        Replaced call to TOGRID with call to ZZTOGRID.
C        Replaced call to PLTREC with call to PLTNRM.
C        Now obtains plate expansion fraction from DSKGTL.
C
C     25-FEB-2015 (NJB) 
C
C        Bug fix: now ray-voxel grid intercept is displaced toward
C        the input ray's vertex only when the vertex is outside
C        the target body's voxel grid.
C
C     10-SEP-2014 (NJB) 
C
C        Bug fix: during an intercept search over the voxel list
C        returned by XDDA, if an intercept outside the current
C        voxel---by more than a small tolerance---is found, the search
C        rejects that intercept and continues until a valid intercept
C        is found and all plates in the voxel containing that intercept
C        have been checked for an intersection. The rejected intercept
C        may later be determined to be a valid solution during a check
C        of plates associated with a voxel that contains that
C        intercept; in fact it is the correct solution if no other
C        plates contain a solution closer to the ray's vertex. (Usually
C        there is a unique voxel containing the intercept, but this is
C        not so if the intercept lies on a voxel boundary not on an
C        edge of the voxel grid.)
C
C        Note that there's no need to look for intersections in voxels
C        further out in the voxel list than the first one that contains
C        an intercept.
C
C        The previous version of the routine terminated the search
C        after checking all plates in the current voxel, after an
C        intercept was found on any plate associated with that voxel.
C        The intercept was not required to be contained in the voxel.
C
C        See the Revisions section for details.
C
C     30-JUN-2014 (NJB) 
C
C        Bug fix: renamed "found" flag returned by ZZRAYBOX.
C
C        Added code to test for empty voxel list after
C        voxel list compression.
C
C        Last update was 15-JUN-2014 (NJB).
C
C        Made some minor edits to in-line comments, and removed
C        comments that had become inapplicable due to code changes.
C
C        Last update was 06-JUN-2014.
C
C        Now expands plates slightly before performing ray-plate
C        intersection computations.
C
C        Bug fix: now calls ZZRAYBOX to find the ray-box intercept.
C        This reduces round-off error in the variable COORD.
C
C-       Last update was 02-MAY-2014 (NJB) 
C
C        Bug fix: added FAILED checks after each DSKI02 and DSKD02
C        call.
C
C        Some precautionary measures were added: a backstop check for
C        an empty voxel list was added after the XDDA call. A backstop
C        initialization of PNTR was added before the plate collection
C        loop. This initialization is needed only if the voxel list
C        returned by XDDA is empty. The list should never be empty.
C       
C        Last update was 2.2.1 25-MAR-2014 (NJB) 
C
C        Bug fix: duplicate plates are now marked so that the
C        unmarked instance is the one in the closest voxel to
C        the ray's origin.
C
C        Bug fix: corrected buffer overflow error detection for
C        insertion of plate IDs into plate ID array.
C
C     20-JUL-2011 (NJB) 
C
C        Bug fix: this routine now tests FAILED after its
C        call to XDDA.
C
C        Header correction: the detailed input section
C        now says that the ray's vertex *is* required to
C        be outside the target body.
C
C     09-JUN-2011 (NJB) 
C
C        All large local arrays are now saved in order to support
C        calling a C translation of this routine from Java.
C
C        The buffer VIDXS is now initialized prior to its
C        first use in an argument list. This was done to 
C        to suppress compiler warnings. The original code was
C        correct, since along with the buffer, an array size
C        of zero was passed to the called function.
C
C        The example program was updated for compatibility with the
C        final DSK descriptor layout. The output format was adjusted.
C        Sample output from the program is now shown in the header.
C
C     13-MAY-2010 (NJB) 
C
C        No longer uses plate records to weed out 
C        plates prior to ray-plate intercept tests.
C        Now uses local vertex buffer. Logic for choosing
C        plate when intercept is on edge or vertex has
C        been simplified.
C
C     06-MAY-2010 (NJB) 
C
C        Now calls DSKB02 rather than DSKP02.
C
C     20-APR-2010 (NJB) 
C
C        Updated header section order.
C
C     26-SEP-2008 (NJB) 
C
C        Moved OBSMAT computation out of loop.
C
C     27-DEC-2006 (NJB) (EDW)
C
C        Header example was updated to show maximum radius
C        being obtained from DSK descriptor rather than via
C        all to DSKD02.
C
C     31-OCT-2006 (NJB) (EDW)
C
C        Modified to work with DLA-based kernels.  Many
C        changes were made to the algorithm to improve
C        execution speed.
C
C     19-AUG-2004 (EDW)
C
C        Implemented "Fast Ray-Box Intersection" algorithm.
C        Renamed routine DSKX02 from PLBORE_3.
C
C     25-FEB-1999 (JAB)
C
C        Based on PLBORE and PLBORE_2.
C
C-&
 
C$ Index_Entries
C
C     plate and plate model point intersected by ray
C
C-&

C$ Revisions
C
C-    10-SEP-2014 (NJB) 
C
C        Bug fix: during an intercept search over the voxel list
C        returned by XDDA, if an intercept outside the current
C        voxel---by more than a small tolerance---is found, the search
C        rejects that intercept and continues until a valid intercept
C        is found and all plates in the voxel containing that intercept
C        have been checked for an intersection. The rejected intercept
C        may later be determined to be a valid solution during a check
C        of plates associated with a voxel that contains that
C        intercept; in fact it is the correct solution if no other
C        plates contain a solution closer to the ray's vertex. (Usually
C        there is a unique voxel containing the intercept, but this is
C        not so if the intercept lies on a voxel boundary not on an
C        edge of the voxel grid.)
C
C        Note that there's no need to look for intersections in voxels
C        further out in the voxel list than the first one that contains
C        an intercept.
C
C        The previous version of the routine terminated the search
C        after checking all plates in the current voxel, after an
C        intercept was found on any plate associated with that voxel.
C        The intercept was not required to be contained in the voxel.
C 
C        In the previous version of the routine, an intercept found
C        outside of the current voxel could effectively mask an
C        intercept closer to the ray's vertex, as shown in the diagram
C        below.
C
C        In this diagram, "V" represents the vertex of the ray.
C        The letter sequences "Q*" and "P*" represent plates. 
C        Here the ray hits voxel 1 and finds an intercept on plate
C        P* at the point marked by "@." No other intercepts on
C        plates in voxel 1 exist, so the search terminates. The
C        intercept marked by "$" is closer to the vertex but is not
C        seen. 
C
C                                V
C                               /
C                              /
C                             /
C           +--------------+-/------------+
C           |    voxel 2   |/    voxel 1  |
C           |              /              |
C           |        QQQQQ$|              |
C           |            / |              |
C           |           /  |              |
C           |          /   |              |
C           |       PP@PPPPPPPPPPPPPPP    |
C           +--------------+--------------+
C
C
C       The updated algorithm, when presented with the situation shown
C       above, will check all plates in voxel 2 before terminating.
C
C       Note that the problem could occur in cases where voxels 1 and 2
C       are not adjacent.
C
C-& 






C
C     SPICELIB Functions
C 
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      VDOT          

      LOGICAL               RETURN
      LOGICAL               FAILED
      LOGICAL               VZERO
 
C
C     Statement function type declarations
C     
      INTEGER               VOX2ID 
      

C
C     Local parameters
C


C
C     Tolerance used for vertex-voxel grid distance test:
C
      DOUBLE PRECISION      VTXTOL
      PARAMETER           ( VTXTOL = 1.D-12 )

C
C     Maximum number of voxels we can accept for
C     one XDDA call.
C
      INTEGER               NVXLST
      PARAMETER           ( NVXLST = 50000 )
          
C
C     Maximum number of plates we work with
C     at any time.
C
      INTEGER               MAXPL
      PARAMETER           ( MAXPL  =  MAXPLT/125 )

C
C     Parameter indicating no coordinates are to be excluded
C     in the test for a point being within segment boundaries.
C
      INTEGER               NONE
      PARAMETER           ( NONE   = 0 )

C
C     Local Variables
C
      INTEGER               ISRCHI

      INTEGER               BUFSIZ
      PARAMETER           ( BUFSIZ = 200 )


      DOUBLE PRECISION      COORD  ( 3 )
      DOUBLE PRECISION      EDGES  ( 3, 3 )
      DOUBLE PRECISION      DSKDSC ( DSKDSZ )
      DOUBLE PRECISION      GRDEXT ( 3 )
      DOUBLE PRECISION      GRDTOL
      DOUBLE PRECISION      GREEDM
      DOUBLE PRECISION      HITCOR ( 3 )
      DOUBLE PRECISION      NEAR
      DOUBLE PRECISION      NORMAL ( 3 )
      DOUBLE PRECISION      OBSMAT ( 3, 3 )
      DOUBLE PRECISION      POINTS ( 3, 3 )
      DOUBLE PRECISION      SCALE
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      VBUFF  ( 3, BUFSIZ )
      DOUBLE PRECISION      VOXORI ( 3 )
      DOUBLE PRECISION      VOXSIZ
      DOUBLE PRECISION      VTX2   ( 3 )
      DOUBLE PRECISION      VTXBDS ( 2, 3 )
      DOUBLE PRECISION      VTXOFF ( 3 )
      DOUBLE PRECISION      XPDFRC
      DOUBLE PRECISION      XPNTS  ( 3, 3 )
      DOUBLE PRECISION      XTOL

      INTEGER               CGREXT ( 3 )
      INTEGER               CGRPTR ( MAXCGR )
      INTEGER               CGSCAL
      INTEGER               CGSCL2
      INTEGER               CGXYZ  ( 3 )
      INTEGER               CORSYS
      INTEGER               CVID
      INTEGER               DIM
      INTEGER               DIM1
      INTEGER               DIM2
      INTEGER               FINAL
      INTEGER               FX
      INTEGER               FY
      INTEGER               FZ
      INTEGER               GROUP
      INTEGER               GRPBEG
      INTEGER               GRPEND
      INTEGER               GRPSIZ
      INTEGER               I
      INTEGER               I1
      INTEGER               I2
      INTEGER               I3
      INTEGER               J
      INTEGER               K
      INTEGER               MINIDX
      INTEGER               NCGR
      INTEGER               NGROUP
      INTEGER               NP
      INTEGER               NPLATE
      INTEGER               NV
      INTEGER               NVBUF
      INTEGER               NVXOUT
      INTEGER               NVXTOT
      INTEGER               OFFSET
      INTEGER               ORDVEC ( MAXPL )
      INTEGER               PLATID ( MAXPL )
      INTEGER               PLROOM
      INTEGER               PNTR
      INTEGER               PRVDSC ( DLADSZ )
      INTEGER               PRVHAN
      INTEGER               SOURCE ( MAXPL )
      INTEGER               START
      INTEGER               TO
      INTEGER               TOTPLT
      INTEGER               VGREXT ( 3 )
      INTEGER               VI
      INTEGER               VIDS   ( 3 )
      INTEGER               VIDXS  (    BUFSIZ )
      INTEGER               VLOC
      INTEGER               VOXLST ( 3, NVXLST )
      INTEGER               VOXNPL
      INTEGER               VOXNPT
      INTEGER               VOXPTR
      INTEGER               VTXNPL
      INTEGER               VXC1
      INTEGER               VXC2
      INTEGER               VXC3
      INTEGER               VXLCG  ( 3, NVXLST )
      INTEGER               VXLOUT ( NVXLST )
      INTEGER               VXLSTR ( NVXLST )
      INTEGER               W

      LOGICAL               BOXHIT
      LOGICAL               EXTRA
      LOGICAL               HAVE
      LOGICAL               HITS
      LOGICAL               INSEG
      LOGICAL               INVOX
      LOGICAL               NEWSEG

C
C     Saved variables
C
      SAVE                  CGREXT
      SAVE                  CGRPTR
      SAVE                  CGSCAL
      SAVE                  CGSCL2
      SAVE                  CORSYS
      SAVE                  DSKDSC
      SAVE                  GRDEXT
      SAVE                  GRDTOL
      SAVE                  NCGR
      SAVE                  NP
      SAVE                  NV
      SAVE                  NVXTOT
      SAVE                  ORDVEC
      SAVE                  PLATID
      SAVE                  PRVDSC
      SAVE                  PRVHAN
      SAVE                  SOURCE
      SAVE                  VBUFF
      SAVE                  VIDXS
      SAVE                  VGREXT
      SAVE                  VOXLST
      SAVE                  VOXNPL
      SAVE                  VOXNPT
      SAVE                  VOXORI
      SAVE                  VOXSIZ
      SAVE                  VTXBDS
      SAVE                  VTXNPL
      SAVE                  VXLCG
      SAVE                  VXLOUT
      SAVE                  VXLSTR
      SAVE                  XTOL

C
C     Initial values
C     
      DATA                  GRDEXT  / 3 * -1.D0 /
      DATA                  PRVHAN  / 0 /
      DATA                  PRVDSC  / DLADSZ * 0 /

C
C     Statement functions
C
      VOX2ID ( I1, I2, I3, DIM1, DIM2 ) =  I1 + DIM1*( I2 +  I3* DIM2 
     .                                                    - (1 + DIM2) )

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF
 
      CALL CHKIN ( 'DSKX02' )

C
C     Until we have better knowledge we assume there is no intersection.
C
      PLID   =  0
      FOUND  = .FALSE.
      HAVE   = .FALSE.
      NEAR   = DPMAX()

C
C     Initialize the vertex buffer.
C
      CALL CLEARI( BUFSIZ, VIDXS )
 
C
C     Check whether the ray direction vector is the zero vector.
C
      IF ( VZERO( RAYDIR ) ) THEN

         CALL SETMSG ( 'Ray direction is the zero vector.' )
         CALL SIGERR ( 'SPICE(RAYISZEROVECTOR)'            )
         CALL CHKOUT ( 'DSKX02'                            )
         RETURN

      END IF
 
C
C     Obtain the unit vector of the ray from the observer.
C
      CALL VHAT ( RAYDIR, UDIR )

C
C     Decide whether we're looking at a new segment.
C
      NEWSEG = .TRUE.

      IF ( HANDLE .EQ. PRVHAN ) THEN
C
C        The input handle matches the previous handle.  Note that the
C        initial value of PRVHAN is 0, which is never a valid handle,
C        so on the first pass, this test will fail.
C
         IF (       ( DLADSC(IBSIDX) .EQ. PRVDSC(IBSIDX) )
     .        .AND. ( DLADSC(DBSIDX) .EQ. PRVDSC(DBSIDX) )
     .        .AND. ( DLADSC(CBSIDX) .EQ. PRVDSC(CBSIDX) )  ) THEN
C
C           All of the DLA segment base addresses match.
C
            NEWSEG = .FALSE.

         END IF

      END IF


      IF ( NEWSEG ) THEN
C
C        Make sure we can't have a match with an earlier
C        segment on a subsequent call, if we exit this
C        block due to an error. 
C
         PRVHAN = 0

C
C        Retrieve the voxel grid origin in model
C        units and calculate the farthest extent of the
C        voxel grid in voxel space.
C
         CALL DSKB02 ( HANDLE, DLADSC, NV,     NP,     NVXTOT,  
     .                 VTXBDS, VOXSIZ, VOXORI, VGREXT,
     .                 CGSCAL, VTXNPL, VOXNPT, VOXNPL         )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKX02' )
            RETURN
         END IF
        
C
C        Compute the grid dimensions in units of km. First check
C        the voxel size.
C
         IF ( VOXSIZ .EQ. 0 ) THEN
            CALL SETMSG ( 'Voxel size is zero. This is an error '
     .      //            'in the DSK file attached to handle #.' )
            CALL ERRINT ( '#', HANDLE                             )
            CALL SIGERR ( 'SPICE(INVALIDVALUE)'                   )
            CALL CHKOUT ( 'DSKX02'                                )
            RETURN

         END IF

         DO I = 1, 3
            GRDEXT(I) = VGREXT(I) * VOXSIZ
         END DO

C
C        Set the margin used for checking whether the ray's vertex
C        is close to the voxel grid.
C
         GRDTOL = VTXTOL * MAX( GRDEXT(1), GRDEXT(2), GRDEXT(3) )

C
C        Check the coarse grid voxel scale.
C        
         IF ( CGSCAL .LT. 1 ) THEN

            CALL SETMSG ( 'Coarse grid scale = #; should be >= 1.' )
            CALL ERRINT ( '#', CGSCAL                              )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
            CALL CHKOUT ( 'DSKX02'                                 )
            RETURN

         END IF

C
C        Get the coarse voxel grid dimensions and the coarse voxel
C        occupancy flags.
C
         DO I = 1, 3
            CGREXT(I) = VGREXT(I) / CGSCAL
         END DO

         CGSCL2 = CGSCAL * CGSCAL
      
         NCGR   = NVXTOT / CGSCAL**3

         IF ( NCGR .GT. MAXCGR ) THEN

            CALL SETMSG ( 'Coarse grid size NCGR = #. Buffer size = #' )
            CALL ERRINT ( '#', NCGR                                    )
            CALL ERRINT ( '#', MAXCGR                                  )
            CALL SIGERR ( 'SPICE(GRIDTOOLARGE)'                        )
            CALL CHKOUT ( 'DSKX02'                                     )
            RETURN

         END IF

         CALL DSKI02 ( HANDLE, DLADSC, KWCGPT, 1, MAXCGR, DIM, CGRPTR )

         CALL DSKGD  ( HANDLE, DLADSC, DSKDSC )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKX02' )
            RETURN
         END IF

         CORSYS = NINT( DSKDSC(SYSIDX) )

         
         PRVHAN = HANDLE

         CALL MOVEI ( DLADSC, DLADSZ, PRVDSC )

      END IF

C
C     Compute tolerance used for determining whether an intercept 
C     is inside a voxel. The expansion fraction must be fetched
C     on every call to DSKX02.
C        
      CALL DSKGTL ( KEYXFR, XPDFRC )

      XTOL   = XPDFRC * MAX( ABS( GRDEXT(1) ), 
     .                       ABS( GRDEXT(2) ), 
     .                       ABS( GRDEXT(3) )  )

C
C     Find the ray intercept on the surface of the fine voxel grid,
C     if the intercept exists.
C
      CALL ZZRAYBOX ( VERTEX, RAYDIR, VOXORI, GRDEXT, VTX2, BOXHIT )

      IF ( .NOT. BOXHIT ) THEN         
         CALL CHKOUT ( 'DSKX02' )
         RETURN
      END IF

C
C     Convert the grid intercept to voxel space coordinates.
C     The result COORD will be used as the ray's vertex in XDDA.
C     
      CALL ZZTOGRID ( VTX2, VOXORI, VOXSIZ, COORD )

C
C     Determine the voxels hit by the ray.
C
      CALL XDDA ( COORD, UDIR, VGREXT, NVXLST, NVXOUT, VOXLST )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'DSKX02' )
         RETURN
      END IF
                 
C
C     We don't expect the voxel list to be empty, but leave now
C     if it is.
C     
      IF ( NVXOUT .EQ. 0 ) THEN
         CALL CHKOUT ( 'DSKX02' )
         RETURN
      END IF

C
C     Rather than using the original observer's location, we use a
C     location derived from COORD, which is the intercept of the ray
C     and the surface of the voxel grid.  We start with COORD, convert
C     this location to the model coordinate system, and back outward a
C     bit to make sure we obtain a location outside the grid (we don't
C     want to miss any plates that might be located right on the grid's
C     surface).
C
C     This vertex change is not performed if the vertex is already
C     inside, or within a small margin away from, the voxel grid.
C
      CALL VSUB ( VERTEX, VOXORI, VTXOFF )

      IF (      (  VTXOFF(1) .LT.   -GRDTOL               )
     .     .OR. (  VTXOFF(1) .GT.  ( GRDTOL + GRDEXT(1) ) )
     .     .OR. (  VTXOFF(2) .LT.   -GRDTOL               )  
     .     .OR. (  VTXOFF(2) .GT.  ( GRDTOL + GRDEXT(2) ) )
     .     .OR. (  VTXOFF(3) .LT.   -GRDTOL               )  
     .     .OR. (  VTXOFF(3) .GT.  ( GRDTOL + GRDEXT(3) ) )   ) THEN

C
C        The vertex is outside of the voxel grid by more than the
C        margin. Move the ray-grid intercept away from the grid to
C        improve numeric performance.
C
         CALL VLCOM3 (  1.D0, VOXORI, VOXSIZ, COORD, 
     .                 -1.D0, UDIR,   VTX2           )

      END IF

C
C     We are going to need to subtract the location of the observer
C     from vertices of a plate. To speed things up a tiny bit, we'll
C     make 3 copies of the observer's location so that we make a single
C     subroutine call to handle the 3 subtractions.
C                                        
      CALL VEQU ( VTX2, OBSMAT(1,1) )
      CALL VEQU ( VTX2, OBSMAT(1,2) )
      CALL VEQU ( VTX2, OBSMAT(1,3) ) 

C
C     Use the coarse voxel grid to compress the voxel list. We
C     remove all voxels belonging to empty coarse voxels.
C
      TO = 0
  
      DO I = 1, NVXOUT
C
C        Find the coordinates, then the ID, of the coarse voxel
C        containing this voxel.
C         
         DO J = 1, 3
            CGXYZ(J) = ( ( VOXLST(J,I) - 1 ) / CGSCAL ) + 1
         END DO


         CVID = VOX2ID ( CGXYZ(1),  CGXYZ(2), CGXYZ(3), 
     .                   CGREXT(1), CGREXT(2)           )

         IF ( CGRPTR(CVID) .GT. 0 ) THEN
C
C           This coarse voxel is non-empty; add the index of the
C           current voxel to the output list.  Save the coordinates of
C           the parent coarse voxel as well.
C           
            TO         = TO + 1
            VXLOUT(TO) = I

            DO J = 1, 3
               VXLCG (J,TO) = CGXYZ (J)
            END DO
C
C           Save the coarse voxel start pointer as well.
C
            VXLSTR(TO) = CGRPTR(CVID)

         END IF

      END DO

C
C     Update NVXOUT to be the number of voxels in the compressed list.
C
      NVXOUT = TO
C
C     If the voxel list became empty after compression, we're
C     done.
C
      IF ( NVXOUT .EQ. 0 ) THEN
         CALL CHKOUT( 'DSKX02' )
         RETURN
      END IF

C
C     Initialize PNTR in case the voxel list is empty.
C     (This is a backstop precaution: the voxel list
C     should never be empty at this point.) PNTR will
C     be referenced after the end of the outer loop below.
C
      PNTR   = 1

C
C     The vertex buffer is empty so far.
C
      NVBUF  = 0

C
C     Break up the list of voxels into groups; process each 
C     group in turn until we find an intersection or run out
C     of voxels.
C
      GRPSIZ =  MAX ( 1,  (NVXOUT+1)/2)

      NGROUP =  (NVXOUT-1)/GRPSIZ + 1
          
      GROUP  =  1


      DO WHILE (  ( GROUP .LE. NGROUP )  .AND.  ( .NOT. HAVE )  )


         PNTR   = 1

         GRPBEG = (GROUP-1)*GRPSIZ + 1

         GRPEND = MIN ( GRPBEG + GRPSIZ - 1,  NVXOUT )

         PLROOM = MAXPL

         DO VI = GRPBEG, GRPEND
C
C           Look up the plate list pointer for this voxel.
C
C
C           We begin by finding the offset of the voxel from
C           the base of its parent coarse voxel.
C
            J      =  VXLOUT(VI)

            FX     =  VOXLST(1,J)  -  CGSCAL * ( VXLCG(1,VI) - 1 )
            FY     =  VOXLST(2,J)  -  CGSCAL * ( VXLCG(2,VI) - 1 )
            FZ     =  VOXLST(3,J)  -  CGSCAL * ( VXLCG(3,VI) - 1 )

            OFFSET = FX + CGSCAL*(FY-1) + CGSCL2*(FZ-1)

C
C           Now compute the index of voxel-plate list pointer in
C           the pointer array, and look up the pointer.
C
            J = VXLSTR(VI) + OFFSET - 1

            CALL DSKI02 ( HANDLE, DLADSC, KWVXPT, J, 1, DIM, VOXPTR )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'DSKX02' )
               RETURN
            END IF

            IF ( VOXPTR .EQ. -1 ) THEN
               
               NPLATE = 0

            ELSE
C
C              Get the plate count for this voxel.
C
               CALL DSKI02 ( HANDLE, DLADSC, KWVXPL, VOXPTR, 
     .                       1,      DIM,    NPLATE          )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'DSKX02' )
                  RETURN
               END IF

            END IF


            IF ( NPLATE .GT. 0 ) THEN

               IF ( NPLATE .LE. PLROOM ) THEN
C
C                 Get the plate list for this voxel.
C
                  CALL DSKI02 ( HANDLE,   DLADSC, KWVXPL,  
     .                          VOXPTR+1, NPLATE, DIM,    PLATID(PNTR) )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'DSKX02' )
                     RETURN
                  END IF

                  PLROOM = PLROOM - NPLATE

               ELSE

                  CALL SETMSG ( 'NPLATE = #. Available room in '
     .            //            'PLATID array = #. Array size = #.' )
                  CALL ERRINT ( '#', NPLATE                         )
                  CALL ERRINT ( '#', PLROOM                         )
                  CALL ERRINT ( '#', MAXPL                          )
                  CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'              )
                  CALL CHKOUT ( 'DSKX02'                            )
                  RETURN

               END IF

C
C              Fill in the corresponding elements of the parallel
C              "source" array with the current voxel loop index.
C              XDDA lists these voxels in the order the ray hits
C              them, so the lowest indexed voxels are hit first.
C
               CALL FILLI ( VI, NPLATE, SOURCE(PNTR) )

            END IF

C
C           NPLATE returns zero or greater.
C              
            PNTR = PNTR + NPLATE

         END DO
C
C        We've collected all plate info for the current voxel group.
C
         TOTPLT = PNTR-1 

C
C        We want to sort the plate ID array and remove duplicates.
C        However, we want to keep the plates ordered according to the
C        sequence in which their containing voxels were hit by the ray.
C        So we find the order vector for the plate ID array, then use
C        this vector to mark duplicates.
C
         CALL ORDERI ( PLATID, TOTPLT, ORDVEC )

C
C        Negate the plate ID of every duplicate we find, leaving
C        the instance in the voxel closest to the ray's origin
C        unmarked. For every pair of plates with the same ID, 
C        we'll mark the one with the greater index in the plate
C        ID array.
C         
C        We use MINDIX to identify the index, in the plate ID array,
C        of the current unmarked plate ID. MINIDX is re-used for
C        each set of distinct plate IDs.
C
C        The following loop considers plate IDs in increasing order.
C
         MINIDX = ORDVEC(1)

         DO I = 2, TOTPLT
C
C           The condition below uses absolute value because the plate
C           ID at index I-1 may have been "marked" via negation.
C 
            IF (              PLATID( ORDVEC(I  ) )  
     .           .EQ.  ABS (  PLATID( ORDVEC(I-1) )  )    ) THEN
C
C              The plates having indices ORDVEC(I-1) and ORDVEC(I) are
C              duplicates. 
C
C              At this point MINIDX is the lowest index in the plate ID
C              array of any plate seen so far having an ID equal to
C              PLATID(ORDVEC(I)).
C
C              If ORDVEC(I) is the new "minimum," negate the plate ID
C              at the old "minimum"; otherwise negate the plate ID at
C              index ORDVEC(I).
C
               IF ( ORDVEC(I) .LT. MINIDX ) THEN
C
C                 The plate that was previously at the minimum index is
C                 now considered a duplicate. The new minimum index for
C                 the current plate ID value is ORDVEC(I).
C              
                  PLATID(MINIDX)  =  - PLATID(MINIDX) 
                  MINIDX          =    ORDVEC(I)

               ELSE
C
C                 The current plate is a duplicate; mark it.
C
                  PLATID( ORDVEC(I) ) =  - PLATID( ORDVEC(I) )

               END IF

            ELSE
C
C              We're looking at a new plate ID. For the moment, this
C              ID has no duplicates.
C 
               MINIDX = ORDVEC(I)

            END IF
 
         END DO
 
C
C        If something went wrong up above there is no point in
C        going on from here.
C 
         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKX02' )
            RETURN
         END IF

C
C        Now examine each plate in the PLATID list for this voxel group.
C        PNTR has the value of the index available for data in
C        PLATID, so the final location of data is at index PNTR - 1.
C
         EXTRA  = .FALSE.
         FINAL  = 0
         NEAR   = 0.D0

         I      = 1

         DO WHILE (  ( I .LE. TOTPLT ) .AND. ( .NOT. EXTRA )  )
C
C           Retrieve the current plate ID.
C
            J = PLATID(I)


            IF ( J .GT. 0 ) THEN
C
C              This is not a duplicate plate; consider it.
C           
               IF ( HAVE ) THEN
C
C                 We already have a hit. See whether this plate resides
C                 in the voxel in which the last hit was found, or in a
C                 later voxel.
C
                  IF ( SOURCE(I) .GT. FINAL ) THEN
C
C                    This is a "late plate": it occurs in a voxel later
C                    than that in which the first valid hit was found.
C
                     EXTRA = .TRUE.

                  END IF

               END IF
 

               IF ( .NOT. EXTRA ) THEN
C
C                 Fetch the vertex IDs of this plate.
C
                  START = 3*(J-1) + 1
      
                  CALL DSKI02 ( HANDLE, DLADSC, KWPLAT, 
     .                          START,  3,      DIM,    VIDS )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'DSKX02' )
                     RETURN
                  END IF

C
C                 Fetch the vertices themselves.
C
                  DO K = 1, 3
C
C                    Any vertex may be buffered already. Look in
C                    the vertex buffer before reading the vertex.
C
                     VLOC = ISRCHI ( VIDS(K), NVBUF, VIDXS )

                     IF ( VLOC .GT. 0 ) THEN
C
C                       The vertex was buffered; just copy it.
C
                        CALL VEQU ( VBUFF(1,VLOC), POINTS(1,K) )

                     ELSE
C
C                       Read in the vertex.
C
                        START = 3*( VIDS(K) - 1 ) + 1

                        CALL DSKD02 ( HANDLE, DLADSC, KWVERT,  
     .                                START,  3,      DIM,    
     .                                POINTS(1,K)            )

                        IF ( FAILED() ) THEN
                           CALL CHKOUT ( 'DSKX02' )
                           RETURN
                        END IF

C
C                       If there's room, buffer this vertex.
C
                        IF ( NVBUF .LT. BUFSIZ ) THEN

                           NVBUF = NVBUF + 1

                           CALL VEQU ( POINTS(1,K), VBUFF(1,NVBUF) )

                           VIDXS( NVBUF ) = VIDS(K)

                        END IF

                     END IF

                  END DO               

               END IF


               IF ( .NOT. EXTRA ) THEN
C                 
C                 The current plate qualifies for testing using INSANG.
C
C                 Retrieve the model coordinates of the J'th plate's
C                 three vertices. Expand the plate slightly to prevent
C                 round-off error from causing the ray to miss the
C                 plate. Compute the edges of the tetrahedral angle
C                 with the observer as the apex and the vertices as
C                 members of the edge rays. Finally see if the
C                 boresight lies inside the tetrahedron.
C
                  CALL VSUBG  ( POINTS,     OBSMAT, 9,  EDGES )

                  CALL PLTEXP ( EDGES,      XPDFRC,     XPNTS )

                  CALL INSANG ( UDIR,
     .                          XPNTS(1,1), XPNTS(1,2), XPNTS(1,3),
     .                          HITS,       SCALE )    

                  IF ( HITS ) THEN
C
C                    Reject intersections with plates that face away
C                    from the ray. Accept intersections with plates
C                    that face toward the ray.
C
                     CALL PLTNRM ( POINTS(1,1), POINTS(1,2), 
     .                             POINTS(1,3), NORMAL      )
                                         
                     HITS = VDOT( UDIR, NORMAL ) .LE. 0.D0 

                  END IF

                  IF ( HITS ) THEN
C
C                    The boresight intersects this plate. 
C
                     IF ( ( .NOT. HAVE ) .OR. ( SCALE .LT. NEAR ) ) THEN
C
C                       Either this is the first intersection we've
C                       found, or this is the closest intersection to
C                       the vertex we've found. Compute the intercept
C                       coordinates and see whether the intercept is
C                       within the current voxel. Use a small tolerance
C                       for the comparison.

C                       If this intersection point is closer to the
C                       ray's vertex than the last one, pick this point
C                       and the plate it's on.
C                          ___   ____   __________
C                          XPT = VTX2 + SCALE*UDIR
C
                        CALL VLCOM ( 1.D0, VTX2, SCALE, UDIR, XPT )

C
C                       Compute the voxel grid coordinates of the
C                       intercept. HITCOR is a double precision vector
C                       having units of voxels (voxel edge length, to
C                       be precise). Note that the components of HITCOR
C                       are zero-based.
C
                        CALL ZZTOGRID ( XPT, VOXORI, VOXSIZ, HITCOR )
C
C                       Look up the voxel grid coordinates (integer,
C                       1-based) of the current voxel.
C
                        K     = VXLOUT ( SOURCE(I) )

                        VXC1  = VOXLST(1,K) 
                        VXC2  = VOXLST(2,K) 
                        VXC3  = VOXLST(3,K) 

                        INVOX =       ( HITCOR(1) .GT. (VXC1-XTOL-1) )
     .                          .AND. ( HITCOR(1) .LT. (VXC1+XTOL  ) )
     .                          .AND. ( HITCOR(2) .GT. (VXC2-XTOL-1) )
     .                          .AND. ( HITCOR(2) .LT. (VXC2+XTOL  ) )
     .                          .AND. ( HITCOR(3) .GT. (VXC3-XTOL-1) )
     .                          .AND. ( HITCOR(3) .LT. (VXC3+XTOL  ) )


                        IF ( INVOX ) THEN
C
C                          Reject solutions that are outside of the
C                          segment's boundaries, where the boundaries
C                          are extended using the "greedy" margin.
C
                           CALL DSKGTL ( KEYSGR, GREEDM )

                           CALL ZZINVELT ( XPT,            
     .                                     CORSYS, 
     .                                     DSKDSC(PARIDX),  
     .                                     DSKDSC(MN1IDX), 
     .                                     GREEDM,
     .                                     NONE,
     .                                     INSEG          )
     
                           IF ( INSEG ) THEN
C
C                             We have a viable intercept. Record the
C                             scale, plate ID, and source voxel index
C                             in the compressed voxel list pointer
C                             array VXLOUT. We won't look for
C                             intercepts beyond the voxel designated by
C                             FINAL.
C
                              HAVE   = .TRUE.
                              NEAR   = SCALE
                              PLID   = J
                              FINAL  = SOURCE(I)
C
C                             Indicate that a solution was found. We'll
C                             keep looking for a better one if PLID is
C                             not the last plate of the current voxel.
C
                              FOUND  = .TRUE.

                           END IF

                        ELSE
C
C                          We must re-consider this plate if we
C                          encounter it in a voxel later in the voxel
C                          list. Remove all duplication markings for
C                          this plate.
C
                           W = ABS(J)

                           DO K = 1, TOTPLT

                              IF ( ABS(PLATID(K)) .EQ. W ) THEN
                                 PLATID(K) = W
                              END IF

                           END DO

                        END IF

                     END IF
C
C                    End of case of checking possible solution
C                    intercept.
C
                  END IF
C
C                 We're done with processing the current hit.
C 
               END IF
C
C              We're done with processing the current qualifying plate.
C
            END IF
C
C           We're done with the current plate.              
C
C           Fetch the next plate for this voxel group.
C
            I = I + 1

         END DO
C
C        We're done with the current voxel group.
C
         GROUP = GROUP + 1

      END DO
C
C     We've either found an intercept or have run out of voxel groups
C     to check.
C
C     That's all folks.
C
      CALL CHKOUT ( 'DSKX02' )
      RETURN
      END
 
