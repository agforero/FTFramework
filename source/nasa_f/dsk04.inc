C
C     Include file dsk04.inc
C
C     This include file declares parameters for DSK data type 4
C     (double precision raster values).
C
C        Version 4.1.0 07-OCT-2016 (NJB)
C     
C           Corrected some long lines for delivery to DSKBRIEF.
C
C        Version 4.0.0 12-NOV-2012 (NJB)
C
C           Added parameter MARGIN, which is used to handle minor
C           round-off discrepancies.
C
C           Changed name of parameter BILCPT to NVBCPT.
C
C        Version 3.0.0 18-SEP-2012 (NJB)
C
C     Each type 4 DSK segment has only a double precision DLA component.
C     The segment layout in DAS address space is as follows:
C
C        D.p. layout:
C            
C
C   +---------------------+
C   | DSK descriptor      |  DSKDSZ elements
C   +---------------------+
C   | Format version ID   |  Format version ID code
C   +---------------------+
C   | NR                  |  number of rows
C   +---------------------+
C   | NC                  |  number of pixel columns
C   +---------------------+
C   | Coord 1 pixel size  |  size in units of coordinate 1
C   +---------------------+
C   | Coord 2 pixel size  |  size in units of coordinate 2
C   +---------------------+
C   | Center coord 1      |  Coordinates of center of
C   +---------------------+  pixel grid
C   | Center coord 2      |
C   +---------------------+
C   | Nulls ok flag       |  0/1 for nulls not allowed/allowed
C   +---------------------+
C   | Null value          |  Integer value indicating "no data"
C   +---------------------+
C   | Height scale        |  Height units in km.
C   +---------------------+
C   | Numeric format      |  Numeric format, e.g. packed 16-bit integer
C   +---------------------+  
C   |Interpolation dsc ptr|  Pointer to interpolation descriptor
C   +---------------------+
C   | Intercept descr ptr |  Pointer to intercept descriptor
C   +---------------------+
C   | Normal descr ptr    |  Pointer to normal vector descriptor
C   +---------------------+
C   | Accel descr ptr     |  Pointer to acceleration descriptor
C   +---------------------+
C   | Ref surface dsc ptr |  Reference surface descriptor pointer
C   +---------------------+
C   | Projection dsc ptr  |  Map projection descriptor pointer
C   +---------------------+
C   | Coarse grid ptr     |  Pointer to coarse grid structure
C   +---------------------+
C   | Data start ptr      |  Pointer to start of pixel data
C   +---------------------+
C   | Interp descr        |  Interpolation method descriptor
C   +---------------------+
C   | Intercept descr     |  Intercept method descriptor
C   +---------------------+
C   | Normal vec descr    |  Normal vector calc descriptor
C   +---------------------+
C   | Accel descr         |  Intercept acceleration descriptor
C   +---------------------+
C   | Ref surface descr   |  Reference surface descriptor
C   +---------------------+
C   | Map projection descr|  Map projection descriptor
C   +---------------------+
C   | Coarse grid struct  |  Coarse grid structure
C   +---------------------+
C
C
C        Integer layout:
C
C   +---------------------+
C   | Grid dimension count|  Depth of nested grids
C   +---------------------+
C   | Grid dimensions     |  Row and column sizes of nested grids
C   +---------------------+
C   | Pixel data          |  Pixels containing height values
C   +---------------------+   
C
C
C
C     Indices of descriptor items at fixed locations
C     ==============================================
C
C     The following indices are relative to the base address
C     of the double precision DAS address range of the DLA
C     segment containing the described type 4 DSK segment.
C     Indices are 1-based: index 1 is the first DAS double
C     precision number of the segment.
C

C     
C     Index of DSK descriptor:
C
      INTEGER               IXDSCR
      PARAMETER           ( IXDSCR = 1 )

C
C     DSK descriptor size:
C
C     This local parameter MUST be kept consistent with
C     the parameter DSKDSZ which is declared in dskdsc.inc.
C
      INTEGER               DSCSIZ
      PARAMETER           ( DSCSIZ = 24 )

C
C     Index of format version:
C
      INTEGER               IXFMTV
      PARAMETER           ( IXFMTV = IXDSCR + DSCSIZ )
      
C
C     Index of row count:
C     
      INTEGER               IXNR
      PARAMETER           ( IXNR   = IXFMTV + 1 )

C
C     Index of column count:
C     
      INTEGER               IXNC
      PARAMETER           ( IXNC   = IXNR   + 1 )

C
C     Index of pixel size in direction of coordinate 1:
C     
      INTEGER               IXPSZ1
      PARAMETER           ( IXPSZ1 = IXNC   + 1 )

C
C     Index of pixel size in direction of coordinate 2:
C     
      INTEGER               IXPSZ2
      PARAMETER           ( IXPSZ2 = IXPSZ1 + 1 )

C
C     Index of coordinate 1 of pixel grid center:
C
      INTEGER               IXCTR1
      PARAMETER           ( IXCTR1 = IXPSZ2 + 1 )

C
C     Index of coordinate 2 of pixel grid center:
C
      INTEGER               IXCTR2
      PARAMETER           ( IXCTR2 = IXCTR1 + 1 )

C
C     Index of flag indicating whether null values are allowed:
C     
      INTEGER               IXNLOK
      PARAMETER           ( IXNLOK = IXCTR2 + 1 )

C
C     Index of null value indicator:
C     
      INTEGER               IXNULL
      PARAMETER           ( IXNULL = IXNLOK + 1 )

C
C     Index of height scale:
C
      INTEGER               IXHSCL
      PARAMETER           ( IXHSCL = IXNULL + 1 )

C
C     Index of numeric format:
C
      INTEGER               IXNUMF
      PARAMETER           ( IXNUMF = IXHSCL + 1 )

C
C     Index of pointer to interpolation algorithm descriptor:
C     
      INTEGER               IXITPP
      PARAMETER           ( IXITPP = IXNUMF + 1 )

C
C     Index of pointer to intercept algorithm descriptor:
C     
      INTEGER               IXXP
      PARAMETER           ( IXXP   = IXITPP + 1 )

C
C     Index of pointer to normal vector computation descriptor:
C     
      INTEGER               IXNVP
      PARAMETER           ( IXNVP  = IXXP   + 1 )

C
C     Index of pointer to intercept acceleration descriptor:
C     
      INTEGER               IXACCP
      PARAMETER           ( IXACCP = IXNVP  + 1 )

C
C     Index of reference surface decriptor pointer:
C
      INTEGER               IXREFP
      PARAMETER           ( IXREFP = IXACCP + 1 )

C
C     Index of projection decriptor pointer:
C
      INTEGER               IXPRJP
      PARAMETER           ( IXPRJP = IXREFP + 1 )

C
C     Index of pointer to coarse grid structure:
C     
      INTEGER               IXCGP
      PARAMETER           ( IXCGP  = IXPRJP + 1 )

C
C     Index of pointer to pixel data:
C     
      INTEGER               IXPIXP
      PARAMETER           ( IXPIXP = IXCGP  + 1 )

 
C
C     Indices of integer items
C     ========================

C
C     Index of grid dimension count:
C
      INTEGER               IXNDIM
      PARAMETER           ( IXNDIM = 1 )


C
C     Index of grid dimensions:
C
      INTEGER               IXGDIM
      PARAMETER           ( IXGDIM = IXNDIM + 1 )

C
C     Maximum number of grid dimensions:
C
      INTEGER               DMAX
      PARAMETER           ( DMAX   = 10 )

C
C     The grid dimension array has fixed size of 2*DMAX integers.
C

C
C     Double precision item keyword parameters used by fetch routines
C     ===============================================================
C
C     DSK descriptor:
C
      INTEGER               KWDSC
      PARAMETER           ( KWDSC  = 1 )

C
C     Format version:
C
      INTEGER               KWFMTV
      PARAMETER           ( KWFMTV = KWDSC  + 1 )

C
C     Keyword for row count:
C     
      INTEGER               KWNR
      PARAMETER           ( KWNR   = KWFMTV + 1 )

C
C     Keyword for column count:
C     
      INTEGER               KWNC
      PARAMETER           ( KWNC   = KWNR   + 1 )

C
C     Keyword for pixel size in direction of coordinate 1:
C     
      INTEGER               KWPSZ1
      PARAMETER           ( KWPSZ1 = KWNC + 1 )

C
C     Keyword for pixel size in direction of coordinate 2:
C     
      INTEGER               KWPSZ2
      PARAMETER           ( KWPSZ2 = KWPSZ1 + 1 )

C
C     Keyword for coordinate 1 of pixel grid center:
C
      INTEGER               KWCTR1
      PARAMETER           ( KWCTR1 = KWPSZ2 + 1 )

C
C     Keyword for coordinate 2 of pixel grid center:
C
      INTEGER               KWCTR2
      PARAMETER           ( KWCTR2 = KWCTR1 + 1 )

C
C     Keyword for flag indicating whether null values are allowed:
C     
      INTEGER               KWNLOK
      PARAMETER           ( KWNLOK = KWCTR2 + 1 )

C
C     Keyword for null value indicator:
C     
      INTEGER               KWNULL
      PARAMETER           ( KWNULL = KWNLOK + 1 )

C
C     Keyword for height scale:
C
      INTEGER               KWHSCL
      PARAMETER           ( KWHSCL = KWNULL + 1 )

C
C     Keyword for numeric format:
C
      INTEGER               KWNUMF
      PARAMETER           ( KWNUMF = KWHSCL + 1 )

C
C     Keyword for pointer to reference surface descriptor:
C     
      INTEGER               KWREFP
      PARAMETER           ( KWREFP = KWNUMF + 1 )

C
C     Keyword for pointer to map projection descriptor:
C     
      INTEGER               KWPRJP
      PARAMETER           ( KWPRJP = KWREFP + 1 )

C
C     Keyword for pointer to interpolation algorithm descriptor:
C     
      INTEGER               KWITPP
      PARAMETER           ( KWITPP = KWPRJP + 1 )

C
C     Keyword for pointer to intercept algorithm descriptor:
C     
      INTEGER               KWXP
      PARAMETER           ( KWXP   = KWITPP + 1 )

C
C     Keyword for pointer to normal vector computation descriptor:
C     
      INTEGER               KWNVP
      PARAMETER           ( KWNVP  = KWXP   + 1 )

C
C     Keyword for pointer to intercept acceleration descriptor:
C     
      INTEGER               KWACCP
      PARAMETER           ( KWACCP = KWNVP  + 1 )

C
C     Keyword for pointer to coarse grid structure:
C     
      INTEGER               KWCGP
      PARAMETER           ( KWCGP  = KWACCP  + 1 )

C
C     Keyword for pointer to pixel data:
C     
      INTEGER               KWPIXP
      PARAMETER           ( KWPIXP = KWCGP  + 1 )

C
C     Keyword for interpolation algorithm descriptor:
C     
      INTEGER               KWITPD
      PARAMETER           ( KWITPD = KWPIXP + 1 )

C
C     Keyword for intercept algorithm descriptor:
C     
      INTEGER               KWXD
      PARAMETER           ( KWXD   = KWITPD + 1 )

C
C     Keyword for normal vector algorithm descriptor:
C     
      INTEGER               KWNVD
      PARAMETER           ( KWNVD  = KWXD   + 1 )

C
C     Keyword for coarse grid descriptor:
C     
      INTEGER               KWCGD
      PARAMETER           ( KWCGD  = KWNVD  + 1 )

C
C     Keyword for reference surface descriptor:
C
      INTEGER               KWREFD
      PARAMETER           ( KWREFD = KWCGD + 1 )

C
C     Keyword for map projection descriptor:
C
      INTEGER               KWPRJD
      PARAMETER           ( KWPRJD = KWREFD + 1 )

C
C     Keyword for intercept acceleration descriptor:
C     
      INTEGER               KWACCD
      PARAMETER           ( KWACCD = KWPRJD  + 1 )

C
C     Integer item keyword parameters used by fetch routines
C     ===============================================================

C
C     Note: keyword codes for integer items occupy a range
C     distinct from that for keywords for d.p. items. The
C     first keyword code for an integer item follows the
C     last d.p. code.
C
C
C     Keyword for nested grid dimension count:
C
      INTEGER               KWNDIM
      PARAMETER           ( KWNDIM = 1 )

C
C     Keyword for grid dimensions:
C
      INTEGER               KWGDIM
      PARAMETER           ( KWGDIM = KWNDIM + 1 )


C     Keyword for pixel data:
C
C     This keyword refers to unpacked data: each 16-bit value is
C     returned in an individual integer. The pixel grid has NR rows and
C     NC columns. There are two 16-bit pixels per stored integer. The
C     data are stored in row-major order.
C     
      INTEGER               KWPIX
      PARAMETER           ( KWPIX   = KWGDIM + 1 )

C
C     Keyword for raw data: 
C
C     This keyword refers to stored integers, each containing two
C     16-bit data values.  The raw grid has NR rows and
C     NC/2 columns. 
C
      INTEGER               KWRAW
      PARAMETER           ( KWRAW   = KWPIX + 1 )

C
C
C     Constants
C     =========
C
C     The constant NULL is used to indicate a null pointer or
C     undefined value. Since this value is stored in double precision
C     DAS words, its type is double precision.
C
      DOUBLE PRECISION      NULL
      PARAMETER           ( NULL   = -1.D0 )

C
C     The constants TRUE and FALSE are numeric codes having
C     the meanings associated with their names.
C
      DOUBLE PRECISION      TRUE
      PARAMETER           ( TRUE   = 1.D0 )

      DOUBLE PRECISION      FALSE
      PARAMETER           ( FALSE  = 0.D0 )

C
C     The constant MAXDSZ is the maximum size of any algorithm
C     descriptor used in the type 4 DSK segment structure.
C
      INTEGER               MAXDSZ
      PARAMETER           ( MAXDSZ = 7 )


C
C     Numeric format
C     ==============
C
C        Codes for 32-bit, packed 16-bit, and packed 8-bit integer
C        data:
C
      INTEGER               NUMI32
      PARAMETER           ( NUMI32 = 1 )

      INTEGER               NUMI16
      PARAMETER           ( NUMI16 = NUMI32 + 1 )

      INTEGER               NUMI08
      PARAMETER           ( NUMI08 = NUMI16 + 1 )

C
C     Code for equivalenced 32-bit real data:
C
      INTEGER               NUMR32
      PARAMETER           ( NUMR32 = NUMI08 + 1 )
      

C
C     Reference surface descriptor
C     ============================
C
C
C         Reference surface codes
C         -----------------------
C
C            Code REFSPH indicates that the reference surface is a
C            sphere.
C
C            For spheres, a single radius is stored in the parameter
C            section of the descriptor.
C
      INTEGER               REFSPH
      PARAMETER           ( REFSPH = 1 )

C
C            Code REFELL indicates that the reference surface is a
C            triaxial ellipsoid.
C
C            For ellipsoids, three radii are stored in the parameter
C            section of the descriptor.
C
      INTEGER               REFELL
      PARAMETER           ( REFELL = REFSPH + 1 )

C
C            Code REFINH indicates that the reference surface is a
C            triaxial ellipsoid having radii inherited from the segment
C            descriptor. Any radii supplied in the reference surface
C            descriptor are ignored.
C
C            The reference surface subtype is also ignored when the
C            reference surface code is REFINH.
C
      INTEGER               REFINH
      PARAMETER           ( REFINH = REFELL + 1 )

C
C         Reference surface subtype codes
C         -------------------------------
C
C            Code REFRAD indicates that surface heights are
C            measured relative to the reference surface
C            in the radial direction.
C
C            REFRAD is the usual choice for segments 
C            using planetocentric coordinates.
C
      INTEGER               REFRAD
      PARAMETER           ( REFRAD = 1 )

C
C            Code REFNRM indicates that surface heights are
C            measured relative to the reference surface
C            in the outward normal direction.
C
C            REFNRM is the usual choice for segments
C            using planetodetic coordinates.
C
      INTEGER               REFNRM
      PARAMETER           ( REFNRM = REFRAD + 1 )


C
C         Reference surface index parameters  
C         ----------------------------------
C
C         Reference surface code within descriptor:
C
      INTEGER               REFXCD
      PARAMETER           ( REFXCD = 1 )

C
C         Index of descriptor size within descriptor:
C
      INTEGER               REFXSZ
      PARAMETER           ( REFXSZ = REFXCD + 1 )

C
C         Index of subtype within descriptor:
C
      INTEGER               REFXST
      PARAMETER           ( REFXST = REFXSZ + 1 )

C
C         Index of parameters within descriptor:
C
      INTEGER               REFXPA
      PARAMETER           ( REFXPA = REFXST + 1 )



C
C     Map projection descriptor
C     =========================
C
C
C         Map projection codes
C         --------------------
C
C            Code PRJEQR indicates that the map projection is
C            equirectangular.
C
      INTEGER               PRJEQR
      PARAMETER           ( PRJEQR = 1 )

C
C            Code PRJSTE indicates that the map projection is
C            stereographic.
C
      INTEGER               PRJSTE
      PARAMETER           ( PRJSTE = PRJEQR + 1 )

C
C            Code PRJMOL indicates that the map projection is
C            MGS MOLA-style stereographic. In this projection,
C            radial distance on the tangent plane is proportional
C            to twice the sine of half of the colatitude of the
C            projected point. In a conventional stereographic
C            projection, the sine function is replaced by the
C            tangent function.
C
      INTEGER               PRJMOL
      PARAMETER           ( PRJMOL = PRJSTE + 1 )


C
C         Map projection descriptor index parameters  
C         ------------------------------------------
C
C         Index of projection code within descriptor:
C
      INTEGER               PRJXCD
      PARAMETER           ( PRJXCD = 1 )

C
C         Index of descriptor size within descriptor:
C
      INTEGER               PRJXSZ
      PARAMETER           ( PRJXSZ = PRJXCD + 1 )

C
C         For stereographic projections (including the MOLA-style
C         projection), the Euler angle matrix product
C
C            [TWIST]  [Pi/2 - COLAT]  [LON]
C                   3               2      3
C
C         describes the orientation of the projection axes relative to
C         the base frame. Angular units are radians.
C
C         The tangent point of the projection plane corresponds to the
C         intersection of the transformed Z-axis with the unit sphere.
C         X and Y pixel scales are used to transform direction vectors
C         in the base frame to points on the tangent plane.
C
C         Indices of longitude, latitude, and twist within the
C         descriptor:
C
      INTEGER               PRJXLO
      PARAMETER           ( PRJXLO = PRJXSZ + 1)

      INTEGER               PRJXLA
      PARAMETER           ( PRJXLA = PRJXLO + 1)

      INTEGER               PRJXTW
      PARAMETER           ( PRJXTW = PRJXLA + 1 )

C
C        Indices of projection scales in the X and Y directions:
C
      INTEGER               PRJXXS
      PARAMETER           ( PRJXXS = PRJXTW + 1 )

      INTEGER               PRJXYS
      PARAMETER           ( PRJXYS = PRJXXS + 1 )


C
C     Interpolation algorithm descriptor
C     ==================================
C
C     The interpolation algorithm descriptor is a variable-size set of
C     parameters that specify the interpolation algorithm used to
C     provide values of the third coordinate for arbitrary pairs of
C     first and second coordinate values.
C
C     The descriptor starts with a code identifying the class of
C     algorithm. This is followed by the size of the descriptor. If
C     applicable, additional parameters follow the size.
C
C         Interpolation algorithm codes
C         -----------------------------
C
C            Code ITPNO indicates that no interpolation
C            is performed.
C
      DOUBLE PRECISION      ITPNO
      PARAMETER           ( ITPNO  = 1 )

C
C            Code ITPBIL denotes bilinear interpolation.
C
      DOUBLE PRECISION      ITPBIL
      PARAMETER           ( ITPBIL = 2 )

C
C            Code ITPNN3 denotes planar interpolation using
C            the three nearest neighboring data values.
C
      DOUBLE PRECISION      ITPNN3
      PARAMETER           ( ITPNN3 = 3 )

C
C         Interpolation descriptor index parameters  
C         -----------------------------------------
C
C         Index of interpolation code within descriptor:
C
      INTEGER               ITPXCD
      PARAMETER           ( ITPXCD = 1 )

C
C         Index of descriptor size within descriptor:
C
      INTEGER               ITPXSZ
      PARAMETER           ( ITPXSZ = ITPXCD + 1 )

C
C         Index of first parameter within descriptor:
C
      INTEGER               ITPXPR
      PARAMETER           ( ITPXPR = ITPXSZ + 1 )
     

C
C     Intercept algorithm descriptor
C     ==============================
C
C     The ray-surface intercept algorithm descriptor is a variable-size
C     set of parameters that specify the intercept algorithm used for
C     the associated segment.
C
C     The descriptor starts with a code identifying the class of
C     algorithm. This is followed by the size of the descriptor. If
C     applicable, additional parameters follow the size.
C
C         Intercept algorithm codes
C         -------------------------
C
C            Code XCASTP indicates that constant angular
C            size stepping is performed. This algorithm is
C            applicable only to class 1 surfaces.
C
      DOUBLE PRECISION      XCASTP
      PARAMETER           ( XCASTP = 1 )

C
C            Code XCPSTP indicates that constant orthogonal projection
C            length stepping is performed. This algorithm is applicable
C            only to surfaces associated with rectangular coordinates.
C
      DOUBLE PRECISION      XCPSTP
      PARAMETER           ( XCPSTP = 2 )

C
C            Code XCELBD denotes stepping from one cell boundary
C            to the next, where a "cell" is region over which
C            the same set of data are used for interpolation. 
C
C            Examples of cells:
C
C               1) When no interpolation is used, cell boundaries
C                  coincide with pixel boundaries.
C
C               2) When bilinear interpolation is used, cell boundaries
C                  are curves, with one of the first two coordinates
C                  held constant, connecting centers of adjacent pixels.
C                  For example, if the coordinate system is latitudinal,
C                  such curves have either constant longitude or 
C                  latitude.
C
C               3) When 3-nearest-neighbor interpolation is used,
C                  cells are triangular regions connecting centers
C                  of the 3 pixels used for interpolation.
C
C
      DOUBLE PRECISION      XCELBD
      PARAMETER           ( XCELBD = 3 )


C
C         Intercept descriptor index parameters  
C         -------------------------------------
C
C         Index of intercept algorithm code within descriptor:
C
      INTEGER               XXCD
      PARAMETER           ( XXCD   = 1 )

C
C         Index of descriptor size within descriptor:
C
      INTEGER               XXSZ
      PARAMETER           ( XXSZ   = XXCD + 1 )

C
C         Index of step parameter within descriptor:
C
      INTEGER               XXSTEP
      PARAMETER           ( XXSTEP  = XXSZ + 1 )

C
C         Index of tolerance parameter within descriptor:
C
      INTEGER               XXTOL
      PARAMETER           ( XXTOL  = XXSTEP + 1 )


C
C     Intercept acceleration algorithm descriptor
C     ===========================================
C

C
C         Intercept acceleration algorithm codes
C         --------------------------------------
C
C
C         Code ACCNO indicates no acceleration algorithm
C         is used.
C
      INTEGER               ACCNO
      PARAMETER           ( ACCNO  = 1 )
C
C         Code ACCBIG indicatest that bounds on height
C         differences in the COORD1 and COORD2 directions
C         are used to calculate a "big" initial step 
C         during the root bracketing process. These bounds
C         are determined during the segment writing process.
C        
C
      INTEGER               ACCBIG
      PARAMETER           ( ACCBIG  = ACCNO + 1 )

C
C         Intercept acceleration descriptor index parameters  
C         --------------------------------------------------

C
C         Index of acceleration algorithm code within descriptor:
C
      INTEGER               AXCD
      PARAMETER           ( AXCD   = 1 )

C
C         Index of descriptor size within descriptor:
C
      INTEGER               AXSZ
      PARAMETER           ( AXSZ   = AXCD + 1 )

C
C         Index of maximum magnitude of height difference
C         in the COORD1 direction:
C
      INTEGER               AXMXD1
      PARAMETER           ( AXMXD1 = AXSZ + 1 )

C
C         Index of maximum magnitude of height difference
C         in the COORD2 direction:
C
      INTEGER               AXMXD2
      PARAMETER           ( AXMXD2 = AXMXD1 + 1 )


C
C     Normal vector calculation algorithm descriptor
C     ==============================================
C
C
C     The normal vector algorithm descriptor is a variable-size set of
C     parameters that specify the normal vector computation algorithm
C     used for the associated segment.
C
C     The descriptor starts with a code identifying the class of
C     algorithm. This is followed by the size of the descriptor. If
C     applicable, additional parameters follow the size.
C
C         Normal vector algorithm codes
C         -----------------------------
C
C            Code NVBCPT indicates that normal vectors are compatible
C            with binlinear interpolation. This algorithm may be used
C            only when the interpolation method is bilinear.
C
C            Normal vector directions at boundaries between
C            interpolation regions---that is, where at least one pixel
C            coordinate is integral---are guaranteed to be compatible
C            with one of the adjacent regions, but are otherwise
C            unspecified. Note that discontinuity of normal vectors
C            should be expected at these boundaries.
C
      DOUBLE PRECISION      NVBCPT
      PARAMETER           ( NVBCPT = 1 )
C
C
C            Code NV1PAR indicates that estimates of the 1st order
C            partial derivatives of the third coordinate in the coord 1
C            and coord 2 directions are used to estimate the normal
C            vector direction.
C
      DOUBLE PRECISION      NV1PAR
      PARAMETER           ( NV1PAR = NVBCPT + 1 )

C
C         Normal vector descriptor index parameters  
C         -------------------------------------
C
C         Index of normal vector algorithm code within descriptor:
C
      INTEGER               NVXCD
      PARAMETER           ( NVXCD   = 1 )

C
C         Index of descriptor size within descriptor:
C
      INTEGER               NVXSZ
      PARAMETER           ( NVXSZ   = NVXCD + 1 )

C
C         Index of first parameter within descriptor:
C
      INTEGER               NVXPR
      PARAMETER           ( NVXPR   = NVXSZ + 1 )



C
C     Additional DSK type 4 parameters
C     ================================
C     
C
C     MARGIN is used to determine when computed coordinates are
C     close enough to the boundaries of a coverage region to be
C     considered to be located within that region.
C
      DOUBLE PRECISION      MARGIN
      PARAMETER           ( MARGIN = 1.D-12 )

C 
C     End of include file dsk04.inc
C

