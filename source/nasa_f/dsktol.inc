
C
C     File: dsktol.inc
C
C
C     This file contains declarations of tolerance and margin values
C     used by the DSK subsystem.
C
C     It is recommended that the default values defined in this file be
C     changed only by expert SPICE users.
C
C     The values declared in this file are accessible at run time
C     through the routines
C
C        DSKGTL  {DSK, get tolerance value}
C        DSKSTL  {DSK, set tolerance value}
C
C     These are entry points of the routine DSKTOL.
C
C        Version 1.0.0 27-FEB-2016 (NJB)
C
C
C
C
C     Parameter declarations
C     ======================
C
C     DSK type 2 plate expansion factor
C     ---------------------------------
C
C     The factor XFRACT is used to slightly expand plates read from DSK
C     type 2 segments in order to perform ray-plate intercept
C     computations.
C
C     This expansion is performed to prevent rays from passing through
C     a target object without any intersection being detected. Such
C     "false miss" conditions can occur due to round-off errors.
C
C     Plate expansion is done by computing the difference vectors
C     between a plate's vertices and the plate's centroid, scaling
C     those differences by (1 + XFRACT), then producing new vertices by
C     adding the scaled differences to the centroid. This process
C     doesn't affect the stored DSK data.
C
C     Plate expansion is also performed when surface points are mapped
C     to plates on which they lie, as is done for illumination angle
C     computations.
C
C     This parameter is user-adjustable. 
C
      DOUBLE PRECISION      XFRACT
      PARAMETER           ( XFRACT = 1.D-10 )

C
C     The keyword for setting or retrieving this factor is
C
      INTEGER               KEYXFR
      PARAMETER           ( KEYXFR = 1 )

C
C     Greedy segment selection factor
C     -------------------------------
C
C     The factor SGREED is used to slightly expand DSK segment
C     boundaries in order to select segments to consider for
C     ray-surface intercept computations. The effect of this factor is
C     to make the multi-segment intercept algorithm consider all
C     segments that are sufficiently close to the ray of interest, even
C     if the ray misses those segments. 
C
C     This expansion is performed to prevent rays from passing through
C     a target object without any intersection being detected. Such
C     "false miss" conditions can occur due to round-off errors.
C
C     The exact way this parameter is used is dependent on the
C     coordinate system of the segment to which it applies, and the DSK
C     software implementation. This parameter may be changed in a
C     future version of SPICE.
C     
      DOUBLE PRECISION      SGREED
      PARAMETER           ( SGREED = 1.D-8 )

C
C     The keyword for setting or retrieving this factor is
C
      INTEGER               KEYSGR
      PARAMETER           ( KEYSGR = KEYXFR + 1 )

C
C     Segment pad margin
C     ------------------
C
C     The segment pad margin is a scale factor used to determine when a
C     point resulting from a ray-surface intercept computation, if
C     outside the segment's boundaries, is close enough to the segment
C     to be considered a valid result.
C
C     This margin is required in order to make DSK segment padding
C     (surface data extending slightly beyond the segment's coordinate
C     boundaries) usable: if a ray intersects the pad surface outside
C     the segment boundaries, the pad is useless if the intercept is
C     automatically rejected.
C
C     However, an excessively large value for this parameter is
C     detrimental, since a ray-surface intercept solution found "in" a
C     segment can supersede solutions in segments farther from the
C     ray's vertex. Solutions found outside of a segment thus can mask
C     solutions that are closer to the ray's vertex by as much as the
C     value of this margin, when applied to a segment's boundary
C     dimensions.

      DOUBLE PRECISION      SGPADM
      PARAMETER           ( SGPADM = 1.D-10 )

C
C     The keyword for setting or retrieving this factor is
C
      INTEGER               KEYSPM
      PARAMETER           ( KEYSPM = KEYSGR + 1 )
      
C
C     Surface-point membership margin
C     -------------------------------
C
C     The surface-point membership margin limits the distance 
C     between a point and a surface to which the point is 
C     considered to belong. The margin is a scale factor applied
C     to the size of the segment containing the surface.
C
C     This margin is used to map surface points to outward
C     normal vectors at those points. 
C
C     If this margin is set to an excessively small value, 
C     routines that make use of the surface-point mapping won't
C     work properly.
C
      DOUBLE PRECISION      PTMEMM
      PARAMETER           ( PTMEMM = 1.D-7 )

C
C     The keyword for setting or retrieving this factor is
C
      INTEGER               KEYPTM
      PARAMETER           ( KEYPTM = KEYSPM + 1 )

C
C     Angular rounding margin
C     -----------------------
C
C     This margin specifies an amount by which angular values
C     may deviate from their proper ranges without a SPICE error
C     condition being signaled.
C
C     For example, if an input latitude exceeds pi/2 radians by a
C     positive amount less than this margin, the value is treated as
C     though it were pi/2 radians.
C
C     Units are radians.
C     
      DOUBLE PRECISION      ANGMRG
      PARAMETER           ( ANGMRG = 1.D-12 )

C
C     This parameter is not user-adjustable. 
C
C     The keyword for retrieving this parameter is
C
      INTEGER               KEYAMG
      PARAMETER           ( KEYAMG = KEYPTM + 1 )


C
C     Longitude alias margin
C     ----------------------
C
C     This margin specifies an amount by which a longitude
C     value can be outside a given longitude range without
C     being considered eligible for transformation by 
C     addition or subtraction of 2*pi radians. 
C
C     A longitude value, when compared to the endpoints of
C     a longitude interval, will be considered to be equal
C     to an endpoint if the value is outside the interval
C     differs from that endpoint by a magnitude less than
C     the alias margin.
C
C
C     Units are radians.
C     
      DOUBLE PRECISION      LONALI
      PARAMETER           ( LONALI = 1.D-12 )

C
C     This parameter is not user-adjustable. 
C
C     The keyword for retrieving this parameter is
C
      INTEGER               KEYLAL
      PARAMETER           ( KEYLAL = KEYAMG + 1 )

C
C     End of include file dsktol.inc
C
