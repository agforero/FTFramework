C$Procedure ZZSFXCOR ( Ray-surface intercept core algorithm )
 
      SUBROUTINE ZZSFXCOR ( UDNEAR, UDMAXR, UDRAYX, TRGCDE,  
     .                      ET,     ABCORR, USELT,  USECN,        
     .                      USESTL, XMIT,   FIXREF, OBSCDE,   
     .                      DFRCDE, DCLASS, DCENTR, DVEC,    
     .                      SPOINT, TRGEPC, SRFVEC, FOUND  )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find the intersection of a ray and the surface described by
C     an ellipsoid or DSK data, where the ray's vertex and the
C     target body's center are associated with ephemeris objects.
C     Use specified aberration corrections. Use callback routines 
C     to carry out low-level ray-surface intercept and ray-surface
C     near point computations.
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
C     GEOMETRY
C     INTERCEPT
C     INTERSECTION
C     RAY
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'frmtyp.inc'
      INCLUDE 'zzabcorr.inc'
      

      EXTERNAL              UDNEAR
      EXTERNAL              UDMAXR
      EXTERNAL              UDRAYX
      INTEGER               TRGCDE
      DOUBLE PRECISION      ET
      CHARACTER*(*)         ABCORR
      LOGICAL               USELT
      LOGICAL               USECN
      LOGICAL               USESTL
      LOGICAL               XMIT
      CHARACTER*(*)         FIXREF
      INTEGER               OBSCDE
      INTEGER               DFRCDE
      INTEGER               DCLASS
      INTEGER               DCENTR
      DOUBLE PRECISION      DVEC    ( 3 )
      DOUBLE PRECISION      SPOINT  ( 3 )
      DOUBLE PRECISION      TRGEPC
      DOUBLE PRECISION      SRFVEC  ( 3 )
      LOGICAL               FOUND

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     UDNEAR     I   Callback routine for ray-surface near point.
C     UDMAXR     I   Callback routine for target outer bounding surface.
C     UDRAYX     I   Callback routine for ray-surface intercept.
C     TRGCDE     I   Target body ID code.
C     ET         I   Epoch at observer.
C     ABCORR     I   Aberration correction string.
C     USELT      I   "Use light time" flag.
C     USECN      I   "Use converged Newtonian correction" flag.
C     USESTL     I   "Use stellar aberration" flag.
C     XMIT       I   "Perform transmission corrections" flag.
C     FIXREF     I   Name of target body-fixed frame.
C     OBSCDE     I   Observer body ID code.
C     DFRCDE     I   Ray direction frame ID code.
C     DCLASS     I   Ray direction frame class.
C     DCENTR     I   Ray direction frame center body ID code.
C     DVEC       I   Ray direction vector.
C     SPOINT     O   Surface intercept, if found.
C     TRGEPC     O   Target epoch, if intercept was found.
C     SRFVEC     O   Observer-to-intercept vector, if found.
C     FOUND      O   Found flag. 
C
C$ Detailed_Input
C
C     UDNEAR     is an external routine that computes the nearest point
C                on a bounding surface to the input ray. The caller
C                must perform any initialization required by UDNEAR.
C
C                The calling sequence of UDNEAR is
C
C                    CALL UDNEAR ( VERTEX, RAYDIR, ET, PNEAR, DIST )
C
C                    DOUBLE PRECISION      VERTEX(3)
C                    DOUBLE PRECISION      RAYDIR(3)
C                    DOUBLE PRECISION      ET
C                    DOUBLE PRECISION      PNEAR (3)
C                    DOUBLE PRECISION      DIST
C
C
C     UDMAXR     is an external routine that computes the radius of an
C                outer bounding surface for the target body. The caller
C                must perform any initialization required by UDMAXR.
C
C                The calling sequence of UDMAXR is
C
C                   CALL UDMAXR ( MAXRAD )
C
C                   DOUBLE PRECISION      MAXRAD
C
C
C     UDRAYX     external routine that computes the intercept of the
C                input ray on the surface associated with the target
C                body. The surface may be represented by an ellipsoid
C                or by DSK data. The caller must perform any
C                initialization required by UDRAYX.
C
C                The calling sequence of UDRAYX is
C
C                    CALL UDRAYX ( VERTEX, RAYDIR, ET, SPOINT, FOUND )
C
C                    DOUBLE PRECISION      VERTEX(3)
C                    DOUBLE PRECISION      RAYDIR(3)
C                    DOUBLE PRECISION      ET
C                    DOUBLE PRECISION      SPOINT(3)
C                    LOGICAL               FOUND
C
C
C     TRGCDE     is the body ID of the target.
C
C     ET         is the epoch of the intersection computation,
C                expressed as seconds past J2000 TDB. This epoch
C                applies to the observer.
C
C     ABCORR     is an aberration correction string. This string
C                is present in addition to the following flags
C                because it's used as an input to SPKEZP.
C
C     USELT      is a logical flag that is .TRUE. if and only if
C                light time corrections are to be performed.
C
C     USECN      is a logical flag that is .TRUE. if and only if
C                converged Newtonian light time corrections are to be
C                performed.
C
C     USESTL     is a logical flag that is .TRUE. if and only if
C                stellar aberration corrections are to be performed.
C
C     XMIT       is a logical flag that is .TRUE. if and only if
C                the aberration corrections are to be performed for
C                the transmission case.
C
C     FIXREF     is the name of a reference frame fixed to, and
C                centered on, the target body.
C
C     OBSCDE     is the body ID code of the observer.
C
C     DFRCDE     is the frame ID code of the frame in which the ray's
C                direction vector is expressed. This may be any frame
C                supported by the SPICE system, including built-in
C                frames (documented in the Frames Required Reading) and
C                frames defined by a loaded frame kernel (FK). The
C                string DREF is case-insensitive, and leading and
C                trailing blanks in DREF are not significant.
C
C                When DRFCDE designates a non-inertial frame, the
C                orientation of the frame is evaluated at an epoch
C                dependent on the frame's center and, if the center is
C                not the observer, on the selected aberration
C                correction. See the description of the direction
C                vector DVEC for details.
C
C     DCLASS     is the frame class of the frame designated by DFRCDE.
C
C     DCENTR     is the body ID code of the center of the frame
C                designated by DFRCDE.
C
C     DVEC       Ray direction vector emanating from the observer. The
C                intercept with the target body's surface of the ray
C                defined by the observer and DVEC is sought.
C
C                DVEC is specified relative to the reference frame
C                designated by DREF.
C
C                Non-inertial reference frames are treated as follows:
C                if the center of the frame is at the observer's
C                location, the frame is evaluated at ET. If the frame's
C                center is located elsewhere, then letting LTCENT be
C                the one-way light time between the observer and the
C                central body associated with the frame, the
C                orientation of the frame is evaluated at ET-LTCENT,
C                ET+LTCENT, or ET depending on whether the requested
C                aberration correction is, respectively, for received
C                radiation, transmitted radiation, or is omitted.
C                LTCENT is computed using the method indicated by
C                ABCORR.
C
C 
C$ Detailed_Output
C
C     SPOINT      is the surface intercept point on the target body of
C                 the ray defined by the observer and the direction
C                 vector. If the ray intersects the target body in
C                 multiple points, the selected intersection point is
C                 the one closest to the observer. The output argument
C                 FOUND (see below) indicates whether an intercept was
C                 found.
C
C                 SPOINT is expressed in Cartesian coordinates,
C                 relative to the target body-fixed frame designated by
C                 FIXREF. The body-fixed target frame is evaluated at
C                 the intercept epoch TRGEPC (see description below).
C
C                 When light time correction is used, the duration of
C                 light travel between SPOINT to the observer is
C                 considered to be the one way light time. When both
C                 light time and stellar aberration corrections are
C                 used, SPOINT is compute such that, when the vector
C                 from the observer to SPOINT is corrected for light
C                 time and stellar aberration, the resulting vector
C                 lies on the ray defined by the observer's location
C                 and DVEC.
C
C                 The components of SPOINT are given in units of km.
C
C
C     TRGEPC      is the "intercept epoch." TRGEPC is defined as
C                 follows: letting LT be the one-way light time between
C                 the observer and the intercept point, TRGEPC is the
C                 epoch ET-LT, ET+LT, or ET depending on whether the
C                 requested aberration correction is, respectively, for
C                 received radiation, transmitted radiation, or
C                 omitted. LT is computed using the method indicated by
C                 ABCORR.
C
C                 TRGEPC is expressed as seconds past J2000 TDB.
C
C
C     SRFVEC      is the vector from the observer's position at ET to
C                 the aberration-corrected (or optionally, geometric)
C                 position of SPOINT, where the aberration corrections
C                 are specified by ABCORR. SRFVEC is expressed in the
C                 target body-fixed reference frame designated by
C                 FIXREF, evaluated at TRGEPC.
C  
C                 The components of SRFVEC are given in units of km.
C
C                 One can use the SPICELIB function VNORM to obtain the
C                 distance between the observer and SPOINT:
C
C                    DIST = VNORM ( SRFVEC )
C
C                 The observer's position OBSPOS, relative to the
C                 target body's center, where the center's position is
C                 corrected for aberration effects as indicated by
C                 ABCORR, can be computed via the call:
C
C                    CALL VSUB ( SPOINT, SRFVEC, OBSPOS )
C
C                 To transform the vector SRFVEC from a reference frame
C                 FIXREF at time TRGEPC to a time-dependent reference
C                 frame REF at time ET, the routine PXFRM2 should be
C                 called. Let XFORM be the 3x3 matrix representing the
C                 rotation from the reference frame FIXREF at time
C                 TRGEPC to the reference frame REF at time ET. Then
C                 SRFVEC can be transformed to the result REFVEC as
C                 follows:
C
C                     CALL PXFRM2 ( FIXREF, REF,    TRGEPC, ET, XFORM )
C                     CALL MXV    ( XFORM,  SRFVEC, REFVEC )
C
C                 The second example in the Examples header section
C                 below presents a complete program that demonstrates
C                 this procedure.
C     
C
C     FOUND       A logical flag indicating whether or not the ray
C                 intersects the target. If an intersection exists
C                 FOUND will be returned as .TRUE. If the ray misses
C                 the target, FOUND will be returned as .FALSE.
C
C
C$ Parameters
C
C     See the include file frmtyp.inc for frame class parameter
C     values.
C
C$ Exceptions
C
C     1)  If the observer-target distance is zero, the error
C         SPICE(NOSEPARATION) is signaled.
C
C     2)  If an error occurs while performing a SPICE kernel lookup, or
C         if an error occurs while performing a geometric computation,
C         the error SHOULD be signaled by a routine in the call tree of
C         this routine. However, see exception (3) below.
C 
C     3)  This routine assumes robust input checking has been
C         performed by the caller. This routine may fail in
C         an unspecified manner if such checking has not 
C         been performed.
C
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - SPK data: ephemeris data for target and observer must be
C          loaded. If aberration corrections are used, the states of
C          target and observer relative to the solar system barycenter
C          must be calculable from the available ephemeris data.
C          Typically ephemeris data are made available by loading one
C          or more SPK files via FURNSH.
C
C        - PCK data: if the computation method is specified as
C          "Ellipsoid," triaxial radii for the target body must be 
C          loaded into the kernel pool. Typically this is done by
C          loading a text PCK file via FURNSH.
C
C        - Further PCK data: rotation data for the target body must
C          be loaded. These may be provided in a text or binary PCK
C          file. 
C
C     The following data may be required:
C
C        - DSK data: if METHOD indicates that DSK data are to be used,
C          DSK files containing topographic data for the target body
C          must be loaded. If a surface list is specified, data for
C          at least one of the listed surfaces must be loaded.
C
C        - Surface name-ID associations: if surface names are specified
C          in METHOD, the association of these names with their
C          corresponding surface ID codes must be established by 
C          assignments of the kernel variables
C
C             NAIF_SURFACE_NAME
C             NAIF_SURFACE_CODE
C             NAIF_SURFACE_BODY
C
C          Normally these associations are made by loading a text
C          kernel containing the necessary assignments. An example
C          of such an assignment is
C
C             NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
C             NAIF_SURFACE_CODE += 1                    
C             NAIF_SURFACE_BODY += 499
C
C        - Frame data: if a frame definition is required to convert
C          the observer and target states to the body-fixed frame of
C          the target, that definition must be available in the kernel
C          pool. Similarly, the frame definition required to map
C          between the frame designated by DREF and the target
C          body-fixed frame must be available. Typically the
C          definitions of frames not already built-in to SPICE are
C          supplied by loading a frame kernel.
C
C        - CK data: if the frame to which DREF refers is fixed to a
C          spacecraft instrument or structure, at least one CK file
C          will be needed to permit transformation of vectors between
C          that frame and both the J2000 and the target body-fixed
C          frames.
C
C        - SCLK data: if a CK file is needed, an associated SCLK
C          kernel is required to enable conversion between encoded SCLK
C          (used to time-tag CK data) and barycentric dynamical time
C          (TDB).
C
C     In all cases, kernel data are normally loaded once per program
C     run, NOT every time this routine is called.
C
C     
C$ Particulars
C
C     This routine implements the core ray-surface intercept algorithm
C     of the SPICELIB routine SINCPT. Parsing and analysis of inputs
C     is largely performed by the caller of this routine. 
C
C     Certain operations performed by this routine, which formerly 
C     were performed in-line by SINCPT, are now handled by callback
C     routines. These are:
C
C         - Computing a ray-surface intercept, where the ray's vertex
C           is expressed in a known target body-fixed frame and 
C           represents an offset from the target body center (NOT
C           from the frame's center, as is the case in lower-level
C           routines.
C
C         - Computing the radius of an outer bounding surface for
C           the target body.
C
C         - Computing the nearest point to a ray on a specified
C           bounding surface for the target body. This computation
C           enables the intercept-observer light time to be estimated
C           for cases where the first iteration of the solution results
C           in non-intersection geometry.
C
C     By using these callback routines, ZZSFXCOR is able to use a
C     single logic case to compute ray-surface intercepts for targets
C     modeled as ellipsoids or as surfaces represented by topographic 
C     data.
C
C$ Examples
C
C     See usage in SINCPT.
C
C$ Restrictions
C
C     1)  This is a private routine. It is meant to be used only by the
C         DSK subsystem.
C
C     2)  This routine relies extensively on the calling routine
C         to check the inputs passed to this routine.
C      
C     3)  This routine relies on the calling routine to perform any
C         initialization functions required by the input 
C         callback routines.
C
C     4)  If the direction vector DVEC is the zero vector, the error
C         SPICE(ZEROVECTOR) will be signaled.
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
C-    SPICELIB Version 1.0.0, 21-FEB-2017 (NJB) 
C
C        Added FAILED tests.
C
C        04-OCT-2016 (NJB) 
C
C        Bug fix: initializes FOUND to .FALSE.
C
C        Bug fix: routine now allows observer to be inside
C        outer bounding sphere of target.
C
C        Re-named callback routines again.
C
C        01-JUN-2016 (NJB) 
C
C           Updated names of callback routines to avoid overlap with
C           names of actual SPICELIB routines.
C     
C        08-FEB-2016 (NJB) 
C
C           Based on version 29-JAN-2015 (NJB)
C
C        Version 1.0.1 29-JAN-2015 (NJB)
C
C           Cleaned up debugging comments.
C
C        Version 1.0.0 15-OCT-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     generalized ray-surface intercept
C     ray-surface intercept core algorithm
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      CLIGHT
      DOUBLE PRECISION      DASINE
      DOUBLE PRECISION      TOUCHD
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VNORM
      DOUBLE PRECISION      VSEP

      LOGICAL               FAILED
      LOGICAL               VZERO
      LOGICAL               RETURN

C
C     Local parameters
C
      CHARACTER*(*)         RNAME
      PARAMETER           ( RNAME  = 'ZZSFXCOR' )


C
C     This value will become system-dependent when systems
C     using 128-bit d.p. numbers are supported by SPICELIB.
C     CNVLIM, when added to 1.0D0, should yield 1.0D0. 
C
      DOUBLE PRECISION      CNVLIM
      PARAMETER           ( CNVLIM = 1.D-17 )

C
C     Round-off error limit for arc sine input:
C
      DOUBLE PRECISION      RNDTOL
      PARAMETER           ( RNDTOL = 1.D-14 )

C
C     Fraction of target body angular radius used to define
C     region outside of which rays are immediately rejected
C     as non-intersecting.
C
      DOUBLE PRECISION      MARGIN
      PARAMETER           ( MARGIN = 1.01D0 )


      INTEGER               MAXITR
      PARAMETER           ( MAXITR =  10 )

C
C     Code for the frame J2000.
C
      INTEGER               J2CODE
      PARAMETER           ( J2CODE = 1 )


C
C     Local variables
C
      CHARACTER*(CORLEN)    LOCCOR
      CHARACTER*(CORLEN)    PRVCOR

      DOUBLE PRECISION      DIST
      DOUBLE PRECISION      ETDIFF
      DOUBLE PRECISION      J2DIR  ( 3 )
      DOUBLE PRECISION      J2EST  ( 3 )
      DOUBLE PRECISION      J2GEOM ( 3 )
      DOUBLE PRECISION      J2POS  ( 3 )
      DOUBLE PRECISION      J2TMAT ( 3, 3 )
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      LTCENT
      DOUBLE PRECISION      LTDIFF
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      NEGPOS ( 3 )
      DOUBLE PRECISION      OBSPOS ( 3 )
      DOUBLE PRECISION      PNEAR  ( 3 )
      DOUBLE PRECISION      PREVET
      DOUBLE PRECISION      PREVLT
      DOUBLE PRECISION      R2JMAT ( 3, 3 )
      DOUBLE PRECISION      RANGE
      DOUBLE PRECISION      RAYALT
      DOUBLE PRECISION      REFEPC
      DOUBLE PRECISION      REJECT
      DOUBLE PRECISION      RELERR
      DOUBLE PRECISION      RPOS   ( 3 )
      DOUBLE PRECISION      S
      DOUBLE PRECISION      SRFLEN
      DOUBLE PRECISION      SSBOST ( 6 )
      DOUBLE PRECISION      SSBTST ( 6 )
      DOUBLE PRECISION      STLDIR ( 3 )
      DOUBLE PRECISION      STLERR ( 3 )
      DOUBLE PRECISION      STLTMP ( 3 )
      DOUBLE PRECISION      TPOS   ( 3 )
      DOUBLE PRECISION      TRGDIR ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      XFORM  ( 3, 3 )

      INTEGER               I
      INTEGER               NITR

      LOGICAL               FIRST

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  LOCCOR
      SAVE                  PRVCOR

C
C     Initial values
C
      DATA                  FIRST  / .TRUE. /
      DATA                  PRVCOR / ' '    /


C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( RNAME )

C
C     Nothing found yet.
C
      FOUND = .FALSE.

C
C     Check for a zero ray direction vector.
C     
      IF ( VZERO(DVEC) ) THEN

         CALL SETMSG ( 'Input ray direction was the zero vector; this '
     .   //            'vector must be non-zero.'                      )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                             )
         CALL CHKOUT ( RNAME                                           )
         RETURN

      END IF

C
C     Get the sign S prefixing LT in the expression for TRGEPC.
C     When light time correction is not used, setting S = 0
C     allows us to seamlessly set TRGEPC equal to ET.
C
      IF ( USELT ) THEN

         IF ( XMIT ) THEN
            S   =  1.D0
         ELSE
            S   = -1.D0
         END IF
         
      ELSE
         S = 0.D0
      END IF
 

      IF ( FIRST .OR. ( ABCORR .NE. PRVCOR ) ) THEN
C
C        Construct aberration correction string without stellar
C        aberration specification.
C
         IF ( USELT ) THEN

            IF ( XMIT ) THEN
               LOCCOR = 'X'
            ELSE
               LOCCOR = ' '
            END IF

            IF ( USECN ) THEN

               CALL SUFFIX ( 'CN', 0, LOCCOR )
            ELSE
               CALL SUFFIX ( 'LT', 0, LOCCOR )            
            END IF

         ELSE
            LOCCOR = 'NONE'
         END IF


         PRVCOR = ABCORR
         FIRST  = .FALSE.

      END IF

C
C     Determine the position of the observer in target
C     body-fixed coordinates.
C
C         -  Call SPKEZP to compute the position of the target body as
C            seen from the observing body and the light time (LT)
C            between them. We request that the coordinates of POS be
C            returned relative to the body fixed reference frame
C            associated with the target body, using aberration
C            corrections specified by LOCCOR; these are the corrections
C            the input argument ABCORR, minus the stellar aberration
C            correction if it was called for.
C
C         -  Call VMINUS to negate the direction of the vector (OBSPOS)
C            so it will be the position of the observer as seen from
C            the target body in target body fixed coordinates.
C
C            Note that this result is not the same as the result of
C            calling SPKEZP with the target and observer switched. We
C            computed the vector FROM the observer TO the target in
C            order to get the proper light time and stellar aberration
C            corrections (if requested). Now we need the inverse of
C            that corrected vector in order to compute the intercept
C            point.
C

      CALL SPKEZP ( TRGCDE, ET, FIXREF, LOCCOR, OBSCDE, TPOS, LT )
 
      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     Negate the target's position to obtain the position of the
C     observer relative to the target.
C
      CALL VMINUS ( TPOS, OBSPOS )

C
C     We now need to convert the direction vector into the
C     body fixed frame associated with the target. The target
C     epoch is dependent on the aberration correction. The
C     coefficient S has been set to give us the correct answer
C     for each case.
C
      TRGEPC = ET  +  S*LT

C
C     Transform the direction vector from frame DREF to the body-fixed
C     frame associated with the target. The epoch TRGEPC associated
C     with the body-fixed frame has been set already.
C     
C     We'll compute the transformation in two parts: first
C     from frame DREF to J2000, then from J2000 to the target
C     frame.
C
      IF ( DCLASS .EQ. INERTL ) THEN
C
C        Inertial frames can be evaluated at any epoch.
C
         REFEPC = ET


      ELSE IF ( .NOT. USELT ) THEN
C
C        We're not using light time corrections (converged or
C        otherwise), so there's no time offset.
C
         REFEPC = ET


      ELSE IF ( DCENTR .EQ. OBSCDE ) THEN
C
C        If the center of frame DREF is the observer (which is
C        usually the case if the observer is a spacecraft), then
C        the epoch of frame DREF is simply ET.
C
C        There's no offset between the center for frame DREF
C        and the observer.
C
         REFEPC = ET

      ELSE
C
C        Find the light time from the observer to the center of 
C        frame DREF.
C           
         CALL SPKEZP ( DCENTR, ET, 'J2000', ABCORR, OBSCDE,
     .                 RPOS,   LTCENT                       )
  
         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         REFEPC  =  ET  +  S * LTCENT

      END IF

C
C     The epoch REFEPC associated with frame DREF has been set.
C
C     Compute the transformation from frame DREF to J2000 and the
C     transformation from J2000 to the target body-fixed frame.
C
C     Map DVEC to both the J2000 and target body-fixed frames. We'll
C     store DVEC, expressed relative to the J2000 frame, in the
C     variable J2DIR. DVEC in the target body-fixed frame will be
C     stored in TRGDIR.
C
C     We may need both versions of DVEC: if we use light time
C     correction, we'll update "intercept epoch", and hence the
C     transformation between J2000 and the target body-fixed frame.
C     The transformation between DREF and J2000 doesn't change, on the
C     other hand, so we don't have to recompute J2DIR. We need TRGDIR
C     in all cases.
C
      CALL REFCHG ( DFRCDE, J2CODE, REFEPC, R2JMAT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      CALL MXV ( R2JMAT, DVEC, J2DIR )

C
C     Save this version of J2DIR as J2GEOM. Later we'll
C     modify J2DIR, if necessary, to account for stellar
C     aberration. 
C
      CALL VEQU ( J2DIR, J2GEOM )

C
C     Map J2DIR (in the J2000 frame) to the target body-fixed
C     frame.
C
      CALL PXFORM ( 'J2000', FIXREF,  TRGEPC,  J2TMAT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      CALL MXV ( J2TMAT, J2DIR, TRGDIR )

C
C     At this point,
C
C        TRGEPC is set.
C        TRGDIR is set.
C        J2DIR is set.
C
C
C     Get the J2000-relative state of the observer relative to
C     the solar system barycenter at ET. We'll use this in
C     several places later.
C
      CALL SPKSSB ( OBSCDE, ET, 'J2000', SSBOST )
 
      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     If we're using stellar aberration correction, at this point we'll
C     account for it. We're going to find a surface point such that
C     the radiation path from that point to the observer, after
C     correction for stellar aberration, is parallel to the ray. So
C     by applying the inverse of the correction to the ray, we obtain
C     the ray with which we must perform our intercept computation.
C 
      IF ( USESTL ) THEN
C
C        We approximate the inverse stellar aberration correction by
C        using the correction for the reverse transmission direction.
C        If we're in the reception case, we apply the transmission
C        stellar aberration correction to J2DIR and vice versa.
C      
C        We iterate our estimates until we have the desired level
C        of convergence or reach the iteration limit.
C
         NITR = 5

         IF ( XMIT ) THEN
C
C           Use reception stellar aberration correction
C           routine STELAB to generate a first estimate of
C           the direction vector after stellar aberration
C           has been "removed"---that is, apply the inverse
C           of the transmission stellar aberration correction
C           mapping to J2DIR.
C
            CALL STELAB ( J2DIR, SSBOST(4), STLDIR )

C
C           Now improve our estimate.
C
            RELERR = 1.D0
            I      = 1
            
            DO WHILE ( ( I .LE. NITR ) .AND. ( RELERR .GT. CNVLIM ) )
C
C              Estimate the error in our previous approximation 
C              by applying the reception stellar aberration
C              to STLDIR and finding the difference with J2DIR.
C
               CALL STLABX ( STLDIR, SSBOST(4), J2EST  )            
               CALL VSUB   ( J2DIR,  J2EST,     STLERR )

C
C              Adding the error in the reception mapping to STLDIR
C              will give us an improved estimate of the inverse.
C
               CALL VADD ( STLERR, STLDIR, STLTMP )
               CALL VEQU ( STLTMP,         STLDIR )

               RELERR = VNORM(STLERR) / VNORM(STLDIR)
               I      = I + 1

            END DO

C
C           At this point we've found a good estimate of the
C           direction vector under the inverse of the transmission
C           stellar aberration correction mapping.
C
         ELSE
C
C           Use transmission stellar aberration correction
C           routine STLABX to generate a first estimate of
C           the direction vector after stellar aberration
C           has been "removed."  
C
            CALL STLABX ( J2DIR, SSBOST(4), STLDIR )

C
C           Now improve our estimate.
C
            RELERR = 1.D0
            I      = 1
            
            DO WHILE ( ( I .LE. NITR ) .AND. ( RELERR .GT. CNVLIM ) )
C
C              Estimate the error in our previous approximation 
C              by applying the reception stellar aberration
C              to STLDIR and finding the difference with J2DIR.
C
               CALL STELAB ( STLDIR, SSBOST(4), J2EST  )            
               CALL VSUB   ( J2DIR,  J2EST,     STLERR )

C
C              Adding the error in the reception mapping to STLDIR
C              will give us an improved estimate of the inverse.
C
               CALL VADD ( STLERR, STLDIR, STLTMP )
               CALL VEQU ( STLTMP,         STLDIR )

               RELERR = VNORM(STLERR) / VNORM(STLDIR)
               I      = I + 1

            END DO

C
C           At this point we've found a good estimate of the
C           direction vector under the inverse of the reception
C           stellar aberration correction mapping.
C
         END IF

C
C        Replace the J2000-relative ray direction with the corrected
C        direction.
C
         CALL VEQU ( STLDIR,  J2DIR )
         CALL MXV  ( J2TMAT,  J2DIR,   TRGDIR )

      END IF

C
C     Find the surface intercept point and distance from observer to 
C     intercept point using the specified geometric definition.
C
C     Find the surface intercept given the target epoch,
C     observer-target position, and target body orientation we've
C     already computed. If we're not using light time correction, this
C     is all we must do. Otherwise, our result will give us an initial
C     estimate of the target epoch, which we'll then improve.
C
C     Make an easy test to see whether we can quit now because an
C     intercept cannot exist. If the ray is separated from the
C     observer-target center vector by more than (MARGIN * the maximum
C     target radius), we're done. Let REJECT be the angular
C     separation limit.
C        
      CALL UDMAXR ( MAXRAD )

      RANGE = VNORM( OBSPOS )

      IF ( RANGE .EQ. 0.D0 ) THEN
C
C        We've already ensured that observer and target are
C        distinct, so this should be a very unusual occurrence.
C
         CALL SETMSG ( 'Observer-target distance is zero. '
     .   //            'Observer ID is #; target ID is #.' )
         CALL ERRINT ( '#', OBSCDE                         )
         CALL ERRINT ( '#', TRGCDE                         )
         CALL SIGERR ( 'SPICE(NOSEPARATION)'               )
         CALL CHKOUT ( RNAME                               )
         RETURN

      END IF


      IF ( RANGE .GT. MARGIN*MAXRAD ) THEN
C
C        Compute the arc sine with SPICE error checking.
C
         REJECT = DASINE ( MARGIN * MAXRAD / RANGE,  RNDTOL ) 

         CALL VMINUS ( OBSPOS, NEGPOS )
         
         IF (  VSEP( NEGPOS, TRGDIR )  .GT.  REJECT  ) THEN
C
C           The angular separation of ray and target is too great
C           for a solution to exist, even with a better light time
C           estimate.
C
            CALL CHKOUT ( RNAME )
            RETURN

         END IF

      END IF


C
C     Locate the intercept of the ray with the target; if there's no
C     intercept, find the closest point on the target to the ray.
C
      CALL UDRAYX ( OBSPOS, TRGDIR, TRGEPC, SPOINT, FOUND )
   
      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     If we found an intercept, and if we're not using light time
C     corrections, we're almost done now. We still need SRFVEC.
C     SPOINT, TRGEPC, and FOUND have already been set.
C        
      IF (  FOUND  .AND.  ( .NOT. USELT )  ) THEN

         CALL VSUB ( SPOINT, OBSPOS, SRFVEC )

         CALL CHKOUT ( RNAME )
         RETURN

      END IF

C
C     From this point onward, we're dealing with a case calling for
C     light time and possibly stellar aberration corrections.
C
      IF ( .NOT. FOUND ) THEN
C
C        If there's no intercept, we're probably done. However, we need
C        to guard against the possibility that the ray does intersect
C        the ellipsoid but we haven't discovered it because our first
C        light time estimate was too poor.
C
C        We'll make an improved light time estimate as follows: Find
C        the nearest point on the ellipsoid to the ray. Find the light
C        time between the observer and this point.
C
C        If we're using converged Newtonian corrections, we iterate
C        this procedure up to three times.
C
         IF ( USECN ) THEN
            NITR = 3
         ELSE
            NITR = 1
         END IF

         I = 1

         DO WHILE (  ( I .LE. NITR ) .AND. ( .NOT. FOUND )  )

            CALL UDNEAR ( OBSPOS, TRGDIR, ET, PNEAR, RAYALT )
 
            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF


            LT =  VDIST ( OBSPOS, PNEAR ) / CLIGHT()

C
C           Use the new light time estimate to repeat the intercept 
C           computation.
C
            TRGEPC  =  ET  +  S*LT

C
C           Get the J2000-relative state of the target relative to
C           the solar system barycenter at the target epoch.
C
            CALL SPKSSB ( TRGCDE, TRGEPC, 'J2000', SSBTST )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

C
C           Find the position of the observer relative to the target.
C           Convert this vector from the J2000 frame to the target
C           frame at TRGEPC.
C
            CALL VSUB   ( SSBOST,  SSBTST, J2POS         )
            CALL PXFORM ( 'J2000', FIXREF, TRGEPC, XFORM )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF 
           
C
C           Convert the observer's position relative to the target
C           from the J2000 frame to the target frame at the target
C           epoch.
C
            CALL MXV ( XFORM, J2POS, OBSPOS )

C
C           Convert the ray's direction vector from the J2000 frame
C           to the target frame at the target epoch.
C   
            CALL MXV ( XFORM, J2DIR, TRGDIR )

C
C           Repeat the intercept computation.
C
            CALL UDRAYX ( OBSPOS, TRGDIR, TRGEPC, SPOINT, FOUND )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

            I = I + 1

         END DO
C
C        If there's still no intercept, we're done.
C
         IF ( .NOT. FOUND ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

      END IF

C
C     Making it to this point means we've got an intersection. 
C
C     Since we're using light time corrections, we're going to make
C     an estimate of light time to the intercept point, then re-do
C     our computation of the target position and orientation using
C     the new light time value.
C
      IF ( USECN ) THEN
         NITR = MAXITR
      ELSE
         NITR = 1
      END IF

C
C     Compute new light time estimate and new target epoch.
C
      DIST  =   VDIST ( OBSPOS, SPOINT )
      LT     =  DIST / CLIGHT()
      TRGEPC =  ET  +  S * LT

      PREVLT = 0.D0
      PREVET = TRGEPC

      I      = 0
      LTDIFF = 1.D0
      ETDIFF = 1.D0

      DO WHILE (       ( I      .LT.   NITR                ) 
     .           .AND. ( LTDIFF .GT. ( CNVLIM * ABS(LT) )  )
     .           .AND. ( ETDIFF .GT.   0.D0                )  )

C
C        Get the J2000-relative state of the target relative to
C        the solar system barycenter at the target epoch.
C
         CALL SPKSSB ( TRGCDE, TRGEPC, 'J2000', SSBTST )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        Find the position of the observer relative to the target.
C        Convert this vector from the J2000 frame to the target frame
C        at TRGEPC.
C
C        Note that SSBOST contains the J2000-relative state of the
C        observer relative to the solar system barycenter at ET.
C
         CALL VSUB   ( SSBOST,  SSBTST, J2POS         )
         CALL PXFORM ( 'J2000', FIXREF, TRGEPC, XFORM )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF
         
C
C        Convert the observer's position relative to the target from
C        the J2000 frame to the target frame at the target epoch.
C
         CALL MXV    ( XFORM,  J2POS, OBSPOS )
         CALL VMINUS ( OBSPOS, NEGPOS )

C
C        Convert the ray's direction vector from the J2000 frame
C        to the target frame at the target epoch.
C
         CALL MXV ( XFORM, J2DIR, TRGDIR )

C
C        Repeat the intercept computation.
C
         CALL UDRAYX ( OBSPOS, TRGDIR, TRGEPC, SPOINT, FOUND )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        If there's no intercept, we're done.
C
         IF ( .NOT. FOUND ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        Compute the distance between intercept and observer.
C
         DIST    =   VDIST ( OBSPOS, SPOINT )

C
C        Compute a new light time estimate and a new target epoch.
C
         LT      =  DIST / CLIGHT()

         TRGEPC  =  ET  +  S * LT
C
C        We use the d.p. identity function TOUCHD to force the compiler
C        to create double precision arguments from the differences
C        LT-PREVLT and TRGEPC-PREVET. Some compilers will perform
C        extended-precision register arithmetic, which can prevent a
C        difference from rounding to zero. Simply storing the result of
C        the subtraction in a double precision variable doesn't solve
C        the problem, because that variable can be optimized out of
C        existence.
C
         LTDIFF  =   ABS( TOUCHD(LT     - PREVLT) )
         ETDIFF  =   ABS( TOUCHD(TRGEPC - PREVET) )
         PREVLT  =   LT
         PREVET  =   TRGEPC        
         I       =   I + 1

      END DO
                   
C
C     FOUND, SPOINT, TRGEPC, and OBSPOS have been set at this point.
C     We need SRFVEC. Since OBSPOS doesn't take into account stellar
C     aberration, we can' derive SRFVEC from OBSPOS as is done in
C     the related routines SUBPNT and SUBSLR. Here, we derive
C     SRFVEC from J2GEOM, which is the input ray direction expressed in 
C     the J2000 frame. We use XFORM, which is computed in the loop
C     above, to convert J2GEOM to FIXREF, evaluated at TRGEPC.
C     
      CALL MXV    ( XFORM, J2GEOM, UDIR )
      CALL VHATIP ( UDIR )

C
C     Let SRFLEN be the length of SRFVEC; we CAN get this
C     length from OBSPOS and SPOINT, since stellar 
C     aberration correction (as implemented in SPICE)
C     doesn't change the length of the vector SPOINT-OBSPOS.
C     
      SRFLEN = VDIST ( SPOINT, OBSPOS )

C
C     Scale UDIR to obtain the desired value of SRFVEC.
C
      CALL VSCL ( SRFLEN, UDIR, SRFVEC )
      

      CALL CHKOUT ( RNAME )
      RETURN
      END
