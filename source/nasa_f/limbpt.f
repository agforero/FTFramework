C$Procedure LIMBPT ( Limb points on an extended object )
 
       SUBROUTINE LIMBPT ( METHOD, TARGET, ET,     FIXREF, ABCORR, 
     .                     CORLOC, OBSRVR, REFVEC, ROLSTP, NCUTS,  
     .                     SCHSTP, SOLTOL, MAXN,   NPTS,   POINTS,  
     .                     EPOCHS, TANGTS                         )
      
C$ Abstract
C
C     Find limb points on a target body. The limb is the set of points
C     of tangency on the target of rays emanating from the observer.
C     The caller specifies half-planes bounded by the observer-target
C     center vector in which to search for limb points.
C
C     The surface of the target body may be represented either by a
C     triaxial ellipsoid or by topographic data.
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
C     CK
C     DSK
C     FRAMES
C     NAIF_IDS
C     PCK
C     SPK
C     TIME
C
C$ Keywords
C
C     GEOMETRY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE               'dsk.inc'
      INCLUDE               'frmtyp.inc'
      INCLUDE               'gf.inc'
      INCLUDE               'zzabcorr.inc'
      INCLUDE               'zzctr.inc'
      INCLUDE               'zzdsk.inc'

      CHARACTER*(*)         METHOD
      CHARACTER*(*)         TARGET
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      CHARACTER*(*)         ABCORR
      CHARACTER*(*)         CORLOC
      CHARACTER*(*)         OBSRVR
      DOUBLE PRECISION      REFVEC ( 3 )
      DOUBLE PRECISION      ROLSTP
      INTEGER               NCUTS
      DOUBLE PRECISION      SCHSTP
      DOUBLE PRECISION      SOLTOL
      INTEGER               MAXN
      INTEGER               NPTS   ( * )
      DOUBLE PRECISION      POINTS ( 3, * )
      DOUBLE PRECISION      EPOCHS ( * )
      DOUBLE PRECISION      TANGTS ( 3, * )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     METHOD     I   Computation method.
C     TARGET     I   Name of target body.
C     ET         I   Epoch in ephemeris seconds past J2000 TDB.
C     FIXREF     I   Body-fixed, body-centered target body frame.
C     ABCORR     I   Aberration correction.
C     CORLOC     I   Aberration correction locus.
C     OBSRVR     I   Name of observing body.
C     REFVEC     I   Reference vector for cutting half-planes.
C     ROLSTP     I   Roll angular step for cutting half-planes.
C     NCUTS      I   Number of cutting half-planes.
C     SCHSTP     I   Angular step size for searching.
C     SOLTOL     I   Solution convergence tolerance.
C     MAXN       I   Maximum number of entries in output arrays.
C     NPTS       O   Counts of limb points corresponding to cuts.
C     POINTS     O   Limb points.
C     EPOCHS     O   Times associated with limb points.
C     TANGTS     O   Tangent vectors emanating from the observer.
C     
C$ Detailed_Input
C
C     METHOD   is a short string providing parameters defining
C              the computation method to be used. In the syntax
C              descriptions below, items delimited by brackets
C              are optional.
C
C              METHOD may be assigned the following values:
C
C                'TANGENT/DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                    The limb point computation uses topographic data
C                    provided by DSK files (abbreviated as "DSK data"
C                    below) to model the surface of the target body. A
C                    limb point is defined as the point of tangency, on
C                    the surface represented by the DSK data, of a ray
C                    emanating from the observer.
C
C                    Limb points are generated within a specified set
C                    of "cutting" half-planes that have as an edge the
C                    line containing the observer-target vector.
C                    Multiple limb points may be found within a given
C                    half-plane, if the target body shape allows for
C                    this.
C
C                    The surface list specification is optional. The
C                    syntax of the list is
C
C                       <surface 1> [, <surface 2>...]
C
C                    If present, it indicates that data only for the
C                    listed surfaces are to be used; however, data need
C                    not be available for all surfaces in the list. If
C                    the list is absent, loaded DSK data for any
C                    surface associated with the target body are used.
C
C                    The surface list may contain surface names or
C                    surface ID codes. Names containing blanks must
C                    be delimited by double quotes, for example
C
C                       SURFACES = "Mars MEGDR 128 PIXEL/DEG"
C
C                    If multiple surfaces are specified, their names
C                    or IDs must be separated by commas.
C
C                    See the Particulars section below for details
C                    concerning use of DSK data.
C
C                    This is the highest-accuracy method supported by
C                    this subroutine. It generally executes much more
C                    slowly than the GUIDED method described below.
C                    
C                    
C                'GUIDED/DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                    This method uses DSK data as described above, but
C                    limb points generated by this method are "guided"
C                    so as to lie in the limb plane of the target
C                    body's reference ellipsoid, on the target body's
C                    surface. This method produces a unique limb point
C                    for each cutting half-plane. If multiple limb
C                    point candidates lie in a given cutting
C                    half-plane, the outermost one is chosen.
C
C                    This method may be used only with the CENTER
C                    aberration correction locus (see the description
C                    of REFLOC below).
C
C                    Limb points generated by this method are
C                    approximations; they are generally not true
C                    ray-surface tangent points. However, these
C                    approximations can be generated much more quickly
C                    than tangent points.
C
C
C                'TANGENT/ELLIPSOID'
C                'GUIDED/ELLIPSOID'
C
C                    Both of these methods generate limb points on the
C                    target body's reference ellipsoid. The TANGENT
C                    option may be used with any aberration correction
C                    locus, while the GUIDED option may be used only
C                    with the CENTER locus (see the description of
C                    REFLOC below). 
C
C                    When the locus is set to 'CENTER', these methods
C                    produce the same results.
C
C
C                 Neither case nor white space are significant in
C                 METHOD, except within double-quoted strings. For
C                 example, the string ' eLLipsoid/tAnGenT ' is valid.
C
C                 Within double-quoted strings, blank characters are
C                 significant, but multiple consecutive blanks are
C                 considered equivalent to a single blank. Case is 
C                 not significant. So
C
C                    "Mars MEGDR 128 PIXEL/DEG"
C
C                 is equivalent to 
C
C                    " mars megdr  128  pixel/deg "
C
C                 but not to
C
C                    "MARS MEGDR128PIXEL/DEG"
C
C                 
C     TARGET      is the name of the target body. The target body is 
C                 an extended ephemeris object.
C
C                 The string TARGET is case-insensitive, and leading
C                 and trailing blanks in TARGET are not significant.
C                 Optionally, you may supply a string containing the
C                 integer ID code for the object. For example both
C                 'MOON' and '301' are legitimate strings that indicate
C                 the Moon is the target body.
C
C                 When the target body's surface is represented by a
C                 tri-axial ellipsoid, this routine assumes that a
C                 kernel variable representing the ellipsoid's radii is
C                 present in the kernel pool. Normally the kernel
C                 variable would be defined by loading a PCK file.
C
C
C     ET          is the epoch of participation of the observer,
C                 expressed as TDB seconds past J2000 TDB: ET is
C                 the epoch at which the observer's state is computed.
C
C                 When aberration corrections are not used, ET is also
C                 the epoch at which the position and orientation of
C                 the target body are computed.
C
C                 When aberration corrections are used, the position
C                 and orientation of the target body are computed at
C                 ET-LT, where LT is the one-way light time between the
C                 aberration correction locus and the observer. The
C                 locus is specified by the input argument CORLOC.
C                 See the descriptions of ABCORR and CORLOC below for
C                 details.
C
C
C     FIXREF      is the name of a body-fixed reference frame centered
C                 on the target body. FIXREF may be any such frame
C                 supported by the SPICE system, including built-in
C                 frames (documented in the Frames Required Reading)
C                 and frames defined by a loaded frame kernel (FK). The
C                 string FIXREF is case-insensitive, and leading and
C                 trailing blanks in FIXREF are not significant.
C
C                 The output limb points in the array POINTS and the
C                 output observer-target tangent vectors in the array
C                 TANGTS are expressed relative to this reference frame.
C
C
C     ABCORR      indicates the aberration corrections to be applied
C                 when computing the target's position and orientation.
C                 Corrections are applied at the location specified by
C                 the aberration correction locus argument CORLOC,
C                 which is described below.
C
C                 For remote sensing applications, where apparent limb
C                 points seen by the observer are desired, normally
C                 either of the corrections
C              
C                    'LT+S' 
C                    'CN+S'
C     
C                 should be used. The correction 'NONE' may be suitable
C                 for cases in which the target is very small and the
C                 observer is close to, and has small velocity relative
C                 to, the target (e.g. comet Churyumov-Gerasimenko and
C                 the Rosetta Orbiter).
C
C                 These and the other supported options are described
C                 below. ABCORR may be any of the following:
C
C                    'NONE'     Apply no correction. Return the
C                               geometric limb points on the target
C                               body.
C
C                 Let LT represent the one-way light time between the
C                 observer and the aberration correction locus. The
C                 following values of ABCORR apply to the "reception"
C                 case in which photons depart from the locus at the
C                 light-time corrected epoch ET-LT and *arrive* at the
C                 observer's location at ET:
C
C
C                    'LT'       Correct for one-way light time (also
C                               called "planetary aberration") using a
C                               Newtonian formulation. This correction
C                               yields the locus at the moment it
C                               emitted photons arriving at the
C                               observer at ET.
C 
C                               The light time correction uses an
C                               iterative solution of the light time
C                               equation. The solution invoked by the
C                               'LT' option uses one iteration.
C
C                               Both the target position as seen by the
C                               observer, and rotation of the target
C                               body, are corrected for light time.
C
C                    'LT+S'     Correct for one-way light time and
C                               stellar aberration using a Newtonian
C                               formulation. This option modifies the
C                               locus obtained with the 'LT' option to
C                               account for the observer's velocity
C                               relative to the solar system
C                               barycenter. These corrections yield
C                               points on the apparent limb.
C
C                    'CN'       Converged Newtonian light time
C                               correction. In solving the light time
C                               equation, the 'CN' correction iterates
C                               until the solution converges. Both the
C                               position and rotation of the target
C                               body are corrected for light time.
C
C                    'CN+S'     Converged Newtonian light time and
C                               stellar aberration corrections. This
C                               option produces a solution that is at
C                               least as accurate at that obtainable
C                               with the `LT+S' option. Whether the
C                               'CN+S' solution is substantially more
C                               accurate depends on the geometry of the
C                               participating objects and on the
C                               accuracy of the input data. In all
C                               cases this routine will execute more
C                               slowly when a converged solution is
C                               computed.
C
C
C     CORLOC      is a string specifying the aberration correction
C                 locus: the point or set of points for which
C                 aberration corrections are performed. CORLOC may be
C                 assigned the values:
C
C                    'CENTER' 
C
C                        Light time and stellar aberration corrections
C                        are applied to the vector from the observer to
C                        the center of the target body. The one way
C                        light time from the target center to the
C                        observer is used to determine the epoch at
C                        which the target body orientation is computed.
C
C                        This choice is appropriate for small target
C                        objects for which the light time from the
C                        surface to the observer varies little across
C                        the entire target. It may also be appropriate
C                        for large, nearly ellipsoidal targets when the
C                        observer is very far from the target.
C
C                        Computation speed for this option is faster
C                        than for the ELLIPSOID LIMB option.
C
C                    'ELLIPSOID LIMB'
C
C                        Light time and stellar aberration corrections
C                        are applied to individual limb points on the
C                        reference ellipsoid. For a limb point on the
C                        surface described by topographic data, lying
C                        in a specified cutting half-plane, the unique
C                        reference ellipsoid limb point in the same
C                        half-plane is used as the locus of the
C                        aberration corrections.
C
C                        This choice is appropriate for large target
C                        objects for which the light time from the limb
C                        to the observer is significantly different
C                        from the light time from the target center to
C                        the observer.
C
C                        Because aberration corrections are repeated for
C                        individual limb points, computational speed for
C                        this option is relatively slow.
C
C
C     OBSRVR      is the name of the observing body. The observing body
C                 is an ephemeris object: it typically is a spacecraft,
C                 the earth, or a surface point on the earth. OBSRVR is
C                 case-insensitive, and leading and trailing blanks in
C                 OBSRVR are not significant. Optionally, you may
C                 supply a string containing the integer ID code for
C                 the object. For example both 'MOON' and '301' are
C                 legitimate strings that indicate the Moon is the
C                 observer.
C
C
C     REFVEC,
C     ROLSTP,
C     NCUTS       are, respectively, a reference vector, a roll step
C                 angle, and a count of cutting half-planes.
C
C                 REFVEC defines the first of a sequence of cutting
C                 half-planes in which limb points are to be found.
C                 Each cutting half-plane has as its edge the line
C                 containing the observer-target vector; the first
C                 half-plane contains REFVEC.
C
C                 REFVEC is expressed in the body-fixed reference frame
C                 designated by FIXREF.
C
C                 ROLSTP is an angular step by which to roll the
C                 cutting half-planes about the observer-target vector.
C                 The first half-plane is aligned with REFVEC; the Ith
C                 half-plane is rotated from REFVEC about the
C                 observer-target vector in the counter-clockwise
C                 direction by (I-1)*ROLSTP. Units are radians.
C                 ROLSTP should be set to 
C
C                    2*pi/NCUTS 
C
C                 to generate an approximately uniform distribution of
C                 limb points along the limb.
C
C                 NCUTS is the number of cutting half-planes used to
C                 find limb points; the angular positions of
C                 consecutive half-planes increase in the positive
C                 sense (counterclockwise) about the target-observer
C                 vector and are distributed roughly equally about that
C                 vector: each half-plane has angular separation of
C                 approximately
C
C                    ROLSTP radians
C
C                 from each of its neighbors. When the aberration
C                 correction locus is set to 'CENTER', the angular
C                 separation is the value above, up to round-off. When
C                 the locus is 'ELLIPSOID LIMB', the separations are
C                 less uniform due to differences in the aberration
C                 corrections used for the respective limb points.
C
C
C     SCHSTP,
C     SOLTOL      are used only for DSK-based surfaces. These inputs
C                 are, respectively, the search angular step size and
C                 solution convergence tolerance used to find tangent
C                 rays and associated limb points within each cutting
C                 half plane. These values are used when the METHOD
C                 argument includes the TANGENT option. In this case,
C                 limb points are found by a two-step search process:
C
C                    1) Bracketing: starting with the direction
C                       opposite the observer-target vector, rays
C                       emanating from the observer are generated
C                       within the half-plane at successively greater
C                       angular separations from the initial direction,
C                       where the increment of angular separation is
C                       SCHSTP. The rays are tested for intersection
C                       with the target surface. When a transition
C                       between non-intersection to intersection is
C                       found, the angular separation of a tangent ray
C                       has been bracketed.
C
C                    2) Root finding: each time a tangent ray is 
C                       bracketed, a search is done to find the angular
C                       separation from the starting direction at which
C                       a tangent ray exists. The search terminates
C                       when successive rays are separated by no more
C                       than SOLTOL. When the search converges, the
C                       last ray-surface intersection point found in
C                       the convergence process is considered to be a
C                       limb point.
C                     
C     
C                  SCHSTP and SOLTOL have units of radians.
C
C                  Target bodies with simple surfaces---for example,
C                  convex shapes---will have a single limb point within
C                  each cutting half-plane. For such surfaces, SCHSTP
C                  can be set large enough so that only one bracketing
C                  step is taken. A value greater than pi, for example
C                  4.D0, is recommended.
C
C                  Target bodies with complex surfaces can have
C                  multiple limb points within a given cutting
C                  half-plane. To find all limb points, SCHSTP must be
C                  set to a value smaller than the angular separation
C                  of any two limb points in any cutting half-plane,
C                  where the vertex of the angle is the observer.
C                  SCHSTP must not be too small, or the search will be
C                  excessively slow.
C
C                  For both kinds of surfaces, SOLTOL must be chosen so
C                  that the results will have the desired precision.
C                  Note that the choice of SOLTOL required to meet a
C                  specified bound on limb point height errors depends
C                  on the observer-target distance.
C
C
C     MAXN         is the maximum number of limb points that can be
C                  stored in the output array POINTS.
C
C
C$ Detailed_Output
C
C
C     NPTS         is an array of counts of limb points within the
C                  specified set of cutting half-planes. The Ith
C                  element of NPTS is the limb point count in the Ith
C                  half-plane. NPTS should be declared with length
C                  at least NCUTS.
C
C                  For most target bodies, there will be one limb point
C                  per half-plane. For complex target shapes, the limb
C                  point count in a given half-plane can be greater
C                  than one (see example 3 below), and it can be zero.
C
C
C     POINTS       is an array containing the limb points found by this
C                  routine. Sets of limb points associated with
C                  half-planes are ordered by the indices of the
C                  half-planes in which they're found. The limb points
C                  in a given half-plane are ordered by decreasing
C                  angular separation from the observer-target
C                  direction; the outermost limb point in a given
C                  half-plane is the first of that set.
C
C                  The limb points for the half-plane containing REFVEC
C                  occupy array elements
C
C                     POINTS(1,1) through POINTS(3,NPTS(1))
C
C                  Limb points for the second half plane occupy
C                  elements
C
C                     POINTS(1, NPTS(1)+1       ) through 
C                     POINTS(3, NPTS(1)+NPTS(2) )
C
C                  and so on.
C
C                  POINTS should be declared with dimensions
C
C                     ( 3, MAXN )
C
C                  Limb points are expressed in the reference frame
C                  designated by FIXREF. For each limb point, the
C                  orientation of the frame is evaluated at the epoch
C                  corresponding to the limb point; the epoch is
C                  provided in the output array EPOCHS (described
C                  below).
C
C                  Units of the limb points are km.
C
C
C     EPOCHS       is an array of epochs associated with the limb
C                  points, accounting for light time if aberration
C                  corrections are used. EPOCHS contains one element
C                  for each limb point. EPOCHS should be declared
C                  with length
C
C                     MAXN
C
C                  The element
C
C                     EPOCHS(I)
C
C                  is associated with the limb point
C
C                     POINTS(J,I), J = 1 to 3
C                  
C                  If CORLOC is set to 'CENTER', all values of EPOCHS
C                  will be the epoch associated with the target body
C                  center. That is, if aberration corrections are used,
C                  and if LT is the one-way light time from the target
C                  center to the observer, the elements of EPOCHS will
C                  all be set to
C
C                     ET - LT
C
C                  If CORLOC is set to 'ELLIPSOID LIMB', all values of
C                  EPOCHS for the limb points in a given half plane
C                  will be those for the reference ellipsoid limb point
C                  in that half plane. That is, if aberration
C                  corrections are used, and if LT(I) is the one-way
C                  light time to the observer from the reference
C                  ellipsoid limb point in the Ith half plane, the
C                  elements of EPOCHS for that half plane will all be
C                  set to
C
C                     ET - LT(I)
C
C                  When the target shape is given by DSK data, there
C                  normally will be a small difference in the light
C                  time between an actual limb point and that implied
C                  by the corresponding element of EPOCHS. See the
C                  description of TANGTS below.
C
C
C     TANGTS       is an array of tangent vectors connecting the
C                  observer to the limb points. The tangent vectors are
C                  expressed in the frame designated by FIXREF. For the
C                  Ith vector, the orientation of the frame is
C                  evaluated at the Ith epoch provided in the output
C                  array EPOCHS (described above).
C
C                  TANGTS should be declared with dimensions
C
C                     ( 3, MAXN )
C
C                  The elements
C
C                     TANGTS(J,I), J = 1 to 3
C
C                  are associated with the limb point
C
C                     POINTS(J,I), J = 1 to 3
C
C                  Units of the tangent vectors are km.
C
C                  When the target shape is given by DSK data, there
C                  normally will be a small difference in the light
C                  time between an actual limb point and that implied
C                  by the corresponding element of EPOCHS. This
C                  difference will affect the orientation of the target
C                  body-fixed frame and the output tangent vectors
C                  returned in the array TANGTS. All other factors
C                  being equal, the error in the tangent vector due to
C                  the light time error is proportional to the
C                  observer-target distance.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the specified aberration correction is unrecognized, the
C         error will be signaled by a routine in the call tree of this
C         routine. If transmission corrections are commanded, the error
C         SPICE(INVALIDOPTION) will be signaled.
C
C     2)  If either the target or observer input strings cannot be
C         converted to an integer ID code, the error
C         SPICE(IDCODENOTFOUND) is signaled.
C
C     3)  If OBSRVR and TARGET map to the same NAIF integer ID code,
C         the error SPICE(BODIESNOTDISTINCT) is signaled.
C
C     4)  If the input target body-fixed frame FIXREF is not
C         recognized, the error SPICE(NOFRAME) is signaled. A frame
C         name may fail to be recognized because a required frame
C         specification kernel has not been loaded; another cause is a
C         misspelling of the frame name.
C
C     5)  If the input frame FIXREF is not centered at the target body,
C         the error SPICE(INVALIDFRAME) is signaled.
C
C     6)  If the input argument METHOD is not recognized, the error
C         SPICE(INVALIDMETHOD) is signaled by this routine, or the
C         error is signaled by a routine in the call tree of this
C         routine.
C
C     7)  If METHOD contains an invalid limb type, the error 
C         SPICE(INVALIDLIMBTYPE) will be signaled.
C
C     8)  If the target and observer have distinct identities but are
C         at the same location the error SPICE(NOSEPARATION) is
C         signaled.
C
C     9)  If insufficient ephemeris data have been loaded prior to
C         calling LIMBPT, the error will be signaled by a routine in
C         the call tree of this routine. When light time correction is
C         used, sufficient ephemeris data must be available to
C         propagate the states of both observer and target to the solar
C         system barycenter.
C
C    10)  If the computation method requires an ellipsoidal target
C         shape and triaxial radii of the target body have not been
C         loaded into the kernel pool prior to calling LIMBPT, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine. 
C
C         If the radii are available in the kernel pool but the count
C         of radii values is not three, the error SPICE(BADRADIUSCOUNT)
C         will be signaled.
C
C         When the target shape is modeled by topographic data, radii
C         of the reference triaxial ellipsoid are still required if
C         the aberration correction locus is ELLIPSOID LIMB or if
C         the limb point generation method is GUIDED.
C
C    11)  The target must be an extended body. If the target body's
C         shape is modeled as an ellipsoid, and if any of the radii of
C         the target body are non-positive, the error will be diagnosed
C         and signaled by routines in the call tree of this routine.
C
C    12)  If PCK data specifying the target body-fixed frame
C         orientation have not been loaded prior to calling LIMBPT,
C         the error will be diagnosed and signaled by a routine in the
C         call tree of this routine.
C
C    13)  If METHOD specifies that the target surface is represented by
C         DSK data, and no DSK files are loaded for the specified
C         target, the error is signaled by a routine in the call tree
C         of this routine. 
C
C    14)  If the array bound MAXN is less than 1, the error
C         SPICE(INVALIDSIZE) will be signaled.
C
C    15)  If the number of cutting half-planes specified by NCUTS
C         is negative or greater than MAXN, the error
C         SPICE(INVALIDCOUNT) will be signaled.
C
C    16)  If the aberration correction locus is not recognized, the
C         error SPICE(INVALIDLOCUS) will be signaled.
C
C    17)  If the aberration correction locus is 'ELLIPSOID LIMB'
C         but limb type is not 'TANGENT', the error 
C         SPICE(BADLIMBLOCUSMIX) will be signaled.
C
C    18)  If the reference vector REFVEC is the zero vector, the 
C         error SPICE(ZEROVECTOR) will be signaled.
C
C    19)  If the reference vector REFVEC and the observer target
C         vector are linearly dependent, the error 
C         SPICE(DEGENERATECASE) will be signaled.
C
C    20)  If the limb computation uses the target ellipsoid limb 
C         plane, and the limb plane normal and reference vector
C         REFVEC are linearly dependent, the error 
C         SPICE(DEGENERATECASE) will be signaled.
C
C    21)  If the limb points cannot all be stored in the output POINTS
C         array, the error SPICE(OUTOFROOM) will be signaled.    
C
C    22)  If the surface is represented by DSK data, and if the search
C         step is non-positive, the error SPICE(INVALIDSEARCHSTEP) will
C         be signaled.
C
C    23)  If the surface is represented by DSK data, and if the search
C         tolerance is non-positive, the error SPICE(INVALIDTOLERANCE)
C         will be signaled.
C  
C    24)  If the roll step is non-positive and NCUTS is greater
C         than 1, the error SPICE(INVALIDROLLSTEP) will be signaled.
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
C        - Target body orientation data: these may be provided in a text
C          or binary PCK file. In some cases, target body orientation
C          may be provided by one more more CK files. In either case,
C          data are made available by loading the files via FURNSH.
C
C        - Shape data for the target body:
C
C            PCK data:
C
C               If the target body shape is modeled as an ellipsoid,
C               triaxial radii for the target body must be loaded into
C               the kernel pool. Typically this is done by loading a
C               text PCK file via FURNSH.
C
C               Triaxial radii are also needed if the target shape is
C               modeled by DSK data but one or both of the GUIDED limb
C               definition method or the ELLIPSOID LIMB aberration
C               correction locus are selected.
C
C            DSK data:
C
C               If the target shape is modeled by DSK data, DSK files
C               containing topographic data for the target body must be
C               loaded. If a surface list is specified, data for at
C               least one of the listed surfaces must be loaded.
C
C     The following data may be required:
C
C        - Frame data: if a frame definition is required to convert the
C          observer and target states to the body-fixed frame of the
C          target, that definition must be available in the kernel
C          pool. Typically the definition is supplied by loading a
C          frame kernel via FURNSH.
C
C        - Surface name-ID associations: if surface names are specified
C          in `method', the association of these names with their
C          corresponding surface ID codes must be established by
C          assignments of the kernel variables
C
C             NAIF_SURFACE_NAME
C             NAIF_SURFACE_CODE
C             NAIF_SURFACE_BODY
C
C          Normally these associations are made by loading a text
C          kernel containing the necessary assignments. An example
C          of such a set of assignments is
C
C             NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
C             NAIF_SURFACE_CODE += 1
C             NAIF_SURFACE_BODY += 499
C
C        - SCLK data: if the target body's orientation is provided by
C          CK files, an associated SCLK kernel must be loaded.
C
C
C     In all cases, kernel data are normally loaded once per program
C     run, NOT every time this routine is called.
C
C
C$ Particulars
C
C
C     Using DSK data
C     ==============
C
C        DSK loading and unloading
C        -------------------------
C
C        DSK files providing data used by this routine are loaded by
C        calling FURNSH and can be unloaded by calling UNLOAD or
C        KCLEAR. See the documentation of FURNSH for limits on numbers
C        of loaded DSK files.
C
C        For run-time efficiency, it's desirable to avoid frequent
C        loading and unloading of DSK files. When there is a reason to
C        use multiple versions of data for a given target body---for
C        example, if topographic data at varying resolutions are to be
C        used---the surface list can be used to select DSK data to be
C        used for a given computation. It is not necessary to unload
C        the data that are not to be used. This recommendation presumes
C        that DSKs containing different versions of surface data for a
C        given body have different surface ID codes.
C
C
C        DSK data priority
C        -----------------
C
C        A DSK coverage overlap occurs when two segments in loaded DSK
C        files cover part or all of the same domain---for example, a
C        given longitude-latitude rectangle---and when the time
C        intervals of the segments overlap as well.
C
C        When DSK data selection is prioritized, in case of a coverage
C        overlap, if the two competing segments are in different DSK
C        files, the segment in the DSK file loaded last takes
C        precedence. If the two segments are in the same file, the
C        segment located closer to the end of the file takes
C        precedence.
C
C        When DSK data selection is unprioritized, data from competing
C        segments are combined. For example, if two competing segments
C        both represent a surface as sets of triangular plates, the
C        union of those sets of plates is considered to represent the
C        surface. 
C
C        Currently only unprioritized data selection is supported.
C        Because prioritized data selection may be the default behavior
C        in a later version of the routine, the UNPRIORITIZED keyword is
C        required in the METHOD argument.
C
C        
C        Syntax of the METHOD input argument
C        -----------------------------------
C
C        The keywords and surface list in the METHOD argument
C        are called "clauses." The clauses may appear in any
C        order, for example
C
C           TANGENT/DSK/UNPRIORITIZED/<surface list>
C           DSK/TANGENT/<surface list>/UNPRIORITIZED
C           UNPRIORITIZED/<surface list>/DSK/TANGENT
C
C        The simplest form of the METHOD argument specifying use of
C        DSK data is one that lacks a surface list, for example:
C
C           'TANGENT/DSK/UNPRIORITIZED'
C           'GUIDED/DSK/UNPRIORITIZED'
C
C        For applications in which all loaded DSK data for the target
C        body are for a single surface, and there are no competing
C        segments, the above strings suffice. This is expected to be
C        the usual case.
C
C        When, for the specified target body, there are loaded DSK
C        files providing data for multiple surfaces for that body, the
C        surfaces to be used by this routine for a given call must be
C        specified in a surface list, unless data from all of the
C        surfaces are to be used together.
C
C        The surface list consists of the string
C
C           SURFACES =
C
C        followed by a comma-separated list of one or more surface
C        identifiers. The identifiers may be names or integer codes in
C        string format. For example, suppose we have the surface
C        names and corresponding ID codes shown below:
C
C           Surface Name                              ID code
C           ------------                              -------
C           'Mars MEGDR 128 PIXEL/DEG'                1
C           'Mars MEGDR 64 PIXEL/DEG'                 2
C           'Mars_MRO_HIRISE'                         3
C
C        If data for all of the above surfaces are loaded, then
C        data for surface 1 can be specified by either
C
C           'SURFACES = 1'
C
C        or
C
C           'SURFACES = "Mars MEGDR 128 PIXEL/DEG"'
C
C        Double quotes are used to delimit the surface name because
C        it contains blank characters. 
C           
C        To use data for surfaces 2 and 3 together, any
C        of the following surface lists could be used:
C
C           'SURFACES = 2, 3'
C
C           'SURFACES = "Mars MEGDR  64 PIXEL/DEG", 3'
C
C           'SURFACES = 2, Mars_MRO_HIRISE'
C
C           'SURFACES = "Mars MEGDR 64 PIXEL/DEG", Mars_MRO_HIRISE'
C         
C        An example of a METHOD argument that could be constructed
C        using one of the surface lists above is
C
C      'TANGENT/DSK/UNPRIORITIZED/SURFACES= "Mars MEGDR 64 PIXEL/DEG",3'
C
C     
C$ Examples
C
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C
C
C     1) Find apparent limb points on Phobos as seen from Mars. 
C
C        Due to Phobos' irregular shape, the TANGENT limb point
C        definition will used. It suffices to compute light time and
C        stellar aberration corrections for the center of Phobos, so
C        the CENTER aberration correction locus will be used. Use
C        converged Newtonian light time and stellar aberration
C        corrections in order to model the apparent position and 
C        orientation of Phobos.
C        
C        For comparison, compute limb points using both ellipsoid
C        and topographic shape models.
C
C        Use the target body-fixed +Z axis as the reference direction
C        for generating cutting half-planes. This choice enables the
C        user to see whether the first limb point is near the target's
C        north pole.
C
C        For each option, use just three cutting half-planes, in order
C        to keep the volume of output manageable. In most applications,
C        the number of cuts and the number of resulting limb points
C        would be much greater.
C
C        Use the meta-kernel below to load the required SPICE 
C        kernels. 
C
C
C           KPL/MK
C
C           File: limbpt_ex1.tm
C
C           This meta-kernel is intended to support operation of SPICE
C           example programs. The kernels shown here should not be
C           assumed to contain adequate or correct versions of data
C           required by SPICE-based user applications.
C
C           In order for an application to use this meta-kernel, the
C           kernels referenced here must be present in the user's
C           current working directory.
C
C           The names and contents of the kernels referenced
C           by this meta-kernel are as follows:
C
C              File name                        Contents
C              ---------                        --------
C              de430.bsp                        Planetary ephemeris
C              mar097.bsp                       Mars satellite ephemeris
C              pck00010.tpc                     Planet orientation and
C                                               radii
C              naif0011.tls                     Leapseconds
C              phobos512.bds                    DSK based on
C                                               Gaskell ICQ Q=512
C                                               Phobos plate model
C           \begindata
C
C              PATH_SYMBOLS    = 'GEN'
C              PATH_VALUES     = '/ftp/pub/naif/generic_kernels'
C
C              KERNELS_TO_LOAD = ( 'de430.bsp',
C                                  'mar097.bsp',
C                                  'pck00010.tpc',
C                                  'naif0011.tls',
C                                  '$GEN/dsk/phobos/phobos512.bds' )
C           \begintext
C
C
C
C     Example code begins here.
C
C
C     C
C     C     LIMBPT example 1
C     C
C     C        Find limb points on Phobos as seen from Mars.
C     C
C     C        Compute limb points using the tangent definition.
C     C        Perform aberration corrections for the target center.
C     C        Use both ellipsoid and DSK shape models.
C     C
C           PROGRAM EX1
C           IMPLICIT NONE
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      DPR
C           DOUBLE PRECISION      PI
C
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         META
C           PARAMETER           ( META   = 'limbpt_ex1.tm' )
C
C           CHARACTER*(*)         FM1
C           PARAMETER           ( FM1     =  '(A,F21.9)' )
C
C           CHARACTER*(*)         FM2
C           PARAMETER           ( FM2     =  '(1X,3F21.9)' )
C
C           INTEGER               BDNMLN
C           PARAMETER           ( BDNMLN = 36 )
C
C           INTEGER               FRNMLN
C           PARAMETER           ( FRNMLN = 32 )
C
C           INTEGER               CORLEN
C           PARAMETER           ( CORLEN = 20 )
C
C           INTEGER               MTHLEN
C           PARAMETER           ( MTHLEN = 50 )
C
C           INTEGER               NMETH
C           PARAMETER           ( NMETH  = 2 )
C
C           INTEGER               MAXN
C           PARAMETER           ( MAXN = 10000 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(CORLEN)    ABCORR
C           CHARACTER*(CORLEN)    CORLOC
C           CHARACTER*(FRNMLN)    FIXREF
C           CHARACTER*(MTHLEN)    METHOD ( NMETH )
C           CHARACTER*(BDNMLN)    OBSRVR
C           CHARACTER*(BDNMLN)    TARGET
C
C           DOUBLE PRECISION      DELROL
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      POINTS ( 3, MAXN )
C           DOUBLE PRECISION      ROLL
C           DOUBLE PRECISION      SCHSTP
C           DOUBLE PRECISION      SOLTOL
C           DOUBLE PRECISION      TANGTS ( 3, MAXN )
C           DOUBLE PRECISION      TRGEPS ( MAXN )
C           DOUBLE PRECISION      Z      ( 3 )
C
C           INTEGER               I
C           INTEGER               J
C           INTEGER               K
C           INTEGER               M
C           INTEGER               NCUTS
C           INTEGER               NPTS   ( MAXN )
C           INTEGER               START
C
C     C
C     C     Initial values
C     C
C           DATA                  METHOD /
C          .                        'TANGENT/ELLIPSOID',
C          .                        'TANGENT/DSK/UNPRIORITIZED'
C          .                             /
C           DATA                  Z      / 0.D0, 0.D0, 1.D0 /
C     C
C     C     Load kernel files via the meta-kernel.
C     C
C           CALL FURNSH ( META )
C     C
C     C     Set target, observer, and target body-fixed,
C     C     body-centered reference frame.
C     C
C           OBSRVR = 'MARS'
C           TARGET = 'PHOBOS'
C           FIXREF = 'IAU_PHOBOS'
C     C
C     C     Set aberration correction and correction locus.
C     C
C           ABCORR = 'CN+S'
C           CORLOC = 'CENTER'
C     C
C     C     Convert the UTC request time string seconds past
C     C     J2000, TDB.
C     C
C           CALL STR2ET ( '2008 AUG 11 00:00:00', ET )
C     C
C     C     Compute a set of limb points using light time and
C     C     stellar aberration corrections. Use both ellipsoid
C     C     and DSK shape models. Use a step size of 100
C     C     microradians to ensure we don't miss the limb.
C     C     Set the convergence tolerance to 100 nanoradians, 
C     C     which will limit the height error to about 1 meter. 
C     C     Compute 3 limb points for each computation method.
C     C
C           SCHSTP = 1.D-4
C           SOLTOL = 1.D-7
C           NCUTS  = 3
C
C           WRITE (*,*) ' '
C           WRITE (*,*) 'Observer:       '//OBSRVR
C           WRITE (*,*) 'Target:         '//TARGET
C           WRITE (*,*) 'Frame:          '//FIXREF
C           WRITE (*,*) ' '
C           WRITE (*,*) 'Number of cuts: ', NCUTS
C           WRITE (*,*) ' '
C
C           DELROL = 2*PI() / NCUTS
C
C           DO I = 1, NMETH
C
C              CALL LIMBPT ( METHOD(I), TARGET, ET,     FIXREF,
C          .                 ABCORR,    CORLOC, OBSRVR, Z,
C          .                 DELROL,    NCUTS,  SCHSTP, SOLTOL, 
C          .                 MAXN,      NPTS,   POINTS, TRGEPS, 
C          .                 TANGTS                            )
C     C
C     C        Write the results.
C     C
C              WRITE(*,*) ' '
C              WRITE(*,*) 'Computation method = ', METHOD(I)
C              WRITE(*,*) 'Locus              = ', CORLOC
C              WRITE(*,*) ' '
C
C
C              START  = 0
C
C              DO J = 1, NCUTS
C
C                 ROLL = (J-1) * DELROL
C
C                 WRITE(*,*)   ' '
C                 WRITE(*,FM1) '  Roll angle (deg) = ', ROLL * DPR()
C                 WRITE(*,FM1) '     Target epoch  = ', TRGEPS(J)
C                 WRITE(*,*)   '    Number of limb points at this '
C          .      //           'roll angle: ',
C          .                   NPTS(J)
C
C                 WRITE (*,*) '      Limb points'
C
C                 DO K = 1, NPTS(J)
C                    WRITE (*,FM2) ( POINTS(M,K+START), M = 1, 3 )
C                 END DO
C
C                 START = START + NPTS(J)
C
C              END DO
C
C              WRITE (*,*) ' '
C
C           END DO
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was:
C
C
C        Observer:       MARS
C        Target:         PHOBOS
C        Frame:          IAU_PHOBOS
C
C        Number of cuts:            3
C
C
C        Computation method = TANGENT/ELLIPSOID
C        Locus              = CENTER
C
C
C         Roll angle (deg) =           0.000000000
C            Target epoch  =   271684865.152078211
C            Number of limb points at this roll angle:            1
C              Limb points
C                  0.016445326         -0.000306114          9.099992715
C
C         Roll angle (deg) =         120.000000000
C            Target epoch  =   271684865.152078211
C            Number of limb points at this roll angle:            1
C              Limb points
C                 -0.204288375         -9.235230829         -5.333237706
C
C         Roll angle (deg) =         240.000000000
C            Target epoch  =   271684865.152078211
C            Number of limb points at this roll angle:            1
C              Limb points
C                  0.242785221          9.234520095         -5.333231253
C
C
C        Computation method = TANGENT/DSK/UNPRIORITIZED
C        Locus              = CENTER
C
C
C         Roll angle (deg) =           0.000000000
C            Target epoch  =   271684865.152078211
C            Number of limb points at this roll angle:            1
C              Limb points
C                 -0.398901673          0.007425178          9.973720555
C
C         Roll angle (deg) =         120.000000000
C            Target epoch  =   271684865.152078211
C            Number of limb points at this roll angle:            1
C              Limb points
C                 -0.959300281         -8.537573427         -4.938700447
C
C         Roll angle (deg) =         240.000000000
C            Target epoch  =   271684865.152078211
C            Number of limb points at this roll angle:            1
C              Limb points
C                 -1.380536729          9.714334047         -5.592916790
C
C
C
C     2) Find apparent limb points on Mars as seen from the earth.
C        Compare results using different computation options.
C
C        Use both the TANGENT and GUIDED limb point definitions. For
C        the tangent limb points, use the ELLIPSOID LIMB aberration
C        correction locus; for the guided limb points, use the CENTER
C        locus. For the GUIDED limb points, also compute the distance
C        of each point from the corresponding point computed using the
C        TANGENT definition.
C
C        For comparison, compute limb points using both ellipsoid and
C        topographic shape models.
C
C        Check the limb points by computing the apparent emission
C        angles at each limb point.
C
C        For the ellipsoid shape model, we expect emission angles very
C        close to 90 degrees, since each illumination angle calculation
C        is done using aberration corrections for the limb point at
C        which the angles are measured.
C
C        Use the target body-fixed +Z axis as the reference direction
C        for generating cutting half-planes. This choice enables the
C        user to see whether the first limb point is near the target's
C        north pole.
C        
C        For each option, use just three cutting half-planes, in order
C        to keep the volume of output manageable. In most applications,
C        the number of cuts and the number of resulting limb points
C        would be much greater.
C
C        Use the meta-kernel shown below.
C
C
C           KPL/MK
C
C           File: limbpt_ex2.tm
C
C           This meta-kernel is intended to support operation of SPICE
C           example programs. The kernels shown here should not be
C           assumed to contain adequate or correct versions of data
C           required by SPICE-based user applications.
C
C           In order for an application to use this meta-kernel, the
C           kernels referenced here must be present in the user's
C           current working directory.
C
C           The names and contents of the kernels referenced
C           by this meta-kernel are as follows:
C
C              File name                        Contents
C              ---------                        --------
C              de430.bsp                        Planetary ephemeris
C              mar097.bsp                       Mars satellite ephemeris
C              pck00010.tpc                     Planet orientation and
C                           radii
C              naif0011.tls                     Leapseconds
C              megr90n000cb_plate.bds           DSK plate model based on
C                                               MGS MOLAR MEGDR DEM, 
C                                               resolution 4 
C                                               pixels/degree.
C
C           \begindata
C
C              KERNELS_TO_LOAD = ( 'de430.bsp',
C                                  'mar097.bsp',
C                                  'pck00010.tpc',
C                                  'naif0011.tls',
C                                  'megr90n000cb_plate.bds' )
C           \begintext
C
C
C
C     Example code begins here. 
C
C
C     C
C     C     LIMBPT example 2
C     C
C     C        Find limb points on Mars as seen from the earth.
C     C
C     C        Compute limb points using both the tangent and
C     C        "guided" definitions.
C     C
C     C        For the tangent limb points, perform aberration
C     C        corrections for the reference ellipsoid limb.
C     C
C     C        Check limb points by computing emission angles at
C     C        each point.
C     C
C     C        Use both ellipsoid and DSK shape models.
C     C
C           PROGRAM EX2
C           IMPLICIT NONE
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      DPR
C           DOUBLE PRECISION      PI
C           DOUBLE PRECISION      VDIST
C           DOUBLE PRECISION      VNORM
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         META
C           PARAMETER           ( META    = 'limbpt_ex2.tm' )
C
C           CHARACTER*(*)         FM1
C           PARAMETER           ( FM1     =  '(A,F21.9)' )
C
C           INTEGER               BDNMLN
C           PARAMETER           ( BDNMLN = 36 )
C
C           INTEGER               FRNMLN
C           PARAMETER           ( FRNMLN = 32 )
C
C           INTEGER               CORLEN
C           PARAMETER           ( CORLEN = 20 )
C
C           INTEGER               MTHLEN
C           PARAMETER           ( MTHLEN = 50 )
C
C           INTEGER               NMETH
C           PARAMETER           ( NMETH  = 3 )
C
C           INTEGER               MAXN
C           PARAMETER           ( MAXN   = 100 )
C     C
C     C     Local variables
C     C
C           CHARACTER*(CORLEN)    ABCORR
C           CHARACTER*(CORLEN)    CORLOC ( NMETH )
C           CHARACTER*(FRNMLN)    FIXREF
C           CHARACTER*(MTHLEN)    ILUMTH ( NMETH )
C           CHARACTER*(BDNMLN)    OBSRVR
C           CHARACTER*(BDNMLN)    TARGET
C           CHARACTER*(MTHLEN)    METHOD ( NMETH )
C
C           DOUBLE PRECISION      ALT
C           DOUBLE PRECISION      DELROL
C           DOUBLE PRECISION      DIST
C           DOUBLE PRECISION      EMISSN
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      F
C           DOUBLE PRECISION      LAT
C           DOUBLE PRECISION      LON
C           DOUBLE PRECISION      LT
C           DOUBLE PRECISION      PHASE
C           DOUBLE PRECISION      POINTS ( 3, MAXN )
C           DOUBLE PRECISION      SVPNTS ( 3, MAXN )
C           DOUBLE PRECISION      POS    ( 3 )
C           DOUBLE PRECISION      RADII  ( 3 )
C           DOUBLE PRECISION      RE
C           DOUBLE PRECISION      ROLL
C           DOUBLE PRECISION      RP
C           DOUBLE PRECISION      SCHSTP
C           DOUBLE PRECISION      SOLAR
C           DOUBLE PRECISION      SOLTOL
C           DOUBLE PRECISION      SRFVEC ( 3 )
C           DOUBLE PRECISION      TANGTS ( 3, MAXN )
C           DOUBLE PRECISION      TRGEPC
C           DOUBLE PRECISION      TRGEPS ( MAXN )
C           DOUBLE PRECISION      Z      ( 3 )
C
C           INTEGER               I
C           INTEGER               J
C           INTEGER               K
C           INTEGER               M
C           INTEGER               N
C           INTEGER               NCUTS
C           INTEGER               NPTS   ( MAXN )
C           INTEGER               START
C
C     C
C     C     Initial values
C     C
C           DATA                  CORLOC /
C          .                        'ELLIPSOID LIMB',
C          .                        'ELLIPSOID LIMB',
C          .                        'CENTER'
C          .                             /
C
C           DATA                  ILUMTH /
C          .                        'ELLIPSOID',
C          .                        'DSK/UNPRIORITIZED',
C          .                        'DSK/UNPRIORITIZED'
C          .                             /
C
C           DATA                  METHOD /
C          .                        'TANGENT/ELLIPSOID',
C          .                        'TANGENT/DSK/UNPRIORITIZED',
C          .                        'GUIDED/DSK/UNPRIORITIZED'
C          .                             /
C
C           DATA                  Z      / 0.D0, 0.D0, 1.D0 /
C     C
C     C     Load kernel files via the meta-kernel.
C     C
C           CALL FURNSH ( META )
C     C
C     C     Set target, observer, and target body-fixed, body-centered
C     C     reference frame.
C     C
C           OBSRVR = 'EARTH'
C           TARGET = 'MARS'
C           FIXREF = 'IAU_MARS'
C     C
C     C     Set the aberration correction. We'll set the
C     C     correction locus below.
C     C
C           ABCORR = 'CN+S'
C     C
C     C     Convert the UTC request time string seconds past
C     C     J2000, TDB.
C     C
C           CALL STR2ET ( '2008 AUG 11 00:00:00', ET )
C     C
C     C     Look up the target body's radii. We'll use these to
C     C     convert Cartesian to planetographic coordinates. Use
C     C     the radii to compute the flattening coefficient of
C     C     the reference ellipsoid.
C     C
C           CALL BODVRD ( TARGET, 'RADII', 3, N, RADII )
C     C
C     C     Compute the flattening coefficient for planetodetic
C     C     coordinates
C     C
C           RE = RADII(1)
C           RP = RADII(3)
C           F  = ( RE - RP ) / RE
C     C
C     C     Compute a set of limb points using light time and
C     C     stellar aberration corrections. Use both ellipsoid
C     C     and DSK shape models.
C     C
C     C     Obtain the observer-target distance at ET.
C     C
C           CALL SPKPOS ( TARGET, ET,  'J2000', ABCORR,
C          .              OBSRVR, POS, LT              )
C           DIST = VNORM( POS )
C     C
C     C     Set the angular step size so that a single step will
C     C     be taken in the root bracketing process; that's all
C     C     that is needed since we don't expect to have multiple
C     C     limb points in any cutting half-plane.
C     C
C           SCHSTP = 4.D0
C     C
C     C     Set the convergence tolerance to minimize the height
C     C     error. We can't achieve the 1 millimeter precision
C     C     suggested by the formula because the earth-Mars
C     C     distance is about 3.5e8 km. Compute 3 limb points
C     C     for each computation method.
C     C
C           SOLTOL = 1.D-6/DIST
C     C
C     C     Set the number of cutting half-planes and roll step.
C     C
C           NCUTS  = 3
C           DELROL = 2*PI() / NCUTS
C
C           WRITE (*,*) ' '
C           WRITE (*,*) 'Observer:       '//OBSRVR
C           WRITE (*,*) 'Target:         '//TARGET
C           WRITE (*,*) 'Frame:          '//FIXREF
C           WRITE (*,*) ' '
C           WRITE (*,*) 'Number of cuts: ', NCUTS
C
C
C           DO I = 1, NMETH
C
C              CALL LIMBPT ( METHOD(I), TARGET,    ET,     FIXREF,
C          .                 ABCORR,    CORLOC(I), OBSRVR, Z,
C          .                 DELROL,    NCUTS,     SCHSTP, SOLTOL, 
C          .                 MAXN,      NPTS,      POINTS, TRGEPS, 
C          .                 TANGTS                                )
C     C
C     C        Write the results.
C     C
C              WRITE(*,*) ' '
C              WRITE(*,*) 'Computation method = ', METHOD(I)
C              WRITE(*,*) 'Locus              = ', CORLOC(I)
C
C
C              START  = 0
C
C              DO J = 1, NCUTS
C
C                 ROLL = (J-1) * DELROL
C
C                 WRITE(*,*)   ' '
C                 WRITE(*,FM1) '   Roll angle (deg) = ', ROLL * DPR()
C                 WRITE(*,FM1) '     Target epoch   = ', TRGEPS(J)
C                 WRITE(*,*)   '    Number of limb points at this '
C          .      //           'roll angle: ',
C          .                   NPTS(J)
C
C                 DO K = 1, NPTS(J)
C
C                    WRITE (*,*) '    Limb point planetodetic '
C          .         //          'coordinates:'
C
C                    CALL RECGEO ( POINTS(1,K+START), RE,  F,
C          .                       LON,               LAT, ALT )
C
C                    WRITE (*,FM1) '      Longitude      (deg): ',
C          .                       LON*DPR()
C                    WRITE (*,FM1) '      Latitude       (deg): ',
C          .                       LAT*DPR()
C                    WRITE (*,FM1) '      Altitude        (km): ',
C          .                       ALT
C
C     C
C     C              Get illumination angles for this limb point.
C     C
C                    M = K+START
C
C                    CALL ILUMIN ( ILUMTH,      TARGET, ET,
C          .                       FIXREF,      ABCORR, OBSRVR,
C          .                       POINTS(1,M), TRGEPC, SRFVEC,
C          .                       PHASE,       SOLAR,  EMISSN  )
C
C                    WRITE (*,FM1) '      Emission angle (deg): ',
C          .                     EMISSN * DPR()
C
C                    IF ( I .EQ. 2 ) THEN
C
C                       CALL VEQU ( POINTS(1,M), SVPNTS(1,M) )
C
C                    ELSE IF ( I .EQ. 3  ) THEN
C
C                       DIST = VDIST( POINTS(1,M), SVPNTS(1,M) )
C
C                       WRITE (*,FM1)
C          .            '      Distance error  (km): ', DIST
C                    END IF
C
C
C                 END DO
C
C                 START = START + NPTS(J)
C
C              END DO
C
C              WRITE (*,*) ' '
C
C           END DO
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was: 
C 
C
C        Observer:       EARTH
C        Target:         MARS
C        Frame:          IAU_MARS
C
C        Number of cuts:            3
C
C        Computation method = TANGENT/ELLIPSOID
C        Locus              = ELLIPSOID LIMB
C
C          Roll angle (deg) =           0.000000000
C            Target epoch   =   271683700.368869901
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):         -19.302258950
C             Latitude       (deg):          64.005620446
C             Altitude        (km):          -0.000000000
C             Emission angle (deg):          90.000000000
C
C          Roll angle (deg) =         120.000000000
C            Target epoch   =   271683700.368948162
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):          85.029135674
C             Latitude       (deg):         -26.912378799
C             Altitude        (km):           0.000000000
C             Emission angle (deg):          90.000000000
C
C          Roll angle (deg) =         240.000000000
C            Target epoch   =   271683700.368949771
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):        -123.633654215
C             Latitude       (deg):         -26.912378799
C             Altitude        (km):          -0.000000000
C             Emission angle (deg):          90.000000000
C
C
C        Computation method = TANGENT/DSK/UNPRIORITIZED
C        Locus              = ELLIPSOID LIMB
C
C          Roll angle (deg) =           0.000000000
C            Target epoch   =   271683700.368869901
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):         -19.302258950
C             Latitude       (deg):          63.893637269
C             Altitude        (km):          -3.667553936
C             Emission angle (deg):          90.112271887
C
C          Roll angle (deg) =         120.000000000
C            Target epoch   =   271683700.368948162
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):          85.434644179
C             Latitude       (deg):         -26.705411231
C             Altitude        (km):          -0.044832377
C             Emission angle (deg):          89.583080113
C
C          Roll angle (deg) =         240.000000000
C            Target epoch   =   271683700.368949771
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):        -123.375003954
C             Latitude       (deg):         -27.043096556
C             Altitude        (km):           3.695628339
C             Emission angle (deg):          90.265135303
C
C
C        Computation method = GUIDED/DSK/UNPRIORITIZED
C        Locus              = CENTER
C
C          Roll angle (deg) =           0.000000000
C            Target epoch   =   271683700.368922532
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):         -19.302259163
C             Latitude       (deg):          64.005910146
C             Altitude        (km):          -3.676424552
C             Emission angle (deg):          89.999998824
C             Distance error  (km):           6.664218206
C
C          Roll angle (deg) =         120.000000000
C            Target epoch   =   271683700.368922532
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):          85.029135793
C             Latitude       (deg):         -26.912405352
C             Altitude        (km):          -0.328988915
C             Emission angle (deg):          89.999999843
C             Distance error  (km):          24.686472808
C
C          Roll angle (deg) =         240.000000000
C            Target epoch   =   271683700.368922532
C            Number of limb points at this roll angle:            1
C            Limb point planetodetic coordinates:
C             Longitude      (deg):        -123.633653487
C             Latitude       (deg):         -26.912086524
C             Altitude        (km):           3.626058850
C             Emission angle (deg):          90.000001307
C             Distance error  (km):          15.716034625
C
C
C     3) Find apparent limb points on comet Churyumov-Gerasimenko
C        as seen from the Rosetta orbiter.
C
C        This computation is an example of a case for which some
C        of the cutting half-planes contain multiple limb points.
C
C        Use the TANGENT limb definition, since the target shape
C        is not well approximated by its reference ellipsoid.
C        Use the CENTER aberration correction locus since the 
C        light time difference across the object is small.
C
C        Use the meta-kernel shown below.
C       
C 
C          KPL/MK
C
C          File: limbpt_ex3.tm
C
C          This meta-kernel is intended to support operation of SPICE
C          example programs. The kernels shown here should not be
C          assumed to contain adequate or correct versions of data
C          required by SPICE-based user applications.
C
C          In order for an application to use this meta-kernel, the
C          paths of the kernels referenced here must be adjusted to 
C          be compatible with the user's host computer directory
C          structure.
C
C          The names and contents of the kernels referenced
C          by this meta-kernel are as follows:
C
C            File name                          Contents
C            ---------                          --------
C            DE405.BSP                          Planetary ephemeris
C            NAIF0011.TLS                       Leapseconds
C            ROS_CG_M004_NSPCESA_N_V1.BDS       DSK plate model based on
C                                               Rosetta NAVCAM data
C            RORB_DV_145_01_______00216.BSP     Rosetta orbiter
C                                               ephemeris
C            CORB_DV_145_01_______00216.BSP     Comet Churyumov-
C                                               Gerasimenko ephemeris
C            ROS_CG_RAD_V10.TPC                 Comet Churyumov-
C                                               Gerasimenko radii
C            ROS_V25.TF                         Comet C-G frame kernel
C                                               (includes SCLK 
C                                               parameters)
C            CATT_DV_145_01_______00216.BC      Comet C-G C-kernel
C
C
C          \begindata
C
C             PATH_VALUES     = (
C
C                '/ftp/pub/naif/pds/data/+'
C                'ro_rl-e_m_a_c-spice-6-v1.0/rossp_1000/DATA'
C
C                               )
C
C             PATH_SYMBOLS    = (
C
C                'KERNELS'
C                               )
C
C             KERNELS_TO_LOAD = (
C
C                '$KERNELS/SPK/DE405.BSP'
C                '$KERNELS/LSK/NAIF0011.TLS'
C                '$KERNELS/SPK/RORB_DV_145_01_______00216.BSP'
C                '$KERNELS/SPK/CORB_DV_145_01_______00216.BSP'
C                '$KERNELS/PCK/ROS_CG_RAD_V10.TPC'
C                '$KERNELS/FK/ROS_V25.TF'
C                '$KERNELS/CK/CATT_DV_145_01_______00216.BC'
C                '$KERNELS/DSK/ROS_CG_M004_NSPCESA_N_V1.BDS'
C
C                               )
C          \begintext
C
C
C     Example code begins here. 
C
C
C       C
C       C     LIMBPT example 3
C       C
C       C        Find limb points on comet Churyumov-Gerasimenko
C       C        as seen from the Rosetta orbiter.
C       C
C       C        Compute limb points using the tangent definition.
C       C        Perform aberration corrections for the target center.
C       C        Use both ellipsoid and DSK shape models.
C       C
C       C        Display only limb points lying in half-planes that
C       C        contain multiple limb points.
C       C
C             PROGRAM EX3
C             IMPLICIT NONE
C       C
C       C     SPICELIB functions
C       C
C             DOUBLE PRECISION      DPR
C             DOUBLE PRECISION      PI
C             DOUBLE PRECISION      RPD
C             DOUBLE PRECISION      VNORM
C       C
C       C     Local parameters
C       C
C             CHARACTER*(*)         META
C             PARAMETER           ( META   = 'limbpt_ex3.tm' )
C
C             CHARACTER*(*)         FM1
C             PARAMETER           ( FM1     =  '(A,F21.9)' )
C
C             CHARACTER*(*)         FM2
C             PARAMETER           ( FM2     =  '(1X,3F21.9)' )
C
C             INTEGER               BDNMLN
C             PARAMETER           ( BDNMLN = 36 )
C
C             INTEGER               FRNMLN
C             PARAMETER           ( FRNMLN = 32 )
C
C             INTEGER               CORLEN
C             PARAMETER           ( CORLEN = 20 )
C
C             INTEGER               MTHLEN
C             PARAMETER           ( MTHLEN = 50 )
C
C             INTEGER               MAXN
C             PARAMETER           ( MAXN = 1000 )
C       C
C       C     Local variables
C       C
C             CHARACTER*(CORLEN)    ABCORR
C             CHARACTER*(CORLEN)    CORLOC
C             CHARACTER*(FRNMLN)    FIXREF
C             CHARACTER*(MTHLEN)    METHOD
C             CHARACTER*(BDNMLN)    OBSRVR
C             CHARACTER*(BDNMLN)    TARGET
C
C             DOUBLE PRECISION      ANGLE
C             DOUBLE PRECISION      AXIS   ( 3 )
C             DOUBLE PRECISION      DELROL
C             DOUBLE PRECISION      ET
C             DOUBLE PRECISION      LT
C             DOUBLE PRECISION      POINTS ( 3, MAXN )
C             DOUBLE PRECISION      REFVEC ( 3 )
C             DOUBLE PRECISION      ROLL
C             DOUBLE PRECISION      SCHSTP
C             DOUBLE PRECISION      SOLTOL
C             DOUBLE PRECISION      TANGTS ( 3, MAXN )
C             DOUBLE PRECISION      TRGEPS ( MAXN )
C             DOUBLE PRECISION      TRGPOS ( 3 )
C             DOUBLE PRECISION      XVEC   ( 3 )
C
C             INTEGER               I
C             INTEGER               J
C             INTEGER               K
C             INTEGER               NCUTS
C             INTEGER               NPTS   ( MAXN )
C             INTEGER               START
C       C
C       C     Initial values
C       C
C             DATA                  METHOD /
C            .                        'TANGENT/DSK/UNPRIORITIZED'
C            .                             /
C             DATA                  XVEC   / 1.D0, 0.D0, 0.D0 /
C       C
C       C     Load kernel files via the meta-kernel.
C       C
C             CALL FURNSH ( META )
C       C
C       C     Set target, observer, and target body-fixed,
C       C     body-centered reference frame.
C       C
C             OBSRVR = 'ROSETTA'
C             TARGET = 'CHURYUMOV-GERASIMENKO'
C             FIXREF = '67P/C-G_CK'
C       C
C       C     Set aberration correction and correction locus.
C       C
C             ABCORR = 'CN+S'
C             CORLOC = 'CENTER'
C       C
C       C     Convert the UTC request time string seconds past
C       C     J2000, TDB.
C       C
C             CALL STR2ET ( '2015 MAY 10 00:00:00', ET )
C       C
C       C     Compute a set of limb points using light time and
C       C     stellar aberration corrections. Use a step size
C       C     corresponding to a 10 meter height error to ensure
C       C     we don't miss the limb. Set the convergence tolerance
C       C     to 1/100 of this amount, which will limit the height
C       C     convergence error to about 10 cm.
C       C
C             CALL SPKPOS ( TARGET, ET,     FIXREF, ABCORR,
C            .              OBSRVR, TRGPOS, LT             )
C
C
C             SCHSTP = 1.D-2  / VNORM(TRGPOS)
C             SOLTOL = SCHSTP / 100.D0
C
C       C
C       C     Set the reference vector to the start of a
C       C     region of the roll domain on which we know
C       C     (from an external computation) that we'll
C       C     find multiple limb points in some half planes.
C       C     Compute 30 limb points, starting with the
C       C     half-plane containing the reference vector.
C       C
C             CALL VMINUS ( TRGPOS, AXIS )
C
C             ANGLE = 310.0D0 * RPD()
C
C             CALL VROTV  ( XVEC, AXIS, ANGLE, REFVEC )
C
C             NCUTS  = 30
C             DELROL = 2*PI() / 1000
C
C             WRITE (*,*) ' '
C             WRITE (*,*) 'Observer:       '//OBSRVR
C             WRITE (*,*) 'Target:         '//TARGET
C             WRITE (*,*) 'Frame:          '//FIXREF
C             WRITE (*,*) ' '
C             WRITE (*,*) 'Number of cuts: ', NCUTS
C             WRITE (*,*) ' '
C
C             CALL LIMBPT ( METHOD, TARGET, ET,     FIXREF,
C            .              ABCORR, CORLOC, OBSRVR, REFVEC,
C            .              DELROL, NCUTS,  SCHSTP, SOLTOL,
C            .              MAXN,   NPTS,   POINTS, TRGEPS,
C            .              TANGTS                          )
C       C
C       C     Write the results.
C       C
C             WRITE(*,*) ' '
C             WRITE(*,*) 'Computation method = ', METHOD
C             WRITE(*,*) 'Locus              = ', CORLOC
C             WRITE(*,*) ' '
C
C             START  = 0
C
C             DO I = 1, NCUTS
C
C                ROLL = (I-1) * DELROL
C
C                IF ( NPTS(I) .GT. 1 ) THEN
C
C                   WRITE(*,*)   ' '
C                   WRITE(*,FM1) '  Roll angle (deg) = ', ROLL * DPR()
C                   WRITE(*,FM1) '     Target epoch  = ', TRGEPS(I)
C                   WRITE(*,*)   '    Number of limb points at this '
C            .      //           'roll angle: ',
C            .                   NPTS(I)
C
C                   WRITE (*,*) '      Limb points'
C
C                   DO J = 1, NPTS(I)
C                      WRITE (*,FM2) ( POINTS(K,J+START), K = 1, 3 )
C                   END DO
C
C                END IF
C
C                START = START + NPTS(I)
C
C             END DO
C             WRITE (*,*) ' '
C
C             END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was (only the first three and last three
C     limb points are shown here): 
C
C
C      Observer:       ROSETTA
C      Target:         CHURYUMOV-GERASIMENKO
C      Frame:          67P/C-G_CK
C
C      Number of cuts:           30
C
C
C      Computation method = TANGENT/DSK/UNPRIORITIZED
C      Locus              = CENTER
C
C
C       Roll angle (deg) =           0.000000000
C          Target epoch  =   484488067.184933782
C          Number of limb points at this roll angle:            3
C            Limb points
C                1.320416231         -0.347379011          1.445260615
C                0.970350318          0.201685071          0.961996205
C                0.436720618          0.048224590          0.442280714
C
C       Roll angle (deg) =           0.360000000
C          Target epoch  =   484488067.184933782
C          Number of limb points at this roll angle:            3
C            Limb points
C                1.330290293         -0.352340416          1.438802587
C                0.965481808          0.202131806          0.946190003
C                0.453917030          0.082062880          0.447624224
C
C       Roll angle (deg) =           0.720000000
C          Target epoch  =   484488067.184933782
C          Number of limb points at this roll angle:            3
C            Limb points
C                1.339037339         -0.357848188          1.431256926
C                0.962159098          0.192370269          0.934342086
C                0.459160821          0.082273840          0.447880429
C
C        ...
C
C
C       Roll angle (deg) =           9.720000000
C          Target epoch  =   484488067.184933782
C          Number of limb points at this roll angle:            3
C            Limb points
C                1.568112677         -0.674947784          1.254880628
C                0.709857306         -0.111495634          0.547778706
C                0.491633785         -0.142729847          0.386229224
C
C       Roll angle (deg) =          10.080000000
C          Target epoch  =   484488067.184933782
C          Number of limb points at this roll angle:            3
C            Limb points
C                1.585230837         -0.663993935          1.249957484
C                0.633077981         -0.300058272          0.502702168
C                0.254736344         -0.760250955          0.266785439
C
C       Roll angle (deg) =          10.440000000
C          Target epoch  =   484488067.184933782
C          Number of limb points at this roll angle:            3
C            Limb points
C                1.599387477         -0.661757808          1.243621216
C                0.633255406         -0.293319746          0.495438969
C                0.271959251         -0.761967204          0.274619198
C
C
C
C$ Restrictions
C
C     1)  The light time approximations made by this routine may be
C         unsuitable for some observation geometries. For example, when
C         computing the limb of Mars as seen from the Earth, the
C         tangent vectors returned by this routine may be in error by
C         several km due to the light time error.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0 08-MAR-2017 (NJB)
C
C        Based on original version 14-NOV-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find limb points on target body
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      CLIGHT
      DOUBLE PRECISION      TOUCHD
      DOUBLE PRECISION      VDOT
      DOUBLE PRECISION      VNORM

      INTEGER               CARDD 

      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local parameters
C
      CHARACTER*(*)         RNAME
      PARAMETER           ( RNAME  = 'LIMBPT' )

      CHARACTER*(*)         IREF
      PARAMETER           ( IREF   = 'J2000' )

C
C     Convergence limit:
C
      DOUBLE PRECISION      CNVLIM
      PARAMETER           ( CNVLIM = 1.D-18 )
C
C     This limit was chosen to achieve convergence on
C     platforms providing extended precision arithmetic.
C

C
C     Maximum number of light time iterations for any
C     aberration correction:
C
      INTEGER               MAXITR
      PARAMETER           ( MAXITR = 5 )
      
      INTEGER               SSB
      PARAMETER           ( SSB   = 0 )

C
C     Saved body name length.
C
      INTEGER               MAXL
      PARAMETER           ( MAXL   = 36 )

C
C     Saved frame name length.
C
      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )

      INTEGER               MAXIVL
      PARAMETER           ( MAXIVL = 1000 )

      INTEGER               MAXWIN
      PARAMETER           ( MAXWIN = 2 * MAXIVL )

C
C     SPICELIB ellipse upper bound:
C
      INTEGER               UBEL
      PARAMETER           ( UBEL   = 9 )

C
C     Local variables
C     
      CHARACTER*(CVTLEN)    LMBSTR
      CHARACTER*(ACLLEN)    NRMLOC
      CHARACTER*(CORLEN)    PRVCOR
      CHARACTER*(ACLLEN)    PRVLOC
      CHARACTER*(MTHLEN)    PRVMTH
      CHARACTER*(SHPLEN)    SHPSTR
      CHARACTER*(SUBLEN)    SUBTYP
      CHARACTER*(CVTLEN)    SVLSTR
      CHARACTER*(TMTLEN)    TRMSTR

      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      CENTER ( 3 )
      DOUBLE PRECISION      CORTRG ( 3 )
      DOUBLE PRECISION      CP     ( 3 )
      DOUBLE PRECISION      CUTNML ( 3 )
      DOUBLE PRECISION      EDIR   ( 3 )
      DOUBLE PRECISION      ENORML ( 3 )
      DOUBLE PRECISION      EPOCH
      DOUBLE PRECISION      EPOINT ( 3 )
      DOUBLE PRECISION      IPOINT ( 3 )
      DOUBLE PRECISION      ISRFVC ( 3 )
      DOUBLE PRECISION      LIMB   ( UBEL )
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      LTERR
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      PLNVEC ( 3 )
      DOUBLE PRECISION      PNTBUF ( 3, MAXWIN )
      DOUBLE PRECISION      POS    ( 3 )
      DOUBLE PRECISION      PRVLT
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      RAYVTX ( 3 )
      DOUBLE PRECISION      RESULT ( LBCELL : MAXWIN )
      DOUBLE PRECISION      ROLL
      DOUBLE PRECISION      SMAJOR ( 3 )
      DOUBLE PRECISION      SMINOR ( 3 )
      DOUBLE PRECISION      PTARG  ( 3 )
      DOUBLE PRECISION      STLOFF ( 3 )
      DOUBLE PRECISION      STLPOS ( 3 )
      DOUBLE PRECISION      STOBS  ( 6 )
      DOUBLE PRECISION      SSBLT
      DOUBLE PRECISION      SSBTRG ( 3 )
      DOUBLE PRECISION      TMPVEC ( 3 )
      DOUBLE PRECISION      TRGEPC
      DOUBLE PRECISION      XFORM  ( 3, 3 )

      INTEGER               FXCENT
      INTEGER               FXCLSS
      INTEGER               FXFCDE
      INTEGER               FXTYID
      INTEGER               I
      INTEGER               J
      INTEGER               LMBTYP
      INTEGER               LOCCDE
      INTEGER               NSURF
      INTEGER               NUMITR
      INTEGER               OBSCDE
      INTEGER               PRVTRG
      INTEGER               ROOM
      INTEGER               SHAPE
      INTEGER               SRFLST ( MAXSRF )
      INTEGER               SVNRAD
      INTEGER               TO
      INTEGER               TOTAL
      INTEGER               TRGCDE
      
      LOGICAL               ATTBLK ( NABCOR )
      LOGICAL               FIRST
      LOGICAL               FND
      LOGICAL               PRI
      LOGICAL               SURFUP
      LOGICAL               USECN
      LOGICAL               USELT
      LOGICAL               USESTL

C
C     Saved name/ID item declarations.
C
      INTEGER               SVCTR1 ( CTRSIZ )
      CHARACTER*(MAXL)      SVTARG
      INTEGER               SVTCDE
      LOGICAL               SVFND1

      INTEGER               SVCTR2 ( CTRSIZ )
      CHARACTER*(MAXL)      SVOBSR
      INTEGER               SVOBSC
      LOGICAL               SVFND2

C
C     Saved frame name/ID item declarations.
C
      INTEGER               SVCTR3 ( CTRSIZ )
      CHARACTER*(FRNMLN)    SVFREF
      INTEGER               SVFXFC


C
C     Saved surface name/ID item declarations.
C
      INTEGER               SVCTR4 ( CTRSIZ )

C
C     Saved target radii declarations.
C
      INTEGER               SVCTR5 ( CTRSIZ )
      DOUBLE PRECISION      SVRADI ( 3 )

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  LMBTYP
      SAVE                  LOCCDE
      SAVE                  NSURF
      SAVE                  PNTBUF
      SAVE                  PRI
      SAVE                  PRVCOR
      SAVE                  PRVLOC
      SAVE                  PRVMTH
      SAVE                  PRVTRG
      SAVE                  SHAPE
      SAVE                  SRFLST
      SAVE                  SUBTYP
      SAVE                  SVNRAD
      SAVE                  USECN
      SAVE                  USELT
      SAVE                  USESTL

C
C     Saved name/ID items.
C
      SAVE                  SVCTR1
      SAVE                  SVTARG
      SAVE                  SVTCDE
      SAVE                  SVFND1

      SAVE                  SVCTR2
      SAVE                  SVOBSR
      SAVE                  SVOBSC
      SAVE                  SVFND2

C
C     Saved frame name/ID items.
C
      SAVE                  SVCTR3
      SAVE                  SVFREF
      SAVE                  SVFXFC


C
C     Saved surface name/ID items.
C
      SAVE                  SVCTR4
 
C
C     Saved reference ellipsoid items.
C
      SAVE                  SVCTR5
      SAVE                  SVRADI

C
C     Initial values
C
      DATA                  FIRST  / .TRUE.  /
      DATA                  PRVCOR / ' '     /
      DATA                  PRVLOC / ' '     /
      DATA                  PRVMTH / ' '     /
      DATA                  PRVTRG / 0       /
      DATA                  SVNRAD / 0       /
      DATA                  USECN  / .FALSE. / 
      DATA                  USELT  / .FALSE. / 
      DATA                  USESTL / .FALSE. / 


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( RNAME )

C
C     Counter initialization is done separately.
C
      IF ( FIRST ) THEN
C
C        Initialize counters.
C
         CALL ZZCTRUIN( SVCTR1 )
         CALL ZZCTRUIN( SVCTR2 )
         CALL ZZCTRUIN( SVCTR3 )
         CALL ZZCTRUIN( SVCTR4 )
         CALL ZZCTRUIN( SVCTR5 )

      END IF


      IF (  FIRST  .OR.  ( ABCORR .NE. PRVCOR )  ) THEN
C
C        Make sure the results of this block won't be reused
C        if we bail out due to an error.
C
         PRVCOR = ' '

C
C        The aberration correction flag differs from the value it
C        had on the previous call, if any. Analyze the new flag.
C
         CALL ZZVALCOR ( ABCORR, ATTBLK )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        In this routine, we don't allow transmission corrections.
C
         IF ( ATTBLK(XMTIDX) ) THEN

            CALL SETMSG ( 'Aberration correction # calls for '
     .      //            'transmission-style corrections. These '
     .      //            'are not supported for limb finding.'   )
            CALL ERRCH  ( '#',  ABCORR                            )
            CALL SIGERR ( 'SPICE(INVALIDOPTION)'                  )
            CALL CHKOUT ( RNAME                                   )
            RETURN

         END IF

C
C        Set logical flags indicating the attributes of the requested
C        correction:
C
C           USELT is .TRUE. when any type of light time correction
C           (normal or converged Newtonian) is specified.
C
C           USECN indicates converged Newtonian light time correction.
C
C           USESTL indicates stellar aberration corrections.
C
C
C        The above definitions are consistent with those used by
C        ZZPRSCOR.
C 
         USELT   =  ATTBLK ( LTIDX  )
         USECN   =  ATTBLK ( CNVIDX )
         USESTL  =  ATTBLK ( STLIDX )

C
C        The aberration correction flag is valid; save it.
C
         PRVCOR = ABCORR

      END IF

C
C     Obtain integer codes for the target and observer.
C 
      CALL ZZBODS2C ( SVCTR1, SVTARG, SVTCDE, SVFND1,
     .                TARGET, TRGCDE, FND    )
            
      IF ( .NOT. FND ) THEN
      
         CALL SETMSG ( 'The target, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit, or that you '
     .   //            'failed to load a kernel containing a '
     .   //            'name-ID mapping for this body.'           )
         CALL ERRCH  ( '#', TARGET                                )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( RNAME                                      )
         RETURN
      
      END IF
      
      
      CALL ZZBODS2C ( SVCTR2, SVOBSR, SVOBSC, SVFND2,
     .                OBSRVR, OBSCDE, FND    )
      
      IF ( .NOT. FND ) THEN

         CALL SETMSG ( 'The observer, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit, or that you '
     .   //            'failed to load a kernel containing a '
     .   //            'name-ID mapping for this body.'           )
         CALL ERRCH  ( '#', OBSRVR                                )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( RNAME                                      )
         RETURN
      
      END IF
      
 
C
C     Check the input body codes. If they are equal, signal
C     an error.
C
      IF ( OBSCDE .EQ. TRGCDE ) THEN
 
         CALL SETMSG ( 'In computing the surface intercept point, ' //
     .                 'the observing body and target body are the '//
     .                 'same. Both are #.'                          )
         CALL ERRCH  ( '#',  OBSRVR                                 )
         CALL SIGERR ( 'SPICE(BODIESNOTDISTINCT)'                   )
         CALL CHKOUT ( RNAME                                        )
         RETURN
 
      END IF

C
C     Determine the attributes of the frame designated by FIXREF.
C
      CALL ZZNAMFRM ( SVCTR3, SVFREF, SVFXFC, FIXREF, FXFCDE )

      CALL FRINFO ( FXFCDE, FXCENT, FXCLSS, FXTYID, FND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      IF ( .NOT. FND ) THEN

         CALL SETMSG ( 'Reference frame # is not recognized by ' //
     .                 'the SPICE frame subsystem. Possibly '    //
     .                 'a required frame definition kernel has ' //
     .                 'not been loaded.'                        )
         CALL ERRCH  ( '#',  FIXREF                              )
         CALL SIGERR ( 'SPICE(NOFRAME)'                          )
         CALL CHKOUT ( RNAME                                     )
         RETURN

      END IF

C
C     Make sure that FIXREF is centered at the target body's center.
C
      IF ( FXCENT .NE. TRGCDE ) THEN

         CALL SETMSG ( 'Reference frame # is not centered at the ' 
     .   //            'target body #. The ID code of the frame '
     .   //            'center is #.'                             )
         CALL ERRCH  ( '#',  FIXREF                               )
         CALL ERRCH  ( '#',  TARGET                               )
         CALL ERRINT ( '#',  FXCENT                               )
         CALL SIGERR ( 'SPICE(INVALIDFRAME)'                      )
         CALL CHKOUT ( RNAME                                      )
         RETURN

      END IF

C
C     Check whether the surface name/ID mapping has been updated.
C
      CALL ZZSRFTRK ( SVCTR4, SURFUP )

C
C     Initialize the SINCPT utility package for the next computation.
C     The choice of initialization routine depends on the target
C     surface type.
C     
      IF ( FIRST  .OR.  SURFUP  .OR.  ( METHOD .NE. PRVMTH ) ) THEN
C
C        Set the previous method string to an invalid value, so it
C        cannot match any future, valid input. This will force this
C        routine to parse the input method on the next call if any
C        failure occurs in this branch. Once success is assured, we can
C        record the current method in the previous method string.
C
         PRVMTH = ' '

C
C        Parse the method string. If the string is valid, the
C        outputs SHAPE and SUBTYP will always be be set. However,
C        SUBTYP is not used in this routine. 
C
C        For DSK shapes, the surface list array and count will be set
C        if the method string contains a surface list.
C
         CALL ZZPRSMET ( TRGCDE, METHOD, MAXSRF, SHPSTR, SUBTYP,
     .                   PRI,    NSURF,  SRFLST, LMBSTR, TRMSTR )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF


         IF (  EQSTR( SHPSTR, 'ELLIPSOID' )  ) THEN

            SHAPE = ELLSHP

         ELSE IF (  EQSTR( SHPSTR, 'DSK' )  ) THEN

            SHAPE = DSKSHP

         ELSE
C
C           This is a backstop check.
C
            CALL SETMSG ( '[1] Returned shape value from method '
     .      //            'string was <#>.'                      )
            CALL ERRCH  ( '#', SHPSTR                            )
            CALL SIGERR ( 'SPICE(BUG)'                           )
            CALL CHKOUT ( RNAME                                  )
            RETURN

         END IF


         IF (  EQSTR( LMBSTR, 'TANGENT' )  ) THEN

            LMBTYP = TANGNT

         ELSE IF (  EQSTR( LMBSTR, 'GUIDED' )  ) THEN

            LMBTYP = GUIDED

         ELSE

            CALL SETMSG ( 'Returned limb type from method '
     .      //            'string was <#>. Value must be '
     .      //            'TANGENT or GUIDED.'              )
            CALL ERRCH  ( '#', LMBSTR                       )
            CALL SIGERR ( 'SPICE(INVALIDLIMBTYPE)'          )
            CALL CHKOUT ( RNAME                             )
            RETURN

         END IF

C
C        There should be no subtype specification in the method 
C        string.
C
         IF ( SUBTYP .NE. ' ' ) THEN
            
            CALL SETMSG ( 'Spurious sub-observer point type <#> '
     .      //            'was present in the method string #. ' 
     .      //            'The sub-observer type is valid in '
     .      //            'the method strings for SUBPNT and '
     .      //            'SUBSLR, but is not applicable for LIMBPT.' )
            CALL ERRCH  ( '#', SUBTYP                                 )
            CALL ERRCH  ( '#', METHOD                                 )
            CALL SIGERR ( 'SPICE(INVALIDMETHOD)'                      )
            CALL CHKOUT ( RNAME                                       )
            RETURN

         END IF

C
C        There should be no terminator specification in the method 
C        string.
C
         IF ( TRMSTR .NE. ' ' ) THEN
            
            CALL SETMSG ( 'Spurious terminator shadow type <#> '
     .      //            'was present in the method string #. ' 
     .      //            'The terminator shadow type is valid in '
     .      //            'the method string for TERMPT, but is not '
     .      //            'applicable for LIMBPT.'                    )
            CALL ERRCH  ( '#', TRMSTR                                 )
            CALL ERRCH  ( '#', METHOD                                 )
            CALL SIGERR ( 'SPICE(INVALIDMETHOD)'                      )
            CALL CHKOUT ( RNAME                                       )
            RETURN

         END IF

         PRVMTH = METHOD

      END IF


C
C     Identify the aberration correction locus.
C
      IF (  FIRST  .OR.  ( CORLOC .NE. PRVLOC )  ) THEN

         CALL LJUCRS ( 1, CORLOC, NRMLOC )

         IF ( NRMLOC .EQ. 'CENTER' ) THEN

            LOCCDE = CTRCOR

         ELSE IF ( NRMLOC .EQ. 'ELLIPSOID LIMB' ) THEN

            LOCCDE = ELLCOR

         ELSE
            
            CALL SETMSG ( 'Aberration correction locus <#> '
     .      //            'was not recognized.'              ) 
            CALL ERRCH  ( '#', CORLOC                        )
            CALL SIGERR ( 'SPICE(INVALIDLOCUS)'              )
            CALL CHKOUT ( RNAME                              )
            RETURN

         END IF
C
C        At this point we have a valid locus. LOCCDE is set.
C        Save the input locus string so we can check for 
C        a change on the next call.
C
         PRVLOC = CORLOC

      END IF

C
C     Check the reference vector.
C     
      IF ( VZERO(REFVEC) ) THEN

         CALL SETMSG ( 'The reference vector was the zero vector.' )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                         )
         CALL CHKOUT ( RNAME                                       )
         RETURN
      
      END IF

C
C     At this point, the first pass actions were successful.
C
      FIRST  = .FALSE.

C
C     Check MAXN.
C
      IF ( MAXN .LT. 1 ) THEN
         
         CALL SETMSG ( 'MAXN = #; MAXN is required '
     .   //            'to be at least 1.'           )
         CALL ERRINT ( '#',  MAXN                    )
         CALL SIGERR ( 'SPICE(INVALIDSIZE)'          )
         CALL CHKOUT ( RNAME                         )
         RETURN

      END IF

C
C     Check NCUTS; there must be room for at least one limb point
C     for each cut. NCUTS may not be negative.
C
      IF (  ( NCUTS .LT. 1 ) .OR. ( NCUTS .GT. MAXN )  ) THEN
         
         CALL SETMSG ( 'NCUTS = #; MAXN = #; NCUTS is required '
     .   //            'to be non-negative and no larger than '
     .   //            'MAXN.'                                 )
         CALL ERRINT ( '#',  NCUTS                             )
         CALL ERRINT ( '#',  MAXN                              )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                   )
         CALL CHKOUT ( RNAME                                   )
         RETURN

      END IF

C
C     Check the angular search step size and convergence
C     tolerance. These checks apply only to DSK shapes.
C
      IF ( SHAPE .EQ. DSKSHP ) THEN

         IF ( SCHSTP .LE. 0.D0 ) THEN

            CALL SETMSG ( 'The angular search step '
     .      //            'SCHSTP = #; SCHSTP is required '
     .      //            'to be positive.'                )
            CALL ERRDP  ( '#',  SCHSTP                     )
            CALL SIGERR ( 'SPICE(INVALIDSEARCHSTEP)'       )
            CALL CHKOUT ( RNAME                            )
            RETURN

         END IF

         IF ( SOLTOL .LE. 0.D0 ) THEN

            CALL SETMSG ( 'The angular search tolerance '
     .      //            'SOLTOL = #; SOLTOL is required '
     .      //            'to be positive.'                )
            CALL ERRDP  ( '#',  SCHSTP                     )
            CALL SIGERR ( 'SPICE(INVALIDTOLERANCE)'        )
            CALL CHKOUT ( RNAME                            )
            RETURN

         END IF

      END IF

C
C     Check the roll step. This value applies only if
C     there are multiple cutting half-planes.
C
      IF (  ( NCUTS .GT. 1 ) .AND. ( ROLSTP .EQ. 0.D0 )  ) THEN

         CALL SETMSG ( 'The angular roll step is 0.D0. '
     .   //            'NCUTS = #. ROLSTP is required '
     .   //            'to be non-zero when NCUTS is ' 
     .   //            'greater than 1.'                )
         CALL ERRINT ( '#',  NCUTS                      )
         CALL SIGERR ( 'SPICE(INVALIDROLLSTEP)'         )
         CALL CHKOUT ( RNAME                            )
         RETURN

      END IF


      IF ( SHAPE .EQ. DSKSHP ) THEN
C
C        This is the DSK case.
C
C        Initialize the intercept algorithm to use a DSK
C        model for the surface of the target body. 
C        
         CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FXFCDE )

      ELSE IF ( SHAPE .NE. ELLSHP ) THEN

         CALL SETMSG ( 'Computation method argument was <#>; this ' 
     .   //            'string must specify a supported shape '
     .   //            'model and computation type. See the '
     .   //            'description of METHOD in the header of '
     .   //            'SUBPNT for details.'                      )
         CALL ERRCH  ( '#',  METHOD                               )
         CALL SIGERR ( 'SPICE(INVALIDMETHOD)'                     )
         CALL CHKOUT ( RNAME                                      )
         RETURN

      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF


C
C     Check MAXN.
C
      IF ( MAXN .LT. 1 ) THEN
         
         CALL SETMSG ( 'MAXN = #; MAXN is required '
     .   //            'to be at least 1.'           )
         CALL ERRINT ( '#',  MAXN                    )
         CALL SIGERR ( 'SPICE(INVALIDSIZE)'          )
         CALL CHKOUT ( RNAME                         )
         RETURN

      END IF

C
C     Check NCUTS; there must be room for at least one limb point
C     for each cut. NCUTS may not be negative.
C
      IF (  ( NCUTS .LT. 1 ) .OR. ( NCUTS .GT. MAXN )  ) THEN
         
         CALL SETMSG ( 'NCUTS = #; MAXN = #; NCUTS is required '
     .   //            'to be non-negative and no larger than '
     .   //            'MAXN.'                                 )
         CALL ERRINT ( '#',  NCUTS                             )
         CALL ERRINT ( '#',  MAXN                              )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                   )
         CALL CHKOUT ( RNAME                                   )
         RETURN

      END IF

C
C     Get target body radii if necessary.
C
      IF (      ( SHAPE  .EQ. ELLSHP ) 
     .     .OR. ( LOCCDE .EQ. ELLCOR ) 
     .     .OR. ( LMBTYP .EQ. GUIDED )  ) THEN

         IF ( TRGCDE .NE. PRVTRG ) THEN
C
C           Reset counter to force lookup.
C           
            CALL ZZCTRUIN ( SVCTR5 )

         END IF
C
C        Look up target radii using counter.
C
         CALL ZZBODVCD ( TRGCDE, 'RADII', 3, SVCTR5, SVNRAD, SVRADI )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         IF ( SVNRAD .NE. 3 ) THEN

            CALL SETMSG ( 'Number of target radii must be 3 '
     .      //            'but was #.'                       )
            CALL ERRINT ( '#',  SVNRAD                       )
            CALL SIGERR ( 'SPICE(BADRADIUSCOUNT)'            )
            CALL CHKOUT ( RNAME                              )
            RETURN

         END IF

         PRVTRG = TRGCDE

      END IF


C
C     Set up activities are complete at this point.
C 
C
C     Find limb points on the target.
C      
      CALL CLEARI ( NCUTS,  NPTS   )
      CALL SSIZED ( MAXWIN, RESULT )

C
C     Get initial observer-target vector, expressed in the target
C     body-fixed frame, evaluated at the target epoch. This vector
C     will be used for all option combinations.
C
      CALL SPKPOS ( TARGET, ET, FIXREF, ABCORR, OBSRVR, POS, LT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      IF ( VZERO(POS) ) THEN

         CALL SETMSG ( 'The distance between the observer and '
     .   //            'target at ET # is zero.'               )
         CALL ERRDP  ( '#',  ET                                )
         CALL SIGERR ( 'SPICE(NOSEPARATION)'                   )
         CALL CHKOUT ( RNAME                                   )
         RETURN

      END IF

C
C     The limb-finding technique depends on the aberration correction
C     locus. Start with the 'CENTER' version, since this is the
C     simpler case.
C

      IF ( LOCCDE .EQ. CTRCOR ) THEN
C
C        Aberration corrections are those applicable at the target
C        center.
C
C        Compute the epoch associated with the target center.
C
         CALL ZZCOREPC ( ABCORR, ET, LT, TRGEPC )

C
C        Compute the central axis, which is also the common ray vertex.
C        The axis points from the target to the observer.
C
         CALL VMINUS ( POS, AXIS )

C
C        Make sure the reference vector and axis are linearly 
C        independent.
C
         CALL VCRSS ( AXIS, REFVEC, CP )

         IF ( VZERO(CP) ) THEN
            
            CALL SETMSG ( 'Input reference vector and observer-target '
     .      //            'vector are linearly dependent.'             )
            CALL SIGERR ( 'SPICE(DEGENERATECASE)'                      )
            CALL CHKOUT ( RNAME                                        )
            RETURN

         END IF

C
C        If we're using an ellipsoidal shape model, or if
C        we're using the "guided" limb option. find the 
C        limb parameters of the reference ellipsoid.
C
         IF (  ( SHAPE .EQ. ELLSHP ) .OR. ( LMBTYP .EQ. GUIDED ) ) THEN

            CALL EDLIMB ( SVRADI(1), SVRADI(2), SVRADI(3), 
     .                    AXIS,      LIMB                 )

            CALL EL2CGV ( LIMB, CENTER, SMAJOR, SMINOR )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

            CALL UCRSS ( SMAJOR, SMINOR, ENORML )

C
C           Make sure ENORML points into the same half-space as
C           AXIS.
C
            IF ( VDOT( ENORML, AXIS ) .LT. 0.D0 ) THEN

               CALL VSCLIP ( -1.D0, ENORML )

            END IF


            IF ( SHAPE .EQ. DSKSHP ) THEN
C
C              Caution: this requires that ZZSUDSKI has been
c              called first.
C           
               CALL ZZMAXRAD ( MAXRAD )

            END IF

         END IF


         TO    = 1
         ROOM  = MAXN
         TOTAL = 0

C
C        Loop over the half planes, collecting limb points for
C        each one.
C
         DO I = 1, NCUTS

            ROLL = ( I - 1 ) * ROLSTP
C
C           Rotation of the half-planes is in the positive
C           sense about AXIS.
C
            CALL VROTV ( REFVEC, AXIS, ROLL, PLNVEC )

C
C           Let CUTNML be a vector normal to the current cutting
C           half-plane. We'll use this vector later.
C
            CALL UCRSS ( AXIS, PLNVEC, CUTNML )


            IF ( SHAPE .EQ. DSKSHP ) THEN
C
C              This is the DSK case.
C
               IF ( LMBTYP .EQ. TANGNT ) THEN
C
C                 This type of solution finds actual tangent rays on
C                 the target.
C
C                 Find the limb points that lie in the current
C                 half-plane.
C
C                 Note that RESULT is a cell, not a window.
C
                  CALL SCARDD ( 0, RESULT )
C
C                 Note that the evaluation epoch for the surface is
C                 optionally corrected for light time.
C
                  CALL ZZTANGNT ( LMBCRV, 0.D0,   SHAPE,  TRGCDE, 
     .                            NSURF,  SRFLST, FXFCDE, TRGEPC,
     .                            PLNVEC, AXIS,   SCHSTP, SOLTOL,
     .                            RESULT, PNTBUF                 )
            
                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  NPTS(I) = CARDD( RESULT )


               ELSE IF ( LMBTYP .EQ. GUIDED ) THEN
C
C                 This option uses the target's reference ellipsoid for
C                 guidance. For DSK shapes, the limb points are
C                 generated by finding surface intercepts of rays
C                 emanating from the center of the limb on the
C                 reference ellipsoid.
C
C                 The limb point we seek must lie in both the limb 
C                 plane and the cutting half-plane. Let EDIR the 
C                 unit direction vector satisfying these constraints.
C                 
                  CALL UCRSS ( CUTNML, ENORML, EDIR )

                  IF ( VZERO(EDIR) ) THEN

                     CALL SETMSG ( 'Vector defining cutting half plane '
     .               //            'and ellipsoid limb normal vector '
     .               //            'are linearly dependent.'           )
                     CALL SIGERR ( 'SPICE(DEGENERATECASE)'             )
                     CALL CHKOUT ( RNAME                               )
                     RETURN

                  END IF

C
C                 Find the intercept on the target surface of the ray
C                 emanating from CENTER in the direction EDIR. We must
C                 use a ray pointed in the opposite direction to
C                 perform this computation, since the surface may be
C                 invisible from the interior of the target.
C
                  CALL VLCOM ( 1.D0, CENTER, 3.D0*MAXRAD, EDIR, RAYVTX )
                  CALL VMINUS( EDIR, RAYDIR )

                  CALL ZZRAYSFX ( RAYVTX,      RAYDIR, TRGEPC, 
     .                            PNTBUF(1,1), FND             )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  IF ( FND ) THEN

                     NPTS(I) = 1
                  ELSE
                     NPTS(I) = 0
                  END IF
            

               ELSE
C
C                 This is a backstop case; it should never be reached.
C
                  CALL SETMSG ( 'Invalid limb type code: #' )
                  CALL ERRINT ( '#', LMBTYP                 )
                  CALL SIGERR ( 'SPICE(BUG)'                )
                  CALL CHKOUT ( RNAME                       )
                  RETURN

               END IF


            ELSE IF ( SHAPE .EQ. ELLSHP ) THEN               
C
C              This is the ellipsoid case.              
C
C              The limb point we seek must lie in both the limb plane
C              and the cutting half-plane. Let EDIR be the unit
C              direction vector satisfying these constraints.
C 
               CALL UCRSS ( CUTNML, ENORML, EDIR )

               IF ( VZERO(EDIR) ) THEN

                  CALL SETMSG ( 'Vector defining cutting half plane '
     .            //            'and ellipsoid limb normal vector '
     .            //            'are linearly dependent.'           )
                  CALL SIGERR ( 'SPICE(DEGENERATECASE)'             )
                  CALL CHKOUT ( RNAME                               )
                  RETURN

               END IF

C
C              Find the intercept on the target surface of the
C              the ray emanating from CENTER in the direction EDIR. 
C
               CALL SURFPT ( CENTER,    EDIR,      SVRADI(1),
     .                       SVRADI(2), SVRADI(3), PNTBUF(1,1), FND )


               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               IF ( .NOT. FND ) THEN

                  CALL SETMSG ( 'Limb point not found on reference '
     .            //            'ellipsoid for cutting half plane '
     .            //            'at index #. The point should always '
     .            //            'be found.'                           )
                  CALL ERRINT ( '#', I                                )
                  CALL SIGERR ( 'SPICE(BUG)'                          )
                  CALL CHKOUT ( RNAME                                 )
                  RETURN

               END IF

               NPTS(I) = 1

            ELSE
C
C              This is a backstop case; it should never be reached.
C
               CALL SETMSG ( 'Invalid shape code: #' )
               CALL ERRINT ( '#', SHAPE             )
               CALL SIGERR ( 'SPICE(BUG)'            )
               CALL CHKOUT ( RNAME                   )
               RETURN

            END IF


            TOTAL = TOTAL + NPTS(I) 

            IF ( NPTS(I) .GT. ROOM ) THEN

               CALL SETMSG ( 'Out of room in output arrays. Index of '
     .         //            'cutting half-plane is # out of #. '
     .         //            'Number of limb points collected so far '
     .         //            'is #. Available room is #.'             )
               CALL ERRINT ( '#', I                                   )
               CALL ERRINT ( '#', NCUTS                               )
               CALL ERRINT ( '#', TOTAL                               )
               CALL ERRINT ( '#', ROOM                                )
               CALL SIGERR ( 'SPICE(OUTOFROOM)'                       )
               CALL CHKOUT ( RNAME                                    )
               RETURN

            END IF

C
C           Transfer the limb points we found to the output limb point
C           array. Set the elements of the surface vector array as we
C           go. Store in each element of the output array the epoch
C           associated with the target center.
C
            DO J = 1, NPTS(I)

               CALL VEQU ( PNTBUF(1,J),       POINTS(1,TO) )
               CALL VSUB ( PNTBUF(1,J), AXIS, TANGTS(1,TO) )

               EPOCHS(TO) = TRGEPC

               TO         = TO + 1

            END DO

         END DO

 

      ELSE IF ( LOCCDE .EQ. ELLCOR ) THEN
C
C        Aberration corrections are done for each cutting half plane.
C        Corrections are performed for the intersections of the 
C        half plane with the reference ellipsoid's limb.
C
C        This locus is supported only for the "tangent" limb point
C        method.
C
         IF ( LMBTYP .NE. TANGNT ) THEN

            CALL SETMSG ( 'Limb type <#> is not supported for the '
     .      //            '# aberration correction locus.'        )
            CALL ERRCH  ( '#',  SVLSTR                            )
            CALL ERRCH  ( '#',  CORLOC                            )
            CALL SIGERR ( 'SPICE(BADLIMBLOCUSMIX)'                )
            CALL CHKOUT ( RNAME                                   )
            RETURN

         END IF

C
C        We need the state of the observer relative to the solar
C        system barycenter. This state is expressed relative to
C        an inertial reference frame. This state is computed once.
C
         CALL SPKSSB ( OBSCDE, ET, IREF, STOBS )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         TO    = 1
         ROOM  = MAXN
         TOTAL = 0
C
C        Loop over the half planes, collecting limb points for
C        each one.
C
         DO I = 1, NCUTS

            ROLL = ( I - 1 ) * ROLSTP


            IF ( USELT ) THEN
C 
C              We'll do an independent light time and stellar
C              aberration correction for each half plane.
C
C              Let NUMITR be the number of iterations we'll perform to
C              compute the light time.
C
               IF ( USECN ) THEN
                  NUMITR = MAXITR
               ELSE
                  NUMITR = 2
               END IF
 
               J     = 0
               LTERR = 1.D0

               DO WHILE (       ( J     .LT. NUMITR )
     .                    .AND. ( LTERR .GT. CNVLIM )  )
C
C                 LT was set either prior to this loop or
C                 during the previous loop iteration.
C
                  EPOCH  =  TOUCHD( ET - LT )

                  CALL SPKGPS ( TRGCDE, EPOCH,  IREF,
     .                          SSB,    SSBTRG, SSBLT )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF
C
C                 Compute the position of the target center relative to
C                 the observer in the inertial frame.
C
                  CALL VSUB ( SSBTRG, STOBS, PTARG )

                  IF ( USESTL ) THEN
C
C                    Apply a stellar aberration correction to the
C                    observer-target center vector. 
C
                     IF ( J .EQ. 0 ) THEN
C
C                       On the first pass, we approximate the
C                       correction by using the correction applicable
C                       to the target center.
C
                        CALL STELAB ( PTARG,  STOBS(4), STLPOS )

                     ELSE
C
C                       We apply the correction found for the previous
C                       limb point estimate.
C
                        CALL VADD ( PTARG, STLOFF, STLPOS )

                     END IF
C
C                    Set CORTRG with the vector corrected for 
C                    stellar aberration.
C
                     CALL VEQU ( STLPOS, CORTRG )

                  ELSE

                     CALL VEQU ( PTARG,  CORTRG )
                     
                  END IF             
C
C                 CORTRG is inertially referenced and includes the
C                 stellar aberration correction, if there is one. PTARG
C                 is inertially referenced and does not include the
C                 stellar aberration correction.
C
C                 Transform the aberration-corrected position vector to
C                 the target body-fixed frame; negate the result. This
C                 gives us the axis for the limb computation.
C                 
                  CALL PXFORM ( IREF, FIXREF, EPOCH, XFORM )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  CALL MXV    ( XFORM,  CORTRG, TMPVEC )
                  CALL VMINUS ( TMPVEC,         AXIS   )     
C
C                 Rotate the reference vector about the axis by 
C                 the current angle to obtain the plane vector.
C          
                  CALL VROTV ( REFVEC, AXIS, ROLL, PLNVEC )
C
C                 Find the limb, the limb center and semi-axes, and
C                 limb plane's normal vector for the current viewing
C                 geometry.
C
                  CALL EDLIMB ( SVRADI(1), SVRADI(2), SVRADI(3), 
     .                          AXIS,      LIMB                 )

                  CALL EL2CGV ( LIMB, CENTER, SMAJOR, SMINOR )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  CALL UCRSS ( SMAJOR, SMINOR, ENORML )

C
C                 Make sure ENORML points into the same half-space as
C                 AXIS.
C
                  IF ( VDOT( ENORML, AXIS ) .LT. 0.D0 ) THEN

                     CALL VSCLIP ( -1.D0, ENORML )

                  END IF

C
C                 Let CUTNML be a vector normal to the current cutting
C                 half-plane.  
C
                  CALL UCRSS ( AXIS, PLNVEC, CUTNML )

C                 The limb point we seek must lie in both the limb 
C                 plane and the cutting half-plane. Let EDIR be the 
C                 unit direction vector satisfying these constraints.
C                 
                  CALL UCRSS ( CUTNML, ENORML, EDIR )


                  IF ( VZERO(EDIR) ) THEN

                     CALL SETMSG ( 'Vector defining cutting half plane '
     .               //            'and ellipsoid limb normal vector '
     .               //            'are linearly dependent. This error '
     .               //            'occurred while computing the limb '
     .               //            'point on the reference ellipsoid '
     .               //            'in half plane #.'                  )
                     CALL ERRINT ( '#',  I                             )
                     CALL SIGERR ( 'SPICE(DEGENERATECASE)'             )
                     CALL CHKOUT ( RNAME                               )
                     RETURN
  
                  END IF

C
C                 Compute the ellipsoid limb point.
C
                  CALL SURFPT ( CENTER,    EDIR,      SVRADI(1),
     .                          SVRADI(2), SVRADI(3), EPOINT,    FND )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  IF ( .NOT. FND ) THEN

                     CALL SETMSG ( 'Limb point not found on reference '
     .               //            'ellipsoid for cutting half plane '
     .               //            'at index #. The point should '
     .               //            'always be found.'                  )
                     CALL ERRINT ( '#', I                              )
                     CALL SIGERR ( 'SPICE(BUG)'                        )
                     CALL CHKOUT ( RNAME                               )
                     RETURN

                  END IF

C
C                 In order to compute the next light time and stellar
C                 aberration correction, we need the inertially
C                 referenced vector from the observer to the light-time
C                 corrected limb point.
C
                  CALL MTXV ( XFORM,  EPOINT, IPOINT )
                  CALL VADD ( IPOINT, PTARG,  ISRFVC )

                  IF ( USESTL ) THEN
C
C                    We're correcting for stellar aberration. Another
C                    loop iteration may occur. Prepare the stellar 
C                    aberration offset for the next loop iteration.
C
C                    Convert the observer-limb vector to the inertial
C                    frame and compute the stellar aberration
C                    correction that applies to this vector.
C
                     CALL STELAB ( ISRFVC, STOBS(4), STLPOS )
                     CALL VSUB   ( STLPOS, ISRFVC,   STLOFF )

                     CALL MXV ( XFORM, STLOFF, TMPVEC )

                  END IF
C
C                 Compute the light time to the limb point.
C
                  PRVLT = LT
                  LT    = TOUCHD(  VNORM(ISRFVC) / CLIGHT()  )

C
C                 LTERR is the magnitude of the change between the
C                 current estimate of light time and the previous
C                 estimate, relative to the previous light time
C                 corrected epoch.
C 
                  LTERR = TOUCHD(    ABS( LT  - PRVLT      ) 
     .                             / MAX( 1.D0, ABS(EPOCH) )  )

                  J     = J + 1
                 
               END DO
C
C              We now have the light time and the stellar aberration
C              offset applicable to the limb point on the ellipsoid for
C              the current half plane. Compute the axis for the DSK
C              limb point computation.
C                   
C              Compute the axis in the body-fixed frame.
C
               CALL MXV    ( XFORM,  CORTRG, TMPVEC )
               CALL VMINUS ( TMPVEC, AXIS )

               EPOCH = ET - LT

            ELSE
C
C              This is the geometric case.
C
C              We'll use the observer target position vector
C              computed above the IF block that branches based
C              on CORLOC.
C
C              Compute the central axis, which is the common ray
C              vertex.
C
               CALL VMINUS ( POS, AXIS )

C
C              The target epoch matches the observer epoch.
C
               EPOCH = ET

C
C              EPOCH and AXIS are set. Reset the plane definition
C              vector PLNVEC based on the new value of AXIS.
C
               CALL VROTV ( REFVEC, AXIS, ROLL, PLNVEC )

C
C              We're ready to compute the limb point in the current
C              half-plane.
C
C
C              Find the limb, the limb center and semi-axes, and
C              limb plane's normal vector for the current viewing
C              geometry.
C
               CALL EDLIMB ( SVRADI(1), SVRADI(2), SVRADI(3), 
     .                       AXIS,      LIMB                 )

               CALL EL2CGV ( LIMB, CENTER, SMAJOR, SMINOR )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               CALL UCRSS ( SMAJOR, SMINOR, ENORML )
C
C              Make sure ENORML points into the same half-space as
C              AXIS.
C
               IF ( VDOT( ENORML, AXIS ) .LT. 0.D0 ) THEN

                  CALL VSCLIP ( -1.D0, ENORML )

               END IF

C
C              Let CUTNML be a vector normal to the current cutting
C              half-plane.  
C
               CALL UCRSS ( AXIS, PLNVEC, CUTNML )
C
C              The limb point we seek must lie in both the limb 
C              plane and the cutting half-plane. Let EDIR be the 
C              unit direction vector satisfying these constraints.
C                 
               CALL UCRSS ( CUTNML, ENORML, EDIR )

               IF ( VZERO(EDIR) ) THEN

                  CALL SETMSG ( 'Vector defining cutting half plane '
     .            //            'and ellipsoid limb normal vector '
     .            //            'are linearly dependent. This '
     .            //            'occurred while computing the limb '
     .            //            'point on the reference ellipsoid '
     .            //            'in half plane #.'                  )
                  CALL ERRINT ( '#',  I                             )
                  CALL SIGERR ( 'SPICE(DEGENERATECASE)'             )
                  CALL CHKOUT ( RNAME                               )
                  RETURN
  
               END IF
C
C              Compute the ellipsoid limb point.
C
               CALL SURFPT ( CENTER,    EDIR,      SVRADI(1),
     .                       SVRADI(2), SVRADI(3), EPOINT,    FND )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               IF ( .NOT. FND ) THEN

                  CALL SETMSG ( 'Limb point not found on reference '
     .            //            'ellipsoid for cutting half plane '
     .            //            'at index #. The point should '
     .            //            'always be found.'                  )
                  CALL ERRINT ( '#', I                              )
                  CALL SIGERR ( 'SPICE(BUG)'                        )
                  CALL CHKOUT ( RNAME                               )
                  RETURN

               END IF

            END IF
C
C           Set the output point (there's exactly 1 in all cases) and
C           the point count here. These values apply to the ellipsoid
C           case. In the DSK case, we'll update the values when we
C           know them.
C
            CALL VEQU ( EPOINT, PNTBUF(1,1) )

            NPTS(I) = 1


            IF ( SHAPE .EQ. DSKSHP ) THEN
C
C              Find the limb points on the target surface as modeled
C              by DSK data. We'll use the axis and epoch we've
C              determined from the ellipsoid approximation.
C
               CALL SCARDD ( 0, RESULT )
C
C              Note that the evaluation epoch for the surface is
C              corrected for light time.
C               
               CALL ZZTANGNT ( LMBCRV, 0.D0,   SHAPE, TRGCDE, 
     .                         NSURF,  SRFLST, FXFCDE, EPOCH,  
     .                         PLNVEC, AXIS,   SCHSTP, SOLTOL,
     .                         RESULT, PNTBUF                 )
            
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF
C
C              Update the limb point count for this cutting
C              half-plane.
C
               NPTS(I) = CARDD( RESULT )


            ELSE IF ( SHAPE .NE. ELLSHP ) THEN

               CALL SETMSG ( 'Backstop error: SHAPE = #.' )
               CALL ERRINT ( '#', SHAPE                   )
               CALL SIGERR ( 'SPICE(BUG)'                  )
               CALL CHKOUT ( RNAME                         )
               RETURN

            END IF


            TOTAL = TOTAL + NPTS(I) 

            IF ( NPTS(I) .GT. ROOM ) THEN

               CALL SETMSG ( 'Out of room in output arrays. Index of '
     .         //            'cutting half-plane is # out of #. '
     .         //            'Number of limb points collected so far '
     .         //            'is #. Available room is #.'             )
               CALL ERRINT ( '#', I                                   )
               CALL ERRINT ( '#', NCUTS                               )
               CALL ERRINT ( '#', TOTAL                               )
               CALL ERRINT ( '#', ROOM                                )
               CALL SIGERR ( 'SPICE(OUTOFROOM)'                       )
               CALL CHKOUT ( RNAME                                    )
               RETURN

            END IF
                 
C
C           Transfer the limb points we found to the output limb
C           point array. Set the elements of the surface vector
C           array as we go. In this case, we set the elements of
C           the output target epoch array as well.
C
            DO J = 1, NPTS(I)

               CALL VEQU ( PNTBUF(1,J),       POINTS(1,TO) )
               CALL VSUB ( PNTBUF(1,J), AXIS, TANGTS(1,TO) )

               EPOCHS(TO) = EPOCH

               TO         = TO + 1

            END DO
C
C           We've found the limb points and tangent vectors
C           for the Ith half-plane.
C
         END DO

      ELSE 
         
         CALL SETMSG ( 'Aberration correction locus # is not '
     .   //            'recognized.'                          )
         CALL ERRCH  ( '#',  CORLOC                           )
         CALL SIGERR ( 'SPICE(INVALIDLOCUS)'                  )
         CALL CHKOUT ( RNAME                                  )
         RETURN

      END IF


      CALL CHKOUT ( RNAME )
      RETURN
      END
