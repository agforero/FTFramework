C$Procedure TERMPT ( Terminator points on an extended object )
 
       SUBROUTINE TERMPT ( METHOD, ILUSRC, TARGET, ET,     FIXREF, 
     .                     ABCORR, CORLOC, OBSRVR, REFVEC, ROLSTP,
     .                     NCUTS,  SCHSTP, SOLTOL, MAXN,   NPTS,   
     .                     POINTS, EPOCHS, TRMVCS                 )
      
C$ Abstract
C
C     Find terminator points on a target body. The caller specifies
C     half-planes, bounded by the illumination source center-target
C     center vector, in which to search for terminator points.
C
C     The terminator can be either umbral or penumbral. The umbral
C     terminator is the boundary of the region on the target surface
C     where no light from the source is visible. The penumbral
C     terminator is the boundary of the region on the target surface
C     where none of the light from the source is blocked by the target
C     itself.
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
C     DSK
C     GEOMETRY
C     SHADOW
C     TERMINATOR
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
      CHARACTER*(*)         ILUSRC
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
      DOUBLE PRECISION      TRMVCS ( 3, * )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     METHOD     I   Computation method.
C     ILUSRC     I   Illumination source.
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
C     NPTS       O   Counts of terminator points corresponding to cuts.
C     POINTS     O   Terminator points.
C     EPOCHS     O   Times associated with terminator points.
C     TRMVCS     O   Terminator vectors emanating from the observer.
C     
C$ Detailed_Input
C
C     METHOD   is a short string providing parameters defining
C              the computation method to be used. In the syntax
C              descriptions below, items delimited by angle brackets
C              '<>' are to be replaced by actual values. Items
C              delimited by brackets '[]' are optional.
C
C              METHOD may be assigned the following values:
C
C                 '<shadow>/<curve type>/<shape specification>'
C
C              An example of such a string is 
C
C                 'UMBRAL/TANGENT/DSK/UNPRIORITIZED'
C
C              In the METHOD string
C
C                 <shadow> may be either of the strings
C
C                    'UMBRAL'    indicates the terminator is the
C                                boundary of the portion of the surface
C                                that receives no light from the
C                                illumination source. The shape of the
C                                source is modeled as a sphere. See the
C                                Particulars section below for details.
C
C                    'PENUMBRAL' indicates the terminator is the
C                                boundary of the portion of the surface
C                                that receives all possible light from
C                                the illumination source. The shape of
C                                the source is modeled as a sphere.

C                                The penumbral terminator bounds the
C                                portion of the surface that is not
C                                subject to self-occultation of light
C                                from the illumination source. Given
C                                that the light source is modeled as a
C                                sphere, from any target surface point
C                                nearer to the source than the
C                                penumbral terminator, the source
C                                appears to be a lit disc. See the
C                                Particulars section below for details.
C
C                                   
C                 <curve type> may be either of the strings 
C
C                    'TANGENT'   for topographic (DSK) target models
C                                indicates that a terminator point is
C                                defined as the point of tangency, on
C                                the surface represented by the
C                                specified data, of a line also tangent
C                                to the illumination source. 
C
C                                For ellipsoidal target models, a
C                                terminator point is a point of
C                                tangency of a plane that is also
C                                tangent to the illumination source.
C                                See the Particulars section below for
C                                details.
C
C                                Terminator points are generated within
C                                a specified set of "cutting"
C                                half-planes that have as an edge the
C                                line containing the illumination
C                                source center-target center vector.
C                                Multiple terminator points may be
C                                found within a given half-plane, if
C                                the target body shape allows for this.
C
C                                This is the highest-accuracy method
C                                supported by this subroutine. It
C                                generally executes much more slowly
C                                than the GUIDED method described
C                                below.
C
C                    'GUIDED'    indicates that terminator points are
C                                "guided" so as to lie on rays
C                                emanating from the target body's
C                                center and passing through the
C                                terminator on the target body's
C                                reference ellipsoid. The terminator
C                                points are constrained to lie on the
C                                target body's surface. As with the
C                                'TANGENT' method (see above), cutting
C                                half-planes are used to generate
C                                terminator points.
C
C                                The GUIDED method produces a unique
C                                terminator point for each cutting
C                                half-plane. If multiple terminator
C                                point candidates lie in a given
C                                cutting half-plane, the outermost one
C                                is chosen.
C
C                                This method may be used only with the
C                                CENTER aberration correction locus
C                                (see the description of CORLOC below).
C
C                                Terminator points generated by this
C                                method are approximations; they are
C                                generally not true ray-surface tangent
C                                points. However, these approximations
C                                can be generated much more quickly
C                                than tangent points.
C
C
C                 <shape specification> may be either of the strings
C
C                    'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       The DSK option indicates that terminator point
C                       computation is to use topographic data provided
C                       by DSK files (abbreviated as "DSK data" below)
C                       to model the surface of the target body.
C 
C                       The surface list specification is optional. The
C                       syntax of the list is
C
C                          <surface 1> [, <surface 2>...]
C
C                       If present, it indicates that data only for the
C                       listed surfaces are to be used; however, data
C                       need not be available for all surfaces in the
C                       list. If the list is absent, loaded DSK data
C                       for any surface associated with the target body
C                       are used.
C
C                       The surface list may contain surface names or
C                       surface ID codes. Names containing blanks must
C                       be delimited by double quotes, for example
C
C                          SURFACES = "Mars MEGDR 128 PIXEL/DEG"
C
C                       If multiple surfaces are specified, their names
C                       or IDs must be separated by commas.
C
C                       See the Particulars section below for details
C                       concerning use of DSK data.
C
C                    
C                    'ELLIPSOID'
C
C                       The ELLIPSOID shape option generates terminator
C                       points on the target body's reference
C                       ellipsoid. When the ELLIPSOID shape is
C                       selected, The TANGENT curve option may be used
C                       with any aberration correction locus, while the
C                       GUIDED option may be used only with the CENTER
C                       locus (see the description of CORLOC below).
C
C                       When the locus is set to 'CENTER', the
C                       'TANGENT' and 'GUIDED' curve options produce
C                       the same results.
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
C     ILUSRC      is the name of the illumination source. This source
C                 may be any ephemeris object. Case, blanks, and
C                 numeric values are treated in the same way as for the
C                 input TARGET.
C
C                 The shape of the illumination source is considered
C                 to be spherical. The radius of the sphere is the
C                 largest radius of the source's reference ellipsoid.
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
C                 The output terminator points in the array POINTS and
C                 the output observer-terminator vectors in the array
C                 TRMVCS are expressed relative to this reference
C                 frame.
C
C
C     ABCORR      indicates the aberration corrections to be applied
C                 when computing the target's position and orientation.
C                 Corrections are applied at the location specified by
C                 the aberration correction locus argument CORLOC,
C                 which is described below.
C
C                 For remote sensing applications, where apparent
C                 terminator points seen by the observer are desired,
C                 normally either of the corrections
C              
C                    'LT+S' 
C                    'CN+S'
C     
C                 should be used. These and the other supported options
C                 are described below. ABCORR may be any of the 
C                 following:
C
C                    'NONE'     Apply no correction. Return the
C                               geometric terminator points on the
C                               target body.
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
C                               body, are corrected for light time. The
C                               position of the illumination source as
C                               seen from the target is corrected as
C                               well.
C
C                    'LT+S'     Correct for one-way light time and
C                               stellar aberration using a Newtonian
C                               formulation. This option modifies the
C                               locus obtained with the 'LT' option to
C                               account for the observer's velocity
C                               relative to the solar system
C                               barycenter. These corrections yield
C                               points on the apparent terminator.
C
C                    'CN'       Converged Newtonian light time
C                               correction. In solving the light time
C                               equation, the 'CN' correction iterates
C                               until the solution converges. Both the
C                               position and rotation of the target
C                               body are corrected for light time. The
C                               position of the illumination source as
C                               seen from the target is corrected as
C                               well.
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
C                        than for the ELLIPSOID TERMINATOR option.
C
C                    'ELLIPSOID TERMINATOR'
C
C                        Light time and stellar aberration corrections
C                        are applied to individual terminator points on
C                        the reference ellipsoid. For a terminator
C                        point on the surface described by topographic
C                        data, lying in a specified cutting half-plane,
C                        the unique reference ellipsoid terminator
C                        point in the same half-plane is used as the
C                        locus of the aberration corrections.
C
C                        This choice is appropriate for large target
C                        objects for which the light time from the
C                        terminator to the observer is significantly
C                        different from the light time from the target
C                        center to the observer.
C
C                        Because aberration corrections are repeated
C                        for individual terminator points,
C                        computational speed for this option is
C                        relatively slow.
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
C                 half-planes in which terminator points are to be
C                 found. Each cutting half-plane has as its edge the
C                 line containing the illumination source center-target
C                 center vector; the first half-plane contains REFVEC.
C
C                 REFVEC is expressed in the body-fixed reference frame
C                 designated by FIXREF.
C
C                 ROLSTP is an angular step by which to roll the
C                 cutting half-planes about the target-illumination
C                 source vector, which we'll call the "axis." The Ith
C                 half-plane is rotated from REFVEC about the axis in
C                 the counter-clockwise direction by (I-1)*ROLSTP.
C                 Units are radians. ROLSTP should be set to
C
C                    2*pi/NCUTS 
C
C                 to generate an approximately uniform distribution of
C                 points along the terminator.
C
C                 NCUTS is the number of cutting half-planes used to
C                 find terminator points; the angular positions of
C                 consecutive half-planes increase in the positive
C                 (counterclockwise) sense about the axis and are
C                 distributed roughly equally about that vector: each
C                 half-plane has angular separation of approximately
C
C                    ROLSTP radians
C
C                 from each of its neighbors. When the aberration
C                 correction locus is set to 'CENTER', the angular
C                 separation is the value above, up to round-off.
C                 When the locus is 'TANGENT', the separations are
C                 less uniform due to differences in the aberration
C                 corrections used for the respective terminator points.
C
C
C     SCHSTP,
C     SOLTOL      are used only for DSK-based surfaces. These inputs
C                 are, respectively, the search angular step size and
C                 solution convergence tolerance used to find tangent
C                 rays and associated terminator points within each
C                 cutting half plane.  These values are used when the
C                 METHOD argument includes the TANGENT option. In this
C                 case, terminator points are found by a two-step
C                 search process:
C
C                    1) Bracketing: starting with a direction having
C                       sufficiently small angular separation from the
C                       axis, rays emanating from the surface of the
C                       illumination source are generated within the
C                       half-plane at successively greater angular
C                       separations from the axis, where the increment
C                       of angular separation is SCHSTP. The rays are
C                       tested for intersection with the target
C                       surface. When a transition from
C                       non-intersection to intersection is found, the
C                       angular separation of a tangent ray has been
C                       bracketed.
C
C                    2) Root finding: each time a tangent ray is 
C                       bracketed, a search is done to find the angular
C                       separation from the starting direction at which
C                       a tangent ray exists. The search terminates
C                       when successive rays are separated by no more
C                       than SOLTOL. When the search converges, the
C                       last ray-surface intersection point found in
C                       the convergence process is considered to be a
C                       terminator point.
C                     
C     
C                  SCHSTP and SOLTOL have units of radians.
C
C                  Target bodies with simple surfaces---for example,
C                  convex shapes---will have a single terminator point
C                  within each cutting half-plane. For such surfaces,
C                  SCHSTP can be set large enough so that only one
C                  bracketing step is taken. A value greater than pi,
C                  for example 4.D0, is recommended.
C
C                  Target bodies with complex surfaces can have
C                  multiple terminator points within a given cutting
C                  half-plane. To find all terminator points, SCHSTP
C                  must be set to a value smaller than the angular
C                  separation of any two terminator points in any
C                  cutting half-plane, where the vertex of the angle is
C                  near a point on the surface of the illumination
C                  source. SCHSTP must not be too small, or the search
C                  will be excessively slow.
C
C                  For both kinds of surfaces, SOLTOL must be chosen so
C                  that the results will have the desired precision.
C                  Note that the choice of SOLTOL required to meet a
C                  specified bound on terminator point height errors
C                  depends on the illumination source-target distance.
C
C
C     MAXN         is the maximum number of terminator points that can
C                  be stored in the output array POINTS.
C
C
C$ Detailed_Output
C
C
C     NPTS         is an array of counts of terminator points within
C                  the specified set of cutting half-planes. The Ith
C                  element of NPTS is the terminator point count in the
C                  Ith half-plane. NPTS should be declared with length
C                  at least NCUTS.
C
C
C     POINTS       is an array containing the terminator points found
C                  by this routine. Terminator points are ordered by
C                  the indices of the half-planes in which they're
C                  found. The terminator points in a given half-plane
C                  are ordered by decreasing angular separation from
C                  the illumination source-target direction; the
C                  outermost terminator point in a given half-plane is
C                  the first of that set.
C
C                  The terminator points for the half-plane containing
C                  REFVEC occupy array elements
C
C                     POINTS(1,1) through POINTS(3,NPTS(1))
C
C                  Terminator points for the second half plane occupy
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
C                  Terminator points are expressed in the reference
C                  frame designated by FIXREF. For each terminator
C                  point, the orientation of the frame is evaluated at
C                  the epoch corresponding to the terminator point; the
C                  epoch is provided in the output array EPOCHS
C                  (described below).
C
C                  Units of the terminator points are km.
C
C
C     EPOCHS       is an array of epochs associated with the terminator
C                  points, accounting for light time if aberration
C                  corrections are used. EPOCHS contains one element
C                  for each terminator point. EPOCHS should be declared
C                  with length
C
C                     MAXN
C
C                  The element
C
C                     EPOCHS(I)
C
C                  is associated with the terminator point
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
C                  If CORLOC is set to 'ELLIPSOID TERMINATOR', all
C                  values of EPOCHS for the terminator points in a
C                  given half plane will be those for the reference
C                  ellipsoid terminator point in that half plane. That
C                  is, if aberration corrections are used, and if LT(I)
C                  is the one-way light time to the observer from the
C                  reference ellipsoid terminator point in the Ith half
C                  plane, the elements of EPOCHS for that half plane
C                  will all be set to
C
C                     ET - LT(I)
C
C
C     TRMVCS       is an array of vectors connecting the observer to
C                  the terminator points. The terminator vectors are
C                  expressed in the frame designated by FIXREF. For the
C                  Ith vector, the orientation of the frame is
C                  evaluated at the Ith epoch provided in the output
C                  array EPOCHS (described above).
C
C                  TRMVCS should be declared with dimensions
C
C                     ( 3, MAXN )
C
C                  The elements
C
C                     TRMVCS(J,I), J = 1 to 3
C
C                  are associated with the terminator point
C
C                     POINTS(J,I), J = 1 to 3
C
C                  Units of the terminator vectors are km.
C
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
C     7)  If METHOD contains an invalid terminator type, the error 
C         SPICE(INVALIDTERMTYPE) will be signaled.
C
C     8)  If the target and observer have distinct identities but are
C         at the same location the error SPICE(NOSEPARATION) is
C         signaled.
C
C     9)  If insufficient ephemeris data have been loaded prior to
C         calling TERMPT, the error will be signaled by a routine in
C         the call tree of this routine. When light time correction is
C         used, sufficient ephemeris data must be available to
C         propagate the states of both observer and target to the solar
C         system barycenter.
C
C    10)  If the computation method requires an ellipsoidal target
C         shape and triaxial radii of the target body have not been
C         loaded into the kernel pool prior to calling TERMPT, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine.
C
C         When the target shape is modeled by topographic data, radii
C         of the reference triaxial ellipsoid are still required if
C         the aberration correction locus is ELLIPSOID TERMINATOR or if
C         the terminator point generation method is GUIDED.
C
C    11)  The target must be an extended body. If the target body's
C         shape is modeled as an ellipsoid, and if any of the radii of
C         the target body are non-positive, the error will be diagnosed
C         and signaled by routines in the call tree of this routine.
C
C    12)  If PCK data specifying the target body-fixed frame
C         orientation have not been loaded prior to calling TERMPT,
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
C    17)  If the GUIDED terminator type is used with the 
C         ELLIPSOID TERMINATOR aberration correction locus, the 
C         error SPICE(BADTERMLOCUSMIX) will be signaled.
C
C    18)  If the reference vector REFVEC is the zero vector, the 
C         error SPICE(ZEROVECTOR) will be signaled.
C
C    19)  If the reference vector REFVEC and the observer target
C         vector are linearly dependent, the error 
C         SPICE(DEGENERATECASE) will be signaled.
C
C    20)  If the terminator points cannot all be stored in the output
C         POINTS array, the error SPICE(OUTOFROOM) will be signaled.
C
C    21)  If NCUTS is greater than 1, the roll step ROLSTP must be
C         positive. Otherwise, the error SPICE(INVALIDROLLSTEP) will
C         be signaled.
C
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - SPK data: ephemeris data for the target, observer, and
C          illumination source must be loaded. If aberration
C          corrections are used, the states of target and observer
C          relative to the solar system barycenter must be calculable
C          from the available ephemeris data. Typically ephemeris data
C          are made available by loading one or more SPK files via
C          FURNSH.
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
C               modeled by DSK data but one or both of the GUIDED
C               terminator definition method or the ELLIPSOID
C               TERMINATOR aberration correction locus are selected.
C
C            DSK data:
C
C               If the target shape is modeled by DSK data, DSK files
C               containing topographic data for the target body must be
C               loaded. If a surface list is specified, data for at
C               least one of the listed surfaces must be loaded.
C
C        - Shape data for the illumination source:
C
C            PCK data:
C
C               Triaxial radii for the illumination source must be
C               loaded into the kernel pool. Typically this is done by
C               loading a text PCK file via FURNSH.
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
C     Terminator definition
C     =====================
C
C     The definitions of terminators used by this routine vary
C     depending on the target surface model.
C     
C     In all cases, the surface of the illumination source is
C     modeled as a sphere. 
C
C
C     Ellipsoidal target surface model
C     --------------------------------
C
C     The umbral terminator is the boundary of the set of target
C     surface points at which the illumination source is completely
C     below the local tangent plane: the entire illumination source is
C     below the horizon as seen from any surface point on the far side,
C     relative to the source, of the umbral terminator. At an umbral
C     terminator point, the target surface tangent plane containing
C     that point is tangent to the surface of the light source as well,
C     and the outward normal vectors at the two points of tangency are
C     parallel.
C
C     The penumbral terminator is the boundary of the set of target
C     surface points at which the illumination source is completely
C     above the local tangent plane: the entire illumination source is
C     above the horizon as seen from any surface point on the near
C     side, relative to the source, of the penumbral terminator. At a
C     penumbral terminator point, the target surface tangent plane
C     containing that point is tangent to the surface of the light
C     source as well, and the outward normal vectors at the two points
C     of tangency are anti-parallel.
C
C
C     Topographic target surface model (DSK case)
C     -------------------------------------------
C
C     The concept of a plane tangent to both a topographic target
C     surface and an illumination source is problematic. If the target
C     tangent point is required to lie in a given cutting half-plane
C     bounded by the line containing the target-source vector, the
C     desired plane may not exist. In general, planes tangent to both
C     the illumination source and the target will rest upon the high
C     points of the target surface.
C
C     For topographic target surface models, this routine uses a
C     modified terminator definition: terminator points are target
C     surface points at which a line is tangent to both the target and
C     the illumination source. The line is constrained to lie in the
C     plane containing the specified cutting half-plane. The concepts
C     of umbral and penumbral terminators still apply. For umbral
C     terminator points, the common tangent line does not cross the
C     target-source line; for penumbral points, it does.
C
C     Note that for ellipsoids, the terminator definitions based on
C     tangent lines are not equivalent to the definitions based on
C     tangent planes. Typically, a plane tangent to the target
C     ellipsoid at a point found by the method described above will not
C     be tangent to the illumination source: it will be rotated about
C     the common tangent line and "cut into" the sphere representing
C     the light source. This implies that some of the source will be
C     visible at umbral terminator points and some will be blocked at
C     penumbral terminator points: both umbral and penumbral terminator
C     points found by this method will lie in a region bounded by the
C     true terminators.
C
C     The two definitions are equivalent for spherical targets.
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
C           UMBRAL/TANGENT/DSK/UNPRIORITIZED/<surface list>
C           DSK/UMBRAL/TANGENT/<surface list>/UNPRIORITIZED
C           UNPRIORITIZED/<surface list>/DSK/TANGENT/UMBRAL
C
C        The simplest form of the METHOD argument specifying use of
C        DSK data is one that lacks a surface list, for example:
C
C           'PENUMBRAL/TANGENT/DSK/UNPRIORITIZED'
C           'UMBRAL/GUIDED/DSK/UNPRIORITIZED'
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
C           'UMBRAL/TANGENT/DSK/UNPRIORITIZED/
C            SURFACES= "Mars MEGDR 64 PIXEL/DEG",3'
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
C     1) Find apparent terminator points on Phobos as seen from Mars. 
C        Use the "umbral" shadow definition.
C
C        Due to Phobos' irregular shape, the TANGENT terminator point
C        definition will be used. It suffices to compute light time and
C        stellar aberration corrections for the center of Phobos, so
C        the CENTER aberration correction locus will be used. Use
C        converged Newtonian light time and stellar aberration
C        corrections in order to model the apparent position and
C        orientation of Phobos.
C        
C        For comparison, compute terminator points using both ellipsoid
C        and topographic shape models.
C
C        Use the target body-fixed +Z axis as the reference direction
C        for generating cutting half-planes. This choice enables the
C        user to see whether the first terminator point is near the
C        target's north pole.
C
C        For each option, use just three cutting half-planes in order
C        to keep the volume of output manageable. In most applications,
C        the number of cuts and the number of resulting terminator
C        points would be much greater.
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
C        C
C        C     TERMPT example 1
C        C
C        C        Find terminator points on Phobos as seen from Mars.
C        C
C        C        Compute terminator points using the tangent
C        C        definition, using the "umbral" shadow type.
C        C        The sun is the illumination source. Perform
C        C        aberration corrections for the target center.
C        C        Use both ellipsoid and DSK shape models.
C        C
C              PROGRAM EX1
C              IMPLICIT NONE
C        C
C        C     SPICELIB functions
C        C
C              DOUBLE PRECISION      DPR
C              DOUBLE PRECISION      PI
C              DOUBLE PRECISION      VNORM
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         META
C              PARAMETER           ( META   = 'termpt_ex1.tm' )
C
C              CHARACTER*(*)         FM1
C              PARAMETER           ( FM1    =  '(A,F21.9)' )
C
C              CHARACTER*(*)         FM2
C              PARAMETER           ( FM2    =  '(1X,3F21.9)' )
C
C              INTEGER               BDNMLN
C              PARAMETER           ( BDNMLN = 36 )
C
C              INTEGER               FRNMLN
C              PARAMETER           ( FRNMLN = 32 )
C
C              INTEGER               CORLEN
C              PARAMETER           ( CORLEN = 20 )
C
C              INTEGER               MTHLEN
C              PARAMETER           ( MTHLEN = 50 )
C
C              INTEGER               NMETH
C              PARAMETER           ( NMETH  = 2 )
C
C              INTEGER               MAXN
C              PARAMETER           ( MAXN = 10000 )
C        C
C        C     Local variables
C        C
C              CHARACTER*(CORLEN)    ABCORR
C              CHARACTER*(CORLEN)    CORLOC
C              CHARACTER*(FRNMLN)    FIXREF
C              CHARACTER*(BDNMLN)    ILUSRC
C              CHARACTER*(MTHLEN)    METHOD ( NMETH )
C              CHARACTER*(BDNMLN)    OBSRVR
C              CHARACTER*(BDNMLN)    TARGET
C
C              DOUBLE PRECISION      DELROL
C              DOUBLE PRECISION      DIST
C              DOUBLE PRECISION      ET
C              DOUBLE PRECISION      LT
C              DOUBLE PRECISION      POINTS ( 3, MAXN )
C              DOUBLE PRECISION      POS    ( 3 )
C              DOUBLE PRECISION      ROLL
C              DOUBLE PRECISION      SCHSTP
C              DOUBLE PRECISION      SOLTOL
C              DOUBLE PRECISION      TRMVCS ( 3, MAXN )
C              DOUBLE PRECISION      TRGEPS ( MAXN )
C              DOUBLE PRECISION      Z      ( 3 )
C
C              INTEGER               I
C              INTEGER               J
C              INTEGER               K
C              INTEGER               M
C              INTEGER               NCUTS
C              INTEGER               NPTS   ( MAXN )
C              INTEGER               START
C        C
C        C     Initial values
C        C
C              DATA                  METHOD /
C             .               'UMBRAL/TANGENT/ELLIPSOID',
C             .               'UMBRAL/TANGENT/DSK/UNPRIORITIZED'
C             .                             /
C              DATA                  Z      / 0.D0, 0.D0, 1.D0 /
C        C
C        C     Load kernel files via the meta-kernel.
C        C
C              CALL FURNSH ( META )
C        C
C        C     Set target, observer, and target body-fixed,
C        C     body-centered reference frame.
C        C
C              ILUSRC = 'SUN'
C              OBSRVR = 'MARS'
C              TARGET = 'PHOBOS'
C              FIXREF = 'IAU_PHOBOS'
C        C
C        C     Set aberration correction and correction locus.
C        C
C              ABCORR = 'CN+S'
C              CORLOC = 'CENTER'
C        C
C        C     Convert the UTC request time string seconds past
C        C     J2000, TDB.
C        C
C              CALL STR2ET ( '2008 AUG 11 00:00:00', ET )
C        C
C        C     Compute a set of terminator points using light
C        C     time and stellar aberration corrections. Use
C        C     both ellipsoid and DSK shape models. Use an
C        C     angular step size corresponding to a height of
C        C     about 100 meters to ensure we don't miss the
C        C     terminator. Set the convergence tolerance to limit
C        C     the height convergence error to about 1 meter.
C        C     Compute 3 terminator points for each computation
C        C     method.
C        C
C        C     Get the approximate light source-target distance
C        C     at ET. We'll ignore the observer-target light
C        C     time for this approximation.
C        C
C              CALL SPKPOS ( ILUSRC, ET,  'J2000', ABCORR,
C             .              TARGET, POS, LT               )
C
C              DIST   = VNORM(POS)
C
C              SCHSTP = 1.D-1 / DIST
C              SOLTOL = 1.D-3 / DIST
C              NCUTS  = 3
C
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Light source:   '//ILUSRC
C              WRITE (*,*) 'Observer:       '//OBSRVR
C              WRITE (*,*) 'Target:         '//TARGET
C              WRITE (*,*) 'Frame:          '//FIXREF
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Number of cuts: ', NCUTS
C              WRITE (*,*) ' '
C
C              DELROL = 2*PI() / NCUTS
C
C              DO I = 1, NMETH
C
C                 CALL TERMPT ( METHOD(I), ILUSRC, TARGET, ET,
C             .                 FIXREF,    ABCORR, CORLOC, OBSRVR,
C             .                 Z,         DELROL, NCUTS,  SCHSTP,
C             .                 SOLTOL,    MAXN,   NPTS,   POINTS,
C             .                 TRGEPS,    TRMVCS                 )
C        C
C        C        Write the results.
C        C
C                 WRITE(*,*) ' '
C                 WRITE(*,*) 'Computation method = ', METHOD(I)
C                 WRITE(*,*) 'Locus              = ', CORLOC
C                 WRITE(*,*) ' '
C
C
C                 START  = 0
C
C                 DO J = 1, NCUTS
C
C                    ROLL = (J-1) * DELROL
C
C                    WRITE(*,*)   ' '
C                    WRITE(*,FM1) '  Roll angle (deg) = ', ROLL * DPR()
C                    WRITE(*,FM1) '     Target epoch  = ', TRGEPS(J)
C                    WRITE(*,*)   '    Number of terminator points  '
C             .      //           'at this roll angle: ',
C             .                   NPTS(J)
C
C                    WRITE (*,*) '      Terminator points'
C
C                    DO K = 1, NPTS(J)
C                       WRITE (*,FM2) ( POINTS(M,K+START), M = 1, 3 )
C                    END DO
C
C                    START = START + NPTS(J)
C
C                 END DO
C
C                 WRITE (*,*) ' '
C
C              END DO
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was:
C
C
C      Light source:   SUN
C      Observer:       MARS
C      Target:         PHOBOS
C      Frame:          IAU_PHOBOS
C
C      Number of cuts:            3
C
C
C      Computation method = UMBRAL/TANGENT/ELLIPSOID
C      Locus              = CENTER
C
C
C       Roll angle (deg) =           0.000000000
C          Target epoch  =   271684865.152078211
C          Number of terminator points  at this roll angle:            1
C            Terminator points:
C                2.040498332          5.012722925          8.047281838
C
C       Roll angle (deg) =         120.000000000
C          Target epoch  =   271684865.152078211
C          Number of terminator points  at this roll angle:            1
C            Terminator points:
C              -11.058054707          0.167672089         -4.782740292
C
C       Roll angle (deg) =         240.000000000
C          Target epoch  =   271684865.152078211
C          Number of terminator points  at this roll angle:            1
C            Terminator points:
C                8.195238564         -6.093889437         -5.122310498
C
C
C      Computation method = UMBRAL/TANGENT/DSK/UNPRIORITIZED
C      Locus              = CENTER
C
C
C       Roll angle (deg) =           0.000000000
C          Target epoch  =   271684865.152078211
C          Number of terminator points  at this roll angle:            1
C            Terminator points:
C                1.626396028          3.995432180          8.853689595
C
C       Roll angle (deg) =         120.000000000
C          Target epoch  =   271684865.152078211
C          Number of terminator points  at this roll angle:            1
C            Terminator points:
C              -11.186659928         -0.142366793         -4.646136984
C
C       Roll angle (deg) =         240.000000000
C          Target epoch  =   271684865.152078211
C          Number of terminator points  at this roll angle:            1
C            Terminator points:
C                9.338447202         -6.091352186         -5.960849442
C
C
C
C     2) Find apparent terminator points on Mars as seen from the
C        earth.
C
C        Use both the "umbral" and "penumbral" shadow definitions. Use
C        only ellipsoid shape models for easier comparison. Find
C        distances between corresponding terminator points on the
C        umbral and penumbral terminators.
C
C        Use the ELLIPSOID TERMINATOR aberration correction locus
C        in order to perform separate aberration corrections for
C        each terminator point. Because of the large size of Mars,
C        corrections for the target center are less accurate.
C
C        For each option, use just three cutting half-planes, in order
C        to keep the volume of output manageable. In most applications,
C        the number of cuts and the number of resulting terminator
C        points would be much greater.
C
C        Use the meta-kernel below to load the required SPICE 
C        kernels. 
C
C
C           KPL/MK
C
C           File: termpt_ex2.tm
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
C              megr90n000cb_plate.bds           Plate model based on
C                                               MEGDR DEM, resolution
C                                               4 pixels/degree.
C
C           \begindata
C
C              KERNELS_TO_LOAD = ( 'de430.bsp',
C                                  'mar097.bsp',
C                                  'pck00010.tpc',
C                                  'naif0011.tls',
C                                  'megr90n000cb_plate.bds' )
C
C           \begintext
C
C
C
C     Example code begins here.
C
C
C      C
C      C     TERMPT example 2
C      C
C      C        Find terminator points on Mars as seen from the
C      C        earth.
C      C
C      C        Use only ellipsoid shape models. Use the
C      C        ELLIPSOID TERMINATOR aberration correction
C      C        locus.
C      C
C      C        Use both UMBRAL and PENUMBRAL shadow definitions.
C      C        Compute the distances between corresponding
C      C        umbral and penumbral terminator points.
C      C
C      C        Check terminator points by computing solar
C      C        incidence angles at each point.
C      C
C      C
C            PROGRAM EX2
C            IMPLICIT NONE
C      C
C      C     SPICELIB functions
C      C
C            DOUBLE PRECISION      DPR
C            DOUBLE PRECISION      PI
C            DOUBLE PRECISION      VDIST
C            DOUBLE PRECISION      VNORM
C      C
C      C     Local parameters
C      C
C            CHARACTER*(*)         META
C            PARAMETER           ( META    = 'termpt_ex2.tm' )
C
C            CHARACTER*(*)         FM1
C            PARAMETER           ( FM1     =  '(A,F21.9)' )
C
C            CHARACTER*(*)         FM2
C            PARAMETER           ( FM2     =  '(A,I2)' )
C
C            INTEGER               BDNMLN
C            PARAMETER           ( BDNMLN = 36 )
C
C            INTEGER               FRNMLN
C            PARAMETER           ( FRNMLN = 32 )
C
C            INTEGER               CORLEN
C            PARAMETER           ( CORLEN = 20 )
C
C            INTEGER               MTHLEN
C            PARAMETER           ( MTHLEN = 50 )
C
C            INTEGER               NMETH
C            PARAMETER           ( NMETH  = 2 )
C
C            INTEGER               MAXN
C            PARAMETER           ( MAXN   = 100 )
C      C
C      C     Local variables
C      C
C            CHARACTER*(CORLEN)    ABCORR
C            CHARACTER*(CORLEN)    CORLOC ( NMETH )
C            CHARACTER*(FRNMLN)    FIXREF
C            CHARACTER*(MTHLEN)    ILUMTH ( NMETH )
C            CHARACTER*(BDNMLN)    ILUSRC
C            CHARACTER*(BDNMLN)    OBSRVR
C            CHARACTER*(BDNMLN)    TARGET
C            CHARACTER*(MTHLEN)    METHOD ( NMETH )
C
C            DOUBLE PRECISION      ADJANG
C            DOUBLE PRECISION      ALT
C            DOUBLE PRECISION      ANGSRC
C            DOUBLE PRECISION      DELROL
C            DOUBLE PRECISION      DIST
C            DOUBLE PRECISION      EMISSN
C            DOUBLE PRECISION      ET
C            DOUBLE PRECISION      F
C            DOUBLE PRECISION      ILUPOS ( 3 )
C            DOUBLE PRECISION      LAT
C            DOUBLE PRECISION      LON
C            DOUBLE PRECISION      LT
C            DOUBLE PRECISION      PHASE
C            DOUBLE PRECISION      POINTS ( 3, MAXN )
C            DOUBLE PRECISION      SVPNTS ( 3, MAXN )
C            DOUBLE PRECISION      TPTILU ( 3 )
C            DOUBLE PRECISION      RADII  ( 3 )
C            DOUBLE PRECISION      RE
C            DOUBLE PRECISION      ROLL
C            DOUBLE PRECISION      RP
C            DOUBLE PRECISION      SCHSTP
C            DOUBLE PRECISION      SOLAR
C            DOUBLE PRECISION      SOLTOL
C            DOUBLE PRECISION      SRCRAD ( 3 )
C            DOUBLE PRECISION      SRFVEC ( 3 )
C            DOUBLE PRECISION      TRMVCS ( 3, MAXN )
C            DOUBLE PRECISION      TRGEPC
C            DOUBLE PRECISION      TRGEPS ( MAXN )
C            DOUBLE PRECISION      Z      ( 3 )
C
C            INTEGER               I
C            INTEGER               J
C            INTEGER               K
C            INTEGER               M
C            INTEGER               N
C            INTEGER               NCUTS
C            INTEGER               NPTS   ( MAXN )
C            INTEGER               START
C
C      C
C      C     Saved variables
C      C
C            SAVE                  METHOD
C      C
C      C     Initial values
C      C
C            DATA                  CORLOC /
C           .                        'ELLIPSOID TERMINATOR',
C           .                        'ELLIPSOID TERMINATOR'
C           .                             /
C
C            DATA                  ILUMTH /
C           .                        'ELLIPSOID',
C           .                        'ELLIPSOID'
C           .                             /
C
C            DATA                  METHOD /
C           .                      'UMBRAL/ELLIPSOID/TANGENT',
C           .                      'PENUMBRAL/ELLIPSOID/TANGENT'
C           .                             /
C
C            DATA                  Z      / 0.D0, 0.D0, 1.D0 /
C      C
C      C     Load kernel files via the meta-kernel.
C      C
C            CALL FURNSH ( META )
C      C
C      C     Set target, observer, and target body-fixed,
C      C     body-centered reference frame.
C      C
C            ILUSRC = 'SUN'
C            OBSRVR = 'EARTH'
C            TARGET = 'MARS'
C            FIXREF = 'IAU_MARS'
C      C
C      C     Set the aberration correction. We'll set the
C      C     correction locus below.
C      C
C            ABCORR = 'CN+S'
C      C
C      C     Convert the UTC request time string seconds past
C      C     J2000, TDB.
C      C
C            CALL STR2ET ( '2008 AUG 11 00:00:00', ET )
C      C
C      C     Look up the target body's radii. We'll use these to
C      C     convert Cartesian to planetographic coordinates. Use
C      C     the radii to compute the flattening coefficient of
C      C     the reference ellipsoid.
C      C
C            CALL BODVRD ( TARGET, 'RADII', 3, N, RADII )
C      C
C      C     Compute the flattening coefficient for planetodetic
C      C     coordinates
C      C
C            RE = RADII(1)
C            RP = RADII(3)
C            F  = ( RE - RP ) / RE
C
C      C
C      C     Get the radii of the illumination source as well.
C      C     We'll use these radii to compute the angular radius
C      C     of the source as seen from the terminator points.
C      C
C            CALL BODVRD ( ILUSRC, 'RADII', 3, N, SRCRAD )
C      C
C      C     Compute a set of terminator points using light time and
C      C     stellar aberration corrections. Use both ellipsoid
C      C     and DSK shape models.
C      C
C      C     Get the approximate light source-target distance
C      C     at ET. We'll ignore the observer-target light
C      C     time for this approximation.
C      C
C            CALL SPKPOS ( ILUSRC, ET,     FIXREF, ABCORR,
C           .              TARGET, ILUPOS, LT             )
C
C            DIST = VNORM( ILUPOS )
C      C
C      C     Set the angular step size so that a single step will
C      C     be taken in the root bracketing process; that's all
C      C     that is needed since we don't expect to have multiple
C      C     terminator points in any cutting half-plane.
C      C
C            SCHSTP = 4.D0
C      C
C      C     Set the convergence tolerance to minimize the
C      C     height error. We can't achieve the precision
C      C     suggested by the formula because the sun-Mars
C      C     distance is about 2.4e8 km. Compute 3 terminator
C      C     points for each computation method.
C      C
C            SOLTOL = 1.D-7/DIST
C      C
C      C     Set the number of cutting half-planes and roll step.
C      C
C            NCUTS  = 3
C            DELROL = 2*PI() / NCUTS
C
C            WRITE (*,*) ' '
C            WRITE (*,*) 'Light source:          '//ILUSRC
C            WRITE (*,*) 'Observer:              '//OBSRVR
C            WRITE (*,*) 'Target:                '//TARGET
C            WRITE (*,*) 'Frame:                 '//FIXREF
C            WRITE (*,*) 'Aberration Correction: '//ABCORR
C            WRITE (*,*) ' '
C            WRITE (*,*) 'Number of cuts: ', NCUTS
C
C            DO I = 1, NMETH
C
C               CALL TERMPT ( METHOD(I), ILUSRC, TARGET,    ET,
C           .                 FIXREF,    ABCORR, CORLOC(I), OBSRVR,
C           .                 Z,         DELROL, NCUTS,     SCHSTP,
C           .                 SOLTOL,    MAXN,   NPTS,      POINTS,
C           .                 TRGEPS,    TRMVCS                    )
C      C
C      C        Write the results.
C      C
C               WRITE(*,*) ' '
C               WRITE(*,*) 'Computation method = ', METHOD(I)
C               WRITE(*,*) 'Locus              = ', CORLOC(I)
C
C
C               START  = 0
C
C               DO J = 1, NCUTS
C
C                  ROLL = (J-1) * DELROL
C
C                  WRITE(*,*)   ' '
C                  WRITE(*,FM1) '   Roll angle (deg) = ', ROLL * DPR()
C                  WRITE(*,FM1) '    Target epoch    = ', TRGEPS(J)
C                  WRITE(*,FM2) '    Number of terminator points at '
C           .      //           'this roll angle: ',
C           .                   NPTS(J)
C
C                  DO K = 1, NPTS(J)
C
C                     WRITE (*,*) '    Terminator point planetodetic '
C           .         //          'coordinates:'
C
C                     CALL RECGEO ( POINTS(1,K+START), RE,  F,
C           .                       LON,               LAT, ALT )
C
C                     WRITE (*,FM1) '      Longitude       (deg): ',
C           .                       LON*DPR()
C                     WRITE (*,FM1) '      Latitude        (deg): ',
C           .                       LAT*DPR()
C                     WRITE (*,FM1) '      Altitude         (km): ',
C           .                       ALT
C
C      C
C      C              Get illumination angles for this terminator point.
C      C
C                     M = K+START
C
C                     CALL ILLUMG ( ILUMTH,      TARGET, ILUSRC, ET,
C           .                       FIXREF,      ABCORR, OBSRVR,
C           .                       POINTS(1,M), TRGEPC, SRFVEC,
C           .                       PHASE,       SOLAR,  EMISSN )
C
C                     WRITE (*,FM1) '      Incidence angle '
C           .         //            '(deg): ', SOLAR * DPR()
C
C
C      C
C      C              Adjust the incidence angle for the angular
C      C              radius of the illumination source. Use the
C      C              epoch associated with the terminator point
C      C              for this lookup.
C      C
C                     CALL SPKPOS ( ILUSRC, TRGEPS(M), FIXREF,
C           .                       ABCORR, TARGET,    TPTILU, LT )
C
C                     DIST   = VNORM( TPTILU )
C
C                     ANGSRC = ASIN (  MAX( SRCRAD(1),
C           .                               SRCRAD(2),
C           .                               SRCRAD(3) )  / DIST  )
C
C                     IF ( I .EQ. 1 ) THEN
C      C
C      C                 For points on the umbral terminator,
C      C                 the ellipsoid outward normal is tilted
C      C                 away from the terminator-source center
C      C                 direction by the angular radius of the
C      C                 source. Subtract this radius from the
C      C                 illumination incidence angle to get the
C      C                 angle between the local normal and the
C      C                 direction to the corresponding tangent
C      C                 point on the source.
C      C
C                        ADJANG = SOLAR - ANGSRC
C
C                     ELSE
C      C
C      C                 For the penumbral case, the outward
C      C                 normal is tilted toward the illumination
C      C                 source by the angular radius of the
C      C                 source. Adjust the illumination
C      C                 incidence angle for this.
C      C
C                        ADJANG = SOLAR + ANGSRC
C
C                     END IF
C
C                     WRITE (*,FM1)  '      Adjusted angle  '
C           .         //             '(deg): ', ADJANG * DPR()
C
C
C                     IF ( I .EQ. 1 ) THEN
C      C
C      C                 Save terminator points for comparison.
C      C
C                        CALL VEQU ( POINTS(1,M), SVPNTS(1,M) )
C
C                     ELSE
C      C
C      C                 Compare terminator points with last
C      C                 saved values.
C      C
C                        DIST = VDIST( POINTS(1,M), SVPNTS(1,M) )
C
C                        WRITE (*,FM1)
C           .            '      Distance offset  (km): ', DIST
C                     END IF
C
C
C                  END DO
C
C                  START = START + NPTS(J)
C
C               END DO
C
C               WRITE (*,*) ' '
C
C            END DO
C            END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was: 
C
C
C        Light source:          SUN
C        Observer:              EARTH
C        Target:                MARS
C        Frame:                 IAU_MARS
C        Aberration Correction: CN+S
C
C        Number of cuts:            3
C
C        Computation method = UMBRAL/ELLIPSOID/TANGENT
C        Locus              = ELLIPSOID TERMINATOR
C
C          Roll angle (deg) =           0.000000000
C           Target epoch    =   271683700.369686902
C           Number of terminator points at this roll angle:  1
C            Terminator point planetodetic coordinates:
C             Longitude       (deg):           4.189318082
C             Latitude        (deg):          66.416132677
C             Altitude         (km):           0.000000000
C             Incidence angle (deg):          90.163842885
C             Adjusted angle  (deg):          89.999999980
C
C          Roll angle (deg) =         120.000000000
C           Target epoch    =   271683700.372003794
C           Number of terminator points at this roll angle:  1
C            Terminator point planetodetic coordinates:
C             Longitude       (deg):         107.074551917
C             Latitude        (deg):         -27.604435701
C             Altitude         (km):           0.000000000
C             Incidence angle (deg):          90.163842793
C             Adjusted angle  (deg):          89.999999888
C
C          Roll angle (deg) =         240.000000000
C           Target epoch    =   271683700.364983618
C           Number of terminator points at this roll angle:  1
C            Terminator point planetodetic coordinates:
C             Longitude       (deg):         -98.695906077
C             Latitude        (deg):         -27.604435700
C             Altitude         (km):          -0.000000000
C             Incidence angle (deg):          90.163843001
C             Adjusted angle  (deg):          90.000000096
C
C
C        Computation method = PENUMBRAL/ELLIPSOID/TANGENT
C        Locus              = ELLIPSOID TERMINATOR
C
C          Roll angle (deg) =           0.000000000
C           Target epoch    =   271683700.369747400
C           Number of terminator points at this roll angle:  1
C            Terminator point planetodetic coordinates:
C             Longitude       (deg):           4.189317837
C             Latitude        (deg):          66.743818467
C             Altitude         (km):           0.000000000
C             Incidence angle (deg):          89.836157094
C             Adjusted angle  (deg):          89.999999999
C             Distance offset  (km):          19.483590936
C
C          Roll angle (deg) =         120.000000000
C           Target epoch    =   271683700.372064054
C           Number of terminator points at this roll angle:  1
C            Terminator point planetodetic coordinates:
C             Longitude       (deg):         107.404259674
C             Latitude        (deg):         -27.456458359
C             Altitude         (km):           0.000000000
C             Incidence angle (deg):          89.836157182
C             Adjusted angle  (deg):          90.000000087
C             Distance offset  (km):          19.411414247
C
C          Roll angle (deg) =         240.000000000
C           Target epoch    =   271683700.365043879
C           Number of terminator points at this roll angle:  1
C            Terminator point planetodetic coordinates:
C             Longitude       (deg):         -99.025614323
C             Latitude        (deg):         -27.456458357
C             Altitude         (km):           0.000000000
C             Incidence angle (deg):          89.836156972
C             Adjusted angle  (deg):          89.999999877
C             Distance offset  (km):          19.411437239
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
C     N.J. Bachman   (JPL)
C
C$ Version
C
C     SPICELIB Version 1.0.0 04-APR-2017 (NJB)
C
C        11-MAR-2016 (NJB)
C
C        Changed ellipsoid algorithm to use ZZEDTMPT. Added ROLSTP
C        argument. Updated calls to ZZTANGNT to accommodate argument
C        list change. Added code examples. Updated Detailed_Input. Made
C        various header corrections.
C
C        Original version 18-NOV-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     find terminator points on target body
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      CLIGHT
      DOUBLE PRECISION      TOUCHD
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
      PARAMETER           ( RNAME  = 'TERMPT' )

      CHARACTER*(*)         IREF
      PARAMETER           ( IREF   = 'J2000' )

C
C     Convergence limit:
C
      DOUBLE PRECISION      CNVLIM
      PARAMETER           ( CNVLIM = 1.D-17 )

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
C     Local variables
C     
      CHARACTER*(CVTLEN)    PNTDEF
      CHARACTER*(ACLLEN)    NRMLOC
      CHARACTER*(CORLEN)    PRVCOR
      CHARACTER*(ACLLEN)    PRVLOC
      CHARACTER*(MTHLEN)    PRVMTH
      CHARACTER*(SHPLEN)    SHPSTR
      CHARACTER*(SUBLEN)    SUBTYP
      CHARACTER*(TMTLEN)    TRMSTR

      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      CP     ( 3 )
      DOUBLE PRECISION      EDIR   ( 3 )
      DOUBLE PRECISION      EPOCH
      DOUBLE PRECISION      ICORVC ( 3 )
      DOUBLE PRECISION      ILUMLT
      DOUBLE PRECISION      ILURAD
      DOUBLE PRECISION      ITRMVC ( 3 )
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      LTERR
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      PLNVEC ( 3 )
      DOUBLE PRECISION      PNTBUF ( 3, MAXWIN )
      DOUBLE PRECISION      PRVLT
      DOUBLE PRECISION      PTARG  ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      RAYVTX ( 3 )
      DOUBLE PRECISION      RESULT ( LBCELL : MAXWIN )
      DOUBLE PRECISION      ROLL
      DOUBLE PRECISION      SRFVEC ( 3 )
      DOUBLE PRECISION      SSBLT
      DOUBLE PRECISION      SSBTRG ( 3 )
      DOUBLE PRECISION      STLOFF ( 3 )
      DOUBLE PRECISION      STOBS  ( 6 )
      DOUBLE PRECISION      TMPVEC ( 3 )
      DOUBLE PRECISION      TRGEPC
      DOUBLE PRECISION      TRGPOS ( 3 )
      DOUBLE PRECISION      XFORM  ( 3, 3 )

      INTEGER               FXCENT
      INTEGER               FXCLSS
      INTEGER               FXFCDE
      INTEGER               FXTYID
      INTEGER               I
      INTEGER               ILUCDE
      INTEGER               J
      INTEGER               LOCCDE
      INTEGER               NRAD
      INTEGER               NSURF
      INTEGER               NUMITR
      INTEGER               OBSCDE
      INTEGER               PRVILU
      INTEGER               PRVTRG
      INTEGER               ROOM
      INTEGER               SHADOW
      INTEGER               SHAPE
      INTEGER               SRFLST ( MAXSRF )
      INTEGER               SVLCOD
      INTEGER               TO
      INTEGER               TOTAL
      INTEGER               TRGCDE
      INTEGER               TRMTYP
      
      LOGICAL               ATTBLK ( NABCOR )
      LOGICAL               FIRST
      LOGICAL               FND
      LOGICAL               PRI
      LOGICAL               SURFUP
      LOGICAL               UFLAG
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

      INTEGER               SVCTR3 ( CTRSIZ )
      CHARACTER*(MAXL)      SVILUM
      INTEGER               SVICDE
      LOGICAL               SVFND3

C
C     Saved frame name/ID item declarations.
C
      INTEGER               SVCTR4 ( CTRSIZ )
      CHARACTER*(FRNMLN)    SVFREF
      INTEGER               SVFXFC


C
C     Saved surface name/ID item declarations.
C
      INTEGER               SVCTR5 ( CTRSIZ )

C
C     Saved target radii declarations.
C
      INTEGER               SVCTR6 ( CTRSIZ )
      DOUBLE PRECISION      SVTRAD ( 3 )

      INTEGER               SVCTR7 ( CTRSIZ )
      DOUBLE PRECISION      SVSRAD ( 3 )

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  ILUCDE
      SAVE                  LOCCDE
      SAVE                  MAXRAD
      SAVE                  NRAD
      SAVE                  NSURF
      SAVE                  PNTBUF
      SAVE                  PRI
      SAVE                  PRVCOR
      SAVE                  PRVILU
      SAVE                  PRVLOC
      SAVE                  PRVMTH
      SAVE                  PRVTRG
      SAVE                  SHADOW
      SAVE                  SHAPE
      SAVE                  SRFLST
      SAVE                  SUBTYP
      SAVE                  SVLCOD
      SAVE                  TRGCDE
      SAVE                  TRMTYP
      SAVE                  UFLAG
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

      SAVE                  SVCTR3
      SAVE                  SVILUM
      SAVE                  SVICDE
      SAVE                  SVFND3

C
C     Saved frame name/ID items.
C
      SAVE                  SVCTR4
      SAVE                  SVFREF
      SAVE                  SVFXFC


C
C     Saved surface name/ID items.
C
      SAVE                  SVCTR5
 
C
C     Saved reference ellipsoid items.
C
      SAVE                  SVCTR6
      SAVE                  SVTRAD
      SAVE                  SVCTR7
      SAVE                  SVSRAD

C
C     Initial values
C
      DATA                  FIRST  / .TRUE.  /
      DATA                  NRAD   / 0       /
      DATA                  PRVCOR / ' '     /
      DATA                  PRVILU / 0       /
      DATA                  PRVLOC / ' '     /
      DATA                  PRVMTH / ' '     /
      DATA                  PRVTRG / 0       /
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
         CALL ZZCTRUIN( SVCTR6 )

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
     .      //            'are not supported for terminator '
     .      //            'finding.'                              )
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
C     Obtain an integer code for the illumination source.
C 
      CALL ZZBODS2C ( SVCTR3, SVILUM, SVICDE, SVFND3,
     .                ILUSRC, ILUCDE, FND    )
            
      IF ( .NOT. FND ) THEN
      
         CALL SETMSG ( 'The illumination source, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit, or that you '
     .   //            'failed to load a kernel containing a '
     .   //            'name-ID mapping for this body.'           )
         CALL ERRCH  ( '#', ILUSRC                                )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( RNAME                                      )
         RETURN
      
      END IF
      
C
C     Check the observer and target body codes. If they are equal,
C     signal an error. The illumination source must be distinct
C     from the target as well.
C
      IF ( OBSCDE .EQ. TRGCDE ) THEN
 
         CALL SETMSG ( 'In computing the terminator, '      
     .   //            'the observing body and target body are the '
     .   //            'same. Both are #.'                          )
         CALL ERRCH  ( '#',  OBSRVR                                 )
         CALL SIGERR ( 'SPICE(BODIESNOTDISTINCT)'                   )
         CALL CHKOUT ( RNAME                                        )
         RETURN
 
      END IF

      IF ( ILUCDE .EQ. TRGCDE ) THEN
 
         CALL SETMSG ( 'In computing the terminator, the ' 
     .   //            'observing body and illumination source '
     .   //            'are the same. Both are #.'               )
         CALL ERRCH  ( '#',  OBSRVR                              )
         CALL SIGERR ( 'SPICE(BODIESNOTDISTINCT)'                )
         CALL CHKOUT ( RNAME                                     )
         RETURN
 
      END IF

C
C     Determine the attributes of the frame designated by FIXREF.
C
      CALL ZZNAMFRM ( SVCTR4, SVFREF, SVFXFC, FIXREF, FXFCDE )

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
      CALL ZZSRFTRK ( SVCTR5, SURFUP )

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
     .                   PRI,    NSURF,  SRFLST, PNTDEF, TRMSTR )

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


         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF


         IF (  EQSTR( PNTDEF, 'TANGENT' )  ) THEN

            TRMTYP = TANGNT

         ELSE IF (  EQSTR( PNTDEF, 'GUIDED' )  ) THEN

            TRMTYP = GUIDED

         ELSE

            CALL SETMSG ( 'Returned point definition from '
     .      //            'method string was <#>. Value '
     .      //            'must be TANGENT or GUIDED.'      )
            CALL ERRCH  ( '#', PNTDEF                       )
            CALL SIGERR ( 'SPICE(INVALIDTERMTYPE)'          )
            CALL CHKOUT ( RNAME                             )
            RETURN

         END IF


         IF (  EQSTR( TRMSTR, 'UMBRAL' )  ) THEN

            SHADOW = UMBRAL
            UFLAG  = .TRUE.

         ELSE IF (  EQSTR( TRMSTR, 'PENUMBRAL' )  ) THEN

            SHADOW = PNMBRL
            UFLAG  = .FALSE.

         ELSE

            CALL SETMSG ( 'Returned shadow type from method '
     .      //            'string was <#>. Value must be '
     .      //            'UMBRAL or PENUMBRAL.'            )
            CALL ERRCH  ( '#', TRMSTR                       )
            CALL SIGERR ( 'SPICE(INVALIDSHADOW)'            )
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
     .      //            'SUBSLR, but is not applicable for TERMPT.' )
            CALL ERRCH  ( '#', SUBTYP                                 )
            CALL ERRCH  ( '#', METHOD                                 )
            CALL SIGERR ( 'SPICE(INVALIDMETHOD)'                      )
            CALL CHKOUT ( RNAME                                       )
            RETURN

         END IF

C
C        Save the current method as the previous method that we've
C        successfully processed the input method.
C         
         PRVMTH = METHOD

      END IF

C
C     Identify the aberration correction locus.
C
      IF (  FIRST  .OR.  ( CORLOC .NE. PRVLOC )  ) THEN

         CALL LJUCRS ( 1, CORLOC, NRMLOC )

         IF ( NRMLOC .EQ. 'CENTER' ) THEN

            LOCCDE = CTRCOR

         ELSE IF ( NRMLOC .EQ. 'ELLIPSOID TERMINATOR' ) THEN

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
         SVLCOD = LOCCDE
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
C     Check NCUTS; there must be room for at least one terminator point
C     for each cut. NCUTS may not be negative.
C
      IF (  ( NCUTS .LT. 0 ) .OR. ( NCUTS .GT. MAXN )  ) THEN
         
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
C     Check ROLSTP, if this step is needed.
C
      IF ( NCUTS .GT. 1 )  THEN
         
         IF ( ROLSTP .EQ. 0.D0 ) THEN

            CALL SETMSG ( 'ROLSTP is zero; NCUTS = #; the roll '
     .      //            'step is required to be non-zero '
     .      //            'when more than one cutting half-plane '
     .      //            'is used.'                              )
            CALL ERRINT ( '#',  NCUTS                             )
            CALL SIGERR ( 'SPICE(INVALIDROLLSTEP)'                )
            CALL CHKOUT ( RNAME                                   )
            RETURN

         END IF

      END IF



      IF ( SHAPE .EQ. DSKSHP ) THEN
C
C        This is the DSK case.
C
C        Initialize the intercept algorithm to use a DSK
C        model for the surface of the target body. 
C        
         CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FXFCDE )
C
C        Save the radius of the outer bounding sphere of
C        the target.
C
         CALL ZZMAXRAD ( MAXRAD )


      ELSE IF ( SHAPE .NE. ELLSHP ) THEN

         CALL SETMSG ( 'Computation method argument was <#>; this ' 
     .   //            'string must specify a supported shape '
     .   //            'model and computation type. See the '
     .   //            'description of METHOD in the header of '
     .   //            'TERMPT for details.'                      )
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
C     Get illumination source radii.
C
      IF ( ILUCDE .NE. PRVILU ) THEN
C
C        Reset counter to force lookup.
C           
         CALL ZZCTRUIN ( SVCTR7 )

      END IF

C
C     Look up illumination source radii using counter.
C
      CALL ZZBODVCD ( ILUCDE, 'RADII', 3, SVCTR7, NRAD, SVSRAD )

      IF ( FAILED() ) THEN
C
C        Make sure we don't reuse the outputs from ZZBODVCD.
C
         PRVILU = 0

         CALL CHKOUT ( RNAME )
         RETURN

      END IF

      IF ( NRAD .NE. 3 ) THEN

         CALL SETMSG ( 'Number of illumination source radii '
     .   //            'must be 3 but was #.'                )
         CALL ERRINT ( '#',  NRAD                            )
         CALL SIGERR ( 'SPICE(BADRADIUSCOUNT)'               )
         CALL CHKOUT ( RNAME                                 )
         RETURN

      END IF

C
C     Obtain the largest radius of the source.
C
      ILURAD = MAX( SVSRAD(1), SVSRAD(2), SVSRAD(3) )

      PRVILU = ILUCDE

C
C     Get target body radii if necessary.
C
      IF (      ( SHAPE   .EQ. ELLSHP ) 
     .     .OR. ( SVLCOD  .EQ. ELLCOR ) 
     .     .OR. ( TRMTYP  .EQ. GUIDED )  ) THEN

         IF ( TRGCDE .NE. PRVTRG ) THEN
C
C           Reset counter to force lookup.
C           
            CALL ZZCTRUIN ( SVCTR6 )

         END IF
C
C        Look up target radii using counter.
C
         CALL ZZBODVCD ( TRGCDE, 'RADII', 3, SVCTR6, NRAD, SVTRAD )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         IF ( NRAD .NE. 3 ) THEN

            CALL SETMSG ( 'Number of target radii must be 3 '
     .      //            'but was #.'                       )
            CALL ERRINT ( '#',  NRAD                         )
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
C     Find terminator points on the target.
C      
      CALL CLEARI ( NCUTS,  NPTS   )
      CALL SSIZED ( MAXWIN, RESULT )

C
C     Get initial observer-target vector, expressed in the target
C     body-fixed frame, evaluated at the target epoch. This vector
C     will be used for all option combinations.
C
      CALL SPKPOS ( TARGET, ET, FIXREF, ABCORR, OBSRVR, TRGPOS, LT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      IF ( VZERO(TRGPOS) ) THEN

         CALL SETMSG ( 'The distance between the observer and '
     .   //            'target at ET # is zero.'               )
         CALL ERRDP  ( '#',  ET                                )
         CALL SIGERR ( 'SPICE(NOSEPARATION)'                   )
         CALL CHKOUT ( RNAME                                   )
         RETURN

      END IF
C
C     The terminator-finding technique depends on the aberration
C     correction locus. Start with the 'CENTER' version, since this is
C     the simpler case.
C
      IF ( SVLCOD .EQ. CTRCOR ) THEN
C
C        Aberration corrections are those applicable at the target
C        center.
C
C        Compute the epoch associated with the target center.
C
         CALL ZZCOREPC ( ABCORR, ET, LT, TRGEPC )

C
C        Get the vector from the target center to the illumination
C        source.
C
         CALL SPKPOS ( ILUSRC, TRGEPC, FIXREF, ABCORR, 
     .                 TARGET, AXIS,   ILUMLT          )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        The target-source vector is the central axis. 
C
C        Make sure the reference vector and axis are linearly 
C        independent.
C
         CALL VCRSS ( AXIS, REFVEC, CP )

         IF ( VZERO(CP) ) THEN
            
            CALL SETMSG ( 'Input reference vector and illumination '
     .      //            'source-target vector are linearly '
     .      //            ' dependent.'                             )
            CALL SIGERR ( 'SPICE(DEGENERATECASE)'                   )
            CALL CHKOUT ( RNAME                                     )
            RETURN

         END IF

 
         TO    = 1
         ROOM  = MAXN
         TOTAL = 0

C
C        Loop over the half planes, collecting terminator points for
C        each one.
C
         DO I = 1, NCUTS

            ROLL = ( I - 1 ) * ROLSTP

C
C           Rotation of the half-planes is in the positive
C           sense about AXIS.
C
            CALL VROTV ( REFVEC, AXIS, ROLL, PLNVEC )


            IF ( SHAPE .EQ. DSKSHP ) THEN
C
C              This is the DSK case.
C
               IF ( TRMTYP .EQ. TANGNT ) THEN
C
C                 This type of solution finds actual tangent rays on
C                 the target.
C
C                 Find the terminator points that lie in the current
C                 half-plane.
C
C                 Note that RESULT is a cell, not a window.
C
                  CALL SCARDD ( 0, RESULT )
C
C                 Note that the evaluation epoch for the surface is
C                 optionally corrected for light time.
C
C                 For this computation, the ray's vertex is computed
C                 on the fly, since it depends on the ray's direction.
C                 The location of the center of the source is passed
C                 to the tangent utilities instead.
C
                  CALL ZZTANGNT ( SHADOW, ILURAD, SHAPE,  TRGCDE,  
     .                            NSURF,  SRFLST, FXFCDE, TRGEPC,
     .                            PLNVEC, AXIS,   SCHSTP, SOLTOL,
     .                            RESULT, PNTBUF                 )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  NPTS(I) = CARDD( RESULT )


               ELSE IF ( TRMTYP .EQ. GUIDED ) THEN
C
C                 This option uses the target's reference ellipsoid for
C                 guidance. For DSK shapes, the DSK terminator points
C                 are generated by finding terminator points on the
C                 target body's reference ellipsoid, then finding
C                 topographic surface intercepts of rays emanating from
C                 the target body's center and passing through the
C                 terminator points on the ellipsoid. If multiple
C                 intercepts are found for a given ray, the outermost
C                 is selected.
C
C                 Find the terminator point on the ellipsoid in the
C                 current cutting half-plane.
C
                  CALL ZZEDTMPT ( UFLAG,    SVTRAD(1), SVTRAD(2), 
     .                            SVTRAD(3), ILURAD,    AXIS,
     .                            PLNVEC,    PNTBUF               ) 
                 
                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  NPTS(I) = 1
 
                  CALL VHAT ( PNTBUF(1,1), EDIR )
C
C                 Find the intercept on the target surface of the ray
C                 emanating from the target in the direction of EDIR.
C                 We must use a ray pointed in the opposite direction
C                 to perform this computation, since the surface may be
C                 invisible from the interior of the target.
C
                  CALL VSCL  ( 3.D0*MAXRAD, EDIR, RAYVTX )
                  CALL VMINUS( EDIR, RAYDIR )

                  CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FXFCDE )

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
                  CALL SETMSG ( 'Invalid terminator type code: #' )
                  CALL ERRINT ( '#', TRMTYP                       )
                  CALL SIGERR ( 'SPICE(BUG)'                      )
                  CALL CHKOUT ( RNAME                             )
                  RETURN

               END IF


            ELSE IF ( SHAPE .EQ. ELLSHP ) THEN               
C
C              This is the ellipsoid case.
C              
C              Find the terminator point in the current cutting
C              half-plane.
C
               CALL ZZEDTMPT ( UFLAG,    SVTRAD(1), SVTRAD(2), 
     .                         SVTRAD(3), ILURAD,    AXIS,
     .                         PLNVEC,    PNTBUF               ) 

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               CALL SCARDD ( 0, RESULT )

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
     .         //            'Number of terminator points collected '
     .         //            'so far is #. Available room is #.'      )
               CALL ERRINT ( '#', I                                   )
               CALL ERRINT ( '#', NCUTS                               )
               CALL ERRINT ( '#', TOTAL                               )
               CALL ERRINT ( '#', ROOM                                )
               CALL SIGERR ( 'SPICE(OUTOFROOM)'                       )
               CALL CHKOUT ( RNAME                                    )
               RETURN

            END IF

C
C           Transfer the terminator points we found to the output
C           terminator point array. Set the elements of the tangent
C           vector array as we go. Store in each element of the output
C           array the epoch associated with the target center.
C
            DO J = 1, NPTS(I)

               CALL VEQU ( PNTBUF(1,J),         POINTS(1,TO) )
               CALL VADD ( PNTBUF(1,J), TRGPOS, TRMVCS(1,TO) )

               EPOCHS( TO ) = TRGEPC

               TO           = TO + 1

               ROOM         = ROOM - NPTS(I)

            END DO

         END DO

 

      ELSE IF ( SVLCOD .EQ. ELLCOR ) THEN
C
C        Aberration corrections are done for each cutting half plane.
C        Corrections are performed for the intersections of the 
C        half plane with the reference ellipsoid's terminator.
C
C        This locus is supported only for the "tangent" terminator
C        point method.
C
         IF ( TRMTYP .NE. TANGNT ) THEN

            CALL SETMSG ( 'Terminator point definition type <#> '
     .      //            'is not supported for the '
     .      //            '# aberration correction locus.'        )
            CALL ERRCH  ( '#',  PNTDEF                            )
            CALL ERRCH  ( '#',  CORLOC                            )
            CALL SIGERR ( 'SPICE(BADTERMLOCUSMIX)'                )
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
C        Loop over the half planes, collecting terminator points for
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
                  NUMITR = 1
               END IF
 
               J     = 0
               LTERR = 1.D0

               DO WHILE (       ( J     .LT. NUMITR ) 
     .                    .AND. ( LTERR .GT. CNVLIM )  )
C
C                 LT was set either prior to this loop or
C                 during the previous loop iteration.
C
                  EPOCH  =  ET - LT 

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

C
C                 Find the position of the illumination source relative
C                 to the target center at EPOCH.
C
                  CALL SPKEZP ( ILUCDE, EPOCH, FIXREF, ABCORR,
     .                          TRGCDE, AXIS,  ILUMLT          ) 

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF                 
C
C                 The illumination source position vector gives us
C                 the axis for the terminator computation.
C
C                 Let PLNVEC be the secondary vector defining the
C                 cutting half-plane. Rotation of the half-planes is in
C                 the positive sense about AXIS.
C
                  CALL VROTV ( REFVEC, AXIS, ROLL, PLNVEC )

C
C                 Find the terminator point on the reference 
C                 ellipsoid, in the cutting half-plane.
C
                  CALL ZZEDTMPT ( UFLAG,    SVTRAD(1), SVTRAD(2), 
     .                            SVTRAD(3), ILURAD,    AXIS,
     .                            PLNVEC,    PNTBUF               ) 

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  NPTS(I) = 1 
C
C
C                 Compute the vector from the observer to the terminator
C                 point.
C
                  CALL PXFORM ( IREF, FIXREF, EPOCH, XFORM )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( RNAME )
                     RETURN
                  END IF

                  CALL MXV  ( XFORM, PTARG, TRGPOS )
                  CALL VADD ( PNTBUF(1,1),  TRGPOS, SRFVEC )

C
C                 Compute the light time to the terminator point.
C
                  PRVLT = LT
                  LT    = TOUCHD(  VNORM(SRFVEC) / CLIGHT()  )

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
C              We now have the light time offset (but not the stellar
C              aberration correction) applicable to the terminator
C              point on the ellipsoid for the current half plane.
C              Compute the vertex and axis for the terminator point
C              computation.
C                           
               EPOCH = ET - LT
C
C              Compute the position of the target at EPOCH relative
C              to the observer at ET. This vector is computed in
C              an inertial frame.
C
               CALL SPKGPS ( TRGCDE, EPOCH, IREF, SSB, SSBTRG, SSBLT )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF
C
C              Compute the position of the target center relative to
C              the observer in the inertial frame.
C
               CALL VSUB ( SSBTRG, STOBS, PTARG )
C
C              Transform the observer-target position to the body-fixed
C              frame at EPOCH.
C
               CALL PXFORM ( IREF, FIXREF, EPOCH, XFORM )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

C
C              If we're using stellar aberration corrections, find the
C              correction applicable to the ellipsoid terminator point.
C              
               IF ( USESTL ) THEN
C
C                 The vector ICORVC below is the inertially-referenced
C                 stellar aberration correction.
C                  
                  CALL MTXV ( XFORM,  PNTBUF(1,1), TMPVEC )
                  CALL VADD ( TMPVEC, PTARG,       ITRMVC )

                  CALL STELAB ( ITRMVC, STOBS(4),  TMPVEC )
                  CALL VSUB   ( TMPVEC, ITRMVC,    ICORVC )
                  
                  CALL MXV  ( XFORM, ICORVC, STLOFF )                  

               END IF


               CALL MXV  ( XFORM,  PTARG, TMPVEC )
               CALL VEQU ( TMPVEC,        TRGPOS )

C
C              Find the apparent position of the illumination source
C              relative to the target at EPOCH.
C
               CALL SPKEZP ( ILUCDE, EPOCH, FIXREF, ABCORR,
     .                       TRGCDE, AXIS,  ILUMLT          ) 

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF                 

            ELSE
C
C              This is the geometric case.
C
C              Get the position of the illumination source
C              as seen from the target at ET.
C
               CALL SPKEZP ( ILUCDE, ET,   FIXREF, ABCORR,
     .                       TRGCDE, AXIS, ILUMLT          ) 

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF                 
C
C              The target epoch matches the observer epoch.
C
               EPOCH = ET
C
C              The position of the target relative to the observer
C              is already present in POS.
C
            END IF                        

C
C           POS, EPOCH, and AXIS are set.
C
C           Reset the plane definition vector PLNVEC based on the new
C           value of AXIS.
C
            CALL VROTV ( REFVEC, AXIS, ROLL, PLNVEC )
C
C           We're ready to compute the terminator point in the current
C           half-plane.
C
            IF ( SHAPE .EQ. DSKSHP ) THEN
C
C              Find the terminator points on the target surface as
C              modeled by DSK data.
C
               CALL SCARDD ( 0, RESULT )
C
C              Note that the evaluation epoch for the surface is
C              corrected for light time.
C                           
               CALL ZZTANGNT ( SHADOW, ILURAD, SHAPE, TRGCDE, NSURF, 
     .                         SRFLST, FXFCDE, EPOCH,  PLNVEC, AXIS,
     .                         SCHSTP, SOLTOL, RESULT, PNTBUF        )
            
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               NPTS(I) = CARDD( RESULT )

            ELSE

C
C              Find the terminator point on the target surface modeled
C              by an ellipsoid.
C
C              If we performed a light time computation, we already
C              have the answer stored in PNTBUF. If this is a geometric
C              computation, we still need to compute the terminator
C              point.
C
               IF ( .NOT. USELT ) THEN

                  CALL ZZEDTMPT ( UFLAG,    SVTRAD(1), SVTRAD(2), 
     .                            SVTRAD(3), ILURAD,    AXIS,
     .                            PLNVEC,    PNTBUF               ) 
               END IF
               
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF
 
               NPTS(I) = 1

            END IF

            TOTAL   = TOTAL + NPTS(I) 

            IF ( NPTS(I) .GT. ROOM ) THEN

               CALL SETMSG ( 'Out of room in output arrays. Index of '
     .         //            'cutting half-plane is # out of #. '
     .         //            'Number of terminator points collected '
     .         //            'so far is #. Available room is #.'      )
               CALL ERRINT ( '#', I                                   )
               CALL ERRINT ( '#', NCUTS                               )
               CALL ERRINT ( '#', TOTAL                               )
               CALL ERRINT ( '#', ROOM                                )
               CALL SIGERR ( 'SPICE(OUTOFROOM)'                       )
               CALL CHKOUT ( RNAME                                    )
               RETURN

            END IF
                 
C
C           Transfer the terminator points we found to the output
C           terminator point array. Set the elements of the tangent
C           vector array as we go. In this case, we set the elements of
C           the output target epoch array as well.
C
            DO J = 1, NPTS(I)

               CALL VEQU ( PNTBUF(1,J),         POINTS(1,TO) )
               CALL VADD ( PNTBUF(1,J), TRGPOS, TRMVCS(1,TO) )

               IF ( USESTL ) THEN
C
C                 Apply the stellar aberration offset for the current
C                 half-plane to each terminator vector in the output
C                 buffer.
C
                  CALL VADD ( TRMVCS(1,TO), STLOFF, TMPVEC )
                  CALL VEQU ( TMPVEC,               TRMVCS(1,TO) )

               END IF

               EPOCHS(TO) = EPOCH

               TO         = TO + 1

               ROOM       = ROOM - NPTS(I)

            END DO
C
C           We've found the terminator points and observer-terminator
C           vectors for the Ith half-plane.
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
