C$Procedure GFOCLT ( GF, find occultation )

      SUBROUTINE GFOCLT ( OCCTYP,  FRONT,   FSHAPE,  FFRAME,
     .                    BACK,    BSHAPE,  BFRAME,  ABCORR,
     .                    OBSRVR,  STEP,    CNFINE,  RESULT  )

C$ Abstract
C
C     Determine time intervals when an observer sees one target body
C     occulted by, or in transit across, another.
C
C     The surfaces of the target bodies may be represented by triaxial
C     ellipsoids or by topographic data provided by DSK files.
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
C     FRAMES
C     GF
C     KERNEL
C     NAIF_IDS
C     SPK
C     TIME
C     WINDOWS
C
C$ Keywords
C
C     EVENT
C     GEOMETRY
C     OCCULTATION
C     SEARCH
C     WINDOW
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE               'gf.inc'
      INCLUDE               'zzholdd.inc'

      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )

      CHARACTER*(*)         OCCTYP
      CHARACTER*(*)         FRONT
      CHARACTER*(*)         FSHAPE
      CHARACTER*(*)         FFRAME
      CHARACTER*(*)         BACK
      CHARACTER*(*)         BSHAPE
      CHARACTER*(*)         BFRAME
      CHARACTER*(*)         ABCORR
      CHARACTER*(*)         OBSRVR
      DOUBLE PRECISION      STEP
      DOUBLE PRECISION      CNFINE ( LBCELL : * )
      DOUBLE PRECISION      RESULT ( LBCELL : * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     LBCELL     P   SPICE Cell lower bound.
C     CNVTOL     P   Convergence tolerance.
C     ZZGET      P   ZZHOLDD retrieves a stored DP value.
C     GF_TOL     P   ZZHOLDD acts on the GF subsystem tolerance.
C     OCCTYP     I   Type of occultation.
C     FRONT      I   Name of body occulting the other.
C     FSHAPE     I   Type of shape model used for front body.
C     FFRAME     I   Body-fixed, body-centered frame for front body.
C     BACK       I   Name of body occulted by the other.
C     BSHAPE     I   Type of shape model used for back body.
C     BFRAME     I   Body-fixed, body-centered frame for back body.
C     ABCORR     I   Aberration correction flag.
C     OBSRVR     I   Name of the observing body.
C     STEP       I   Step size in seconds for finding occultation
C                    events.
C     CNFINE     I   SPICE window to which the search is restricted.
C     RESULT     O   SPICE window containing results.
C
C$ Detailed_Input
C
C
C     OCCTYP     indicates the type of occultation that is to be found.
C                Note that transits are considered to be a type of
C                occultation.
C
C                Supported values and corresponding definitions are:
C
C                   'FULL'      denotes the full occultation of the
C                               body designated by BACK by the body
C                               designated by FRONT, as seen from the
C                               location of the observer. In other
C                               words, the occulted body is completely
C                               invisible as seen from the observer's
C                               location.
C
C                   'ANNULAR'   denotes an annular occultation: the
C                               body designated by FRONT blocks part
C                               of, but not the limb of, the body
C                               designated by BACK, as seen from the
C                               location of the observer.
C
C                   'PARTIAL'   denotes a partial, non-annular
C                               occultation: the body designated by
C                               FRONT blocks part, but not all, of the
C                               limb of the body designated by BACK, as
C                               seen from the location of the observer.
C
C                   'ANY'       denotes any of the above three types of
C                               occultations: 'PARTIAL', 'ANNULAR', or
C                               'FULL'.
C
C                               'ANY' should be used to search for
C                               times when the body designated by FRONT
C                               blocks any part of the body designated
C                               by BACK.
C
C                               The option 'ANY' must be used if either
C                               the front or back target body is
C                               modeled as a point.
C
C                Case and leading or trailing blanks are not
C                significant in the string OCCTYP.
C
C
C     FRONT      is the name of the target body that occults---that is,
C                passes in front of---the other. Optionally, you may
C                supply the integer NAIF ID code for the body as a
C                string. For example both 'MOON' and '301' are
C                legitimate strings that designate the Moon.
C
C                Case and leading or trailing blanks are not
C                significant in the string FRONT.
C
C
C     FSHAPE     is a string indicating the geometric model used to
C                represent the shape of the front target body. The
C                supported options are:
C
C                   'ELLIPSOID'    
C
C                       Use a triaxial ellipsoid model with radius
C                       values provided via the kernel pool. A kernel
C                       variable having a name of the form
C
C                          'BODYnnn_RADII'
C
C                       where nnn represents the NAIF integer code
C                       associated with the body, must be present in
C                       the kernel pool. This variable must be
C                       associated with three numeric values giving the
C                       lengths of the ellipsoid's X, Y, and Z
C                       semi-axes.
C
C                   'POINT'       
C
C                       Treat the body as a single point. When a point
C                       target is specified, the occultation type must
C                       be set to 'ANY'.
C
C                   'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       Use topographic data provided by DSK files to
C                       model the body's shape. These data must be
C                       provided by loaded DSK files.
C
C                       The surface list specification is optional. The
C                       syntax of the list is
C
C                          <surface 1> [, <surface 2>...]
C
C                       If present, it indicates that data only for the
C                       listed surfaces are to be used; however, data
C                       need not be available for all surfaces in the
C                       list. If absent, loaded DSK data for any surface
C                       associated with the target body are used.
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
C                The combinations of the shapes of the target bodies
C                FRONT and BACK must be one of:
C
C                   One ELLIPSOID, one POINT
C                   Two ELLIPSOIDs
C                   One DSK, one POINT
C
C                Case and leading or trailing blanks are not
C                significant in the string FSHAPE.
C
C
C     FFRAME     is the name of the body-fixed, body-centered reference
C                frame associated with the front target body. Examples
C                of such names are 'IAU_SATURN' (for Saturn) and
C                'ITRF93' (for the Earth).
C
C                If the front target body is modeled as a point, FFRAME
C                should be left blank.
C
C                Case and leading or trailing blanks bracketing a
C                non-blank frame name are not significant in the string
C                FFRAME.
C
C
C     BACK       is the name of the target body that is occulted
C                by---that is, passes in back of---the other.
C                Optionally, you may supply the integer NAIF ID code
C                for the body as a string. For example both 'MOON' and
C                '301' are legitimate strings that designate the Moon.
C
C                Case and leading or trailing blanks are not
C                significant in the string BACK.
C
C
C     BSHAPE     is the shape specification for the body designated
C                by BACK. The supported options are those for
C                FSHAPE. See the description of FSHAPE above for
C                details.
C
C
C     BFRAME     is the name of the body-fixed, body-centered reference
C                frame associated with the ``back'' target body.
C                Examples of such names are 'IAU_SATURN' (for Saturn)
C                and 'ITRF93' (for the Earth).
C
C                If the back target body is modeled as a point, BFRAME
C                should be left blank.
C
C                Case and leading or trailing blanks bracketing a
C                non-blank frame name are not significant in the string
C                BFRAME.
C
C
C     ABCORR     indicates the aberration corrections to be applied to
C                the state of each target body to account for one-way
C                light time.  Stellar aberration corrections are
C                ignored if specified, since these corrections don't
C                improve the accuracy of the occultation determination.
C
C                See the header of the SPICE routine SPKEZR for a
C                detailed description of the aberration correction
C                options. For convenience, the options supported by
C                this routine are listed below:
C
C                   'NONE'     Apply no correction.
C
C                   'LT'       "Reception" case:  correct for
C                              one-way light time using a Newtonian
C                              formulation.
C
C                   'CN'       "Reception" case:  converged
C                              Newtonian light time correction.
C
C                   'XLT'      "Transmission" case:  correct for
C                              one-way light time using a Newtonian
C                              formulation.
C
C                   'XCN'      "Transmission" case:  converged
C                              Newtonian light time correction.
C
C                Case and blanks are not significant in the string
C                ABCORR.
C
C
C     OBSRVR     is the name of the body from which the occultation is
C                observed. Optionally, you may supply the integer NAIF
C                ID code for the body as a string.
C
C                Case and leading or trailing blanks are not
C                significant in the string OBSRVR.
C
C
C     STEP       is the step size to be used in the search. STEP must
C                be shorter than any interval, within the confinement
C                window, over which the specified occultation condition
C                is met. In other words, STEP must be shorter than the
C                shortest occultation event that the user wishes to
C                detect; STEP must also be shorter than the shortest
C                time interval between two occultation events that
C                occur within the confinement window (see below).
C                However, STEP must not be *too* short, or the search
C                will take an unreasonable amount of time.
C
C                The choice of STEP affects the completeness but not
C                the precision of solutions found by this routine; the
C                precision is controlled by the convergence tolerance.
C                See the discussion of the parameter CNVTOL for
C                details.
C
C                STEP has units of TDB seconds.
C
C
C     CNFINE     is a SPICE window that confines the time period over
C                which the specified search is conducted. CNFINE may
C                consist of a single interval or a collection of
C                intervals.
C
C                The endpoints of the time intervals comprising CNFINE
C                are interpreted as seconds past J2000 TDB.
C
C                See the Examples section below for a code example
C                that shows how to create a confinement window.
C
C                CNFINE must be initialized by the caller via the
C                SPICELIB routine SSIZED.
C
C
C$ Detailed_Output
C
C     RESULT     is a SPICE window representing the set of time
C                intervals, within the confinement window, when the
C                specified occultation occurs.
C
C                The endpoints of the time intervals comprising RESULT
C                are interpreted as seconds past J2000 TDB.
C
C                If RESULT is non-empty on input, its contents
C                will be discarded before GFOCLT conducts its
C                search.
C
C$ Parameters
C
C     LBCELL     is the lower bound for SPICE cell arrays.
C
C     CNVTOL     is the convergence tolerance used for finding
C                endpoints of the intervals comprising the result
C                window. CNVTOL is used to determine when binary
C                searches for roots should terminate: when a root is
C                bracketed within an interval of length CNVTOL, the
C                root is considered to have been found.
C
C                The accuracy, as opposed to precision, of roots found
C                by this routine depends on the accuracy of the input
C                data. In most cases, the accuracy of solutions will be
C                inferior to their precision.
C
C
C     See INCLUDE file gf.inc for declarations and descriptions of
C     parameters used throughout the GF system.
C
C$ Exceptions
C
C     1)  In order for this routine to produce correct results,
C         the step size must be appropriate for the problem at hand.
C         Step sizes that are too large may cause this routine to miss
C         roots; step sizes that are too small may cause this routine
C         to run unacceptably slowly and in some cases, find spurious
C         roots.
C
C         This routine does not diagnose invalid step sizes, except
C         that if the step size is non-positive, the error
C         SPICE(INVALIDSTEPSIZE) will be signaled.
C
C     2)  Due to numerical errors, in particular,
C
C            - Truncation error in time values
C            - Finite tolerance value
C            - Errors in computed geometric quantities
C
C         it is *normal* for the condition of interest to not always be
C         satisfied near the endpoints of the intervals comprising the
C         result window.
C
C         The result window may need to be contracted slightly by the
C         caller to achieve desired results. The SPICE window routine
C         WNCOND can be used to contract the result window.
C
C     3)  If name of either target or the observer cannot be translated
C         to a NAIF ID code, the error will be diagnosed by a routine
C         in the call tree of this routine.
C
C     4)  If the radii of a target body modeled as an ellipsoid cannot
C         be determined by searching the kernel pool for a kernel
C         variable having a name of the form
C
C            'BODYnnn_RADII'
C
C         where nnn represents the NAIF integer code associated with
C         the body, the error will be diagnosed by a routine in the
C         call tree of this routine.
C
C     5)  If either of the target bodies FRONT or BACK coincides with
C         the observer body OBSRVR, the error will be diagnosed by a
C         routine in the call tree of this routine.
C
C     6)  If the body designated by FRONT coincides with that
C         designated by BACK, the error will be diagnosed by a routine
C         in the call tree of this routine.
C
C     7)  If either of the body model specifiers FSHAPE or BSHAPE
C         is not recognized, the error will be diagnosed by a routine
C         in the call tree of this routine.
C
C     8)  If both of the body model specifiers FSHAPE and BSHAPE
C         specify point targets, the error will be diagnosed by a
C         routine in the call tree of this routine.
C
C     9)  If a target body-fixed reference frame associated with a
C         non-point target is not recognized, the error will be
C         diagnosed by a routine in the call tree of this routine.
C
C     10) If a target body-fixed reference frame is not centered at
C         the corresponding target body,  the error will be
C         diagnosed by a routine in the call tree of this routine.
C
C     11) If the loaded kernels provide insufficient data to
C         compute any required state vector, the deficiency will
C         be diagnosed by a routine in the call tree of this routine.
C
C     12) If an error occurs while reading an SPK or other kernel file,
C         the error will be diagnosed by a routine in the call tree
C         of this routine.
C
C     13) If a point target is specified and the occultation
C         type is set to a valid value other than 'ANY', the
C         error will be diagnosed by a routine in the call tree
C         of this routine.
C
C     14) If the output SPICE window RESULT has insufficient capacity
C         to contain the number of intervals on which the specified
C         occultation condition is met, the error will be diagnosed
C         by a routine in the call tree of this routine. If the result
C         window has size less than 2, the error SPICE(WINDOWTOOSMALL)
C         will be signaled by this routine.
C
C     15) Invalid occultation types will be diagnosed by a routine in
C         the call tree of this routine.
C
C     16) Invalid aberration correction specifications will be
C         diagnosed by a routine in the call tree of this routine.
C
C     17) If either FSHAPE or BSHAPE specifies that the target surface
C         is represented by DSK data, and no DSK files are loaded for
C         the specified target, the error is signaled by a routine in
C         the call tree of this routine.
C
C     18) If either FSHAPE or BSHAPE specifies that the target surface
C         is represented by DSK data, but the shape specification is
C         invalid, the error is signaled by a routine in the call tree
C         of this routine.
C
C
C$ Files
C
C     Appropriate SPICE kernels must be loaded by the calling program
C     before this routine is called.
C
C     The following data are required:
C
C        - SPK data: the calling application must load ephemeris data
C          for the targets, source and observer that cover the time
C          period specified by the window CNFINE. If aberration
C          corrections are used, the states of the target bodies and of
C          the observer relative to the solar system barycenter must be
C          calculable from the available ephemeris data. Typically
C          ephemeris data are made available by loading one or more SPK
C          files via FURNSH.
C
C        - PCK data: bodies modeled as triaxial ellipsoids must have
C          semi-axis lengths provided by variables in the kernel pool.
C          Typically these data are made available by loading a text
C          PCK file via FURNSH.
C
C        - FK data: if either of the reference frames designated by
C          BFRAME or FFRAME are not built in to the SPICE system,
C          one or more FKs specifying these frames must be loaded.
C
C     The following data may be required:
C
C        - DSK data: if either FSHAPE or BSHAPE indicates that DSK
C          data are to be used, DSK files containing topographic data
C          for the target body must be loaded. If a surface list is
C          specified, data for at least one of the listed surfaces must
C          be loaded.
C
C        - Surface name-ID associations: if surface names are specified
C          in FSHAPE or BSHAPE, the association of these names with
C          their corresponding surface ID codes must be established by
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
C        - CK data: either of the body-fixed frames to which FFRAME or
C          BFRAME refer might be a CK frame. If so, at least one CK
C          file will be needed to permit transformation of vectors
C          between that frame and the J2000 frame.
C
C        - SCLK data: if a CK file is needed, an associated SCLK
C          kernel is required to enable conversion between encoded SCLK
C          (used to time-tag CK data) and barycentric dynamical time
C          (TDB).
C
C     Kernel data are normally loaded once per program run, NOT every
C     time this routine is called.
C
C$ Particulars
C
C     This routine provides a simpler, but less flexible, interface
C     than does the SPICELIB routine GFOCCE for conducting searches for
C     occultation events. Applications that require support for
C     progress reporting, interrupt handling, non-default step or
C     refinement functions, or non-default convergence tolerance should
C     call GFOCCE rather than this routine.
C
C     This routine determines a set of one or more time intervals
C     within the confinement window when a specified type of
C     occultation occurs. The resulting set of intervals is returned as
C     a SPICE window.
C
C     Below we discuss in greater detail aspects of this routine's
C     solution process that are relevant to correct and efficient
C     use of this routine in user applications.
C
C
C     The Search Process
C     ==================
C
C     The search for occultations is treated as a search for state
C     transitions: times are sought when the state of the BACK body
C     changes from "not occulted" to "occulted" or vice versa.
C
C     Step Size
C     =========
C
C     Each interval of the confinement window is searched as follows:
C     first, the input step size is used to determine the time
C     separation at which the occultation state will be sampled.
C     Starting at the left endpoint of the interval, samples of the
C     occultation state will be taken at each step. If a state change
C     is detected, a root has been bracketed; at that point, the
C     "root"--the time at which the state change occurs---is found by a
C     refinement process, for example, via binary search.
C
C     Note that the optimal choice of step size depends on the lengths
C     of the intervals over which the occultation state is constant:
C     the step size should be shorter than the shortest occultation
C     duration and the shortest period between occultations, within
C     the confinement window.
C
C     Having some knowledge of the relative geometry of the targets and
C     observer can be a valuable aid in picking a reasonable step size.
C     In general, the user can compensate for lack of such knowledge by
C     picking a very short step size; the cost is increased computation
C     time.
C
C     Note that the step size is not related to the precision with which
C     the endpoints of the intervals of the result window are computed.
C     That precision level is controlled by the convergence tolerance.
C
C
C     Convergence Tolerance
C     =====================
C
C     Once a root has been bracketed, a refinement process is used to
C     narrow down the time interval within which the root must lie.
C     This refinement process terminates when the location of the root
C     has been determined to within an error margin called the
C     "convergence tolerance." The default convergence tolerance 
C     used by this routine is set by the parameter CNVTOL (defined
C     in gf.inc).
C 
C     The value of CNVTOL is set to a "tight" value so that the
C     tolerance doesn't become the limiting factor in the accuracy of
C     solutions found by this routine. In general the accuracy of input
C     data will be the limiting factor.
C
C     The user may change the convergence tolerance from the default
C     CNVTOL value by calling the routine GFSTOL, e.g.
C
C        CALL GFSTOL( tolerance value )
C
C     Call GFSTOL prior to calling this routine. All subsequent
C     searches will use the updated tolerance value.
C
C     Setting the tolerance tighter than CNVTOL is unlikely to be 
C     useful, since the results are unlikely to be more accurate. 
C     Making the tolerance looser will speed up searches somewhat,
C     since a few convergence steps will be omitted. However, in most
C     cases, the step size is likely to have a much greater effect 
C     on processing time than would the convergence tolerance.
C
C
C     The Confinement Window
C     ======================
C
C     The simplest use of the confinement window is to specify a time
C     interval within which a solution is sought.
C
C     The confinement window also can be used to restrict a search to
C     a time window over which required data (typically ephemeris
C     data, in the case of occultation searches) are known to be
C     available.
C
C     In some cases, the confinement window be used to make searches
C     more efficient. Sometimes it's possible to do an efficient search
C     to reduce the size of the time period over which a relatively
C     slow search of interest must be performed. See the "CASCADE"
C     example program in gf.req for a demonstration.
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
C        required in the FSHAPE and BSHAPE arguments.
C
C        
C        Syntax of the shape input arguments for the DSK case
C        ----------------------------------------------------
C
C        The keywords and surface list in the target shape arguments
C        FSHAPE and BSHAPE, when DSK shape models are specified, are
C        called "clauses." The clauses may appear in any order, for
C        example
C
C           DSK/<surface list>/UNPRIORITIZED
C           DSK/UNPRIORITIZED/<surface list>
C           UNPRIORITIZED/<surface list>/DSK
C
C        The simplest form of a target argument specifying use of
C        DSK data is one that lacks a surface list, for example:
C
C           'DSK/UNPRIORITIZED'
C
C        For applications in which all loaded DSK data for the target
C        body are for a single surface, and there are no competing
C        segments, the above string suffices. This is expected to be
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
C        An example of a shape argument that could be constructed
C        using one of the surface lists above is
C
C          'DSK/UNPRIORITIZED/SURFACES = "Mars MEGDR 64 PIXEL/DEG", 3'
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
C     1) Find occultations of the Sun by the Moon (that is, solar
C        eclipses) as seen from the center of the Earth over the month
C        December, 2001.
C
C        Use light time corrections to model apparent positions of Sun
C        and Moon. Stellar aberration corrections are not specified
C        because they don't affect occultation computations.
C
C        We select a step size of 3 minutes, which means we
C        ignore occultation events lasting less than 3 minutes,
C        if any exist.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C           KPL/MK
C
C           File name: standard.tm
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
C
C           \begindata
C
C              KERNELS_TO_LOAD = ( 'de421.bsp',
C                                  'pck00008.tpc',
C                                  'naif0009.tls'  )
C
C           \begintext
C
C
C       Example code begins here.
C
C
C           PROGRAM EX1
C
C           IMPLICIT NONE
C
C           INTEGER               WNCARD
C
C           CHARACTER*(*)         TIMFMT
C           PARAMETER           ( TIMFMT =
C          .   'YYYY MON DD HR:MN:SC.###### (TDB)::TDB' )
C
C           INTEGER               MAXWIN
C           PARAMETER           ( MAXWIN = 2 * 100 )
C
C           INTEGER               TIMLEN
C           PARAMETER           ( TIMLEN = 40 )
C
C           INTEGER               LBCELL
C           PARAMETER           ( LBCELL = -5 )
C
C           CHARACTER*(TIMLEN)    WIN0
C           CHARACTER*(TIMLEN)    WIN1
C           CHARACTER*(TIMLEN)    BEGSTR
C           CHARACTER*(TIMLEN)    ENDSTR
C
C           DOUBLE PRECISION      CNFINE ( LBCELL : MAXWIN )
C           DOUBLE PRECISION      ET0
C           DOUBLE PRECISION      ET1
C           DOUBLE PRECISION      LEFT
C           DOUBLE PRECISION      RESULT ( LBCELL : MAXWIN )
C           DOUBLE PRECISION      RIGHT
C           DOUBLE PRECISION      STEP
C
C           INTEGER               I
C
C     C
C     C     Load kernels.
C     C
C           CALL FURNSH ( 'standard.tm' )
C
C     C
C     C     Initialize the confinement and result windows.
C     C
C           CALL SSIZED ( MAXWIN, CNFINE )
C           CALL SSIZED ( MAXWIN, RESULT )
C
C     C
C     C     Obtain the TDB time bounds of the confinement
C     C     window, which is a single interval in this case.
C     C
C           WIN0 = '2001 DEC 01 00:00:00 TDB'
C           WIN1 = '2002 JAN 01 00:00:00 TDB'
C
C           CALL STR2ET ( WIN0, ET0 )
C           CALL STR2ET ( WIN1, ET1 )
C
C     C
C     C     Insert the time bounds into the confinement
C     C     window.
C     C
C           CALL WNINSD ( ET0, ET1, CNFINE )
C
C     C
C     C     Select a 3-minute step. We'll ignore any occultations
C     C     lasting less than 3 minutes. Units are TDB seconds.
C     C
C           STEP = 180.D0
C
C     C
C     C     Perform the search.
C     C
C           CALL GFOCLT ( 'ANY',
C          .              'MOON',  'ellipsoid', 'IAU_MOON',
C          .              'SUN',   'ellipsoid', 'IAU_SUN',
C          .              'LT',    'EARTH',     STEP,
C          .              CNFINE,  RESULT                  )
C
C
C           IF ( WNCARD(RESULT) .EQ. 0 ) THEN
C
C              WRITE (*,*) 'No occultation was found.'
C
C           ELSE
C
C              DO I = 1, WNCARD(RESULT)
C     C
C     C           Fetch and display each occultation interval.
C     C
C                 CALL WNFETD ( RESULT, I, LEFT, RIGHT )
C
C                 CALL TIMOUT ( LEFT,  TIMFMT, BEGSTR )
C                 CALL TIMOUT ( RIGHT, TIMFMT, ENDSTR )
C
C                 WRITE (*,*) 'Interval ', I
C                 WRITE (*,*) '   Start time: '//BEGSTR
C                 WRITE (*,*) '   Stop time:  '//ENDSTR
C
C              END DO
C
C           END IF
C
C           END
C
C
C     When this program was executed on a PC/Linux/g77 platform, the
C     output was:
C
C        Interval  1
C           Start time: 2001 DEC 14 20:10:14.195952 (TDB)
C           Stop time:  2001 DEC 14 21:35:50.317994 (TDB)
C
C
C
C     2) Find occultations of Titan by Saturn or of Saturn by
C        Titan as seen from the center of the Earth over the
C        last four months of 2008. Model both target bodies as
C        ellipsoids. Search for every type of occultation.
C
C        Use light time corrections to model apparent positions of
C        Saturn and Titan. Stellar aberration corrections are not
C        specified because they don't affect occultation computations.
C
C        We select a step size of 15 minutes, which means we
C        ignore occultation events lasting less than 15 minutes,
C        if any exist.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C
C           KPL/MK
C
C           File name: gfoclt_ex2.tm
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
C              File name                     Contents
C              ---------                     --------
C              de421.bsp                     Planetary ephemeris
C              sat288.bsp                    Satellite ephemeris for
C                                            Saturn
C              pck00008.tpc                  Planet orientation and
C                                            radii
C              naif0009.tls                  Leapseconds
C
C           \begindata
C
C              KERNELS_TO_LOAD = ( 'de421.bsp',
C                                  'sat286.bsp',
C                                  'pck00008.tpc',
C                                  'naif0009.tls'  )
C
C           \begintext
C
C           End of meta-kernel
C
C
C       Example code begins here.
C
C
C           PROGRAM EX2
C           IMPLICIT NONE
C     C
C     C     SPICELIB functions
C     C
C           INTEGER               WNCARD
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         TIMFMT
C           PARAMETER           ( TIMFMT =
C          .   'YYYY MON DD HR:MN:SC.###### (TDB)::TDB' )
C
C           INTEGER               MAXWIN
C           PARAMETER           ( MAXWIN = 2 * 100 )
C
C           INTEGER               TIMLEN
C           PARAMETER           ( TIMLEN = 40 )
C
C           INTEGER               BDNMLN
C           PARAMETER           ( BDNMLN = 36 )
C
C           INTEGER               FRNMLN
C           PARAMETER           ( FRNMLN = 32 )
C     C
C     C     Number of occultation types:
C     C
C           INTEGER               NTYPES
C           PARAMETER           ( NTYPES = 4 )
C     C
C     C     Occultation type name length:
C     C
C           INTEGER               OCNMLN
C           PARAMETER           ( OCNMLN = 10 )
C     C
C     C     Output line length:
C     C
C           INTEGER               LNSIZE
C           PARAMETER           ( LNSIZE = 80 )
C
C           INTEGER               LBCELL
C           PARAMETER           ( LBCELL = -5 )
C
C           CHARACTER*(BDNMLN)    BACK
C           CHARACTER*(FRNMLN)    BFRAME
C           CHARACTER*(FRNMLN)    FFRAME
C           CHARACTER*(BDNMLN)    FRONT
C           CHARACTER*(LNSIZE)    LINE
C           CHARACTER*(BDNMLN)    OBSRVR
C           CHARACTER*(OCNMLN)    OCCTYP ( NTYPES )
C           CHARACTER*(LNSIZE)    TEMPLT ( NTYPES )
C           CHARACTER*(TIMLEN)    TIMSTR
C           CHARACTER*(LNSIZE)    TITLE
C           CHARACTER*(TIMLEN)    WIN0
C           CHARACTER*(TIMLEN)    WIN1
C
C           DOUBLE PRECISION      CNFINE ( LBCELL : MAXWIN )
C           DOUBLE PRECISION      ET0
C           DOUBLE PRECISION      ET1
C           DOUBLE PRECISION      FINISH
C           DOUBLE PRECISION      RESULT ( LBCELL : MAXWIN )
C           DOUBLE PRECISION      START
C           DOUBLE PRECISION      STEP
C
C           INTEGER               I
C           INTEGER               J
C           INTEGER               K
C     C
C     C     Saved variables
C     C
C     C     The confinement and result windows CNFINE
C     C     and RESULT are saved because this practice
C     C     helps to prevent stack overflow.
C     C
C     C     The variables OCCTYP and TEMPLT are
C     C     saved to facilitate turning this main program into
C     C     a subroutine. In a main program, it's not
C     C     necessary to save these variables.
C     C
C           SAVE                  CNFINE
C           SAVE                  OCCTYP
C           SAVE                  RESULT
C           SAVE                  TEMPLT
C     C
C     C     Initial values
C     C
C           DATA                  OCCTYP / 'FULL',
C          .                               'ANNULAR',
C          .                               'PARTIAL',
C          .                               'ANY'     /
C
C           DATA                  TEMPLT /
C          .      'Condition: # occultation of # by #',
C          .      'Condition: # occultation of # by #',
C          .      'Condition: # occultation of # by #',
C          .      'Condition: # occultation of # by #'      /
C
C     C
C     C     Load kernels.
C     C
C           CALL FURNSH ( 'gfoclt_ex2.tm' )
C
C     C
C     C     Initialize the confinement and result windows.
C     C
C           CALL SSIZED ( MAXWIN, CNFINE )
C           CALL SSIZED ( MAXWIN, RESULT )
C
C     C
C     C     Obtain the TDB time bounds of the confinement
C     C     window, which is a single interval in this case.
C     C
C           WIN0 = '2008 SEP 01 00:00:00 TDB'
C           WIN1 = '2009 JAN 01 00:00:00 TDB'
C
C           CALL STR2ET ( WIN0, ET0 )
C           CALL STR2ET ( WIN1, ET1 )
C     C
C     C     Insert the time bounds into the confinement
C     C     window.
C     C
C           CALL WNINSD ( ET0, ET1, CNFINE )
C     C
C     C     Select a 15-minute step. We'll ignore any occultations
C     C     lasting less than 15 minutes. Units are TDB seconds.
C     C
C           STEP = 900.D0
C     C
C     C     The observation location is the Earth.
C     C
C           OBSRVR = 'EARTH'
C
C     C
C     C     Loop over the occultation types.
C     C
C           DO I = 1, NTYPES
C     C
C     C        For each type, do a search for both transits of
C     C        Titan across Saturn and occultations of Titan by
C     C        Saturn.
C     C
C              DO J = 1, 2
C
C                 IF ( J .EQ. 1 ) THEN
C
C                    FRONT  = 'TITAN'
C                    FFRAME = 'IAU_TITAN'
C                    BACK   = 'SATURN'
C                    BFRAME = 'IAU_SATURN'
C
C                 ELSE
C
C                    FRONT  = 'SATURN'
C                    FFRAME = 'IAU_SATURN'
C                    BACK   = 'TITAN'
C                    BFRAME = 'IAU_TITAN'
C
C                 END IF
C     C
C     C           Perform the search. The target body shapes
C     C           are modeled as ellipsoids.
C     C
C                 CALL GFOCLT ( OCCTYP(I),
C          .                    FRONT,  'ELLIPSOID', FFRAME,
C          .                    BACK,   'ELLIPSOID', BFRAME,
C          .                    'LT',   OBSRVR,      STEP,
C          .                    CNFINE, RESULT              )
C     C
C     C           Display the results.
C     C
C                 WRITE (*,*) ' '
C     C
C     C           Substitute the occultation type and target
C     C           body names into the title string:
C     C
C                 CALL REPMC ( TEMPLT(I), '#', OCCTYP(I), TITLE )
C                 CALL REPMC ( TITLE,     '#', BACK,      TITLE )
C                 CALL REPMC ( TITLE,     '#', FRONT,     TITLE )
C
C                 WRITE (*, '(A)' ) TITLE
C
C                 IF ( WNCARD(RESULT) .EQ. 0 ) THEN
C
C                    WRITE (*, '(A)' ) ' Result window is empty: '
C          .         //                'no occultation was found.'
C
C                 ELSE
C
C                    WRITE (*, '(A)' ) ' Result window start, '
C          .         //                'stop times:'
C
C                    DO K = 1, WNCARD(RESULT)
C     C
C     C                 Fetch the endpoints of the Kth interval
C     C                 of the result window.
C     C
C                       CALL WNFETD ( RESULT, K, START, FINISH )
C
C                       LINE = '  #  #'
C
C                       CALL TIMOUT ( START, TIMFMT, TIMSTR )
C
C                       CALL REPMC  ( LINE, '#', TIMSTR, LINE )
C
C                       CALL TIMOUT ( FINISH, TIMFMT, TIMSTR )
C
C                       CALL REPMC  ( LINE, '#', TIMSTR, LINE )
C
C                       WRITE ( *, '(A)' ) LINE
C
C                    END DO
C
C                 END IF
C     C
C     C           We've finished displaying the results of the
C     C           current search.
C     C
C              END DO
C     C
C     C        We've finished displaying the results of the
C     C        searches using the current occultation type.
C     C
C           END DO
C
C           WRITE (*,*) ' '
C
C           END
C
C     When this program was executed on a PC/Linux/g77 platform, the
C     output was:
C
C
C Condition: FULL occultation of SATURN by TITAN
C Result window is empty: no occultation was found.
C
C Condition: FULL occultation of TITAN by SATURN
C  Result window start, stop times:
C   2008 OCT 27 22:08:01.627053 (TDB)  2008 OCT 28 01:05:03.375236 (TDB)
C   2008 NOV 12 21:21:59.252262 (TDB)  2008 NOV 13 02:06:05.053051 (TDB)
C   2008 NOV 28 20:49:02.402832 (TDB)  2008 NOV 29 02:13:58.986344 (TDB)
C   2008 DEC 14 20:05:09.246177 (TDB)  2008 DEC 15 01:44:53.523002 (TDB)
C   2008 DEC 30 19:00:56.577073 (TDB)  2008 DEC 31 00:42:43.222909 (TDB)
C
C Condition: ANNULAR occultation of SATURN by TITAN
C  Result window start, stop times:
C   2008 OCT 19 21:29:20.599087 (TDB)  2008 OCT 19 22:53:34.518737 (TDB)
C   2008 NOV 04 20:15:38.620368 (TDB)  2008 NOV 05 00:18:59.139978 (TDB)
C   2008 NOV 20 19:38:59.647712 (TDB)  2008 NOV 21 00:35:26.725908 (TDB)
C   2008 DEC 06 18:58:34.073268 (TDB)  2008 DEC 07 00:16:17.647040 (TDB)
C   2008 DEC 22 18:02:46.288289 (TDB)  2008 DEC 22 23:26:52.712459 (TDB)
C
C Condition: ANNULAR occultation of TITAN by SATURN
C  Result window is empty: no occultation was found.
C
C Condition: PARTIAL occultation of SATURN by TITAN
C  Result window start, stop times:
C   2008 OCT 19 20:44:30.326771 (TDB)  2008 OCT 19 21:29:20.599087 (TDB)
C   2008 OCT 19 22:53:34.518737 (TDB)  2008 OCT 19 23:38:26.250580 (TDB)
C   2008 NOV 04 19:54:40.339331 (TDB)  2008 NOV 04 20:15:38.620368 (TDB)
C   2008 NOV 05 00:18:59.139978 (TDB)  2008 NOV 05 00:39:58.612935 (TDB)
C   2008 NOV 20 19:21:46.689523 (TDB)  2008 NOV 20 19:38:59.647712 (TDB)
C   2008 NOV 21 00:35:26.725908 (TDB)  2008 NOV 21 00:52:40.604703 (TDB)
C   2008 DEC 06 18:42:36.100544 (TDB)  2008 DEC 06 18:58:34.073268 (TDB)
C   2008 DEC 07 00:16:17.647040 (TDB)  2008 DEC 07 00:32:16.324244 (TDB)
C   2008 DEC 22 17:47:10.776722 (TDB)  2008 DEC 22 18:02:46.288289 (TDB)
C   2008 DEC 22 23:26:52.712459 (TDB)  2008 DEC 22 23:42:28.850542 (TDB)
C
C Condition: PARTIAL occultation of TITAN by SATURN
C  Result window start, stop times:
C   2008 OCT 27 21:37:16.970175 (TDB)  2008 OCT 27 22:08:01.627053 (TDB)
C   2008 OCT 28 01:05:03.375236 (TDB)  2008 OCT 28 01:35:49.266506 (TDB)
C   2008 NOV 12 21:01:47.105498 (TDB)  2008 NOV 12 21:21:59.252262 (TDB)
C   2008 NOV 13 02:06:05.053051 (TDB)  2008 NOV 13 02:26:18.227357 (TDB)
C   2008 NOV 28 20:31:28.522707 (TDB)  2008 NOV 28 20:49:02.402832 (TDB)
C   2008 NOV 29 02:13:58.986344 (TDB)  2008 NOV 29 02:31:33.691598 (TDB)
C   2008 DEC 14 19:48:27.094229 (TDB)  2008 DEC 14 20:05:09.246177 (TDB)
C   2008 DEC 15 01:44:53.523002 (TDB)  2008 DEC 15 02:01:36.360243 (TDB)
C   2008 DEC 30 18:44:23.485898 (TDB)  2008 DEC 30 19:00:56.577073 (TDB)
C   2008 DEC 31 00:42:43.222909 (TDB)  2008 DEC 31 00:59:17.030568 (TDB)
C
C Condition: ANY occultation of SATURN by TITAN
C  Result window start, stop times:
C   2008 OCT 19 20:44:30.326771 (TDB)  2008 OCT 19 23:38:26.250580 (TDB)
C   2008 NOV 04 19:54:40.339331 (TDB)  2008 NOV 05 00:39:58.612935 (TDB)
C   2008 NOV 20 19:21:46.689523 (TDB)  2008 NOV 21 00:52:40.604703 (TDB)
C   2008 DEC 06 18:42:36.100544 (TDB)  2008 DEC 07 00:32:16.324244 (TDB)
C   2008 DEC 22 17:47:10.776722 (TDB)  2008 DEC 22 23:42:28.850542 (TDB)
C
C Condition: ANY occultation of TITAN by SATURN
C  Result window start, stop times:
C   2008 OCT 27 21:37:16.970175 (TDB)  2008 OCT 28 01:35:49.266506 (TDB)
C   2008 NOV 12 21:01:47.105498 (TDB)  2008 NOV 13 02:26:18.227357 (TDB)
C   2008 NOV 28 20:31:28.522707 (TDB)  2008 NOV 29 02:31:33.691598 (TDB)
C   2008 DEC 14 19:48:27.094229 (TDB)  2008 DEC 15 02:01:36.360243 (TDB)
C   2008 DEC 30 18:44:23.485898 (TDB)  2008 DEC 31 00:59:17.030568 (TDB)
C
C
C
C     3) Find occultations of the Mars Reconaissance Orbiter (MRO)
C        by Mars or transits of the MRO spacecraft across Mars
C        as seen from the DSN station DSS-14 over a period of a
C        few hours on FEB 28 2015.
C
C        Use both ellipsoid and DSK shape models for Mars.
C
C        Use light time corrections to model apparent positions of
C        Mars and MRO. Stellar aberration corrections are not
C        specified because they don't affect occultation computations.
C
C        We select a step size of 3 minutes, which means we
C        ignore occultation events lasting less than 3 minutes,
C        if any exist.
C     
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C
C          KPL/MK
C
C          File: gfoclt_ex3.tm
C
C          This meta-kernel is intended to support operation of SPICE
C          example programs. The kernels shown here should not be
C          assumed to contain adequate or correct versions of data
C          required by SPICE-based user applications.
C
C          In order for an application to use this meta-kernel, the
C          kernels referenced here must be present in the user's
C          current working directory.
C
C          The names and contents of the kernels referenced
C          by this meta-kernel are as follows:
C
C             File name                        Contents
C             ---------                        --------
C             de410.bsp                        Planetary ephemeris
C             mar063.bsp                       Mars satellite ephemeris
C             pck00010.tpc                     Planet orientation and
C                                              radii
C             naif0011.tls                     Leapseconds
C             earthstns_itrf93_050714.bsp      DSN station ephemeris
C             earth_latest_high_prec.bpc       Earth orientation
C             mro_psp34.bsp                    MRO ephemeris
C             megr90n000cb_plate.bds           Plate model based on
C                                              MEGDR DEM, resolution
C                                              4 pixels/degree.
C
C          \begindata
C
C             PATH_SYMBOLS    = ( 'MRO', 'GEN' )
C
C             PATH_VALUES     = (
C                                 '/ftp/pub/naif/pds/data+'
C                                 '/mro-m-spice-6-v1.0/+'
C                                 'mrosp_1000/data/spk',
C                                 '/ftp/pub/naif/generic_kernels'
C                               )
C
C             KERNELS_TO_LOAD = ( '$MRO/de410.bsp',
C                                 '$MRO/mar063.bsp',
C                                 '$MRO/mro_psp34.bsp',
C                                 '$GEN/spk/stations/+'
C                                 'earthstns_itrf93_050714.bsp',
C                                 '$GEN/pck/earth_latest_high_prec.bpc',
C                                 'pck00010.tpc',
C                                 'naif0011.tls',
C                                 'megr90n000cb_plate.bds'
C                               )
C          \begintext
C
C
C       Example code begins here.
C 
C
C              PROGRAM EX3
C
C              IMPLICIT NONE
C        C
C        C     SPICELIB functions
C        C
C              INTEGER               WNCARD
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         META
C              PARAMETER           ( META   = 'gfoclt_ex3.tm' )
C
C              CHARACTER*(*)         TIMFMT
C              PARAMETER           ( TIMFMT =
C             .   'YYYY MON DD HR:MN:SC.###### (TDB)::TDB' )
C
C              INTEGER               MAXWIN
C              PARAMETER           ( MAXWIN = 2 * 100 )
C
C              INTEGER               CORLEN
C              PARAMETER           ( CORLEN = 10 )
C
C              INTEGER               TIMLEN
C              PARAMETER           ( TIMLEN = 40 )
C
C              INTEGER               BDNMLN
C              PARAMETER           ( BDNMLN = 36 )
C
C              INTEGER               FRNMLN
C              PARAMETER           ( FRNMLN = 32 )
C
C              INTEGER               SHPLEN
C              PARAMETER           ( SHPLEN = 100 )
C
C              INTEGER               OTYPLN
C              PARAMETER           ( OTYPLN = 20 )
C
C              INTEGER               LBCELL
C              PARAMETER           ( LBCELL = -5 )
C
C        C
C        C     Local variables
C        C
C              CHARACTER*(CORLEN)    ABCORR
C              CHARACTER*(BDNMLN)    BACK
C              CHARACTER*(FRNMLN)    BFRAME
C              CHARACTER*(SHPLEN)    BSHAPE
C              CHARACTER*(BDNMLN)    FRONT
C              CHARACTER*(SHPLEN)    FSHAPE
C              CHARACTER*(FRNMLN)    FFRAME
C              CHARACTER*(OTYPLN)    OCCTYP
C              CHARACTER*(BDNMLN)    OBSRVR
C              CHARACTER*(TIMLEN)    WIN0
C              CHARACTER*(TIMLEN)    WIN1
C              CHARACTER*(TIMLEN)    BEGSTR
C              CHARACTER*(TIMLEN)    ENDSTR
C
C              DOUBLE PRECISION      CNFINE ( LBCELL : MAXWIN )
C              DOUBLE PRECISION      ET0
C              DOUBLE PRECISION      ET1
C              DOUBLE PRECISION      LEFT
C              DOUBLE PRECISION      RESULT ( LBCELL : MAXWIN )
C              DOUBLE PRECISION      RIGHT
C              DOUBLE PRECISION      STEP
C
C              INTEGER               I
C              INTEGER               J
C              INTEGER               K
C        C
C        C     Load kernels.
C        C
C              CALL FURNSH ( META )
C
C        C
C        C     Initialize the confinement and result windows.
C        C
C              CALL SSIZED ( MAXWIN, CNFINE )
C              CALL SSIZED ( MAXWIN, RESULT )
C        C
C        C     Set the observer and aberration correction.
C        C
C              OBSRVR = 'DSS-14'
C              ABCORR = 'CN'
C        C
C        C     Set the occultation type.
C        C
C              OCCTYP = 'ANY'
C        C
C        C     Set the TDB time bounds of the confinement
C        C     window, which is a single interval in this case.
C        C
C              WIN0 = '2015 FEB 28 07:00:00 TDB'
C              WIN1 = '2015 FEB 28 12:00:00 TDB'
C
C              CALL STR2ET ( WIN0, ET0 )
C              CALL STR2ET ( WIN1, ET1 )
C        C
C        C     Insert the time bounds into the confinement
C        C     window.
C        C
C              CALL WNINSD ( ET0, ET1, CNFINE )
C        C
C        C     Select a 3-minute step. We'll ignore any occultations
C        C     lasting less than 3 minutes. Units are TDB seconds.
C        C
C              STEP = 180.D0
C
C        C
C        C     Perform both spacecraft occultation and spacecraft
C        C     transit searches.
C        C
C              WRITE (*,*) ' '
C
C              DO I = 1, 2
C
C                 IF ( I .EQ. 1 ) THEN
C        C
C        C           Perform a spacecraft occultation search.
C        C
C                    FRONT  = 'MARS'
C                    FFRAME = 'IAU_MARS'
C
C                    BACK   = 'MRO'
C                    BSHAPE = 'POINT'
C                    BFRAME = ' '
C
C                 ELSE
C        C
C        C           Perform a spacecraft transit search.
C        C
C                    FRONT  = 'MRO'
C                    FSHAPE = 'POINT'
C                    FFRAME = ' '
C
C                    BACK   = 'MARS'
C                    BFRAME = 'IAU_MARS'
C
C                 END IF
C
C
C                 DO J = 1, 2
C
C                    IF ( J .EQ. 1 ) THEN
C        C
C        C              Model the planet shape as an ellipsoid.
C        C
C                       IF ( I .EQ. 1 ) THEN
C                          FSHAPE = 'ELLIPSOID'
C                       ELSE
C                          BSHAPE = 'ELLIPSOID'
C                       END IF
C
C                    ELSE
C        C
C        C              Model the planet shape using DSK data.
C        C
C                       IF ( I .EQ. 1 ) THEN
C                          FSHAPE = 'DSK/UNPRIORITIZED'
C                       ELSE
C                          BSHAPE = 'DSK/UNPRIORITIZED'
C                       END IF
C
C                    END IF
C
C        C
C        C           Perform the spacecraft occultation or
C        C           transit search.
C
C                    IF ( I .EQ. 1 ) THEN
C                       CALL TOSTDO ( 'Using shape model '//FSHAPE     )
C                       CALL TOSTDO ( 'Starting occultation search...' )
C                    ELSE
C                       CALL TOSTDO ( 'Using shape model '//BSHAPE )
C                       CALL TOSTDO ( 'Starting transit search...' )
C                    END IF
C
C                    CALL GFOCLT ( OCCTYP,
C             .                    FRONT,  FSHAPE, FFRAME,
C             .                    BACK,   BSHAPE, BFRAME,
C             .                    ABCORR, OBSRVR, STEP,
C             .                    CNFINE, RESULT         )
C
C                    IF ( WNCARD(RESULT) .EQ. 0 ) THEN
C
C                       WRITE (*,*) 'No event was found.'
C
C                    ELSE
C
C                       DO K = 1, WNCARD(RESULT)
C        C
C        C                 Fetch and display each event interval.
C        C
C                          CALL WNFETD ( RESULT, K, LEFT, RIGHT )
C
C                          CALL TIMOUT ( LEFT,  TIMFMT, BEGSTR )
C                          CALL TIMOUT ( RIGHT, TIMFMT, ENDSTR )
C
C                          WRITE (*,*) '   Interval ', K
C                          WRITE (*,*) '      Start time: '//BEGSTR
C                          WRITE (*,*) '      Stop time:  '//ENDSTR
C
C                       END DO
C
C                    END IF
C
C                    WRITE (*,*) ' '
C
C                 END DO
C
C              END DO
C
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit
C     platform, the output was:
C 
C
C        Using shape model ELLIPSOID
C        Starting occultation search...
C            Interval            1
C               Start time: 2015 FEB 28 07:17:35.379879 (TDB)
C               Stop time:  2015 FEB 28 07:50:37.710284 (TDB)
C            Interval            2
C               Start time: 2015 FEB 28 09:09:46.920140 (TDB)
C               Stop time:  2015 FEB 28 09:42:50.497193 (TDB)
C            Interval            3
C               Start time: 2015 FEB 28 11:01:57.845730 (TDB)
C               Stop time:  2015 FEB 28 11:35:01.489716 (TDB)
C
C        Using shape model DSK/UNPRIORITIZED
C        Starting occultation search...
C            Interval            1
C               Start time: 2015 FEB 28 07:17:38.130608 (TDB)
C               Stop time:  2015 FEB 28 07:50:38.310802 (TDB)
C            Interval            2
C               Start time: 2015 FEB 28 09:09:50.314903 (TDB)
C               Stop time:  2015 FEB 28 09:42:55.369626 (TDB)
C            Interval            3
C               Start time: 2015 FEB 28 11:02:01.756296 (TDB)
C               Stop time:  2015 FEB 28 11:35:08.368384 (TDB)
C
C        Using shape model ELLIPSOID
C        Starting transit search...
C            Interval            1
C               Start time: 2015 FEB 28 08:12:21.112018 (TDB)
C               Stop time:  2015 FEB 28 08:45:48.401746 (TDB)
C            Interval            2
C               Start time: 2015 FEB 28 10:04:32.682324 (TDB)
C               Stop time:  2015 FEB 28 10:37:59.920302 (TDB)
C            Interval            3
C               Start time: 2015 FEB 28 11:56:39.757564 (TDB)
C               Stop time:  2015 FEB 28 12:00:00.000000 (TDB)
C
C        Using shape model DSK/UNPRIORITIZED
C        Starting transit search...
C            Interval            1
C               Start time: 2015 FEB 28 08:12:15.750020 (TDB)
C               Stop time:  2015 FEB 28 08:45:43.406870 (TDB)
C            Interval            2
C               Start time: 2015 FEB 28 10:04:29.031706 (TDB)
C               Stop time:  2015 FEB 28 10:37:55.565509 (TDB)
C            Interval            3
C               Start time: 2015 FEB 28 11:56:34.634642 (TDB)
C               Stop time:  2015 FEB 28 12:00:00.000000 (TDB)
C
C
C$ Restrictions
C
C     The kernel files to be used by GFOCLT must be loaded (normally
C     via the SPICELIB routine FURNSH) before GFOCLT is called.
C
C$ Literature_References
C
C    None.
C
C$ Author_and_Institution
C
C    N. J. Bachman  (JPL)
C    L. S. Elson    (JPL)
C    E. D. Wright   (JPL)
C
C$ Version
C
C-    SPICELIB Version 2.0.0 29-FEB-2016 (NJB)
C
C        Header was updated. An example program demonstrating
C        DSK usage was added.
C
C        04-MAR-2015 (NJB)
C
C        Upgraded to support surfaces represented by DSKs. 
C
C-    SPICELIB Version 1.1.0  31-AUG-2010 (EDW)
C
C        Implemented use of ZZHOLDD to allow user to alter convergence
C        tolerance.
C
C        Removed the STEP > 0 error check. The GFSSTP call includes
C        the check.
C
C-    SPICELIB Version 1.0.0  07-APR-2009 (NJB) (LSE) (EDW)
C
C-&

C$ Index_Entries
C
C     GF occultation search
C
C-&

C$ Revisions
C
C     None.
C
C-&


C
C     SPICELIB functions
C
      INTEGER               SIZED

      LOGICAL               RETURN

C
C     Local variables.
C
      DOUBLE PRECISION      TOL
      LOGICAL               OK

C
C     External routines
C
C
C     Interrupt handler:
C
      LOGICAL               GFBAIL
      EXTERNAL              GFBAIL

C
C     Routines to set step size, refine transition times
C     and report work:
C
      EXTERNAL              GFREFN
      EXTERNAL              GFREPF
      EXTERNAL              GFREPI
      EXTERNAL              GFREPU
      EXTERNAL              GFSTEP

C
C     Local parameters
C
C
C     Geometric quantity  bail switch:
C
      LOGICAL               BAIL
      PARAMETER           ( BAIL = .FALSE. )

C
C     Progress report switch:
C
      LOGICAL               RPT
      PARAMETER           ( RPT  = .FALSE. )


C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'GFOCLT' )

C
C     Note to maintenance programmer: input exception checks
C     are delegated to GFOCCE. If the implementation of that
C     routine changes, or if this routine is modified to call
C     a different routine in place of GFOCCE, then the error
C     handling performed by GFOCCE will have to be performed
C     here or in a routine called by this routine.
C
C     Check the result window's size.
C
      IF ( SIZED(RESULT) .LT. 2 ) THEN

         CALL SETMSG ( 'Result window size must be at least 2 '
     .   //            'but was #.'                            )
         CALL ERRINT ( '#',  SIZED(RESULT)                     )
         CALL SIGERR ( 'SPICE(WINDOWTOOSMALL)'                 )
         CALL CHKOUT ( 'GFOCLT'                                )
         RETURN

      END IF

C
C     Check and set the step size.
C
      CALL GFSSTP ( STEP )

C
C     Retrieve the convergence tolerance, if set.
C
      CALL ZZHOLDD ( ZZGET, GF_TOL, OK, TOL )

C
C     Use the default value CNVTOL if no stored tolerance value.
C
      IF ( .NOT. OK ) THEN

         TOL = CNVTOL

      END IF

C
C     Look for solutions.
C
      CALL GFOCCE ( OCCTYP,  FRONT,   FSHAPE,  FFRAME,
     .              BACK,    BSHAPE,  BFRAME,  ABCORR,
     .              OBSRVR,  TOL,     GFSTEP,  GFREFN,
     .              RPT,     GFREPI,  GFREPU,  GFREPF,
     .              BAIL,    GFBAIL,  CNFINE,  RESULT )


      CALL CHKOUT ( 'GFOCLT' )

      RETURN
      END
