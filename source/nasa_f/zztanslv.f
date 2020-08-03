C$Procedure ZZTANSLV ( Private --- tangent point solver )

      SUBROUTINE ZZTANSLV ( UDCOND, UDSTEP, UDREFN, 
     .                      CSTEP,  STEP,   START,  FINISH,
     .                      TOL,    RESULT, POINTS, ENDFLG  )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     This routine finds tangent points of rays on a target surface,
C     where the rays are confined to a specified half-plane. It may
C     be used for limb and terminator computations.
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
C     ROOT
C     SEARCH
C     WINDOWS
C
C$ Declarations

      IMPLICIT NONE

      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )

      EXTERNAL              UDCOND
      EXTERNAL              UDSTEP
      EXTERNAL              UDREFN
      LOGICAL               CSTEP
      DOUBLE PRECISION      STEP
      DOUBLE PRECISION      START
      DOUBLE PRECISION      FINISH
      DOUBLE PRECISION      TOL
      DOUBLE PRECISION      RESULT ( LBCELL : * )
      DOUBLE PRECISION      POINTS ( 3, * )
      LOGICAL               ENDFLG ( 2 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     UDCOND     I   Name of the routine that compares the current
C                    state and ray intercept.
C     UDSTEP     I   Name of the routine that computes a step
C     UDREFN     I   Name of the routine that computes a refined input.
C     CSTEP      I   Logical indicating constant step size.
C     STEP       I   Constant step size for finding geometric events.
C     START      I   Beginning of the search interval.
C     FINISH     I   End of the search interval.
C     TOL        I   Maximum error in detection of state transitions.
C     RESULT    I-O  SPICE window containing results.
C     POINTS     O   Array of points associated with transitions.
C     ENDFLG     O   Endpoint transition flags.  
C
C$ Detailed_Input
C
C     For the purpose of solving for ray tangency points on a target
C     body, the independent variable can be considered to be angular
C     separation of a ray from an axis vector. For a limb computation,
C     the axis vector points from the target body center toward the
C     observer. For a terminator computation, the axis vector points
C     from the target body center towards the center of the
C     illumination source. The "system state" for these computations is
C     the condition of the ray intersecting the target body: if an
C     intersection exists, the state is "true."
C
C     The discussion below is more general; we'll use the terms
C     "abscissa" or "x-value" rather than "angle" to describe the
C     independent variable. The "system state" is simply a boolean
C     function of the independent variable.
C     
C     The first three inputs to this routine are names of subroutines
C     that this routine will call.
C
C     These routines should meet the following specifications.
C
C     UDCOND     the routine that determines if the system state
C                satisfies some constraint at a given independent
C                variable value X.
C
C                The calling sequence:
C
C                   CALL UDCOND ( X, IN_CON, POINT )
C
C                where:
C
C                   X        a double precision value at which to
C                            evaluate the state.
C
C                   IN_CON   a logical value indicating whether
C                            or not the quantity satisfies the 
C                            constraint at X (TRUE) or not (FALSE).
C
C                   POINT    is a 3-vector associated with X. POINT
C                            is defined if and only if IN_CON is .TRUE.
C
C
C     UDSTEP     the routine that computes a step in an attempt to
C                find a transition of the state of the specified 
C                coordinate. In the context of this routine's algorithm,
C                a "state transition" occurs where the geometric state 
C                changes from being in the desired geometric condition 
C                event to not, or vice versa.
C
C                This routine relies on UDSTEP returning step sizes
C                small enough so that state transitions within the
C                interval [START, FINISH] are not overlooked.  There
C                must never be two roots A and B separated by less than
C                STEP, where STEP is the minimum step size returned by
C                UDSTEP for any value of X in the interval [A, B].
C
C                The calling sequence for UDSTEP is:
C
C                   CALL UDSTEP ( X, STEP )
C
C                where:
C
C                   X       a double precision value from which the
C                           algorithm is to search forward for a state
C                           transition.
C
C                   STEP    is the output step size. STEP indicates how
C                           far to advance X so that X and X+STEP may
C                           bracket a state transition and definitely
C                           do not bracket more than one state
C                           transition.
C
C                If a constant step size is desired, the routine
C
C                   GFSTEP
C
C                may be used. This is the default option. If using
C                GFSTEP, the step size must be set by calling
C
C                   GFSSTP(STEP)
C
C                prior to calling this routine.
C
C
C     UDREFN     the routine that computes a refinement in the abscissa
C                values that bracket a transition point. In other
C                words, once a pair of abscissa values have been
C                detected such that the system is in different states
C                at each of the two values, UDREFN selects an
C                intermediate abscissa value which should be closer to
C                the transition state than one of the two known X
C                values. The calling sequence for UDREFN is:
C
C                   CALL UDREFN ( X1, X2, S1, S2, T )
C
C                where the inputs are:
C
C                   X1    an X (abscissa) value at which the system is
C                         in state S1.
C
C                   X2    an X value at which the system is in state
C                         S2. X2 is assumed to be larger than X1.
C
C                   S1    a logical indicating the state of the system
C                         at X1.
C
C                   S2    a logical indicating the state of the system
C                         at X2.
C
C                UDREFN may use or ignore the S1 and S2 values.
C
C                The output is:
C
C                   T     an X value to check for a state transition
C                         between X1 and X2.
C
C                If a simple bisection method is desired, the routine
C                GFREFN may be used. This is the default option.
C
C     CSTEP      is a logical indicating whether or not the step size
C                used in searching is constant.  If it is, the value
C                STEP is used. Note that even if UDSTEP has the value
C                GFSTEP, i.e. the public, constant step routine, CSTEP
C                should still be .FALSE., in which case STEP is ignored.
C
C     STEP       is the step size to be used in the search. STEP must
C                be short enough for a search using this step size to
C                locate the intervals where the geometric event
C                function is monotone increasing or decreasing.
C                However, STEP must not be *too* short, or the search
C                will take an unreasonable amount of time.
C
C                The choice of STEP affects the completeness but not
C                the precision of solutions found by this routine;
C                precision is controlled by the convergence the
C                tolerance, TOL.
C
C     START      is the beginning of the interval over which the state
C                is to be detected.
C
C     FINISH     is the end of the interval over which the state is
C                to be detected.
C
C     TOL        is a tolerance value used to determine convergence of
C                root-finding operations. TOL is measured in the units
C                of the independent variable and is greater than zero.
C
C     RESULT     is an initialized SPICE window. RESULT must be large
C                enough to hold all of the intervals found by the
C                search.
C
C$ Detailed_Output
C
C     RESULT     is a SPICE window containing the results of the
C                search. With the exception of the first and last
C                endpoints of the window, the endpoints of the
C                intervals in RESULT always are abscissa values of
C                state transitions. The first and last endpoints may or
C                may not correspond to state transitions.
C
C                The first and last endpoints are abscissa values of
C                state transitions if and only if the state function is
C                .FALSE. at those points.
C
C                Note that, in the special case where the state function
C                is .TRUE. at the first endpoint and .FALSE. in a
C                half-neighborhood to the right of that endpoint, it is
C                possible for this function to find a transition at that
C                endpoint and assign it to the right endpoint of the
C                degenerate interval
C
C                    [ left endpoint,  left endpoint ]
C
C                Analogously, it is possible to find a state transition
C                at the last endpoint if the state function is .TRUE.
C                at that endpoint and .FALSE. in a half-neighborhood
C                to the left of that point.
C
C                The output ENDFLG indicates whether the first and
C                last endpoints of RESULT are transitions.
C
C
C     POINTS     is an array of 3-vectors associated with the endpoints
C                of the intervals of RESULT. Elements
C
C                   POINTS(J,I), J = 1 .. 3
C
C                constitute a vector associated with 
C
C                   RESULT(I)
C
C                The first and last vectors of POINTS are valid if
C                and only if the corresponding elements of ENDFLG
C                are .TRUE.
C
C                Presuming this routine is used to solve for tangent
C                points on a target body, the vectors contained in
C                POINTS, when valid, are such tangent points.
C               
C                POINTS must be declared with size at least 3 times
C                that of RESULT.
C               
C 
C     ENDFLG     is an array of two logical flags that indicate
C                whether state transitions occur at the initial
C                and final endpoints of the result window. Element
C                1 of this array is .TRUE. if and only if
C                there is a state transition at the first endpoint of
C                RESULT (in element RESULT(1)); element 2 is .TRUE.
C                if and only if there is a state transition at the
C                last element of RESULT.
C               
C
C$ Parameters
C
C     LBCELL     is the SPICELIB cell lower bound.
C
C$ Exceptions
C
C     1)  If TOL is negative, the error SPICE(INVALIDTOLERANCE)
C         will be signaled.
C
C     2)  If START +/- TOL is indistinguishable from START or
C         FINISH +/- TOL is indistinguishable from FINISH, the
C         error SPICE(INVALIDTOLERANCE) will be signaled.
C
C     3)  If START is greater than FINISH, the error
C         SPICE(BOUNDSOUTOFORDER) will be signaled.
C
C     4)  If a constant step is used, the step must be positive and
C         have magnitude large enough so that, when added to a value in
C         the range [START, FINISH], it will yield a distinct value.
C         Otherwise, the error SPICE(INVALIDCONSTSTEP) will be
C         signaled.
C
C     5)  If the step function is used and it returns a value less than
C         the preceding value, the error SPICE(INVALIDSTEP) will
C         be signaled.
C
C     6)  If the inner convergence loop fails to converge to TOL within
C         MXLOOP iterations, the error SPICE(NOCONVERGENCE) will be
C         signaled.
C
C     7)  If the POINTS array doesn't have enough room to store
C         the points associated with state transitions, the error 
C         SPICE(ARRAYTOOSMALL) will be signaled.
C
C     8)  If the result window doesn't have enough room to store
C         the abscissas associated with state transitions, the error 
C         will be signaled by a routine in the call tree of this
C         routine.
C
C$ Files
C
C     Kernels used by this routine are those needed by the input
C     routines
C
C        UDCOND
C        UDGETP
C        UDSTEP
C        UDREFN
C
C$ Particulars
C
C     This routine supports limb and terminator point detection.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     It is important that the user understand how the routines UDCOND,
C     UDSTEP and UDREFN are to be used and that the calling sequences
C     match precisely with the descriptions given here.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C     L.S. Elson     (JPL)
C     W.L. Taber     (JPL)
C     I.M. Underwood (JPL)
C     E.D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 30-JUN-2016 (NJB) (EDW)
C
C        Updated logic for placement of points in output POINTS
C        array. If the first element of RESULT is equal to BEGIN,
C        then space will be reserved in the first element of POINTS,
C        so as to keep the output points synced with the elements
C        of RESULT.
C
C        Updated short error messages.
C
C        Updated header I/O sections.
C
C        12-FEB-2016 (NJB) (EDW)
C
C        Derived from ZZGFSOLV Version 1.1.0, 24-OCT-2010 (EDW)
C
C-&

C$ Index_Entries
C
C     find tangent points on target
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      BRCKTD
      DOUBLE PRECISION      TOUCHD
      
      INTEGER               SIZED

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C
      CHARACTER*(256)       CONTXT

      DOUBLE PRECISION      BEGIN
      DOUBLE PRECISION      CURX
      DOUBLE PRECISION      MAXMAG
      DOUBLE PRECISION      PRVPNT ( 3 )
      DOUBLE PRECISION      SVDX
      DOUBLE PRECISION      T
      DOUBLE PRECISION      X1
      DOUBLE PRECISION      X2
      DOUBLE PRECISION      XSTEP
      DOUBLE PRECISION      TRNSTN
      DOUBLE PRECISION      XPOINT ( 3 )

      INTEGER               NLOOP
      INTEGER               ROOM 
      INTEGER               TO

      LOGICAL               CURSTA
      LOGICAL               INSTAT
      LOGICAL               S
      LOGICAL               STATE1
      LOGICAL               STATE2
      LOGICAL               SAVSTA

      LOGICAL               PRVSET


C
C     The maximum number of search loop iterations to execute.
C     The default refinement method is bisection, a very slow
C     method to convergence. Since 2**1000 ~ 10**301,
C     1000 loop iterations represents enough effort to assume
C     either the search will not converge or that the refinement
C     function operates slower than would bisection, in which
C     case the user should use the default GFREFN function.
C
      INTEGER               MXLOOP
      PARAMETER           ( MXLOOP = 1000 )

C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZTANSLV' )


C
C     Check the convergence tolerance.
C
      IF ( TOL .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Tolerance must be positive but was #.' )
         CALL ERRDP  ( '#',  TOL                               )
         CALL SIGERR ( 'SPICE(INVALIDTOLERANCE)'               )
         CALL CHKOUT ( 'ZZTANSLV'                              )
         RETURN

      END IF

C
C     Make sure that START is not greater than FINISH. Signal an
C     error for START > FINISH.
C
      IF ( START .GT. FINISH ) THEN

         CALL SETMSG ( 'Bad input interval: '     
     .   //            'START = # > FINISH = #.' )
         CALL ERRDP  ( '#', START                )
         CALL ERRDP  ( '#', FINISH               )
         CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)' )
         CALL CHKOUT ( 'ZZTANSLV'                )
         RETURN

      END IF

C
C     Make sure that TOL is not too small, i.e. that neither
C     START + TOL nor START - TOL equals START.
C
      IF (  ( TOUCHD (START - TOL) .EQ. START )
     .     .OR.
     .      ( TOUCHD (START + TOL) .EQ. START ) ) THEN

         CALL SETMSG ('TOL has value #1. '                      //
     .                'This value is too small to distinguish ' //
     .                'START - TOL or START + TOL from '        //
     .                'START, #2.'                               )
         CALL ERRDP  ( '#1', TOL                                 )
         CALL ERRDP  ( '#2', START                               )
         CALL SIGERR ( 'SPICE(INVALIDTOLERANCE)'                 )
         CALL CHKOUT ( 'ZZTANSLV'                                )
         RETURN

      END IF


C
C     Make sure that TOL is not too small, i.e. that neither
C     FINISH + TOL nor FINISH - TOL equals FINISH.
C
      IF (  ( TOUCHD (FINISH - TOL) .EQ. FINISH )
     .     .OR.
     .      ( TOUCHD (FINISH + TOL) .EQ. FINISH ) ) THEN

         CALL SETMSG ('TOL has value #1. '                      //
     .                'This value is too small to distinguish ' //
     .                'FINISH - TOL or FINISH + TOL from '      //
     .                'FINISH, #2.'                              )
         CALL ERRDP  ( '#1', TOL                                 )
         CALL ERRDP  ( '#2', FINISH                              )
         CALL SIGERR ( 'SPICE(INVALIDTOLERANCE)'                 )
         CALL CHKOUT ('ZZTANSLV'                                 )
         RETURN

      END IF

C
C     Make sure that STEP is not too small: it must be greater
C     than TOL.
C
      IF ( CSTEP ) THEN

         IF ( STEP .LE. 0.D0 ) THEN

            CALL SETMSG ( 'STEP has value #1. The search step '
     .      //            'must be positive.'                 )
            CALL ERRDP  ( '#1', STEP                          )
            CALL SIGERR ( 'SPICE(INVALIDCONSTSTEP)'           )
            CALL CHKOUT ('ZZTANSLV'                           )
            RETURN

         END IF

         MAXMAG = MAX( ABS(START), ABS(FINISH) )


         IF ( TOUCHD(MAXMAG+STEP) .EQ. MAXMAG ) THEN

            CALL SETMSG ( 'STEP has value #1. This value '              
     .      //            'is too small to guarantee that '
     .      //            'the search will advance.'        )
            CALL ERRDP  ( '#1', STEP                        )
            CALL SIGERR ( 'SPICE(INVALIDCONSTSTEP)'         )
            CALL CHKOUT ('ZZTANSLV'                         )
            RETURN

         END IF

      END IF

C
C     This algorithm determines those intervals when a given state is
C     observed to occur within a specified search interval.
C
C     Pairs of X values are recorded. The first member of each pair
C     denotes the X value at which the system changes to the state of
C     interest. The second denotes a transition out of that state.
C
C     If the state is .TRUE. at the beginning of the interval, the
C     beginning of the X interval will be recorded. This may or may not
C     be a transition point.
C
C     Similarly if the state is .TRUE. at the end of the interval, the
C     end of the interval will be recorded. Again, this may or may not
C     be a transition point.
C
C     Initially the current X value is the beginning of the search
C     interval.
C
      CURX   = START

      TO     = 1
      ROOM   = SIZED( RESULT )
      PRVSET = .FALSE.

C
C     Determine if the state at the current X value satisfies the
C     constraint.
C
      CALL UDCOND ( CURX, CURSTA, XPOINT )

      IF ( FAILED() ) THEN
         CALL CHKOUT (  'ZZTANSLV' )
         RETURN
      END IF

      IF ( CURSTA ) THEN
         CALL VEQU ( XPOINT, PRVPNT )
         PRVSET = .TRUE.
      END IF

C
C     If the system is in the state of interest, record the initial
C     X value of the search interval. The variable BEGIN will be
C     used to store the starting point of an interval over which
C     the state is .TRUE.
C
      IF ( CURSTA ) THEN

         INSTAT    = .TRUE.
         BEGIN     = CURX
         ENDFLG(1) = .FALSE.

C
C        BEGIN will be the first element of RESULT, presuming
C        a state transition is found later. We'll shift the 
C        pointer for the output point so the Ith point will
C        correspond to the Ith element of RESULT.
C
C        We don't have to check ROOM yet because we're not 
C        inserting anything into POINTS.
C
         TO   = TO   + 1
         ROOM = ROOM - 1

      ELSE

         INSTAT    = .FALSE.
         ENDFLG(1) = .TRUE.

      END IF

C
C     If the step size is constant, use the value supplied.
C
      IF ( CSTEP ) THEN
         XSTEP = STEP
      END IF

C
C     Save the current X value and state.
C
      SVDX   = CURX
      SAVSTA = CURSTA

C
C     Once initializations have been performed keep working
C     until the search interval has been exhausted.
C
C     While the last X value precedes the end of the interval:
C
      DO WHILE ( SVDX .LT. FINISH )
C
C        Attempt to bracket a state change.
C
C        Using the current window and internally stored information
C        about the current state, select a new current X.
C
         IF ( .NOT. CSTEP ) THEN

            CALL UDSTEP ( CURX, XSTEP )

            IF ( FAILED() ) THEN
               CALL CHKOUT (  'ZZTANSLV' )
               RETURN
            END IF

         END IF

C
C        Add the X step to the current X.  Make sure that the
C        X does not move beyond the end of the search interval.
C
         CURX = MIN ( TOUCHD(CURX + XSTEP), FINISH )

C
C        Compute the state at CURX.
C
         CALL UDCOND ( CURX, CURSTA, XPOINT )

         IF ( FAILED() ) THEN
            CALL CHKOUT (  'ZZTANSLV' )
            RETURN
         END IF

         IF ( CURSTA ) THEN
            CALL VEQU ( XPOINT, PRVPNT )
            PRVSET = .TRUE.
         END IF

C
C        While the state remains unchanged and the interval has not
C        been completely searched ...
C
         DO WHILE (       ( SAVSTA .EQV. CURSTA )
     .              .AND. ( SVDX   .LT.  FINISH )  )
C
C           Save the current X and state.
C
            SVDX   = CURX
            SAVSTA = CURSTA

C
C           Compute a new current X so that we will not step
C           past the end of the interval. 
C
            IF ( .NOT. CSTEP ) THEN

               CALL UDSTEP ( CURX, XSTEP )

               IF ( FAILED() ) THEN
                  CALL CHKOUT (  'ZZTANSLV' )
                  RETURN
               END IF

            END IF


            CURX = MIN ( TOUCHD(CURX + XSTEP), FINISH )
C
C           Compute the current state.
C
            CALL UDCOND ( CURX, CURSTA, XPOINT )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZTANSLV' )
               RETURN
            END IF

            IF ( CURSTA ) THEN
C
C              Save the associated vector for the X value CURX. In
C              normal usage, XPOINT is a surface intercept point.
C
               CALL VEQU ( XPOINT, PRVPNT )
               PRVSET = .TRUE.

            END IF
C
C           Loop back to see if the state has changed.
C
         END DO
C
C        At this point, SVDX and CURX are the X-values at the previous
C        and latest steps, respectively. SAVSTA and CURSTA are the
C        states at these X-values, respectively.
C
C        If we have detected a state change and not merely run out
C        of the search interval...
C
         IF ( SAVSTA .NEQV. CURSTA ) THEN
C
C           Call the previous state STATE1.
C           Call the current state STATE2.
C
C           Let X1 be the X value at state STATE1.
C           Let X2 be the X value at state STATE2.
C
C           Save the current X.
C
            STATE1 = SAVSTA
            STATE2 = CURSTA
            X1     = SVDX
            X2     = CURX

C
C           Make sure that X1 is not greater than X2. Signal an
C           error for X1 > X2.
C
            IF ( X1 .GT. X2 ) THEN

               CALL SETMSG ( 'Bad x interval result: X1 = # > X2 = #.' )
               CALL ERRDP  ( '#',  X1                                  )
               CALL ERRDP  ( '#',  X2                                  )
               CALL SIGERR ( 'SPICE(INVALIDSTEP)'                      )
               CALL CHKOUT ( 'ZZTANSLV'                                )
               RETURN

            END IF

C
C           Update the saved X and state values to those on the 
C           right side of the bracketing interval. We'll use these
C           values for the next bracketing step after a root is
C           found.
C
            SVDX   = CURX
            SAVSTA = CURSTA

C
C           X1 and X2 bracket the X value of transition. Squeeze this
C           interval down until it is less than some tolerance in
C           length. Do it as described below...
C
C           Loop while the difference between the X values X1 and X2
C           exceeds a specified tolerance.
C
            NLOOP = 0

            DO WHILE ( TOUCHD (X2 - X1) .GT. TOL )

               NLOOP = NLOOP + 1
C
C              This loop count error exists to catch pathologies
C              in the refinement function. The default bisection
C              refinement will converge before 1000 iterations if
C              a convergence is numerically possible. Any other
C              refinement function should require fewer iterations
C              compared to bisection. If not, the user should
C              probably use bisection.
C
               IF ( NLOOP .GE. MXLOOP ) THEN

                  CALL SETMSG ( 'Loop run exceeds maximum loop count. '
     .         //               'Unable to converge to TOL value #1 '
     .         //               'within MXLOOP value #2 iterations.')
                  CALL ERRDP  ( '#1', TOL              )
                  CALL ERRINT ( '#2', MXLOOP           )
                  CALL SIGERR ( 'SPICE(NOCONVERGENCE)' )
                  CALL CHKOUT ( 'ZZTANSLV'             )
                  RETURN

               END IF

C
C              Select an X value T, between X1 and X2 (possibly based
C              on the state values).
C
               CALL UDREFN ( X1, X2, STATE1, STATE2, T )

C
C              Check for an error signal. The default refinement
C              routine, GFREFN, does not include error checks.
C
               IF ( FAILED() ) THEN
                  CALL CHKOUT (  'ZZTANSLV' )
                  RETURN
               END IF

C
C              Check whether T is between X1 and X2.  If
C              not then assume that we have gone as far as
C              we can in refining our estimate of the transition
C              point. Set X1 and X2 equal to T.
C

               T = BRCKTD ( T, X1, X2 )

               IF ( T .EQ. X1 ) THEN
C
C                 This assignment may break the invariant that
C                 the state at X2 is STATE2. This is allowed
C                 because we'll exit the loop immediately.
C
                  X2 = T

               ELSE IF ( T .EQ. X2 ) THEN
C
C                 This assignment may break the invariant that
C                 the state at X1 is STATE1. This is allowed
C                 because we'll exit the loop immediately.

                  X1 = T

               ELSE
C
C                 Compute the state at X value T. If this state, S,
C                 equals STATE1, set X1 to T, otherwise set X2 to T.
C
                  CALL UDCOND ( T, S, XPOINT )

                  IF ( S ) THEN
C
C                    Save the latest point associated with a 
C                    .TRUE. state.
C
                     CALL VEQU ( XPOINT, PRVPNT )
                     PRVSET = .TRUE.

                  END IF

C
C                 Narrow the interval. Either increase X1 or decrease
C                 X2 by setting one of these endpoints to T. Maintain
C                 the invariant that the state is STATE1 at X1 and
C                 STATE2 at X2.
C
                  IF ( S .EQV. STATE1 ) THEN

                     X1 = T

                  ELSE

                     X2 = T

                  END IF

               END IF

            END DO

C
C           Let TRNSTN be the midpoint of [X1, X2]. Record this
C           abscissa value as marking the transition from STATE1 to
C           STATE2.
C
            TRNSTN = BRCKTD( (X1 + X2)*0.5D0, X1, X2 )

C
C           In state-of-interest or not? INSTAT indicates that STATE1
C           was .TRUE. We record intervals where the state is .TRUE.
C           when we detect the right hand endpoints of these intervals.
C
            IF ( INSTAT ) THEN
C
C              We were in the state of interest. TRNSTN marks the point
C              on the X-axis when the state changed to .FALSE. We need
C              to record the interval from BEGIN to FINISH and note
C              that the state has become .FALSE.
C
C              Add an interval starting at BEGIN and ending at TRNSTN
C              to the result window.
C
               CONTXT = 'Adding interval [BEGIN,TRNSTN] to RESULT. '
     .         //       'TRNSTN represents time of passage out of the '
     .         //       'state-of-interest.'

               CALL ZZWNINSD ( BEGIN, TRNSTN, CONTXT, RESULT )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZTANSLV' )
                  RETURN
               END IF

            ELSE
C
C              The previous state was .FALSE. As a result TRNSTN marks
C              the point where the state becomes .TRUE. Note that we
C              have transitioned to the state of interest and record
C              the X-value at which the transition occurred.
C
               BEGIN = TRNSTN

            END IF

C
C           A transition occurred either from from in-state to
C           out-of-state or the inverse. Reverse the value of the
C           INSTAT flag to signify the transition event.
C
            INSTAT = .NOT. INSTAT

C
C           For all state transitions, record the last point found
C           by the state function.
C           
            IF ( ROOM .GT. 0 ) THEN
C
C              Add the last point found during the transition search to
C              the POINTS array.
C
               IF ( PRVSET ) THEN

                  CALL VEQU ( PRVPNT, POINTS(1,TO) ) 

                  TO     = TO   + 1
                  ROOM   = ROOM - 1
                  PRVSET = .FALSE.

               ELSE

                  CALL SETMSG ( 'PRVPNT should always be set when a '
     .            //            'transition is detected. We found a '
     .            //            'transition at #, but PRVSET '
     .            //            'indicates we don''t have a previous '
     .            //            'point saved.'                        )
                  CALL ERRDP  ( '#', TRNSTN                           )
                  CALL SIGERR ( 'SPICE(BUG)'                          )
                  CALL CHKOUT ( 'ZZTANSLV'                            )
                  RETURN

               END IF


            ELSE
C
C              We ran out of room in the output point array. Note that
C              this error can occur before the result window insertion
C              fails, since that insertion takes place when the state
C              becomes .FALSE.
C
               CALL SETMSG ( 'Out of room in the POINTS array. Room '
     .         //            'is assumed to be adequate for SIZED('
     .         //            'RESULT) 3-vectors; this size is #.'    )
               CALL ERRINT ( '#', SIZED(RESULT)                      )
               CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                  )
               CALL CHKOUT ( 'ZZTANSLV'                              )
               RETURN

            END IF
C
C           That's it for this detection of state change.
C
         END IF
C
C        Continue if the search interval extends to the right
C        of the latest step. 
C
C        SVDX and SAVSTA are already set to the values at the
C        right side of the bracketing interval.
C
      END DO

C
C     Check if in-state at this abscissa value (FINISH). INSTAT is the
C     latest state value. If so record the interval.
C
      IF ( INSTAT ) THEN
C
C        The state is .TRUE. at FINISH.
C
C        Add an interval starting at BEGIN and ending at FINISH to the
C        window.
C
         CONTXT = 'Adding interval [BEGIN,FINISH] to RESULT. FINISH '
     .   //       'represents end of the search interval.'

         CALL ZZWNINSD ( BEGIN, FINISH, CONTXT, RESULT )

         ENDFLG(2) = .FALSE.
      ELSE         
         ENDFLG(2) = .TRUE.
      END IF
 
      CALL CHKOUT (  'ZZTANSLV' )
      RETURN
      END

