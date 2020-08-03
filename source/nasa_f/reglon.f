C$Procedure REGLON ( Regularize longitude intervals )

      SUBROUTINE REGLON ( NIVALS, BOUNDS, MAXN,   NOUT,
     .                    MINLON, MAXLON, OUTBDS, SRCS )

C$ Abstract
C
C     Regularize a set of longitude intervals. The output
C     intervals have their endpoints in ascending order,
C     and all output bounds are in the range MINLON:MAXLON.
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
C     DSKBRIEF
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsktol.inc'

      INTEGER               NIVALS
      DOUBLE PRECISION      BOUNDS ( 2, NIVALS )
      INTEGER               MAXN
      INTEGER               NOUT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      OUTBDS ( 2, * )
      INTEGER               SRCS   ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NIVALS     I   is the number of input longitude intervals.
C     BOUNDS     I   is an array of bounds of the input intervals.
C     MAXN       I   is the maximum number of output intervals.
C     NOUT       O   is the number of output intervals.
C     MINLON     O   is the minimum longitude of all output intervals.
C     MAXLON     O   is the maximum longitude of all output intervals.
C     OUTBDS     O   is an array of output longitude intervals.
C     SRCS       O   is an array mapping output to input intervals.
C
C$ Detailed_Input
C
C     NIVALS     is the number of input longitude intervals.
C
C     BOUNDS     is an array of upper and lower interval bounds.
C                Units are radians. 
C    
C                The elements
C
C                   BOUNDS(1,I)
C                   BOUNDS(2,I)
C
C                are, respectively, the lower and upper bounds of
C                the Ith interval. The upper bound may be less
C                than the lower bound. If the bounds are equal, they
C                are treated as though the upper bound exceeds the
C                lower bound by 2*pi.
C
C                 Bounds must be in the range
C
C                   [ -2*pi,  2*pi ]
C
C     MAXN       is the maximum number of output intervals that
C                can be returned.
C
C$ Detailed_Output
C
C     NOUT       is the number of output intervals. NOUT is greater
C                than or equal to NIVALS.
C
C     MINLON,
C     MAXLON     are, respectively, lower and upper bounds on 
C                the range of the output longitudes.
C
C                   {MINLON, MAXLON} 
C
C                is either 
C
C                   {0, 2*pi}  or {-pi, pi}
C
C     OUTBDS     is an array of output longitude bounds.
C                Units are radians. 
C
C                Each output interval represents a subset of some input
C                interval; the endpoints of an output interval may be
C                shifted by 2*pi relative to the endpoints of the
C                subset.
C                
C                Input intervals for which the lower bound is greater
C                than the upper bound are broken up into two output
C                intervals.
C
C                The elements
C
C                   OUTBDS(1,I)
C                   OUTBDS(2,I)
C
C                are, respectively, the lower and upper bounds of
C                the Ith interval. The upper bound is always 
C                greater than or equal to the lower bound.
C                Bounds are in the range
C
C                   [MINLON, MAXLON]
C
C
C     SRCS       is an array of indices that map the output intervals
C                to the source input intervals from which they were
C                derived.
C     
C$ Parameters
C
C     ANGMRG     See the description in dsktol.inc.
C
C$ Exceptions
C
C     1)  If the output array doesn't have enough room to store
C         the output intervals, the error SPICE(ARRAYTOOSMALL) 
C         is signaled.
C
C     2)  Longitudes outside the range
C
C            -2*pi - ANGMRG : 2*pi + ANGMRG 
C
C         are not accepted: if such a value is encountered, the
C         error will be diagnosed by a routine in the call tree 
C         of this routine.
C
C     3)  If the lower bound of a longitude interval matches the 
C         upper bound, the error SPICE(ZEROBOUNDSEXTENT) is
C         signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine "regularizes" a set of longitude intervals: it maps
C     them to a set of output longitude intervals that has the same
C     coverage on the unit circle as the input set, such that each of
C     the output intervals has its endpoints in increasing order, and
C     all of the output intervals have their endpoints in a common
C     interval of length 2*pi.
C
C     The set of output intervals has the property that order
C     relationships between endpoints are valid.
C
C     This routine supports coverage gap determination.
C 
C$ Examples
C
C     See usage in DSKBRIEF.
C
C$ Restrictions
C
C     1) For use only within program DSKBRIEF.
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
C     DSKBRIEF Version 1.0.0, 21-FEB-2017 (NJB)
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPMAX
      DOUBLE PRECISION      DPMIN
      DOUBLE PRECISION      DPR
      DOUBLE PRECISION      PI
      DOUBLE PRECISION      TWOPI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      LB
      DOUBLE PRECISION      LOCLB
      DOUBLE PRECISION      LOCUB
      DOUBLE PRECISION      UB

      INTEGER               I
      INTEGER               NREQ
      

      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'REGLON' )

C
C     No output intervals have been found yet.
C
      NOUT = 0

      IF ( NIVALS .EQ. 0 ) THEN

         CALL CHKOUT ( 'REGLON')
         RETURN

      END IF

C
C     Get lower and upper bounds of input values. 
C
      MINLON = DPMAX()
      MAXLON = DPMIN()

      DO I = 1, NIVALS   
      
         LB = BOUNDS(1,I)
         UB = BOUNDS(2,I)
C
C        Rectangles of zero longitude extent not allowed.
C
         IF ( LB .EQ. UB ) THEN

            CALL SETMSG ( 'Longitude lower bound # (# degrees) '
     .      //            'equals upper bound.'                  )
            CALL ERRDP  ( '#', LB                                )
            CALL ERRDP  ( '#', LB * DPR()                        )
            CALL SIGERR ( 'SPICE(ZEROBOUNDSEXTENT)'              )
            CALL CHKOUT ( 'REGLON'                               )
            RETURN

         END IF 

C  
C        Adjust UB if necessary before deciding on the output
C        range.
C
         IF ( UB .LT. LB ) THEN
            UB = UB + TWOPI()
         END IF

         MINLON = MIN( LB, UB, MINLON )
         MAXLON = MAX( LB, UB, MAXLON )

      END DO
      
C
C     If MAXLON and MINLON lie within the range 
C
C        0 - ANGMRG : 2*pi + ANGMRG
C
C     we'll set the output longitudes to lie in the range
C
C        0 : 2*pi
C
C
C     If MAXLON and MINLON lie within the range 
C
C        -pi - ANGMRG : pi + ANGMRG
C
C     we'll set the output longitudes to lie in the range
C
C        -pi : pi
C
C
C     We use the latter range if neither of the first two
C     conditions are met.
C
C     
      IF (      ( MINLON .GT.  -ANGMRG              )
     .    .AND. ( MAXLON .LT. ( TWOPI() + ANGMRG )  )  ) THEN

         A = 0.D0
         B = TWOPI()

      ELSE
C
C        We arbitrarily pick the output longitude range
C
C           -pi : pi
C
         A = -PI()
         B =  PI()

      END IF

C
C     Set the output values of MINLON and MAXLON.
C     
      MINLON = A
      MAXLON = B

C
C     Process each input interval.
C
      DO I = 1, NIVALS

         LB = BOUNDS(1,I)
         UB = BOUNDS(2,I)
C
C        We'll adjust the inputs to ensure they're in range.
C
C        First, make sure we're starting with values in 
C        the range [-2*pi, 2*pi]. 
C
         CALL ZZNRMLON ( LB, UB, ANGMRG, LOCLB, LOCUB )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'REGLON' )
            RETURN
         END IF

C
C        Move each output into the range [A, B].
C
         IF ( LOCLB .LT. A ) THEN

            LOCLB = LOCLB + TWOPI()

         ELSE IF ( LOCLB .GT. B ) THEN

            LOCLB = LOCLB - TWOPI()

         END IF


         IF ( LOCUB .LT. A ) THEN

            LOCUB = LOCUB + TWOPI()

         ELSE IF ( LOCUB .GT. B ) THEN

            LOCUB = LOCUB - TWOPI()

         END IF

C
C        Now the bounds are in range, but they may be
C        out of order.
C           
         IF ( LOCLB .LT. LOCUB ) THEN
C
C           The bounds are in order. Add the interval to
C           the list of output intervals.
C
            NOUT = NOUT + 1

            IF ( NOUT .GT. MAXN ) THEN
C
C              We're out of room.
C
               CALL SETMSG ( 'Output arrays have room for # '
     .         //            'intervals we have found # '
     .         //            'output intervals so far.'      )
               CALL ERRINT ( '#', MAXN                       )
               CALL ERRINT ( '#', NOUT                       )
               CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'          )
               CALL CHKOUT ( 'REGLON'                        )
               RETURN
   
            END IF

            OUTBDS(1,NOUT) = LOCLB
            OUTBDS(2,NOUT) = LOCUB
            SRCS  (  NOUT) = I


         ELSE
C
C           The bounds are in range but out of order.
C           We'll split the input interval into two
C           output intervals.
C
            NREQ = 0

            IF ( A .LT. LOCUB ) THEN
               NREQ = NREQ + 1
            END IF
               
            IF ( LOCLB .LT. B ) THEN
               NREQ = NREQ + 1
            END IF

            
            IF ( NOUT+NREQ .GT. MAXN ) THEN
C
C              We're out of room.
C
               CALL SETMSG ( 'Output arrays have room for # '
     .         //            'intervals we have found # '
     .         //            'output intervals so far.'      )
               CALL ERRINT ( '#', MAXN                       )
               CALL ERRINT ( '#', NOUT+NREQ                  )
               CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'          )
               CALL CHKOUT ( 'REGLON'                        )
               RETURN
 
            END IF
C
C           The input interval "wraps around" the output boundaries.
C
C           The output intervals extend from A to the upper bound and
C           from the lower bound to B.
C
            IF ( A .LT. LOCUB ) THEN

               NOUT = NOUT + 1
 
               OUTBDS(1,NOUT) = A
               OUTBDS(2,NOUT) = LOCUB
               SRCS  (  NOUT) = I

            END IF

            IF ( LOCLB .LT. B ) THEN

               NOUT = NOUT + 1
 
               OUTBDS(1,NOUT) = LOCLB
               OUTBDS(2,NOUT) = B
               SRCS  (  NOUT) = I

            END IF

         END IF
 
      END DO

      CALL CHKOUT ( 'REGLON' )
      RETURN
      END 



