C$Procedure  ZZCXBRUT ( Cone-segment intersection by brute force )
 
      SUBROUTINE ZZCXBRUT ( APEX,   AXIS,   ANGLE, 
     .                      ENDPT1, ENDPT2, XPT,   ISBRCK )
 
C$ Abstract
C
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due
C     to the volatile nature of this routine.
C
C     Compute a bracketed point of intersection of a specified nappe of
C     a cone and a line segment by "brute force"---specifically by
C     bisection.
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
C     CONE
C     GEOMETRY
C     INTERSECTION
C     LINE
C     MATH
C     SEGMENT
C
C$ Declarations
 
      IMPLICIT NONE
      
      DOUBLE PRECISION      APEX   ( 3 )
      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      ANGLE
      DOUBLE PRECISION      ENDPT1 ( 3 )
      DOUBLE PRECISION      ENDPT2 ( 3 )
      DOUBLE PRECISION      XPT    ( 3 )
      LOGICAL               ISBRCK
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     APEX       I   Apex of cone.
C     AXIS       I   Axis of cone.
C     ANGLE      I   Angle of cone.
C     ENDPT1,
C     ENDPT2     I   Endpoints of line segment.
C     XPT        O   Intersection point, if it exists.
C     ISBRCK     O   Logical flag indicating whether root is bracketed.
C
C$ Detailed_Input
C
C     APEX       is the apex (tip) of the cone. In this routine's
C                documentation, we'll consider the cone to be a
C                semi-infinite pyramid with circular cross-section. In
C                some contexts, this object is called one "nappe" of
C                the complete cone.
C
C     AXIS       is an axis vector of the cone.                  
C
C     ANGLE      is the angular separation from AXIS of the rays
C                comprising the cone. Let the notation
C
C                   < A, B >
C 
C                denote the dot product of vectors A and B, and let
C
C                   ||A||
C
C                denote the norm of vector A. Then the cone is the set
C                of points
C
C                             X-APEX       AXIS
C                   { X:  < ----------,  -------- >  =  cos(ANGLE) }
C                           ||X-APEX||   ||AXIS||
C
C
C     ENDPT1,
C     ENDPT2     are endpoints of a line segment. These points 
C                must be distinct.
C
C                Exactly one of ENDPT1 and ENDPT2 must be inside
C                the cone. 
C          
C$ Detailed_Output
C
C     NXPTS      is the number of points of intersection of the input
C                line segment and cone.
C
C     XPT        is the unique point of intersection of the segment and
C                cone that lies on the line segment connecting ENDPT1
C                and ENDPT2.
C
C     ISBRCK     is a logical flag that is set to .TRUE. if and only if
C                ENDPT1 and ENDPT2 bracket a root. Equivalently, the
C                endpoints bracket a root when exactly one of ENDPT1
C                and ENDPT2 is inside the cone.
C                 
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If AXIS is the zero vector, the error SPICE(ZEROVECTOR)
C         will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is used by the SPICELIB routine INCNSG to handle
C     cases where solution of a quadratic equation is subject to
C     excessive loss of precision.
C
C$ Examples
C
C     See usage in INCNSG.
C
C$ Restrictions
C
C     This is a private SPICELIB routine.
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
C-    SPICELIB Version 1.0.0 30-SEP-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     brute force intersection of line segment and cone
C
C-&







C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      PI
      DOUBLE PRECISION      VDOT

      LOGICAL               RETURN
      LOGICAL               STATE1
      LOGICAL               STATE2
      LOGICAL               VZERO

C
C     Local parameters
C
      DOUBLE PRECISION      CNVLIM
      PARAMETER           ( CNVLIM = 1.D-15 )

      INTEGER               MAXITR
      PARAMETER           ( MAXITR = 1000 )

C
C     Local variables
C 
      DOUBLE PRECISION      COSANG
      DOUBLE PRECISION      DELTA
      DOUBLE PRECISION      DP
      DOUBLE PRECISION      DP1
      DOUBLE PRECISION      DP2
      DOUBLE PRECISION      HIGH
      DOUBLE PRECISION      LOCANG
      DOUBLE PRECISION      LOCAXI ( 3 )
      DOUBLE PRECISION      LOW
      DOUBLE PRECISION      MIDPT  
      DOUBLE PRECISION      OFF1   ( 3 )
      DOUBLE PRECISION      OFF2   ( 3 )
      DOUBLE PRECISION      PRVDLT
      DOUBLE PRECISION      SEG    ( 3 )
      DOUBLE PRECISION      UOFF1  ( 3 )
      DOUBLE PRECISION      UOFF2  ( 3 )
      DOUBLE PRECISION      UX     ( 3 )
      DOUBLE PRECISION      X      ( 3 )

      INTEGER               NITR

      LOGICAL               STATE

C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Check the axis.
C
      IF ( VZERO(AXIS) ) THEN

         CALL CHKIN  ( 'ZZCXBRUT'                     )
         CALL SETMSG ( 'Cone axis is the zero vector' )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'            )
         CALL CHKOUT ( 'ZZCXBRUT'                     )
         RETURN

      END IF

C
C     Make a local version of the cone's axis and angle. The
C     angle will be less than or equal to pi/2 radians.
C
      IF ( ANGLE .GT. HALFPI() ) THEN
         
         LOCANG = PI() - ANGLE
         CALL VMINUS ( AXIS, LOCAXI )
      ELSE
         LOCANG = ANGLE
         CALL VEQU   ( AXIS, LOCAXI )
      END IF

      CALL VHATIP ( LOCAXI )

      COSANG = COS( LOCANG )

C
C     Calculate the offsets of the endpoints from the apex,
C     and get unit-length versions of these.
C     
      CALL VSUB ( ENDPT1, APEX, OFF1 )
      CALL VSUB ( ENDPT2, APEX, OFF2 )

      CALL VHAT ( OFF1, UOFF1 )
      CALL VHAT ( OFF2, UOFF2 )

C
C     Get the dot products of the unit offsets with the axis.
C     These will serve as proxies for latitude.
C
      DP1 = VDOT( UOFF1, LOCAXI )
      DP2 = VDOT( UOFF2, LOCAXI )

C
C     The "state" variables at the endpoints are .TRUE. if
C     the endpoints are on or inside the cone. 
C
      STATE1 = DP1 .GE. COSANG
      STATE2 = DP2 .GE. COSANG

C
C     The intersection is supposed to be bracketed. Return
C     if not, indicating the situation via ISBRCK.
C
      ISBRCK = STATE1 .NEQV. STATE2

      IF ( .NOT. ISBRCK) THEN
         RETURN
      END IF
 
C
C     Prepare for a solution by bisection.
C
      CALL VSUB ( OFF2, OFF1, SEG )

      LOW    = 0.D0
      HIGH   = 1.D0
      DELTA  = ABS( HIGH - LOW )
      PRVDLT = 2.D0
      NITR   = 0

      DO WHILE (       ( DELTA .GT. CNVLIM )
     .           .AND. ( DELTA .LT. PRVDLT ) 
     .           .AND. ( NITR  .LT. MAXITR )  )

         MIDPT = ( LOW + HIGH ) / 2
         
         CALL VLCOM ( 1.D0, OFF1, MIDPT, SEG, X )

         CALL VHAT ( X, UX )

         DP    = VDOT ( UX, LOCAXI )

         STATE = DP .GE. COSANG

         IF ( STATE .EQV. STATE1 ) THEN
C
C           There has been no state change between OFF1
C           and XPT.
C            
            LOW  = MIDPT
         ELSE
            HIGH = MIDPT
         END IF

         PRVDLT = DELTA
         DELTA  = ABS( HIGH - LOW )         

         NITR   = NITR + 1

      END DO
C
C     X is an offset from APEX. The solution is an offset from the
C     origin.
C
      CALL VADD ( APEX, X, XPT )

      RETURN
      END



