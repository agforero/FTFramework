C$Procedure ZZNRMLON ( Normalize longitude bounds )

      SUBROUTINE ZZNRMLON ( INMIN, INMAX, TOL, OUTMIN, OUTMAX )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Normalize longitude bounds: map bounds into the interval
C
C        [ -2*pi, 2*pi ]
C
C     Put the bounds in order and ensure that the bounds differ by no
C     more than 2*pi.
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
C     GEOMETRY
C     LONGITUDE
C     MATH
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      INMIN
      DOUBLE PRECISION      INMAX
      DOUBLE PRECISION      TOL
      DOUBLE PRECISION      OUTMIN
      DOUBLE PRECISION      OUTMAX

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INMIN      I   Longitude interval lower bound.
C     INMAX      I   Longitude interval upper bound.
C     TOL        I   Round-off tolerance.
C     OUTMIN     O   Adjusted longitude interval lower bound.
C     OUTMAX     O   Adjusted longitude interval upper bound.
C
C$ Detailed_Input
C
C     INMIN,
C     INMAX      are, respectively, the lower and upper
C                bounds of a longitude interval. Units are
C                radians. INMIN and INMAX must lie in the range
C
C                   [-2*pi, 2pi] 
C
C                INMAX is allowed to be less than INMIN, but must
C                not be equal to it.
C                
C
C     TOL        is a non-negative tolerance value. If an input bound
C                lies outside of the range 
C
C                   [-2*pi, 2pi] 
C
C                by less than TOL, it is interpreted as being equal
C                to the nearest interval endpoint.
C
C                If INMAX exceeds INMIN by less than TOL, the bounds
C                are interpreted as being 2*pi radians apart.
C    
C
C$ Detailed_Output
C
C     OUTMIN,
C     OUTMAX     are, respectively, the normalized lower and upper
C                bounds of the longitude interval described by the
C                inputs. Units are radians.
C
C                "Normalization" means the bounds are modified if
C                necessary so that they represent the same interval
C                as [INMIN, INMAX], but also meet the following
C                criteria:
C
C                   1) Both OUTMIN and OUTMAX lie in the interval
C                      [-2*pi, 2*pi].
C
C                   2) OUTMIN is strictly less than OUTMAX.
C
C                   3) OUTMAX does not exceed OUTMIN by more than 2*pi.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If TOL is negative, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     2)  If either input longitude is less than 2*pi-TOL or
C         greater than 2*pi + TOL, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     3)  If INMAX equals INMIN, or if INMAX is less than INMIN
C         by an integer multiple of 2*pi, the error 
C         SPICE(ZEROBOUNDSEXTENT) is signaled.
C
C$ Files
C
C     None.
C     
C$ Particulars
C
C     This routine centralizes an oft-repeated algorithm. It is
C     called by several DSK routines.
C
C$ Examples
C
C     See usage in ZZINLAT.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by SPICE
C     routines.
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
C-    SPICELIB Version 1.0.0, 11-OCT-2016 (NJB) 
C
C        Original version 16-MAY-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     normalize longitude interval
C     normalize longitude boundaries
C     
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPR
      DOUBLE PRECISION      TOUCHD
      DOUBLE PRECISION      TWOPI

C
C     Local variables
C
      DOUBLE PRECISION      DELTA
      DOUBLE PRECISION      PI2

      LOGICAL               FIRST

C
C     Saved values
C
      SAVE                  FIRST
      SAVE                  PI2

C
C     Initial values
C
      DATA                  FIRST / .TRUE. /

C
C     Use discovery check-in. Don't check RETURN.
C
      IF ( FIRST ) THEN

         PI2   = TWOPI()
         FIRST = .FALSE.

      END IF

C
C     TOL cannot be negative.
C
      IF ( TOL .LT. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZNRMLON'                                  )
         CALL SETMSG ( 'Tolerance must be non-negative but was #.' )
         CALL ERRDP  ( '#', TOL                                    )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                    )
         CALL CHKOUT ( 'ZZNRMLON'                                  )
         RETURN

      END IF

C
C     Reject inputs that lie outside of [-2*pi, 2*pi], accounting
C     for the tolerance.
C
      IF (      ( INMIN .LT. -PI2-TOL )
     .     .OR. ( INMIN .GT.  PI2+TOL )  ) THEN

         CALL CHKIN  ( 'ZZNRMLON'                                  )
         CALL SETMSG ( 'Longitude lower bound INMIN = # (radians), '
     .   //            ' = # (deg). The minimum allowed value is '
     .   //            ' -2*pi - TOL = # (radians), = # (deg).'    )
         CALL ERRDP  ( '#', INMIN                                  )
         CALL ERRDP  ( '#', INMIN * DPR()                          )
         CALL ERRDP  ( '#',   -PI2  - TOL                          )
         CALL ERRDP  ( '#', ( -PI2  - TOL ) * DPR()                )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                    )
         CALL CHKOUT ( 'ZZNRMLON'                                  )
         RETURN

      END IF

C
C     The input bounds may not be equal.
C
      IF ( INMIN .EQ. INMAX ) THEN

         CALL CHKIN  ( 'ZZNRMLON'                                  )
         CALL SETMSG ( 'Longitude lower bound INMIN = # (radians), '
     .   //            ' = # (deg), is equal to upper bound.'      )
         CALL ERRDP  ( '#', INMIN                                  )
         CALL ERRDP  ( '#', INMIN * DPR()                          )
         CALL SIGERR ( 'SPICE(ZEROBOUNDSEXTENT)'                   )
         CALL CHKOUT ( 'ZZNRMLON'                                  )
         RETURN

      END IF

C
C     The input longitude is within range or is out of range by at most
C     |TOL| radians. Bracket it.

      OUTMIN = MAX( -PI2,  MIN( INMIN, PI2 ) )

C
C     Same deal for the upper bound.
C
      IF (      ( INMAX .LT. -PI2-TOL )
     .     .OR. ( INMAX .GT.  PI2+TOL )  ) THEN

         CALL CHKIN  ( 'ZZNRMLON'                                  )
         CALL SETMSG ( 'Longitude upper bound INMAX = # (radians), '
     .   //            ' = # (deg). The minimum allowed value is '
     .   //            ' -2*pi - TOL = # (radians), = # (deg).'    )
         CALL ERRDP  ( '#', INMAX                                  )
         CALL ERRDP  ( '#', INMAX * DPR()                          )
         CALL ERRDP  ( '#',   -PI2  - TOL                          )
         CALL ERRDP  ( '#', ( -PI2  - TOL ) * DPR()                )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                    )
         CALL CHKOUT ( 'ZZNRMLON'                                  )
         RETURN

      END IF

      OUTMAX = MAX( -PI2,  MIN( INMAX, PI2 ) )

C
C     If the bounds are out of order, put them in order. 
C     It is assumed that no interval has length zero.
C
C     If the upper bound is greater than the lower bound by
C     less than TOL, the bounds are considered to be "out of
C     order."
C
      IF ( OUTMAX .LE. TOUCHD(OUTMIN + TOL) ) THEN
C
C        Shift one of the bounds by 2*pi, while keeping
C        the bounds in the range [-2pi, 2pi].
C
         IF ( OUTMAX .LE. 0.D0 ) THEN
C
C           OUTMAX is non-positive. Shift it to the right.
C            
            OUTMAX = MIN ( TOUCHD(OUTMAX + PI2),  PI2 )

            IF ( OUTMAX .LT. OUTMIN ) THEN
C
C              If the bounds are still out of order, shift the lower
C              bound left.
C
               OUTMIN = MAX ( TOUCHD(OUTMIN - PI2), -PI2 )

            END IF

         ELSE
C
C           OUTMAX is > 0. Shift the lower bound left.
C
            OUTMIN = MAX ( TOUCHD(OUTMIN - PI2), -PI2 )

         END IF

      END IF

C
C     If the bounds are too far apart, move them together. Note
C     that OUTMIN and OUTMAX are already set at this point.
C
      DELTA = TOUCHD( OUTMAX - OUTMIN )

      IF (  DELTA  .GT.  TOUCHD( PI2 + TOL )  ) THEN
C
C        Shift the upper bound lower by 2*pi.
C
         OUTMAX = TOUCHD( OUTMAX - PI2 )

      END IF


C
C     The output bounds must not be equal. We could end up with
C     equal output bounds if the input maximum is less than
C     the input minimum and the bounds differ by an integer 
C     multiple of 2*pi.
C
      IF ( OUTMIN .EQ. OUTMAX ) THEN

         CALL CHKIN  ( 'ZZNRMLON'                                     )
         CALL SETMSG ( 'After adjustment, input longitude lower '
     .   //            'bound INMIN = # (radians),  = # (deg), '
     .   //            'is equal to adjusted longitude upper bound. '
     .   //            'Input upper bound = # (radians),  = # (deg). '
     .   //            'When the input upper bound is less than the '
     .   //            'input lower bound, the difference must not '
     .   //            'be an integer multiple of 2*pi.'              )
         CALL ERRDP  ( '#', INMIN                                     )
         CALL ERRDP  ( '#', INMIN * DPR()                             )
         CALL ERRDP  ( '#', INMAX                                     )
         CALL ERRDP  ( '#', INMAX * DPR()                             )
         CALL SIGERR ( 'SPICE(ZEROBOUNDSEXTENT)'                      )
         CALL CHKOUT ( 'ZZNRMLON'                                     )
         RETURN

      END IF

      END
