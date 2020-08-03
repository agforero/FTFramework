C$Procedure ZZEDTMPT ( Ellipsoid terminator point in half-plane )
 
      SUBROUTINE ZZEDTMPT ( UMBRAL, A, B, C, R, AXIS, PLNVEC, POINT )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Compute an umbral or penumbral terminator point on an ellipsoidal
C     target. The point is confined to a specified half-plane. The 
C     illumination source is spherical.
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
C     ILLUMINATION
C     TERMINATOR
C
C$ Declarations
 
      IMPLICIT NONE

      LOGICAL               UMBRAL
      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      C
      DOUBLE PRECISION      R
      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      PLNVEC ( 3 )
      DOUBLE PRECISION      POINT  ( 3 )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     UMBRAL     I   Flag indicating whether terminator is umbral.
C     A          I   Semi-axis length in the X direction.
C     B          I   Semi-axis length in the Y direction.
C     C          I   Semi-axis length in the Z direction.
C     R          I   Radius of illumination source.
C     AXIS       I   Axis vector from target to source.
C     PLNVEC     I   Vector in cutting half-plane.
C     POINT      O   Terminator point.
C
C$ Detailed_Input
C
C     UMBRAL     is a logical flag that must be set by the caller
C                to .TRUE. if an umbral terminator point is to be
C                computed, and set to .FALSE. if a penumbral 
C                terminator point is to be computed.
C
C     A,
C     B,
C     C          are, respectively, the semi-axis lengths of a 
C                triaxial ellipsoid representing a target object.
C                The axes of the ellipsoid are aligned with the 
C                axes of the Cartesian coordinate system.
C
C
C     R          is the radius of the spherical model of the
C                illumination source. 
C
C
C     AXIS       is the position of the center of the illumination
C                source relative to the center of the target. AXIS is
C                contained in the line bounding cutting half-planes in
C                which terminator points are found.
C
C
C     PLNVEC     is a vector that, together with AXIS, defines
C                a cutting half-plane. The half-plane contains 
C                PLNVEC and has the line containing AXIS as a
C                boundary.
C
C                The terminator point that is sought lies within
C                the half-plane.
C
C$ Detailed_Output
C
C     POINT      is a terminator point lying within the half-plane
C                defined by AXIS and PLNVEC. 
C
C                When UMBRAL is set to .TRUE., POINT lies on the 
C                umbral terminator: the plane tangent to the target
C                ellipsoid at POINT is also tangent to the illumination
C                source. This plane does not intersect the vector
C                AXIS.
C
C                When UMBRAL is set to .FALSE., POINT lies on the
C                penumbral terminator: the plane tangent to the target
C                ellipsoid at POINT is also tangent to the illumination
C                source. This plane intersects the vector AXIS.
C
C                Unless the target is spherical, the plane tangent
C                to the target at POINT is usually tangent to the
C                illumination source at a point outside the plane
C                containing the cutting half-plane.
C                
C$ Parameters
C
C     None.
C
C$ Exceptions
C 
C     1)  If any of the target ellipsoid's semi-axis lengths is
C         non-positive, the error SPICE(INVALIDAXISLENGTH) is signaled.
C
C     2)  If the radius of the illumination source is non-positive, the
C         error SPICE(INVALIDRADIUS) is signaled.
C
C     3)  If AXIS is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C
C     4)  If the target and illumination source are separated by 
C         less than the sum of the radius of the source sphere 
C         and the maximum radius of the target ellipsoid, the
C         error SPICE(OBJECTSTOOCLOSE) is signaled.
C
C     5)  If PLNVEC is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C
C     6)  If AXIS and PLNVEC are linearly dependent, the error 
C         SPICE(DEGENERATECASE) is signaled.
C
C     7)  If the terminator solution doesn't converge, the error
C         SPICE(NOCONVERGENCE) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The "umbral" terminator on an ellipsoid is the boundary of the
C     region that is in total shadow. At any point P on this boundary,
C     the illumination source is below and tangent to the local
C     horizontal plane at P.
C
C     The "penumbral" terminator on an ellipsoid is the boundary of the
C     region that totally illuminated. At any point P on this boundary,
C     the illumination source is above and tangent to the local
C     horizontal plane at P.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     This is a SPICELIB private routine.
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
C-    SPICELIB Version 2.0.0, 30-JUN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     compute ellipsoid terminator point in half-plane
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DASINE
      DOUBLE PRECISION      HALFPI
      DOUBLE PRECISION      TOUCHD
      DOUBLE PRECISION      VDOT
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VNORM

      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local parameters
C
C
C     Tolerance used for arcsine arguments:
C
      DOUBLE PRECISION      ATOL
      PARAMETER           ( ATOL = 1.D-14 )

C
C     Angular error used to determine convergence:
C
      DOUBLE PRECISION      CNVLIM
      PARAMETER           ( CNVLIM = 1.D-15 )

C
C     Maximum number of iterations allowed for root finding:
C
      INTEGER               MAXITR
      PARAMETER           ( MAXITR = 20 )

C
C     Local variables
C
      DOUBLE PRECISION      ANGERR
      DOUBLE PRECISION      ANGLE
      DOUBLE PRECISION      CONST
      DOUBLE PRECISION      D
      DOUBLE PRECISION      H
      DOUBLE PRECISION      HPLNML ( 3 )
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      NORMAL ( 3 )
      DOUBLE PRECISION      PROJ   ( 3 )
      DOUBLE PRECISION      S
      DOUBLE PRECISION      SGNNML ( 3 )
      DOUBLE PRECISION      SRCPNT ( 3 )
      DOUBLE PRECISION      TA
      DOUBLE PRECISION      TARGPT ( 3 )
      DOUBLE PRECISION      TAXIS  ( 3 )
      DOUBLE PRECISION      TB
      DOUBLE PRECISION      TC
      DOUBLE PRECISION      THETA
      DOUBLE PRECISION      TMPVEC ( 3 )
      DOUBLE PRECISION      TPLNVC ( 3 )
      DOUBLE PRECISION      TRANS  ( 3, 3 )
      DOUBLE PRECISION      UTAXIS ( 3 )
      DOUBLE PRECISION      XA
      DOUBLE PRECISION      XB
      DOUBLE PRECISION      XC

      INTEGER               NITR


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZEDTMPT' )

C
C     Check A, B, C, and R.
C
      IF (      ( A .LE. 0.D0 )
     .     .OR. ( B .LE. 0.D0 )
     .     .OR. ( C .LE. 0.D0 )  ) THEN

         CALL SETMSG ( 'Target radii must be strictly positive '
     .   //            'but were #, #, #.'                      )
         CALL ERRDP  ( '#',  A                                  )
         CALL ERRDP  ( '#',  B                                  )
         CALL ERRDP  ( '#',  C                                  )
         CALL SIGERR ( 'SPICE(INVALIDAXISLENGTH)'               )
         CALL CHKOUT ( 'ZZEDTMPT'                               )
         RETURN

      END IF

      IF ( R .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Source radius must be strictly positive '
     .   //            'but was #.'                             )
         CALL ERRDP  ( '#',  R                                  )
         CALL SIGERR ( 'SPICE(INVALIDRADIUS)'                   )
         CALL CHKOUT ( 'ZZEDTMPT'                               )
         RETURN

      END IF

C
C     Check AXIS and PLNVEC.
C
      IF ( VZERO( AXIS ) ) THEN

         CALL SETMSG ( 'AXIS must be a non-zero vector but is '
     .   //            'in fact zero.'                         )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                     )
         CALL CHKOUT ( 'ZZEDTMPT'                              )
         RETURN

      END IF

      IF (  R + MAX(A,B,C)  .GE.  VNORM(AXIS)  ) THEN

         CALL SETMSG ( 'Centers of source and target are too '
     .   //            'close together; distance is #. Radius '
     .   //            'of source is #; semi-axis lengths are '
     .   //            '#, #, #.'                               )
         CALL ERRDP  ( '#',  VNORM(AXIS)                        )
         CALL ERRDP  ( '#',  R                                  )
         CALL ERRDP  ( '#',  A                                  )
         CALL ERRDP  ( '#',  B                                  )
         CALL ERRDP  ( '#',  C                                  )
         CALL SIGERR ( 'SPICE(OBJECTSTOOCLOSE)'                 )
         CALL CHKOUT ( 'ZZEDTMPT'                               )
         RETURN
         
      END IF


      IF ( VZERO( PLNVEC ) ) THEN

         CALL SETMSG ( 'PLNVEC must be a non-zero vector but is '
     .   //            'in fact zero.'                           )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                       )
         CALL CHKOUT ( 'ZZEDTMPT'                                )
         RETURN

      END IF

C
C     Transform the source, target, axis, and plane vector
C     so that the target becomes a unit sphere.
C
      CALL CLEARD ( 9, TRANS )

      TA         = 1.D0 / A 
      TB         = 1.D0 / B 
      TC         = 1.D0 / C 

      XA         = TA * R
      XB         = TB * R
      XC         = TC * R

      TRANS(1,1) = TA
      TRANS(2,2) = TB
      TRANS(3,3) = TC
      
C
C     TNEGAX is the negative of the transformed axis.
C     UTAXIS is the unit vector in the direction of TNEGAX.
C
      CALL MXV    ( TRANS,  PLNVEC, TPLNVC )
      CALL MXV    ( TRANS,  AXIS,   TAXIS )
      CALL VHAT   ( TAXIS,  UTAXIS )

C
C     Let HPLNML be a normal vector to the plane containing
C     the transformed axis and plane vectors.
C
      CALL VCRSS ( TPLNVC, TAXIS, HPLNML )

      IF ( VZERO(HPLNML) ) THEN

         CALL SETMSG ( 'Plane reference vector and axis are '
     .   //            'linearly dependent.'                 )
         CALL SIGERR ( 'SPICE(DEGENERATECASE)'               )
         CALL CHKOUT ( 'ZZEDTMPT'                            )
         RETURN

      END IF

C
C     Let MAXR be an outer bounding radius for the transformed
C     source sphere.
C
      MAXR = MAX ( XA, XB, XC ) 

      D    = VNORM ( TAXIS )

      IF ( UMBRAL ) THEN
C
C        Find the angle between the negative axis and a ray tangent to
C        both the transformed target and the outer bounding sphere of
C        the transformed source. Here a tangent point on the
C        transformed target is the vertex, and the tangent ray is
C        confined to the half-plane normal to HPLNML, containing
C        TPLNVC, and bounded by the line containing TNEGAX. The
C        tangent ray does not cross the line containing TNEGAX.
C
         ANGLE = DASINE ( (MAXR-1.D0) / D,  ATOL  ) 

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZEDTMPT' )
            RETURN
         END IF

C
C        Create the tangent point on the transformed target. 
C
         THETA = - ( HALFPI() + ANGLE )

         CALL VROTV ( UTAXIS, HPLNML, THETA, TARGPT )

C
C        S is the sign applied to pi/2 - ANGLE.
C
         S = 1.D0

      ELSE
C
C        This is the penumbral case. The tangent ray crosses
C        the line containing TNEGAX.
C
C
C        The tangent line always slopes downward (toward AXIS)
C        toward the light source.
C           
         ANGLE = DASINE ( (MAXR+1.D0) / D,  ATOL  ) 


         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZEDTMPT' )
            RETURN
         END IF

C
C        Create the tangent point on the transformed target. 
C
         THETA = ANGLE - HALFPI()

         CALL VROTV ( UTAXIS, HPLNML, THETA, TARGPT )

         S = -1.D0

      END IF

C
C     The tangent point is also a normal direction for the plane
C     tangent to both objects. Get the corresponding unit normal and
C     the plane constant.
C
      CALL VHAT ( TARGPT, NORMAL )

      CONST = VDOT( NORMAL, TARGPT )
C
C     Find the height of the plane relative to the transformed source.
C     We'll find the unique point on the transformed source where the
C     outward normal is parallel to NORMAL and find the height of this
C     point relative to the plane.
C
C     Let SGNNML be the "signed" normal which is parallel to NORMAL
C     in the umbral case and anti-parallel otherwise.
C
      CALL VSCL   ( S,  NORMAL, SGNNML )

      CALL EDNMPT ( XA, XB, XC, SGNNML, SRCPNT )

C
C     Express the source point as an offset from the transformed
C     target center.
C
      CALL VADD ( SRCPNT, TAXIS, TMPVEC )
      CALL VEQU ( TMPVEC,        SRCPNT )

C
C     H is the height of the surface point on the source, relative
C     to the plane tangent to the target at TARGPT. ANGERR is the
C     corresponding angular error estimate: an estimate of the 
C     amount by which TARGPT needs to be rotated in the positive
C     sense about HPLNML to make the plane contain SRCPNT.
C
      H      = VDOT( SRCPNT, NORMAL ) - CONST       
         
      ANGERR = TOUCHD( -H / D )

      NITR   = 0

C
C     The loop terminates when the angular error magnitude
C     stops decreasing. If the iteration count exceeds the
C     limit, an error will be signaled.
C
      DO WHILE (       ( ABS(ANGERR) .GT. CNVLIM )
     .           .AND. ( NITR        .LE. MAXITR )  )
C
C        Rotate the target point about HPLNML in the positive sense
C        by the angular error. This should make the tangent plane
C        closer to the source point.
C
         CALL VROTV ( TARGPT, HPLNML, ANGERR, TMPVEC )
         CALL VEQU  ( TMPVEC,                 TARGPT )
         CALL VHAT  ( TARGPT, NORMAL )

C
C        Re-compute the normal and constant of the tangent plane.
C
         CONST = VDOT( NORMAL, TARGPT )

         CALL VSCL   ( S,  NORMAL, SGNNML )

C
C        Find the near point on the source to the tangent plane.
C
         CALL EDNMPT ( XA, XB, XC, SGNNML, SRCPNT )

         CALL VADD ( SRCPNT, TAXIS, TMPVEC )
         CALL VEQU ( TMPVEC,        SRCPNT )
            
C
C        Re-compute the height error and angular error.
C
         H      = VDOT  ( SRCPNT, NORMAL ) - CONST

         CALL VPERP ( SRCPNT, HPLNML, PROJ )

         D      = VDIST ( PROJ, TARGPT )
         ANGERR = TOUCHD( -H / D )

         NITR   = NITR + 1

         IF ( NITR .GT. MAXITR ) THEN

            CALL SETMSG ( 'Tangent finding loop failed to '
     .      //            'converge. Iteration count = #.' )
            CALL ERRINT ( '#', NITR                        )
            CALL SIGERR ( 'SPICE(NOCONVERGENCE)'           )
            CALL CHKOUT ( 'ZZEDTMPT'                       )
            RETURN

         END IF

      END DO

C
C     Apply the inverse distortion transformation to TARGPT in order to
C     obtain the tangent point on the original, ellipsoidal target.
C
      POINT(1) = A * TARGPT(1)
      POINT(2) = B * TARGPT(2)
      POINT(3) = C * TARGPT(3)

      CALL CHKOUT ( 'ZZEDTMPT' )
      RETURN
      END
