C$Procedure ZZTANUTL ( DSK, tangent ray utilities )
 
      SUBROUTINE ZZTANUTL ( CURVE,  SRCRAD, SHAPE, TRGCDE, NSURF,  
     .                      SRFLST, FIXFID, ET,    PLNVEC, AXIS,
     .                      ANGLE,  OCULTD, POINT                )
  
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     This is the umbrella routine for utilities supporting the 
C     tangent ray finding capability used by LIMBPT and TERMPT.
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
C     DLA
C     DSK
C     LIMB
C     RAY
C     TANGENT
C     TERMINATOR
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'zzdsk.inc'

      INTEGER               CURVE
      DOUBLE PRECISION      SRCRAD
      INTEGER               SHAPE
      INTEGER               TRGCDE
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      INTEGER               FIXFID
      DOUBLE PRECISION      ET
      DOUBLE PRECISION      PLNVEC ( 3 )
      DOUBLE PRECISION      AXIS   ( 3 )
      DOUBLE PRECISION      ANGLE
      LOGICAL               OCULTD
      DOUBLE PRECISION      POINT  ( 3 )
 
C$ Brief_I/O
C
C     Variable  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     CURVE      I   ZZTANINI
C     SRCRAD     I   ZZTANINI
C     SHAPE      I   ZZTANINI
C     TRGCDE     I   ZZTANINI
C     NSURF      I   ZZTANINI
C     SRFLST     I   ZZTANINI
C     FIXFID     I   ZZTANINI
C     ET         I   ZZTANINI
C     PLNVEC     I   ZZTANINI
C     AXIS       I   ZZTANINI
C     ANGLE      I   ZZTANSTA
C     OCULTD     O   ZZTANSTA
C     POINT      O   ZZTANSTA
C
C$ Detailed_Input
C
C     See the entry points.
C
C$ Detailed_Output
C
C     See the entry points. 
C
C$ Parameters
C
C     See the entry points.
C
C$ Exceptions
C
C     1)  If this routine is called directly, it signals the error 
C         SPICE(BOGUSENTRY).
C        
C     See the entry points for descriptions of errors specific to
C     those routines.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C
C     This routine uses PCK data for target body reference
C     ellipsoids.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine contains the following entry points that support
C     the generalized ray-surface intercept algorithm used by SINCPT:
C
C        ZZTANINI:   Initialize tangent ray utilities for limb or
C                    terminator computation. Terminators may be
C                    umbral or penumbral.
C
C        ZZTANSTA:   Compute state of ray for a given input angle.
C                    State is "occulted" or "not occulted." When
C                    the state is occulted, return the ray-surface 
C                    intercept point.
C
C$ Examples
C
C     See usage in ZZTANGNT.
C
C$ Restrictions
C
C     1) This is a private routine. 
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
C-     SPICELIB Version 1.0.0 16-FEB-2016 (NJB)
C
C        Added error handling and headers.
C
C        Original version 26-OCT-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     umbrella for tangent ray utilities
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      PI

      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local variables
C
      DOUBLE PRECISION      APEX   ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      SVAXIS ( 3 )
      DOUBLE PRECISION      SVET
      DOUBLE PRECISION      SVIRAD
      DOUBLE PRECISION      SVNRML ( 3 )
      DOUBLE PRECISION      SVVRTX ( 3 )
      DOUBLE PRECISION      VRTOFF ( 3 )
 
      INTEGER               SVCURV

C
C     Saved variables
C
      SAVE                  SVAXIS
      SAVE                  SVCURV
      SAVE                  SVET
      SAVE                  SVIRAD
      SAVE                  SVNRML
      SAVE                  SVVRTX
 
      CALL CHKIN  ( 'ZZTANUTL'          )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)' )
      CALL CHKOUT ( 'ZZTANUTL'          )
      RETURN


    

C$Procedure ZZTANINI ( DSK, tangent utility initialization )
 
      ENTRY ZZTANINI ( CURVE,  SRCRAD, SHAPE, TRGCDE, NSURF, 
     .                 SRFLST, FIXFID, ET,    PLNVEC, AXIS   )
  
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Initialize the data used by the tangent ray state routine
C     ZZTANSTA, which is used directly by ZZTANGNT and indirectly by
C     LIMBPT and TERMPT.
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
C     DLA
C     DSK
C     LIMB
C     RAY
C     TANGENT
C     TERMINATOR
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     INTEGER               CURVE
C     DOUBLE PRECISION      SRCRAD
C     INTEGER               SHAPE
C     INTEGER               TRGCDE
C     INTEGER               NSURF
C     INTEGER               SRFLST ( * )
C     INTEGER               FIXFID
C     DOUBLE PRECISION      ET
C     DOUBLE PRECISION      PLNVEC ( 3 )
C     DOUBLE PRECISION      AXIS   ( 3 )
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     CURVE      I   Type of curve: limb or terminator type.
C     SRCRAD     I   Radius of illumination source.
C     SHAPE      I   Target shape.
C     TRGCDE     I   Target body ID code.
C     NSURF      I   Number of surfaces in list.
C     SRFLST     I   Surface ID list.
C     FIXFID     I   Frame ID of target body-fixed frame.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     PLNVEC     I   Reference vector contained in cutting half-plane.
C     AXIS       I   Axis vector: edge of cutting half-plane.
C     ELLSHP     P   Ellipsoid shape code.
C     DSKSHP     P   DSK shape code.
C     LMBCRV     P   Limb code.
C     UMBRAL     P   Umbral terminator code.
C     PNMBRL     P   Penumbral terminator code.
C
C$ Detailed_Input
C
C     CURVE      is an integer code indicating the type of set on the
C                target body on which tangent points are to be found.
C                When the target is an ellipsoid, the set is literally
C                a curve. When the target is a DSK shape, the set is
C                often not connected; so it is not actually a curve;
C                still, we retain the familiar terminology for the
C                ellipsoid case.
C
C                Possible values and meanings are:
C
C                   LMBCRV            limb
C                   UMBRAL            umbral terminator
C                   PNMBRL            penumbral terminator      
C
C                Terminator computations are performed assuming
C                the light source is spherical.
C
C
C     SRCRAD     is the radius of the illumination source. This
C                value is used for terminator computations only.
C                When used, SRCRAD must be strictly positive.
C                Units are km.
C
c
C     SHAPE      is an integer code indicating the target body shape.
C                
C                Possible values and meanings are:
C
C                   ELLSHP            shape is modeled as an ellipsoid
C                   DSKSHP            shape is provided by DSK data
C
C
C     TRGCDE     is the body ID code of the target body.
C
c
C     NSURF,
C     SRFLST     are, respectively, the count of surface IDs and the
C                surface ID list.
C
c
C     FIXFID     is the frame ID code of a body-fixed frame centered
C                on the target body. The output tangent points will
C                be expressed in this frame.
C
C
C     ET         is the computation epoch, expressed as seconds past
C                J2000 TDB. ET is used for DSK segment selection.
C                If the target shape is modeled as an ellipsoid, ET
C                is ignored.
C
C
C     PLNVEC     is a vector used to define a half-plane in which
C                to search for tangent rays. The half-plane contains
C                PLNVEC, the target body center, and AXIS.
C
C                For limb and umbral terminator computations, 
C                tangent rays lie in the half-plane containing PLNVEC.
C        
C                For penumbral terminator computations, tangent rays
C                touch the target in the half-plane containing PLNVEC,
C                but touch the light source in the complementary
C                half-plane bounded by the line containing AXIS.
C     
C
C     AXIS       is a second vector used to define a half-plane in 
C                which to search for tangent vectors. AXIS lies in
C                the line forming the boundary of the half-plane.
C
C                For limb computations, AXIS points from the target
C                toward the observation point.
C
C                For terminator computations, AXIS points from the
C                target toward the center of the illumination source.
C                
C
C$ Detailed_Output
C
C     None. This routine operates by side effects.
C
C$ Parameters
C
C     See the INCLUDE file zzdsk.inc for declarations of these 
C     parameters.
C
C
C     ELLSHP     is a code specifying an ellipsoidal shape model.
C
C     DSKSHP     is a code specifying a DSK-based shape model.
C
C     LMBCRV     is a code specifying a limb computation.
C
C     UMBRAL     is a code specifying an umbral terminator computation.
C               
C                The umbral terminator is the boundary of the portion
C                of the target surface that receives no light from the
C                illumination source.
C
C     PNMBRL     is a code specifying an penumbral terminator
C                computation.
C
C                The penumbral terminator is the boundary of the
C                portion of the target surface that is not subject to
C                self-occultation of light from the illumination
C                source. Given that the light source is modeled as
C                a sphere, from any target surface point nearer to
C                the source than the penumbral terminator, the source
C                appears to be a lit disc. 
C                
C
C$ Exceptions
C
C     1)  If AXIS is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C        
C     2)  If PLNDEF is the zero vector, the error SPICE(ZEROVECTOR) is
C         signaled.
C
C     3)  If the curve type code is unrecognized, the error
C         SPICE(BADCURVETYPE) is signaled.
C
C     4)  If a terminator curve is specified by the source radius
C         is non-positive, the error SPICE(BADSOURCERADIUS) is
C         signaled.
C
C     5)  If AXIS and PLNDEF are linearly dependent, the error
C         SPICE(DEGENERATECASE) is signaled.
C
C     6)  If the target shape code is unrecognized, the error
C         SPICE(BADSHAPE) is signaled.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C
C     This routine uses PCK data for target body reference
C     ellipsoids.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine initializes the tangent ray utilities for limb or
C     terminator computation.
C
C$ Examples
C
C     See usage in ZZTANGNT.
C
C$ Restrictions
C
C     1) This is a private routine. 
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
C-     SPICELIB Version 1.0.0 11-FEB-2016 (NJB)
C
C        Original version 26-OCT-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     initialize tangent ray finding state function
C
C-&
 

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZTANINI' )

C
C     Check for zero vectors on input.
C
      IF ( VZERO(AXIS) ) THEN

         CALL SETMSG ( 'Input axis vector is the zero vector.' )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                     )
         CALL CHKOUT ( 'ZZTANINI'                              )
         RETURN

      END IF

      IF ( VZERO(PLNVEC) ) THEN

         CALL SETMSG ( 'Input reference vector is the zero vector.' )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                          )
         CALL CHKOUT ( 'ZZTANINI'                                   )
         RETURN

      END IF
      
C
C     Save the curve type.
C
      IF (       ( CURVE .NE. LMBCRV ) 
     .    .AND.  ( CURVE .NE. UMBRAL ) 
     .    .AND.  ( CURVE .NE. PNMBRL )  ) THEN

         CALL SETMSG ( 'Curve type code # was not recognized.' )
         CALL ERRINT ( '#', CURVE                              )
         CALL SIGERR ( 'SPICE(BADCURVETYPE)'                   )
         CALL CHKOUT ( 'ZZTANINI'                              )
         RETURN         

      END IF

      SVCURV = CURVE

C
C     Save the illumination source radius.
C
      IF (  ( CURVE .EQ. UMBRAL ) .OR. ( CURVE .EQ. PNMBRL )  ) THEN

         IF ( SRCRAD .LE. 0.D0 ) THEN

            CALL SETMSG ( 'The source radius was #. The radius '
     .      //            'must be positive for a terminator '
     .      //            'computation.'                         )
            CALL ERRDP  ( '#', SRCRAD                            )
            CALL SIGERR ( 'SPICE(BADSOURCERADIUS)'               )
            CALL CHKOUT ( 'ZZTANINI'                             )
            RETURN

         END IF

      END IF

      SVIRAD = SRCRAD

C
C     Compute a normal vector to the plane defined by 
C     AXIS and PLNVEC. The direction of positive rotation
C     about the normal is from AXIS toward PLNVEC.
C
      CALL VCRSS ( AXIS, PLNVEC, SVNRML )

      IF ( VZERO(SVNRML) ) THEN

         CALL SETMSG ( 'Input reference vector and axis vector '
     .   //            'are linearly dependent.'                 )
         CALL SIGERR ( 'SPICE(DEGENERATECASE)'                   )
         CALL CHKOUT ( 'ZZTANINI'                                )
         RETURN

      END IF

C
C     Scale the normal vector to unit length.
C
      CALL VHATIP ( SVNRML )

C
C     Save a unit-length copy of the input axis.
C     Save the original axis as the ray's vertex; this
C     will be used directly in the limb computation and 
C     will be used, after addition of an offset, in the
C     terminator computations. Save the evaluation epoch.
C
      CALL VEQU ( AXIS, SVVRTX )
      CALL VHAT ( AXIS, SVAXIS )

      SVET = ET
C
C     Prepare the DSK SINCPT utilities for a computation with
C     the input surface set.
C
      IF ( SHAPE .EQ. ELLSHP ) THEN
C
C        This is the ellipsoid case.
C
         CALL ZZSUELIN ( TRGCDE )

      ELSE IF ( SHAPE .EQ. DSKSHP ) THEN
C
C        This is the DSK case.
C
         CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FIXFID )

      ELSE

         CALL SETMSG ( 'Target shape code # was not recognized.' )
         CALL ERRINT ( '#', SHAPE                                )
         CALL SIGERR ( 'SPICE(BADSHAPE)'                         )
         CALL CHKOUT ( 'ZZTANINI'                                )
         RETURN         
         
      END IF

      CALL CHKOUT ( 'ZZTANINI' )
      RETURN




C$Procedure ZZTANSTA ( DSK, tangent ray state )
 
      ENTRY ZZTANSTA ( ANGLE, OCULTD, POINT )
  
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     This is the state callback routine used by ZZTANSLV to 
C     find tangent rays on the target body's surface.
C
C     Indicate whether a vector emanating from a given vertex, and
C     rotated about the normal by the specified angle, measured from
C     the AXIS direction, intersects the target. If an intersection
C     exists, return the intercept closest to the ray's vertex.
C    
C     The ray's vertex depends on the curve type set via a call
C     to ZZTANINI.
C
C     This routine is called directly by ZZTANGNT and is used
C     indirectly by LIMBPT and TERMPT.
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
C     DLA
C     DSK
C     LIMB
C     RAY
C     TANGENT
C     TERMINATOR
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
C
C     DOUBLE PRECISION      ANGLE
C     LOGICAL               OCULTD
C     DOUBLE PRECISION      POINT  ( 3 )
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     ANGLE      I   Angle of ray.
C     OCULTD     O   Occultation state flag. True if ray is occulted.
C     POINT      O   Tangent point on target.
C
C$ Detailed_Input
C
C     ANGLE      is the angle between a ray and the AXIS vector 
C                stored by ZZTANINI. AXIS points away from the
C                target.
C
C                The vertex of the vector depends on the curve type,
C                which also is stored by ZZTANINI. 
C
C                For a limb computation, the vertex is the observer's
C                location.
C
C                For an umbral terminator computation, the vertex 
C                is on the surface of the light source, in the half-
C                plane defined by PLNDEF and AXIS. The line containing
C                the ray is tangent to the light source at the vertex.
C
C
C                For a penumbral terminator computation, the vertex is
C                on the surface of the light source, in the half-plane
C                complementary to that defined by PLNDEF and AXIS. The
C                line containing the ray is tangent to the light source
C                at the vertex.
C
C                Units are radians.
C
C$ Detailed_Output
C
C     OCULTD     is a logical flag that is .TRUE. if and only if the
C                ray defined by ANGLE and the values set by ZZTANINI
C                intersects the target.
C
C     POINT      is the ray-surface intercept closest to the ray's
C                vertex, if an intercept exists.
C
C$ Parameters
C
C     See zzdsk.inc for declarations of parameters used internally.
C
C$ Exceptions
C
C     1)  Any errors that occur while looking up DSK data will be
C         signaled by a routine in the call tree of this routine.
C
C     2)  Any errors that occur while computing ray-surface intercepts
C         will be signaled by a routine in the call tree of this
C         routine.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C
C     This routine uses PCK data for target body reference
C     ellipsoids.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine computes the state of ray for a given input angle.
C     The state is "occulted" or "not occulted." When the state is
C     occulted, this routine returns the surface intercept point
C     nearest to the ray's vertex.
C
C     This is the state callback routine used by ZZTANSLV to 
C     find tangent rays on the target body's surface.
C
C$ Examples
C
C     See usage in ZZTANGNT and ZZTANSLV.
C
C$ Restrictions
C
C     1) This is a private routine. 
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
C-     SPICELIB Version 1.0.0 11-FEB-2016 (NJB)
C
C        Original version 26-OCT-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     tangent ray finding occultation state and intercept
C
C-&
 

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZTANSTA' )


      IF ( SVCURV .EQ. LMBCRV ) THEN
C
C        This is the limb case.
C
C        We'll rotate SVAXIS by ANGLE to achieve the desired result.
C
         CALL VROTV ( SVAXIS, SVNRML, ANGLE, RAYDIR )
      
         CALL ZZRAYSFX ( SVVRTX, RAYDIR, SVET, POINT, OCULTD )


      ELSE IF ( SVCURV .EQ. UMBRAL ) THEN
C
C        This is the umbral terminator case.
C
C        Produce the ray's direction vector by rotating
C        the axis about the cutting half-plane normal by 
C        the input angle.
C
         CALL VROTV ( SVAXIS, SVNRML, ANGLE, RAYDIR )
C
C        Produce the offset of the ray's vertex from the
C        center of the source by rotating the axis
C        vector by ANGLE-pi/2 radians. The length
C        of the vector must be SVIRAD. The saved axis
C        has unit length.
C
         CALL VROTV  ( SVAXIS, SVNRML, ANGLE-(PI()/2), VRTOFF )
         CALL VSCLIP ( SVIRAD, VRTOFF )
         CALL VADD   ( SVVRTX, VRTOFF, APEX )

         CALL ZZRAYSFX ( APEX, RAYDIR, SVET, POINT, OCULTD )


      ELSE IF ( SVCURV .EQ. PNMBRL ) THEN
C
C        This is the penumbral terminator case.
C
C        Produce the ray's direction vector by rotating
C        the axis about the cutting half-plane normal by 
C        the *negative* of the input angle.
C
         CALL VROTV ( SVAXIS, SVNRML, -ANGLE, RAYDIR )
C
C        Produce the ray's vertex by rotating the axis
C        vector about the normal, *not its negative,*
C        by 3*pi/2 - ANGLE radians. The length of the vector 
C        must be SRCRAD. The saved axis has unit length.
C
         CALL VROTV  ( SVAXIS, SVNRML, (1.5D0*PI())-ANGLE, VRTOFF )

         CALL VSCLIP ( SVIRAD, VRTOFF )

         CALL VADD   ( SVVRTX, VRTOFF, APEX )

         CALL ZZRAYSFX ( APEX, RAYDIR, SVET, POINT, OCULTD )

      ELSE
C
C        This case should have been ruled out by a check in
C        ZZTANINI. Check again anyway.
C
         CALL SETMSG ( 'Bad curve type code #.' )
         CALL ERRINT ( '#', SVCURV              )
         CALL SIGERR ( 'SPICE(BUG)'             )
         CALL CHKOUT ( 'ZZTANSTA'               )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZTANSTA' )
      RETURN
      END
