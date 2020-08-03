C$Procedure ZZINRYPL ( Simplified intersection of ray and plane )
 
      SUBROUTINE ZZINRYPL ( VERTEX, UDIR, UPLNML, 
     .                      CONST,  MAXD, NXPTS,  XPT )
      IMPLICIT NONE
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Quickly find the intersection of a ray and a plane, without
C     performing normal error handling.
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
C     PLANES
C
C$ Keywords
C
C     GEOMETRY
C
C$ Declarations
 
      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      UPLNML ( 3 )
      DOUBLE PRECISION      CONST 
      DOUBLE PRECISION      MAXD
      INTEGER               NXPTS
      DOUBLE PRECISION      XPT    ( 3 )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     VERTEX,
C     UDIR       I   Vertex and unit length direction vector of ray.
C     UPLNML     I   Unit plane normal vector.
C     CONST      I   Plane constant.
C     MAXD       I   Maximum distance of intersection from vertex.
C     NXPTS      O   Number of intersection points of ray and plane.
C     XPT        O   Intersection point, if NXPTS = 1.
C
C$ Detailed_Input
C
C     VERTEX,
C     UDIR           are a point and unit-length direction vector that
C                    define a ray in three-dimensional space. The ray
C                    is the set of points in 3-dimensional space 
C
C                       { X : X = VERTEX + s*UDIR, s >= 0 }
C
C
C     UPLNML,
C     CONST          are a unit-length plane normal vector and plane
C                    constant. The plane is the set of points in
C                    3-dimensional space
C
C                       { X :  < X, UPLNML > = CONST }
C
C
C     MAXD           is the maximum length of a ray-plane intercept
C                    from the ray's vertex. If an intercept exists but
C                    has distance from VERTEX greater than or equal to
C                    MAXD, the intercept will be considered NOT to
C                    exist.
C
C$ Detailed_Output
C
C     NXPTS          is the number of points of intersection of the
C                    input ray and plane.  Values and meanings of
C                    NXPTS are:
C
C                       0     Either no intersection, or the ray
C                             lies in the plane.
C
C                       1     One point of intersection. Note that
C                             this case may occur when the ray's
C                             vertex is in the plane.
C
C
C     XPT            is the point of intersection of the input ray
C                    and plane, when there is exactly one point of
C                    intersection. Otherwise, XPT is the zero vector.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  The ray's direction vector must have unit length. The
C         outputs of this routine will be invalid otherwise.
C
C     2)  The plane's normal vector must have unit length, or 
C         else the outputs of this routine will be invalid.
C     
C     3)  If an intercept exists but is too far from the ray's
C         vertex, no intersection will be considered to exist.
C
C     4)  If the input ray lies in the input plane, no intersection
C         is considered to exist.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine performs a simplified version of the ray-plane
C     intersection computation performed by INRYPL. This routine
C     dispenses with error handling and nuanced handling of near-
C     singular geometry in return for speed. 
C
C     On a PC/Linux/gfortran/64bit platform on which this routine
C     was tested, it ran about 8 times faster than the public
C     SPICELIB routine INRYPL.
C     
C     This routine is meant to be used only by the DSK subsystem.
C
C$ Examples
C
C     See usage in ZZRYXLAT.
C
C$ Restrictions
C
C     1) All inputs must be checked by the caller; they are not checked
C        here.
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
C-    SPICELIB Version 1.0.0, 11-JAN-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     simplified intersection of ray and plane
C
C-&
 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      VDOT

C
C     Local variables
C
      DOUBLE PRECISION      DIRCON
      DOUBLE PRECISION      H
      DOUBLE PRECISION      LPAR
      DOUBLE PRECISION      S
      DOUBLE PRECISION      VTXCON

C
C     Start out indicating no intersection.
C
      NXPTS = 0

C
C     VTXCON is the plane constant of the ray's vertex.
C
      VTXCON = VDOT( VERTEX, UPLNML )
C
C     DIRCON is the length of the component of the ray's
C     direction vector in the direction of UPLNML.
C
      DIRCON = VDOT( UDIR, UPLNML )

C
C     Dispose of the easy non-intersection cases. (The ray 
C     lying in the plane is considered a non-intersection case,
C     by the way.)
C     
      IF ( ( VTXCON .GT. CONST ) .AND. ( DIRCON .GT. 0.D0 ) ) THEN
         RETURN
      END IF

      IF ( ( VTXCON .LT. CONST ) .AND. ( DIRCON .LT. 0.D0 ) ) THEN
         RETURN
      END IF


      IF ( VTXCON .EQ. CONST ) THEN
C
C        The ray's vertex lies in the plane.
C        
         IF ( DIRCON .NE. 0.D0 ) THEN
C
C           The ray does not lie in the plane. The
C           intercept is the ray's vertex.
C
            NXPTS = 1
            CALL VEQU ( VERTEX, XPT )

         END IF

         RETURN

      END IF

C
C     Let UPAR and UPERP be, respectively, the components of UDIR
c     parallel to and perpendicular to UPLNML. 
C
C     Compute the maximum allowed length of UPERP.
C     
C
      H     = ABS( VTXCON - CONST )
      LPAR  = ABS( DIRCON )

C
C     To prevent overflow, we require 
C
C          H
C        ----  <= MAXD
C        LPAR 
C
C     or equivalently
C
C        H  <=  MAXD * LPAR 
C
      IF ( H .GT. (MAXD * LPAR) ) THEN
         RETURN
      END IF

C
C     For safety, return if we could have a divide-by-zero error.
C
      IF ( LPAR .EQ. 0.D0 ) THEN
         RETURN
      END IF

C
C     Still being here means we can compute XPT, provided
C     the given value of MAXD was reasonable.
C
C     Note that the earlier tests we performed should 
C     rule out the case 
C
C        DIRCON = 0
C
C     We have also ruled out overflow in the computation below.
C
      S      = H / LPAR

      XPT(1) = VERTEX(1) + S*UDIR(1)
      XPT(2) = VERTEX(2) + S*UDIR(2)
      XPT(3) = VERTEX(3) + S*UDIR(3)

      NXPTS  = 1

      END
