C$Procedure ZZRYXSPH ( Intersection of ray and sphere )
 
      SUBROUTINE ZZRYXSPH ( VERTEX, UDIR, R, XPT, FOUND )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Quickly find the intersection of a ray and a sphere, without
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
 
      IMPLICIT NONE

      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      UDIR   ( 3 )
      DOUBLE PRECISION      R
      DOUBLE PRECISION      XPT    ( 3 )
      LOGICAL               FOUND
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     VERTEX,
C     UDIR       I   Vertex and unit length direction vector of ray.
C     R          I   Radius of sphere.
C     XPT        O   Intersection point, if NXPTS = 1.
C     FOUND      O   Flag indicating whether intersection exists.
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
C     R              is the radius of a sphere. The sphere is centered
C                    at the origin.
C
C$ Detailed_Output
C
C     XPT            is the point of intersection nearest to the ray's
C                    vertex of the input ray and sphere, if the
C                    intersection exists. Otherwise, XPT is undefined.
C
C     FOUND          is a local flag that is .TRUE. if and only if 
C                    the input ray intersects the sphere.
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
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine performs a simplified version of the ray-ellipsoid
C     intersection computation performed by SURFPT; this routine works
C     with spheres only. This routine dispenses with error handling and
C     nuanced handling of near- singular geometry in return for speed.
C
C     On a PC/Linux/gfortran/64bit platform on which this routine
C     was tested, it ran about 5 times faster than the public
C     SPICELIB routine SURFPT.
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
C     intersection of ray and sphere
C
C-&
 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      VDOT

C
C     Local variables
C     
      DOUBLE PRECISION      CPAR
      DOUBLE PRECISION      PERP   ( 3 )
      DOUBLE PRECISION      PMAG2 
      DOUBLE PRECISION      R2
      DOUBLE PRECISION      S
      DOUBLE PRECISION      VMAG2

      FOUND = .FALSE.
 
C
C     Find the component of VERTEX orthogonal to UDIR. If the magnitude
C     of this component exceeds R, there's no intercept.
C     
      CPAR = VDOT ( VERTEX, UDIR )

      CALL VLCOM ( 1.D0, VERTEX, -CPAR, UDIR, PERP )

      PMAG2 = VDOT( PERP, PERP )
      R2    = R * R

C
C     Compare squares of magnitudes, rather than magnitudes, for
C     efficiency.
C
      IF ( PMAG2 .GT. R2 ) THEN 
         RETURN
      END IF
         
      S     = SQRT (  MAX( 0.D0,  R2-PMAG2 )  )

      VMAG2 = VDOT(VERTEX,VERTEX)

      IF ( VMAG2 .GT. R2 ) THEN
C
C        If the magnitude of the vertex exceeds R, the vertex is
C        outside the sphere. Above, we have compared squares of
C        magnitudes for efficiency.
C
         IF ( CPAR .GT. 0.D0 ) THEN
C
C           The ray points away from the sphere; there can be no
C           intersection.
C
            RETURN

         END IF
C
C        Given that an intercept exists, we can find it between VERTEX
C        and VPERP by following -UDIR from PERP towards VERTEX.
C     
         XPT(1) = PERP(1) - S*UDIR(1)
         XPT(2) = PERP(2) - S*UDIR(2)
         XPT(3) = PERP(3) - S*UDIR(3)

      ELSE IF ( VMAG2 .LT. R2 ) THEN
C
C        The vertex is inside the sphere. We can calculate the exit
C        point by using PERP as a vertex.
C
         XPT(1) = PERP(1) + S*UDIR(1)
         XPT(2) = PERP(2) + S*UDIR(2)
         XPT(3) = PERP(3) + S*UDIR(3)

      ELSE
C
C        PERP is the sole intercept.
C
         XPT(1) = PERP(1) 
         XPT(2) = PERP(2)
         XPT(3) = PERP(3) 

      END IF

      FOUND = .TRUE.

      END
