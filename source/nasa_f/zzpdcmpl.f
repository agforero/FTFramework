C$Procedure ZZPDCMPL (Planetodetic coordinates, compare latitudes )
 
      SUBROUTINE ZZPDCMPL ( RE, F, P, LAT, REL )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Compare the planetodetic latitude of a point in 3-dimensional
C     space against a specified value, without converting the point to
C     planetodetic coordinates.
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
C     COORDINATES
C     GEOMETRY
C     PLANETODETIC
C     LATITUDE
C     MATH
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      RE
      DOUBLE PRECISION      F
      DOUBLE PRECISION      P   ( 3 )
      DOUBLE PRECISION      LAT
      INTEGER               REL

      INTEGER               LT
      PARAMETER           ( LT = -1 )
      
      INTEGER               EQ
      PARAMETER           ( EQ =  0 )

      INTEGER               GT
      PARAMETER           ( GT =  1 )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     RE         I   Equatorial radius.
C     F          I   Flattening coefficient.
C     P          I   Three-dimensional point.
C     LAT        I   Planetodetic latitude.
C     REL        O   Relation code.
C     LT         P   Code indicating latitude of P < LAT.
C     EQ         P   Code indicating latitude of P = LAT.
C     GT         P   Code indicating latitude of P > LAT.
C
C$ Detailed_Input
C
C     RE,
C     F          are, respectively, the equatorial radius
C                and flattening coefficient of a biaxial 
C                spheroid. 
C
C                The polar radius RP of the spheroid is
C                
C                   RP = RE * ( 1 - F )
C
C                RP may be less than, equal to, or greater than RE.
C
C
C     P          is a point (equivalently, a vector) in
C                three-dimensional space. P is expressed in Cartesian
C                coordinates.
C
C                The units of P must be consistent with those of RE.
C
C
C     LAT        is a planetodetic latitude value to be compared
C                against the planetodetic latitude of P. Units
C                are radians.
C
C$ Detailed_Output
C
C     REL        is an integer code that indicates the order 
C                relation between the planetodetic latitude of P
C                and LAT. The planetodetic coordinate system is
C                defined by the inputs RE and F.
C
C                The code <rel> indicates that the relation
C
C                   <latitude of P>  <rel>  LAT
C
C                is true. See the Parameters section below for
C                the parameter names.
C
C
C$ Parameters
C
C     LT,
C     EQ,
C     GT         are, respectively, codes indicating the relationship
C                of the planetodetic latitude of the input vector to
C                the input latitude value. Let LP represent the 
C                planetodetic latitude of P. 
C
C                   Code LT indicates   LP < LAT
C                   Code EQ indicates   LP = LAT
C                   Code GT indicates   LP > LAT
C
C$ Exceptions
C
C     1)  If either the equatorial radius or flattening coefficient
C         is invalid, the error will be signaled by a routine in the
C         call tree of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine performs a planetodetic latitude comparison more
C     efficiently than can be done using a rectangular-to-planetodetic
C     coordinate conversion.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     There are some cases for which this routine cannot be applied. 
C     See the SPICELIB routine ZZPDPLTC and its usage in ZZRYTPDT.
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
C-    SPICELIB Version 1.0.0, 19-JAN-2017 (NJB)
C
C       Original version 22-AUG-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     compare planetodetic latitude of vector against value
C     compare planetodetic latitude of point against value
C
C-&

C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      APEX   ( 3 )
      DOUBLE PRECISION      LON
      DOUBLE PRECISION      OFFPCL
      DOUBLE PRECISION      OFFSET ( 3 )
      DOUBLE PRECISION      R
      DOUBLE PRECISION      RP
      DOUBLE PRECISION      XINCPT
      DOUBLE PRECISION      YINCPT

C
C     Saved variables
C
      SAVE                  APEX

C
C     Initial values
C
      DATA                  APEX / 3  *  0.D0       /

      
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZPDCMPL' )

C
C     Treat points on the Z axis as a special case. The
C     computations performed in the general case may introduce
C     round-off errors that will lead to false results for
C     this case.
C
      IF (  ( P(1) .EQ. 0.D0 ) .AND. ( P(2) .EQ. 0.D0 )  ) THEN

         IF ( P(3) .GT. 0.D0 ) THEN

            IF ( LAT .EQ. HALFPI() ) THEN

               REL = EQ
            ELSE
               REL = GT
            END IF

         ELSE IF ( P(3) .EQ. 0.D0 ) THEN
C
C           We consider the latitude of P to be zero.
C
            IF ( LAT .GT. 0.D0 ) THEN
               
               REL = LT
          
            ELSE IF ( LAT .EQ. 0.D0 ) THEN

               REL = EQ
            ELSE
               REL = GT
            END IF
           
         ELSE
C
C           P(3) < 0.
C
            IF ( LAT .EQ. -HALFPI() ) THEN

               REL = EQ
            ELSE
               REL = LT
            END IF

         END IF

         CALL CHKOUT ( 'ZZPDCMPL' )
         RETURN

      END IF

C
C     Latitude zero is a special case. The planetodetic latitude of the
C     input point has the same sign as the Z component of the point.
C
      RP = RE * ( 1.D0 - F )
C
C     Get the y-intercept of the latitude cone for LAT. Note that a
C     result is defined for LAT = +/- pi/2.
C
      CALL ZZELNAXX ( RE, RP, LAT, XINCPT, YINCPT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZPDCMPL' )
         RETURN
      END IF

C
C     Ideally YINCPT is zero if and only if LAT is zero.
C     We'll group these conditions together.
C
      IF (  ( LAT .EQ. 0.D0 ) .OR. ( YINCPT .EQ. 0.D0 )  ) THEN

         IF ( P(3) .GT. 0.D0 ) THEN

            REL = GT

         ELSE IF ( P(3) .EQ. 0.D0 ) THEN

            REL = EQ
         ELSE            
            REL = LT
         END IF

         CALL CHKOUT ( 'ZZPDCMPL' )
         RETURN
         
      END IF

C
C     This is the normal case.
C
C     Find the offset of the point from the latitude cone's apex.
C     Create a unit-length copy of the offset vector.
C
      APEX(3) = YINCPT

      CALL VSUB ( P, APEX, OFFSET )
 
C     We'll use the planetocentric [sic] latitude of the offset
C     vector for comparison.
C
      CALL RECLAT ( OFFSET, R, LON, OFFPCL )
 

      IF ( LAT .GT. 0.D0 ) THEN

         IF ( YINCPT .GT. 0 ) THEN
C
C           This is the prolate case. 
C              
            IF ( OFFPCL .GT. LAT ) THEN

               REL = GT

            ELSE IF ( OFFPCL .EQ. LAT ) THEN

               REL = EQ
            ELSE
               REL = LT
            END IF

         ELSE
C
C           YINCPT = 0 was handled previously, so YINCPT < 0.
C
C           This is the oblate case. 
C
C           In addition to the comparison of angles, we need to know
C           the input point is above the X-Y plane in order for the
C           GT or EQ relations to hold.
C
            IF ( P(3) .GT. 0.D0 ) THEN

               IF ( OFFPCL .GT. LAT ) THEN

                  REL = GT

               ELSE IF ( OFFPCL .EQ. LAT ) THEN

                  REL = EQ
               ELSE
                  REL = LT
               END IF

            ELSE
C
C              The input latitude is positive, while the point
C              is on or below the X-Y plane.
C
               REL = LT

            END IF

         END IF 


      ELSE
C
C        LAT < 0, since the case LAT = 0 has already been handled.
C
         IF ( YINCPT .LT. 0.D0 ) THEN
C
C           This is the prolate case.
C
            IF ( OFFPCL .GT. LAT ) THEN

               REL = GT

            ELSE IF ( OFFPCL .EQ. LAT ) THEN
               
               REL = EQ
            ELSE 
               REL = LT
            END IF

         ELSE
C
C           YINCPT > 0, since the case YINCPT = 0 was handled
C           previously.
C
C           This is the oblate case. 
C
            IF ( P(3) .LT. 0.D0 ) THEN

               IF ( OFFPCL .GT. LAT ) THEN
                  
                  REL = GT
               
               ELSE IF ( OFFPCL .EQ. LAT ) THEN
                  
                  REL = EQ
               ELSE
                  REL = LT
               END IF

            ELSE
C
C              The input latitude is negative, while the point
C              is on or above the X-Y plane.
C
               REL = GT

            END IF

         END IF

      END IF

      CALL CHKOUT ( 'ZZPDCMPL' )
      RETURN
      END
