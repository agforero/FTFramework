C$Procedure MAKVTX ( MKDSK: create vertex from coordinates )

      SUBROUTINE MAKVTX ( CORSYS, CORPAR, COORDS, 
     .                    REFVAL, HEIGHT, VERTEX )
      
C$ Abstract
C
C     Create a vertex from domain coordinates and height.
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
C     MKDSK
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'
      
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      COORDS ( 2 )
      DOUBLE PRECISION      REFVAL
      DOUBLE PRECISION      HEIGHT
      DOUBLE PRECISION      VERTEX ( 3 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     CORSYS     I   Coordinate system.
C     CORPAR     I   Coordinate parameters.
C     COORDS     I   Domain coordinates.
C     REFVAL     I   Height reference value.
C     HEIGHT     I   Height.
C     VERTEX     O   Vertex.
C
C$ Detailed_Input
C 
C     CORSYS     is a DSK subsystem code designating the coordinate
C                system of the input coordinates.
C
C     CORPAR     is an array containing parameters associated with
C                the input coordinate system. The contents of the
C                array are as described in the DSK include file
C                dskdsc.inc.
C
C     COORDS     is a pair of domain coordinates: these may be,
C               
C                   - planetocentric longitude and latitude
C
C                   - planetodetic longitude and latitude
C
C                   - X and Y
C
C                For a given coordinate system, the order of the
C                elements of COORDS is that of the coordinate names in
C                the list above.
C
C     REFVAL     is a reference value to be added to the input height.
C                REFVAL is used only for latitudinal and rectangular
C                coordinates.
C
C                REFVAL must be non-negative.
C
C                Units are km.
C
C
C     HEIGHT     is a height datum. Units are km.
C
C
C$ Detailed_Output
C
C     VERTEX     is a 3-vector corresponding to the input coordinates,
C                height, and if applicable, height reference.
C
C                Units are always km.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the coordinate system code is not recognized, the error
C         SPICE(NOTSUPPORTED) is signaled.
C
C     2)  If an error occurs while converting planetodetic coordinates
C         to rectangular coordinates, the error will be diagnosed by a
C         routine in the call tree of this routine.
C
C     3)  If REFVAL is negative, the error SPICE(VALUEOUTOFRANGE) is
C         signaled. REFVAL is checked whether or not it is applicable.
C     
C$ Files
C
C     None.
C
C$ Particulars
C
C     None.
C
C$ Examples
C
C     See usage in the MKDSK routine MKVARR.
C
C$ Restrictions
C
C     1) For use only within program MKDSK.
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
C-    MKDSK Version 1.0.0, 25-FEB-2017 (NJB)
C
C-&





C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C     
      DOUBLE PRECISION      ALT
      DOUBLE PRECISION      F
      DOUBLE PRECISION      RADIUS
      DOUBLE PRECISION      RE

      

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'MAKVTX' )


      IF ( REFVAL .LT. 0.D0 ) THEN

         CALL SETMSG ( 'Reference value # must be '
     .   //            'non-negative.'             )
         CALL ERRDP  ( '#', REFVAL                 )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'    )
         CALL CHKOUT ( 'MAKVTX'                    )
         RETURN

      END IF



      IF ( CORSYS .EQ. LATSYS ) THEN
C
C        REFVAL has the same units as HEIGHT.
C        HSCALE converts these units to km.
C
         RADIUS = REFVAL + HEIGHT

         CALL LATREC ( RADIUS, COORDS(1), COORDS(2), VERTEX )


      ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C        Height is relative to the system's reference spheroid.
C
         RE  = CORPAR(1)
         F   = CORPAR(2)
         ALT = HEIGHT

         CALL GEOREC ( COORDS(1), COORDS(2), ALT, RE, F, VERTEX )

      ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C        Height is relative to the reference Z-value.
C
C        REFVAL has the same units as HEIGHT.
C        HSCALE converts these units to km.
C
         VERTEX(1) = COORDS(1)
         VERTEX(2) = COORDS(2)
         VERTEX(3) = REFVAL + HEIGHT

      ELSE

         CALL SETMSG ( 'Coordinate system code # is not '
     .   //            'recognized.'                     )
         CALL ERRINT ( '#', CORSYS                       )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'             )
         CALL CHKOUT ( 'MAKVTX'                          )
         RETURN

      END IF

      CALL CHKOUT ( 'MAKVTX' )
      RETURN
      END


