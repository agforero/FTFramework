C$Procedure ZZDBRGAP ( DSKBRIEF, compute gaps in coverage )

      SUBROUTINE ZZDBRGAP ( CORSYS, NREC,  BDS1,  BDS2, 
     .                      MAXN,   NCOMP, CBDS1, CBDS2 )

C$ Abstract
C
C     Find spatial coverage gaps for a set of rectangles representing
C     the coverage of DSK segments.
C
C     The gap region is represented as a list of rectangles in a
C     specified coordinate system. The rectangles are disjoint except
C     for possible overlap at their edges.
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

      INCLUDE 'dskdsc.inc'

      INTEGER               MAXC
      PARAMETER           ( MAXC   =  100000 )
      
      INTEGER               MAXGRD
      PARAMETER           ( MAXGRD = 1000000 )

      INTEGER               CORSYS
      INTEGER               NREC
      DOUBLE PRECISION      BDS1    ( 2, NREC )
      DOUBLE PRECISION      BDS2    ( 2, NREC )
      INTEGER               MAXN
      INTEGER               NCOMP
      DOUBLE PRECISION      CBDS1   ( 2, * )
      DOUBLE PRECISION      CBDS2   ( 2, * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     CORSYS     I   Coordinate system code.
C     NREC       I   Number of input rectangles.
C     BDS1       I   Bounds of the rectangles' first coordinates.
C     BDS2       I   Bounds of the rectangles' second coordinates.
C     MAXN       I   Maximum number of output rectangles.
C     NCOMP      O   Number of output components.
C     CBDS1      O   Bounds of components' first coordinates.
C     CBDS2      O   Bounds of components' second coordinates.
C     MAXC       P   Maximum number of components.
C     MAXGRD     P   Maximum size of pixel grid.
C     
C$ Detailed_Input
C
C     CORSYS     is the DSK code for the coordinate system in
C                which the input rectangles are represented.
C
C                The supported coordinate systems are
C
c                   Latitudinal
c                   Planetodetic
c                   Rectangular
C
C                The coordinate system plays little role in
C                determination of gaps, other than that longitude
C                requires special treatment so that longitude values
C                can be compared.
C
C     NREC       is the number of input rectangles.
C
C     BDS1,
C     BDS2       are, respectively, the bounds of the first and
C                second coordinates of the input rectangles.
C                The elements
C
C                   BDS1(J,I), J = 1, 2
C
C                are the lower and upper bounds of the first 
C                coordinate. The elements
C
C                   BDS2(J,I), J = 1, 2
C
C                are the lower and upper bounds of the second
C                coordinate.
C
C                Units are radians for angular coordinates and
C                km for distance coordinates.
C
C     MAXN       is the maximum number of gap components that
C                can be returned.
C
C$ Detailed_Output
C
C     NCOMP      is the number of gap components found.
C
C     CBDS1,
C     CBDS2      are, respectively, the bounds of the first and second
C                coordinates of the rectangular components comprising
C                the gap region.
C
C                The bounds are expressed in the coordinate system
C                designated by CORSYS.
C
C                The components are disjoint except for possible
C                overlap at their edges.
C        
C                Units of the bounds are those of the corresponding
C                input coordinates.
C
C$ Parameters
C
C     MAXC       is the maximum number of components that can be 
C                accommodated by this routine's workspace arrays.
C
C     MAXGRD     is the maximum number of pixels in the workspace
C                grid used by this routine.
C
C$ Exceptions
C
C     1)  If any array used in the gap computation is too small, an
C         error will be signaled by a routine in the call tree of this
C         routine.
C 
C     2)  Longitudes outside the range
C
C            -2*pi - ANGMRG : 2*pi + ANGMRG 
C
C         are not accepted: if such a value is encountered, the
C         error will be diagnosed by a routine in the call tree 
C         of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports determination of spatial coverage gaps of a
C     set of DSK segments.
C
C     The algorithms used in this routine and those it calls attempt
C     to operate efficiently: sorting operations are kept to a minimum.
C
C     The algorithms attempt to achieve efficiency in return for memory
C     consumption.
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
C-    DSKBRIEF Version 1.0.0, 30-JAN-2017 (NJB)
C
C        Updated to accommodate FNDCMP's new argument list order.
C
C        Original version 1.0.0, 07-OCT-2016 (NJB)
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )

      INTEGER               MAXC2
      PARAMETER           ( MAXC2  = 2 * MAXC )

C
C     Local variables
C
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      OUTXBD ( 2, MAXC )
      DOUBLE PRECISION      OUTYBD ( 2, MAXC )
      
      INTEGER               CIVORX ( MAXC2 )
      INTEGER               CIVORY ( MAXC2 )
      INTEGER               CMPORX ( MAXC2 )
      INTEGER               CMPORY ( MAXC2 )
      INTEGER               H
      INTEGER               I
      INTEGER               J
      INTEGER               K     
      INTEGER               MAXPXX ( MAXC2 )
      INTEGER               MAXPXY ( MAXC2 )
      INTEGER               MINPXX ( MAXC2 )
      INTEGER               MINPXY ( MAXC2 )
      INTEGER               MRKSET ( LBCELL : MAXC2 )
      INTEGER               NCOLS
      INTEGER               NR
      INTEGER               NROWS
      INTEGER               ORDX   ( MAXC2 )
      INTEGER               ORDY   ( MAXC2 )
      INTEGER               SRCS   ( MAXC2 )
      INTEGER               TMPSET ( LBCELL : MAXC2 )
      INTEGER               VSET   ( LBCELL : MAXC2 )

      LOGICAL               COVERD
      LOGICAL               GAPVAL
      LOGICAL               GRID   ( MAXGRD )

C
C     Saved variables
C
      SAVE


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDBRGAP' )

C
C     Initialize the workspace sets.
C
      CALL SSIZEI ( MAXC, MRKSET )
      CALL SSIZEI ( MAXC, TMPSET )
      CALL SSIZEI ( MAXC, VSET   )

      IF (  ( CORSYS .EQ. LATSYS ) .OR. ( CORSYS .EQ. PDTSYS )  ) THEN
C
C        Adjust the input longitudes to make them usable by RC2GRD.
C
         CALL REGLON ( NREC,   BDS1,   MAXN,   NR,
     .                 MINLON, MAXLON, OUTXBD, SRCS )
C
C        Since REGLON may create new rectangles, the input set of "Y"
C        bounds (actually latitude) may not match up with the bounds in
C        OUTXBD. Create a new array of "Y" bounds parallel to the array
C        of X bounds. In this array, each rectangle has the Y bounds of
C        the source box from which it was created.
C
         DO I = 1, NR

            OUTYBD(1,I) = BDS2( 1, SRCS(I) )
            OUTYBD(2,I) = BDS2( 2, SRCS(I) )

         END DO

         
      ELSE
C
C        Just transfer the input X and Y bounds.
C        
         CALL MOVED ( BDS1, 2*NREC, OUTXBD )
         CALL MOVED ( BDS2, 2*NREC, OUTYBD )

         NR = NREC

      END IF

C
C     Map the coordinate rectangles to a pixel grid. Mark
C     the coverage with the value .TRUE.
C
      COVERD = .TRUE.
      GAPVAL = .FALSE.

      CALL RC2GRD ( NR,     OUTXBD, OUTYBD, MAXGRD, MAXC, 
     .              COVERD, ORDX,   ORDY,   CIVORX, CIVORY,  
     .              CMPORX, CMPORY, NROWS,  NCOLS,  GRID   )
C
C     Map the gaps in the pixel grid to a set of rectangles in pixel
C     space.
C     
      I = NROWS*NCOLS
      K = 0

      DO J = 1, I
         
         IF ( GRID(J) ) THEN
            K = K + 1
         END IF

      END DO


      CALL FNDCMP ( NROWS,  NCOLS,  GAPVAL, MAXN,  GRID,
     .              VSET,   MRKSET, TMPSET, NCOMP,
     .              MINPXX, MAXPXX, MINPXY, MAXPXY       )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZDBRGAP' )
         RETURN
      END IF

C
C     Map the gap rectangles, which are expressed in pixel coordinates,
C     to a set of rectangles in the input coordinate system.
C
      DO I = 1, NCOMP
C
C        If the range of a pixel coordinate is 
C
C           A:B
C
C        then the range of the indices of the corresponding bounds in
C        the compressed, ordered set of bounds is
C
C           A : B+1
C
C        Map each bound index to the corresponding value using the
C        mappings output by RC2GRD.
C
C        We need to deal with the fact that the arrays CMPORX and
C        CMPORY treat the bounds arrays OUTXBD and OUTYBD as
C        one-dimensional.
C
         J = CMPORX( MINPXX(I) )
         
         K = 1 + ( (J-1)/2 ) 
         H = J - ( (K-1)*2 )

         CBDS1(1,I) = OUTXBD( H, K )


         J = CMPORX( MAXPXX(I)+1 )
         
         K = 1 + ( (J-1)/2 ) 
         H = J - ( (K-1)*2 )

         CBDS1(2,I) = OUTXBD( H, K )


         J = CMPORY( MINPXY(I) )
         
         K = 1 + ( (J-1)/2 ) 
         H = J - ( (K-1)*2 )

         CBDS2(1,I) = OUTYBD( H, K )


         J = CMPORY( MAXPXY(I)+1 )
         
         K = 1 + ( (J-1)/2 ) 
         H = J - ( (K-1)*2 )

         CBDS2(2,I) = OUTYBD( H, K )
        
      END DO

      CALL CHKOUT ( 'ZZDBRGAP' )
      RETURN
      END
