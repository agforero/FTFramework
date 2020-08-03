C$Procedure IOVCMP ( Inverse order vector with compressed range )

      SUBROUTINE IOVCMP ( DARRAY, NDIM, IORDER, INVORD, RNGMAX )

C$ Abstract
C
C     Create an inverse order vector having a compressed range.
C     In this vector, the Ith element is the index of DARRAY(I)
C     in the array that would be produced by sorting DARRAY and
C     removing duplicates.
C
C     Produce this result in O(N) time.
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

      DOUBLE PRECISION      DARRAY ( * )
      INTEGER               NDIM
      INTEGER               IORDER ( * )
      INTEGER               INVORD ( * )
      INTEGER               RNGMAX

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     DARRAY     I   Double precision array.
C     NDIM       I   Size of DARRAY.
C     IORDER     O   Order vector for DARRAY.
C     INVORD     O   Compressed inverse order vector.
C     RNGMAX     O   Maximum value in range of INVORD.
C
C$ Detailed_Input
C
C     DARRAY     is an array of double precision numbers.
C
C     NDIM       is the size of DARRAY.
C
C$ Detailed_Output
C
C     IORDER     is an order vector for DARRAY.
C
C     INVORD     is a compressed inverse order vector for DARRAY.
C
C                INVORD(I) is the position that DARRAY(I) would
C                have in an ordered set created by sorting DARRAY
C                and removing duplicates. 
C
C                INVORD has size NDIM. Its elements belong to the 
C                set {1, ..., RNGMAX}.
C
C     RNGMAX     is the largest element in INVORD. This value is the
C                number of distinct values in DARRAY.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     Error free.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports efficient determination of spatial
C     coverage gaps by DSKBRIEF.
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
C     DSKBRIEF Version 1.0.0, 05-OCT-2016 (NJB)
C
C-&




C
C     Local variables
C
      INTEGER               I
      INTEGER               NUPRED


C
C     First step: create an order vector for DARRAY.
C
      CALL ORDERD ( DARRAY, NDIM, IORDER )

C
C     Produce the corresponding inverse order vector.
C
      DO I = 1, NDIM

         INVORD( IORDER(I) ) = I

      END DO
      
C
C     Step through the order vector, keeping track of the count of
C     unique predecessors, in the array that would be produced by
C     sorting DARRAY, of each element pointed to by an element of the
C     order vector.
C
C     The element of DARRAY at index IORDER(1) has no predecessors,
C     and the element INVORD( IORDER(1) ) is already correct. So
C     we start at the second element of IORDER (if it exists).
C
C     Initialize NUPRED to the number of unique predecessors of
C     the first value.
C
      NUPRED = 0
      
      DO I = 2, NDIM
C
C        At this point, NUPRED is the number of unique predecessors of
C        DARRAY(I). I is greater than or equal to 2.
C
         IF (  DARRAY( IORDER(I) )  .GT.  DARRAY( IORDER(I-1) )  ) THEN
C
C           DARRAY( IORDER(I) ) is strictly greater than, and hence not
C           a copy of, its predecessor. It has NUPRED + 1 unique
C           predecessors in the array produced by sorting DARRAY.
C
            NUPRED = NUPRED + 1

C
C           The position of DARRAY( IORDER(I) ) in the sorted,
C           compressed set derived from DARRAY is one more than the
C           count of unique predecessors.
C
            INVORD( IORDER(I) ) = NUPRED + 1

         ELSE
C
C           DARRAY( IORDER(I) ) is a duplicate. Its position in the
C           sorted, compressed array derived from DARRAY is the same
C           as that of DARRAY( IORDER(I-1) ).
C
            INVORD( IORDER(I) ) = INVORD( IORDER(I-1) )

         END IF

      END DO

C
C     INVORD has been updated so that its elements belong to the
C     set { 1 : NUPRED+1 }.
C     
C     Set the maximum range value of INVORD.
C     
      RNGMAX = NUPRED + 1

      END







