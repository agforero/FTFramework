C$Procedure             PRODAI ( Product of an integer array )
 
      INTEGER FUNCTION  PRODAI ( ARRAY, N )
 
C$ Abstract
C
C      Return the product of the elements of an integer array.
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
C      ARRAY,  MATH,  UTILITY
C
C$ Declarations
 
      INTEGER          ARRAY  ( * )
      INTEGER          N
 
C$ Brief_I/O
C
C      VARIABLE  I/O  DESCRIPTION
C      --------  ---  --------------------------------------------------
C      ARRAY      I   Input array.
C      N          I   Number of elements in ARRAY.
C      PRODAI     O   Product of the elements of ARRAY.
C
C$ Detailed_Input
C
C      ARRAY       is the input array.
C
C      N           is the number of elements in the array.
C
C$ Detailed_Output
C
C      PRODAI      is the product of the elements of the input array.
C                  That is,
C
C                     PRODAI = ARRAY(1) * ARRAY(2) * ... * ARRAY(N)
C
C                  If N is zero or negative, PRODAI is one.
C
C$ Parameters
C
C     None.
C
C$ Particulars
C
C      The value of the function is initially set to one. The elements
C      of the array are then multiplied. If the number of elements is
C      zero or negative, PRODAI is one.
C
C$ Examples
C
C      Let ARRAY contain the following elements.
C
C            ARRAY(1) = 12
C            ARRAY(2) =  2
C            ARRAY(3) =  4
C            ARRAY(4) = 75
C            ARRAY(5) = 18
C
C      Then
C
C            PRODAI ( ARRAY,   -3 )       =      1
C            PRODAI ( ARRAY,    0 )       =      1
C            PRODAI ( ARRAY,    1 )       =     12
C            PRODAI ( ARRAY,    2 )       =     24
C            PRODAI ( ARRAY,    5 )       = 129600
C            PRODAI ( ARRAY(3), 3 )       =   5400
C
C
C$ Restrictions
C
C      PRODAI does not check for overflow. (For integers, this can
C      occur relatively quickly.)
C
C$ Exceptions
C
C      Error free.
C
C$ Files
C
C      None.
C
C$ Author_and_Institution
C
C      I.M. Underwood  (JPL)
C
C$ Literature_References
C
C      None.
C
C$ Version
C
C-     SPICELIB Version 1.0.1, 10-MAR-1992 (WLT)
C
C         Comment section for permuted index source lines was added
C         following the header.
C
C-     SPICELIB Version 1.0.0, 31-JAN-1990 (IMU)
C
C-&
 
C$ Index_Entries
C
C     product of an integer array
C
C-&
 
 
 
 
C
C     Local variables
C
      INTEGER          PROD
      INTEGER          I
 
 
C
C     Begin at one.
C
      PROD = 1
 
C
C     Multiply the elements. If N is zero or negative, nothing happens.
C
      DO I = 1, N
         PROD = PROD * ARRAY(I)
      END DO
 
C
C     Return the product.
C
      PRODAI = PROD
 
      RETURN
      END
