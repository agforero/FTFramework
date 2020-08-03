C$Procedure      VNORMG ( Vector norm, general dimension )
 
      DOUBLE PRECISION FUNCTION VNORMG ( V1, NDIM )
 
C$ Abstract
C
C      Compute the magnitude of a double precision vector of arbitrary
C      dimension.
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
C      VECTOR
C
C$ Declarations
 
      INTEGER            NDIM
      DOUBLE PRECISION   V1 ( NDIM )
 
C$ Brief_I/O
C
C      VARIABLE  I/O  DESCRIPTION
C      --------  ---  --------------------------------------------------
C       V1        I     Vector whose magnitude is to be found.
C       NDIM      I     Dimension of V1.
C
C$ Detailed_Input
C
C      V1      This may be any double precision vector or arbitrary
C              size.
C
C$ Detailed_Output
C
C      VNORMG is the magnitude of V1 calculated in a numerically stable
C      way.
C
C$ Parameters
C
C     None.
C
C$ Particulars
C
C      VNORMG finds the component of V1 whose magnitude is the largest.
C      If the absolute magnitude of that component indicates that a
C      numeric overflow would occur when it is squared, or if it
C      indicates that an underflow would occur when squared (falsely
C      giving a magnitude of zero) then the following expression is
C      used:
C
C      VNORMG = V1MAX * MAGNITUDE OF [ (1/V1MAX)*V1 ]
C
C      Otherwise a simpler expression is used:
C
C      VNORMG = MAGNITUDE OF [ V1 ]
C
C      Beyond the logic described above, no further checking of the
C      validity of the input is performed.
C
C$ Examples
C
C      The following table show the correlation between various input
C      vectors V1 and VNORMG:
C
C      NDIM   V1(NDIM)                               VNORMG
C      -----------------------------------------------------------------
C      1      (-7.0D20)                              7.D20
C      3      (1.D0, 2.D0, 2.D0)                     3.D0
C      4      (3.D0, 3.D0, 3.D0, 3.D0)               6.D0
C      5      (5.D0, 12.D0, 0.D0, 0.D0, 0.D0)       13.D0
C      3      (-5.D-17, 0.0D0, 12.D-17)             13.D-17
C
C$ Restrictions
C
C      None.
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
C      W.M. Owen       (JPL)
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
C-     SPICELIB Version 1.0.0, 31-JAN-1990 (WMO)
C
C-&
 
C$ Index_Entries
C
C     norm of n-dimensional vector
C
C-&
 
      INTEGER               I
      DOUBLE PRECISION      V1MAX
      DOUBLE PRECISION      A
 
C
C  Determine the maximum component of the vector.
C
      V1MAX = 0.D0
 
      DO I=1,NDIM
 
         IF ( DABS(V1(I)) .GT. V1MAX )
     .      V1MAX = DABS(V1(I))
 
      END DO
C
C  If the vector is zero, return zero; otherwise normalize first.
C  Normalizing helps in the cases where squaring would cause overflow
C  or underflow.  In the cases where such is not a problem it not worth
C  it to optimize further.
C
      IF ( V1MAX .EQ. 0.D0 ) THEN
 
         VNORMG = 0.D0
 
      ELSE
 
         VNORMG = 0.D0
 
         DO I = 1, NDIM
            A      = V1(I)  / V1MAX
            VNORMG = VNORMG + A*A
         END DO
 
         VNORMG = V1MAX * DSQRT (VNORMG)
 
      END IF
C
      RETURN
      END
