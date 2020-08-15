      SUBROUTINE iau_D2TF ( NDP, DAYS, SIGN, IHMSF )
*+
*  - - - - - - - - -
*   i a u _ D 2 T F
*  - - - - - - - - -
*
*  Decompose days to hours, minutes, seconds, fraction.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     NDP       i        resolution (Note 1)
*     DAYS      d        interval in days
*
*  Returned:
*     SIGN      c        '+' or '-'
*     IHMSF     i(4)     hours, minutes, seconds, fraction
*
*  Notes:
*
*  1) NDP is interpreted as follows:
*
*     NDP         resolution
*      :      ...0000 00 00
*     -7         1000 00 00
*     -6          100 00 00
*     -5           10 00 00
*     -4            1 00 00
*     -3            0 10 00
*     -2            0 01 00
*     -1            0 00 10
*      0            0 00 01
*      1            0 00 00.1
*      2            0 00 00.01
*      3            0 00 00.001
*      :            0 00 00.000...
*
*  2) The largest positive useful value for NDP is determined by the
*     size of DAYS, the format of DOUBLE PRECISION floating-point
*     numbers on the target platform, and the risk of overflowing
*     IHMSF(4).  On a typical platform, for DAYS up to 1D0, the
*     available floating-point precision might correspond to NDP=12.
*     However, the practical limit is typically NDP=9, set by the
*     capacity of a 32-bit IHMSF(4).
*
*  3) The absolute value of DAYS may exceed 1D0.  In cases where it
*     does not, it is up to the caller to test for and handle the
*     case where DAYS is very nearly 1D0 and rounds up to 24 hours,
*     by testing for IHMSF(1)=24 and setting IHMSF(1-4) to zero.
*
*  This revision:  2005 August 26
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NDP
      DOUBLE PRECISION DAYS
      CHARACTER SIGN*(*)
      INTEGER IHMSF(4)

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

      INTEGER NRS, N
      DOUBLE PRECISION RS, RM, RH, A, AH, AM, AS, AF

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Handle sign.
      IF ( DAYS .GE. 0D0 ) THEN
         SIGN = '+'
      ELSE
         SIGN = '-'
      END IF

*  Interval in seconds.
      A = D2S * ABS(DAYS)

*  Pre-round if resolution coarser than 1 second (then pretend NDP=1).
      IF ( NDP .LT. 0 ) THEN
         NRS = 1
         DO 1 N=1,-NDP
            IF ( N.EQ.2 .OR. N.EQ.4 ) THEN
               NRS = NRS * 6
            ELSE
               NRS = NRS * 10
            END IF
 1       CONTINUE
         RS = DBLE(NRS)
         A = RS * ANINT(A/RS)
      END IF

*  Express the unit of each field in resolution units.
      NRS = 1
      DO 2 N=1,NDP
         NRS = NRS * 10
 2    CONTINUE
      RS = DBLE(NRS)
      RM = RS * 60D0
      RH = RM * 60D0

*  Round the interval and express in resolution units.
      A = ANINT(RS*A)

*  Break into fields.
      AH = AINT(A/RH)
      A = A - AH*RH
      AM = AINT(A/RM)
      A = A - AM*RM
      AS = AINT(A/RS)
      AF = A - AS*RS

*  Return results.
      IHMSF(1) = NINT(AH)
      IHMSF(2) = NINT(AM)
      IHMSF(3) = NINT(AS)
      IHMSF(4) = NINT(AF)

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2020
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
*  CONDITIONS WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The names of all routines in your derived work shall not
*        include the prefix "iau" or "sofa" or trivial modifications
*        thereof such as changes of case.
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  5. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  6. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  In any published work or commercial product which uses the SOFA
*  software directly, acknowledgement (see www.iausofa.org) is
*  appreciated.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*-----------------------------------------------------------------------

      END
