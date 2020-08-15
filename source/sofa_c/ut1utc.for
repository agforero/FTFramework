      SUBROUTINE iau_UT1UTC ( UT11, UT12, DUT1, UTC1, UTC2, J )
*+
*  - - - - - - - - - - -
*   i a u _ U T 1 U T C
*  - - - - - - - - - - -
*
*  Time scale transformation:  Universal Time, UT1, to Coordinated
*  Universal Time, UTC.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical.
*
*  Given:
*     UT11,UT12    d      UT1 as a 2-part Julian Date (Note 1)
*     DUT1         d      Delta UT1: UT1-UTC in seconds (Note 2)
*
*  Returned:
*     UTC1,UTC2    d      UTC as a 2-part quasi Julian Date (Notes 3,4)
*     J            i      status: +1 = dubious year (Note 5)
*                                  0 = OK
*                                 -1 = unacceptable date
*
*  Notes:
*
*  1) UT11+UT12 is Julian Date, apportioned in any convenient way
*     between the two arguments, for example where UT11 is the Julian
*     Day Number and UT12 is the fraction of a day.  The returned UTC1
*     and UTC2 form an analogous pair, except that a special convention
*     is used, to deal with the problem of leap seconds - see Note 3.
*
*  2) Delta UT1 can be obtained from tabulations provided by the
*     International Earth Rotation and Reference Systems Service.  The
*     value changes abruptly by 1s at a leap second;  however, close to
*     a leap second the algorithm used here is tolerant of the "wrong"
*     choice of value being made.
*
*  3) JD cannot unambiguously represent UTC during a leap second unless
*     special measures are taken.  The convention in the present routine
*     is that the returned quasi JD day UTC1+UTC2 represents UTC days
*     whether the length is 86399, 86400 or 86401 SI seconds.
*
*  4) The routine iau_D2DTF can be used to transform the UTC quasi-JD
*     into calendar date and clock time, including UTC leap second
*     handling.
*
*  5) The warning status "dubious year" flags UTCs that predate the
*     introduction of the time scale or that are too far in the future
*     to be trusted.  See iau_DAT for further details.
*
*  Called:
*     iau_JD2CAL   JD to Gregorian calendar
*     iau_DAT      delta(AT) = TAI-UTC
*     iau_CAL2JD   Gregorian calendar to JD
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992)
*
*  This revision:  2019 June 20
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION UT11, UT12, DUT1, UTC1, UTC2
      INTEGER J

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

      LOGICAL BIG1
      INTEGER I, IY, IM, ID, JS
      DOUBLE PRECISION DUTS, U1, U2, D1, DATS1, D2, FD, DATS2, DDATS,
     :                 US1, US2, DU

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  UT1-UTC in seconds.
      DUTS = DUT1

*  Put the two parts of the UT1 into big-first order.
      BIG1 = ABS(UT11) .GE. ABS(UT12)
      IF ( BIG1 ) THEN
         U1 = UT11
         U2 = UT12
      ELSE
         U1 = UT12
         U2 = UT11
      END IF

*  See if the UT1 can possibly be in a leap-second day.
      D1 = U1
      DATS1 = 0D0
      DO 1 I=-1,3
         D2 = U2 + DBLE(I)
         CALL iau_JD2CAL ( D1, D2, IY, IM, ID, FD, JS )
         IF ( JS.NE.0 ) GO TO 9
         CALL iau_DAT ( IY, IM, ID, 0D0, DATS2, JS )
         IF ( JS.LT.0 ) GO TO 9
         IF ( I.EQ.-1 ) DATS1 = DATS2
         DDATS = DATS2 - DATS1
         IF ( ABS(DDATS).GE.0.5D0 ) THEN

*        Yes, leap second nearby: ensure UT1-UTC is "before" value.
            IF ( DDATS*DUTS.GE.0D0 ) DUTS = DUTS-DDATS

*        UT1 for the start of the UTC day that ends in a leap.
            CALL iau_CAL2JD ( IY, IM, ID, D1, D2, JS )
            US1 = D1
            US2 = D2 - 1D0 + DUTS/D2S

*        Is the UT1 after this point?
            DU = U1 - US1
            DU = DU + ( U2 - US2 )
            IF ( DU.GT.0D0 ) THEN

*           Yes:  fraction of the current UTC day that has elapsed.
               FD = DU * D2S / ( D2S + DDATS )

*           Ramp UT1-UTC to bring about SOFA's JD(UTC) convention.
               DUTS = DUTS + DDATS*MIN(FD,1D0)
            END IF

*        Break.
            GO TO 2
         END IF
         DATS1 = DATS2
 1    CONTINUE
 2    CONTINUE

*  Subtract the (possibly adjusted) UT1-UTC from UT1 to give UTC.
      U2 = U2 - DUTS/D2S

*  Result, safeguarding precision.
      IF ( BIG1 ) THEN
         UTC1 = U1
         UTC2 = U2
      ELSE
         UTC1 = U2
         UTC2 = U1
      END IF

*  Return the status.
 9    CONTINUE
      J = JS

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
