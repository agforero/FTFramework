      SUBROUTINE iau_TAIUTC ( TAI1, TAI2, UTC1, UTC2, J )
*+
*  - - - - - - - - - - -
*   i a u _ T A I U T C
*  - - - - - - - - - - -
*
*  Time scale transformation:  International Atomic Time, TAI, to
*  Coordinated Universal Time, UTC.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical.
*
*  Given:
*     TAI1,TAI2    d      TAI as a 2-part Julian Date (Note 1)
*
*  Returned:
*     UTC1,UTC2    d      UTC as a 2-part quasi Julian Date (Notes 1-3)
*     J            i      status: +1 = dubious year (Note 4)
*                                  0 = OK
*                                 -1 = unacceptable date
*
*  Notes:
*
*  1) TAI1+TAI2 is Julian Date, apportioned in any convenient way
*     between the two arguments, for example where TAI1 is the Julian
*     Day Number and TAI2 is the fraction of a day.  The returned UTC1
*     and UTC2 form an analogous pair, except that a special convention
*     is used, to deal with the problem of leap seconds - see the next
*     note.
*
*  2) JD cannot unambiguously represent UTC during a leap second unless
*     special measures are taken.  The convention in the present routine
*     is that the JD day represents UTC days whether the length is
*     86399, 86400 or 86401 SI seconds.  In the 1960-1972 era there were
*     smaller jumps (in either direction) each time the linear UTC(TAI)
*     expression was changed, and these "mini-leaps" are also included
*     in the SOFA convention.
*
*  3) The routine iau_D2DTF can be used to transform the UTC quasi-JD
*     into calendar date and clock time, including UTC leap second
*     handling.
*
*  4) The warning status "dubious year" flags UTCs that predate the
*     introduction of the time scale or that are too far in the future
*     to be trusted.  See iau_DAT for further details.
*
*  Called:
*     iau_UTCTAI   UTC to TAI
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Explanatory Supplement to the Astronomical Almanac,
*     P. Kenneth Seidelmann (ed), University Science Books (1992)
*
*  This revision:  2020 June 22
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION TAI1, TAI2, UTC1, UTC2
      INTEGER J

      LOGICAL BIG1
      INTEGER I, JS
      DOUBLE PRECISION A1, A2, U1, U2, G1, G2

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Put the two parts of the TAI into big-first order.
      BIG1 = ABS(TAI1) .GE. ABS(TAI2)
      IF ( BIG1 ) THEN
         A1 = TAI1
         A2 = TAI2
      ELSE
         A1 = TAI2
         A2 = TAI1
      END IF

*  Initial guess for UTC.
      U1 = A1
      U2 = A2

*  Iterate (though in most cases just once is enough).
      DO 1 I=1,3

*     Guessed UTC to TAI.
         CALL iau_UTCTAI ( U1, U2, G1, G2, JS )
         IF ( JS.LT.0 ) GO TO 9

*     Adjust guessed UTC.
         U2 = U2 + (A1-G1)
         U2 = U2 + (A2-G2)

 1    CONTINUE

*  Return the UTC result, preserving the TAI order.
      IF ( BIG1 ) THEN
         UTC1 = U1
         UTC2 = U2
      ELSE
         UTC1 = U2
         UTC2 = U1
      END IF

*  Status.
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
