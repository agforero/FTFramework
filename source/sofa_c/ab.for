      SUBROUTINE iau_AB ( PNAT, V, S, BM1, PPR )
*+
*  - - - - - - -
*   i a u _ A B
*  - - - - - - -
*
*  Apply aberration to transform natural direction into proper
*  direction.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*    PNAT      d(3)   natural direction to the source (unit vector)
*    V         d(3)   observer barycentric velocity in units of c
*    S         d      distance between the Sun and the observer (au)
*    BM1       d      sqrt(1-|v|^2): reciprocal of Lorenz factor
*
*  Returned:
*    PPR       d(3)   proper direction to source (unit vector)
*
*  Notes:
*
*  1) The algorithm is based on Expr. (7.40) in the Explanatory
*     Supplement (Urban & Seidelmann 2013), but with the following
*     changes:
*
*     o  Rigorous rather than approximate normalization is applied.
*
*     o  The gravitational potential term from Expr. (7) in
*        Klioner (2003) is added, taking into account only the Sun's
*        contribution.  This has a maximum effect of about
*        0.4 microarcsecond.
*
*  2) In almost all cases, the maximum accuracy will be limited by the
*     supplied velocity.  For example, if the SOFA iau_EPV00 routine is
*     used, errors of up to 5 microarcseconds could occur.
*
*  References:
*
*     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
*     the Astronomical Almanac, 3rd ed., University Science Books
*     (2013).
*
*     Klioner, Sergei A., "A practical relativistic model for micro-
*     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
*
*  Called:
*     iau_PDP      scalar product of two p-vectors
*
*  This revision:   2013 August 31
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION PNAT(3), V(3), S, BM1, PPR(3)

*  Schwarzschild radius of the Sun (au)
*  = 2 * 1.32712440041 D20 / (2.99792458 D8)^2 / 1.49597870700 D11
      DOUBLE PRECISION SRS
      PARAMETER ( SRS = 1.97412574336D-08 )

      INTEGER I
      DOUBLE PRECISION PDV, W1, W2, R2, W, P(3), R

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      CALL iau_PDP ( PNAT, V, PDV )
      W1 = 1D0 + PDV/(1D0+BM1)
      W2 = SRS / S
      R2 = 0D0
      DO 1 I=1,3
         W = PNAT(I)*BM1 + W1*V(I) + W2*(V(I)-PDV*PNAT(I))
         P(I) = W
         R2 = R2 + W*W
 1    CONTINUE
      R = SQRT ( R2 )
      DO 2 I=1,3
         PPR(I) = P(I) / R
 2    CONTINUE

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
