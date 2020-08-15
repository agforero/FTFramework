      SUBROUTINE iau_LDN ( N, B, OB, SC, SN )
*+
*  - - - - - - - -
*   i a u _ L D N
*  - - - - - - - -
*
*  For a star, apply light deflection by multiple solar-system bodies,
*  as part of transforming coordinate direction into natural direction.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     N     i        number of bodies (Note 1)
*     B     d(8,N)   data for each of the N bodies (Notes 1,2):
*            (1,I)     mass of the body (solar masses, Note 3)
*            (2,I)     deflection limiter (Note 4)
*            (3-5,I)   barycentric position of the body (au)
*            (6-8,I)   barycentric velocity of the body (au/day)
*     OB    d(3)     barycentric position of the observer (au)
*     SC    d(3)     observer to star coordinate direction (unit vector)
*
*  Returned:
*     SN    d(3)     observer to deflected star (unit vector)
*
*  1) The array B contains N entries, one for each body to be
*     considered.  If N = 0, no gravitational light deflection will be
*     applied, not even for the Sun.
*
*  2) The array B should include an entry for the Sun as well as for any
*     planet or other body to be taken into account.  The entries should
*     be in the order in which the light passes the body.
*
*  3) In the entry in the B array for body I, the mass parameter B(1,I)
*     can, as required, be adjusted in order to allow for such effects
*     as quadrupole field.
*
*  4) The deflection limiter parameter B(2,I) is phi^2/2, where phi is
*     the angular separation (in radians) between star and body at which
*     limiting is applied.  As phi shrinks below the chosen threshold,
*     the deflection is artificially reduced, reaching zero for phi = 0.
*     Example values suitable for a terrestrial observer, together with
*     masses, are as follows:
*
*        body I     B(1,I)         B(2,I)
*
*        Sun        1D0            6D-6
*        Jupiter    0.00095435D0   3D-9
*        Saturn     0.00028574D0   3D-10
*
*  5) For cases where the starlight passes the body before reaching the
*     observer, the body is placed back along its barycentric track by
*     the light time from that point to the observer.  For cases where
*     the body is "behind" the observer no such shift is applied.  If
*     a different treatment is preferred, the user has the option of
*     instead using the iau_LD routine.  Similarly, iau_LD can be used
*     for cases where the source is nearby, not a star.
*
*  6) The returned vector SN is not normalized, but the consequential
*     departure from unit magnitude is always negligible.
*
*  7) For efficiency, validation is omitted.  The supplied masses must
*     be greater than zero, the position and velocity vectors must be
*     right, and the deflection limiter greater than zero.
*
*  Reference:
*
*     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
*     the Astronomical Almanac, 3rd ed., University Science Books
*     (2013), Section 7.2.4.
*
*  Called:
*     iau_CP       copy p-vector
*     iau_PDP      scalar product of two p-vectors
*     iau_PMP      p-vector minus p-vector
*     iau_PPSP     p-vector plus scaled p-vector
*     iau_PN       decompose p-vector into modulus and direction
*     iau_LD       light deflection by a solar-system body
*
*  This revision:   2017 March 16
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION B(8,N), OB(3), SC(3), SN(3)

*  Seconds per day.
      DOUBLE PRECISION DAYSEC
      PARAMETER ( DAYSEC = 86400D0 )

*  Speed of light (m/s)
      DOUBLE PRECISION CMPS
      PARAMETER ( CMPS = 299792458D0 )

*  Astronomical unit (m, IAU 2012)
      DOUBLE PRECISION AUM
      PARAMETER ( AUM = 149597870.7D3 )

*  Light time for 1 au (day)
      DOUBLE PRECISION CR
      PARAMETER ( CR = AUM/CMPS/DAYSEC )

      INTEGER I
      DOUBLE PRECISION S(3), V(3), D , DT, EV(3), EM, E(3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Star direction prior to deflection.
      CALL iau_CP ( SC, S )

*  Body by body.
      DO 1 I=1,N

*     Body to observer vector at epoch of observation (au).
         CALL iau_PMP ( OB, B(3,I), V )

*     Minus the time since the light passed the body (days).
         CALL iau_PDP ( S, V, D )
         DT = D * CR

*     Neutralize if the star is "behind" the observer.
         DT = MIN ( DT, 0D0 )

*     Backtrack the body to the time the light was passing the body.
         CALL iau_PPSP ( V, -DT, B(6,I), EV )

*     Separate the body to observer vector into magnitude and direction.
         CALL iau_PN ( EV, EM, E )

*     Apply light deflection for this body.
         CALL iau_LD ( B(1,I), S, S, E, EM, B(2,I), V )

*     Update the star direction.
         CALL iau_CP ( V, S )

*     Next body.
 1    CONTINUE

*  Return the deflected star direction.
      CALL iau_CP ( S, SN )

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
