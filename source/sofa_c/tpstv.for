      SUBROUTINE iau_TPSTV ( XI, ETA, V0, V )
*
*  - - - - - - - - - -
*   i a u _ T P S T V
*  - - - - - - - - - -
*
*  In the tangent plane projection, given the star's rectangular
*  coordinates and the direction cosines of the tangent point, solve
*  for the direction cosines of the star.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     XI,ETA    d       rectangular coordinates of star image (Note 2)
*     V0        d(3)    tangent point's direction cosines (Note 4)
*
*  Returned:
*     V         d(3)    star's direction cosines
*
*  1) The tangent plane projection is also called the "gnomonic
*     projection" and the "central projection".
*
*  2) The eta axis points due north in the adopted coordinate system.
*     If the direction cosines represent observed (RA,Dec), the tangent
*     plane coordinates (xi,eta) are conventionally called the "standard
*     coordinates".  If the direction cosines are with respect to a
*     right-handed triad, (xi,eta) are also right-handed.  The units of
*     (xi,eta) are, effectively, radians at the tangent point.
*
*  3) The method used is to complete the star vector in the (xi,eta)
*     based triad and normalize it, then rotate the triad to put the
*     tangent point at the pole with the x-axis aligned to zero
*     longitude.  Writing (a0,b0) for the celestial spherical
*     coordinates of the tangent point, the sequence of rotations is
*     (b0-pi/2) around the x-axis followed by (-a0-pi/2) around the
*     z-axis.
*
*  4) If vector V0 is not of unit length, the returned vector V will
*     be wrong.
*
*  5) If vector V0 points at a pole, the returned vector V will be
*     based on the arbitrary assumption that the longitude coordinate
*     of the tangent point is zero.
*
*  6) This routine is a member of the following set:
*
*         spherical       vector       solve for
*
*         iau_TPXES      iau_TPXEV      xi,eta
*         iau_TPSTS    > iau_TPSTV <     star
*         iau_TPORS      iau_TPORV      origin
*
*  References:
*
*     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
*     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
*
*     Green, R.M., "Spherical Astronomy", Cambridge University Press,
*     1987, Chapter 13.
*
*  This revision:   2018 January 2
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION XI, ETA, V0(3), V(3)

      DOUBLE PRECISION X, Y, Z, R, F

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Tangent point.
      X = V0(1)
      Y = V0(2)
      Z = V0(3)

*  Deal with polar case.
      R = SQRT(X*X+Y*Y)
      IF ( R .EQ. 0D0 ) THEN
         R = 1D-20
         X = R
      END IF

*  Star vector length to tangent plane.
      F = SQRT(1D0+XI*XI+ETA*ETA)

*  Apply the transformation and normalize.
      V(1) = ( X - (XI*Y+ETA*X*Z) / R ) / F
      V(2) = ( Y + (XI*X-ETA*Y*Z) / R ) / F
      V(3) = ( Z + ETA*R ) / F

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
