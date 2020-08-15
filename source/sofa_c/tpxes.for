      SUBROUTINE iau_TPXES ( A, B, A0, B0, XI, ETA, J )
*+
*  - - - - - - - - - -
*   i a u _ T P X E S
*  - - - - - - - - - -
*
*  In the tangent plane projection, given celestial spherical
*  coordinates for a star and the tangent point, solve for the star's
*  rectangular coordinates in the tangent plane.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     A,B       d       star's spherical coordinates
*     A0,B0     d       tangent point's spherical coordinates
*
*  Returned:
*     XI,ETA    d       rectangular coordinates of star image (Note 2)
*     J         i       status:  0 = OK
*                                1 = star too far from axis
*                                2 = antistar on tangent plane
*                                3 = antistar too far from axis
*
*  Notes:
*
*  1) The tangent plane projection is also called the "gnomonic
*     projection" and the "central projection".
*
*  2) The eta axis points due north in the adopted coordinate system.
*     If the spherical coordinates are observed (RA,Dec), the tangent
*     plane coordinates (xi,eta) are conventionally called the "standard
*     coordinates".  For right-handed spherical coordinates, (xi,eta)
*     are also right-handed.  The units of (xi,eta) are, effectively,
*     radians at the tangent point.
*
*  3) All angular arguments are in radians.
*
*  4) This routine is a member of the following set:
*
*         spherical       vector       solve for
*
*       > iau_TPXES <    iau_TPXEV       xi,eta
*         iau_TPSTS      iau_TPSTV        star
*         iau_TPORS      iau_TPORV       origin
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

      DOUBLE PRECISION A, B, A0, B0, XI, ETA
      INTEGER J

      DOUBLE PRECISION TINY
      PARAMETER ( TINY = 1D-6 )

      DOUBLE PRECISION SB0, SB, CB0, CB, DA, SDA, CDA, D

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Functions of the spherical coordinates.
      SB0 = SIN(B0)
      SB = SIN(B)
      CB0 = COS(B0)
      CB = COS(B)
      DA = A - A0
      SDA = SIN(DA)
      CDA = COS(DA)

*  Reciprocal of star vector length to tangent plane.
      D = SB*SB0 + CB*CB0*CDA

*  Check for error cases.
      IF ( D .GT. TINY ) THEN
         J = 0
      ELSE IF ( D .GE. 0D0 ) THEN
         J = 1
         D = TINY
      ELSE IF ( D .GT. -TINY ) THEN
         J = 2
         D = -TINY
      ELSE
         J = 3
      END IF

*  Return the tangent plane coordinates (even in dubious cases).
      XI = CB*SDA / D
      ETA = ( SB*CB0 - CB*SB0*CDA ) / D

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
