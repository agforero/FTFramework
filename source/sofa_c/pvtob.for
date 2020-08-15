      SUBROUTINE iau_PVTOB ( ELONG, PHI, HM, XP, YP, SP, THETA, PV )
*+
*  - - - - - - - - - -
*   i a u _ P V T O B
*  - - - - - - - - - -
*
*  Position and velocity of a terrestrial observing station.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     ELONG      d       longitude (radians, east +ve, Note 1)
*     PHI        d       latitude (geodetic, radians, Note 1)
*     HM         d       height above reference ellipsoid (geodetic, m)
*     XP,YP      d       coordinates of the pole (radians, Note 2)
*     SP         d       the TIO locator s' (radians, Note 2)
*     THETA      d       Earth rotation angle (radians, Note 3)
*
*  Returned:
*     PV         d(3,2)  position/velocity vector (m, m/s, CIRS)
*
*  Notes:
*
*  1) The terrestrial coordinates are with respect to the WGS84
*     reference ellipsoid.
*
*  2) XP and YP are the coordinates (in radians) of the Celestial
*     Intermediate Pole with respect to the International Terrestrial
*     Reference System (see IERS Conventions 2003), measured along the
*     meridians 0 and 90 deg west respectively.  SP is the TIO locator
*     s', in radians, which positions the Terrestrial Intermediate
*     Origin on the equator.  For many applications, XP, YP and
*     (especially) SP can be set to zero.
*
*  3) If THETA is Greenwich apparent sidereal time instead of Earth
*     rotation angle, the result is with respect to the true equator
*     and equinox of date, i.e. with the x-axis at the equinox rather
*     than the celestial intermediate origin.
*
*  4) The velocity units are meters per UT1 second, not per SI second.
*     This is unlikely to have any practical consequences in the modern
*     era.
*
*  5) No validation is performed on the arguments.  Error cases that
*     could lead to arithmetic exceptions are trapped by the iau_GD2GC
*     routine, and the result set to zeros.
*
*  References:
*
*     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
*     IERS Technical Note No. 32, BKG (2004)
*
*     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
*     the Astronomical Almanac, 3rd ed., University Science Books
*     (2013), Section 7.4.3.3.
*
*  Called:
*     iau_GD2GC    geodetic to geocentric transformation
*     iau_POM00    polar motion matrix
*     iau_TRXP     product of transpose of r-matrix and p-vector
*
*  This revision:   2013 June 25
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION ELONG, PHI, HM, XP, YP, SP, THETA, PV(3,2)

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Earth rotation rate in radians per UT1 second
      DOUBLE PRECISION OM
      PARAMETER ( OM = 1.00273781191135448D0 * D2PI / D2S )

      INTEGER J
      DOUBLE PRECISION XYZM(3), RPM(3,3), XYZ(3), X, Y, Z, S, C


*  Geodetic to geocentric transformation (WGS84).
      CALL iau_GD2GC ( 1, ELONG, PHI, HM, XYZM, J )

*  Polar motion and TIO position.
      CALL iau_POM00 ( XP, YP, SP, RPM )
      CALL iau_TRXP ( RPM, XYZM, XYZ )
      X = XYZ(1)
      Y = XYZ(2)
      Z = XYZ(3)

*  routines of ERA.
      S = SIN(THETA)
      C = COS(THETA)

*  Position.
      PV(1,1) = C*X - S*Y
      PV(2,1) = S*X + C*Y
      PV(3,1) = Z

*  Velocity.
      PV(1,2) = OM * ( -S*X - C*Y )
      PV(2,2) = OM * (  C*X - S*Y )
      PV(3,2) = 0D0

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
