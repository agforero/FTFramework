      SUBROUTINE iau_PMPX ( RC, DC, PR, PD, PX, RV, PMT, POB, PCO )
*+
*  - - - - - - - - -
*   i a u _ P M P X
*  - - - - - - - - -
*
*  Proper motion and parallax.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     RC,DC    d      ICRS RA,Dec at catalog epoch (radians)
*     PR       d      RA proper motion (radians/year; Note 1)
*     PD       d      Dec proper motion (radians/year)
*     PX       d      parallax (arcsec)
*     RV       d      radial velocity (km/s, +ve if receding)
*     PMT      d      proper motion time interval (SSB, Julian years)
*     POB      d(3)   SSB to observer vector (au)
*
*  Returned:
*     PCO      d(3)   coordinate direction (BCRS unit vector)
*
*  Notes:
*
*  1) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
*
*  2) The proper motion time interval is for when the starlight
*     reaches the solar system barycenter.
*
*  3) To avoid the need for iteration, the Roemer effect (i.e. the
*     small annual modulation of the proper motion coming from the
*     changing light time) is applied approximately, using the
*     direction of the star at the catalog epoch.
*
*  References:
*
*     1984 Astronomical Almanac, pp B39-B41.
*
*     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
*     the Astronomical Almanac, 3rd ed., University Science Books
*     (2013), Section 7.2.
*
*  Called:
*     iau_PDP      scalar product of two p-vectors
*     iau_PN       decompose p-vector into modulus and direction
*
*  This revision:   2017 March 11
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RC, DC, PR, PD, PX, RV, PMT, POB(3), PCO(3)

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

*  Days per Julian year
      DOUBLE PRECISION DJY
      PARAMETER ( DJY = 365.25D0 )

*  Days per Julian millennium
      DOUBLE PRECISION DJM
      PARAMETER ( DJM = 365250D0 )

*  Speed of light (m/s)
      DOUBLE PRECISION CMPS
      PARAMETER ( CMPS = 299792458D0 )

*  Astronomical unit (m, IAU 2012)
      DOUBLE PRECISION AUM
      PARAMETER ( AUM = 149597870.7D3 )

*  Light time for 1 au, Julian years
      DOUBLE PRECISION AULTY
      PARAMETER ( AULTY = AUM/CMPS/D2S/DJY )

*  Km/s to au/year
      DOUBLE PRECISION VF
      PARAMETER ( VF = D2S*DJM/AUM )

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

      INTEGER I
      DOUBLE PRECISION SR, CR, SD, CD, X, Y, Z, P(3), PDB,
     :                 DT, PXR, W, PDZ, PM(3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Catalog spherical coordinates to unit vector (and useful functions).
      SR = SIN(RC)
      CR = COS(RC)
      SD = SIN(DC)
      CD = COS(DC)
      X = CR*CD
      Y = SR*CD
      Z = SD
      P(1) = X
      P(2) = Y
      P(3) = Z

*  Component of observer vector in star direction.
      CALL iau_PDP ( P, POB, PDB )

*  Proper motion time interval (y), including Roemer effect.
      DT = PMT + PDB*AULTY

*  Space motion (radians per year).
      PXR = PX*DAS2R
      W = VF*RV*PXR
      PDZ = PD*Z
      PM(1) = - PR*Y - PDZ*CR + W*X
      PM(2) =   PR*X - PDZ*SR + W*Y
      PM(3) =          PD*CD  + W*Z

*  Coordinate direction of star (unit vector, BCRS).
      DO 1 I=1,3
         P(I) = P(I) + DT*PM(I) - PXR*POB(I)
 1    CONTINUE
      CALL iau_PN ( P, W, PCO )

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
