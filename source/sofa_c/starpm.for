      SUBROUTINE iau_STARPM ( RA1, DEC1, PMR1, PMD1, PX1, RV1,
     :                        EP1A, EP1B, EP2A, EP2B,
     :                        RA2, DEC2, PMR2, PMD2, PX2, RV2, J )
*+
*  - - - - - - - - - - -
*   i a u _ S T A R P M
*  - - - - - - - - - - -
*
*  Star proper motion:  update star catalog data for space motion.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     RA1      d         right ascension (radians), before
*     DEC1     d         declination (radians), before
*     PMR1     d         RA proper motion (radians/year), before
*     PMD1     d         Dec proper motion (radians/year), before
*     PX1      d         parallax (arcseconds), before
*     RV1      d         radial velocity (km/s, +ve = receding), before
*     EP1A     d         "before" epoch, part A (Note 1)
*     EP1B     d         "before" epoch, part B (Note 1)
*     EP2A     d         "after" epoch, part A (Note 1)
*     EP2B     d         "after" epoch, part B (Note 1)
*
*  Returned:
*     RA2      d         right ascension (radians), after
*     DEC2     d         declination (radians), after
*     PMR2     d         RA proper motion (radians/year), after
*     PMD2     d         Dec proper motion (radians/year), after
*     PX2      d         parallax (arcseconds), after
*     RV2      d         radial velocity (km/s, +ve = receding), after
*     J        i         status:
*                          -1 = system error (should not occur)
*                           0 = no warnings or errors
*                           1 = distance overridden (Note 6)
*                           2 = excessive velocity (Note 7)
*                           4 = solution didn't converge (Note 8)
*                        else = binary logical OR of the above warnings
*
*  Notes:
*
*  1) The starting and ending TDB epochs EP1A+EP1B and EP2A+EP2B are
*     Julian Dates, apportioned in any convenient way between the two
*     parts (A and B).  For example, JD(TDB)=2450123.7 could be
*     expressed in any of these ways, among others:
*
*             EPnA          EPnB
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) In accordance with normal star-catalog conventions, the object's
*     right ascension and declination are freed from the effects of
*     secular aberration.  The frame, which is aligned to the catalog
*     equator and equinox, is Lorentzian and centered on the SSB.
*
*     The proper motions are the rate of change of the right ascension
*     and declination at the catalog epoch and are in radians per TDB
*     Julian year.
*
*     The parallax and radial velocity are in the same frame.
*
*  3) Care is needed with units.  The star coordinates are in radians
*     and the proper motions in radians per Julian year, but the
*     parallax is in arcseconds.
*
*  4) The RA proper motion is in terms of coordinate angle, not true
*     angle.  If the catalog uses arcseconds for both RA and Dec proper
*     motions, the RA proper motion will need to be divided by cos(Dec)
*     before use.
*
*  5) Straight-line motion at constant speed, in the inertial frame,
*     is assumed.
*
*  6) An extremely small (or zero or negative) parallax is interpreted
*     to mean that the object is on the "celestial sphere", the radius
*     of which is an arbitrary (large) value (see the iau_STARPV routine
*     for the value used).  When the distance is overridden in this way,
*     the status, initially zero, has 1 added to it.
*
*  7) If the space velocity is a significant fraction of c (see the
*     constant VMAX in the routine iau_STARPV), it is arbitrarily set to
*     zero.  When this action occurs, 2 is added to the status.
*
*  8) The relativistic adjustment carried out in the iau_STARPV routine
*     involves an iterative calculation.  If the process fails to
*     converge within a set number of iterations, 4 is added to the
*     status.
*
*  Called:
*     iau_STARPV   star catalog data to space motion pv-vector
*     iau_PVU      update a pv-vector
*     iau_PDP      scalar product of two p-vectors
*     iau_PVSTAR   space motion pv-vector to star catalog data
*
*  This revision:  2017 March 16
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RA1, DEC1, PMR1, PMD1, PX1, RV1,
     :                 EP1A, EP1B, EP2A, EP2B,
     :                 RA2, DEC2, PMR2, PMD2, PX2, RV2
      INTEGER J

*  Days to seconds
      DOUBLE PRECISION D2S
      PARAMETER ( D2S = 86400D0 )

*  Speed of light (m/s)
      DOUBLE PRECISION CMPS
      PARAMETER ( CMPS = 299792458D0 )

*  Astronomical unit (m, IAU 2012)
      DOUBLE PRECISION AUM
      PARAMETER ( AUM = 149597870.7D3 )

*  Speed of light (au per day)
      DOUBLE PRECISION C
      PARAMETER ( C = D2S*CMPS/AUM )

      DOUBLE PRECISION PV1(3,2), R, TL1, DT, PV(3,2), R2, RDV, V2,
     :                 C2MV2, TL2, PV2(3,2)
      INTEGER J1, J2

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  RA,Dec etc. at the "before" epoch to space motion pv-vector.
      CALL iau_STARPV ( RA1, DEC1, PMR1, PMD1, PX1, RV1, PV1, J1 )

*  Light time when observed (days).
      CALL iau_PM ( PV1, R )
      TL1 = R / C

*  Time interval, "before" to "after" (days).
      DT = ( EP2A-EP1A ) + ( EP2B-EP1B )

*  Move star along track from the "before" observed position to the
*  "after" geometric position.
      CALL iau_PVU ( DT+TL1, PV1, PV )

*  From this geometric position, deduce the observed light time (days)
*  at the "after" epoch (with theoretically unneccessary error check).
      CALL iau_PDP ( PV(1,1), PV(1,1), R2 )
      CALL iau_PDP ( PV(1,1), PV(1,2), RDV )
      CALL iau_PDP ( PV(1,2), PV(1,2), V2 )
      C2MV2 = C*C - V2
      IF ( C2MV2 .LE. 0D0 ) THEN
         J = -1
         GO TO 9
      END IF
      TL2 = ( - RDV + SQRT(RDV*RDV + C2MV2*R2) ) / C2MV2

*  Move the position along track from the observed place at the
*  "before" epoch to the observed place at the "after" epoch.
      CALL iau_PVU ( DT + ( TL1-TL2 ), PV1, PV2 )

*  Space motion pv-vector to RA,Dec etc. at the "after" epoch.
      CALL iau_PVSTAR ( PV2, RA2, DEC2, PMR2, PMD2, PX2, RV2, J2 )

*  Return the status.
      IF ( J2 .NE. 0 ) J1 = -1
      J = J1

*  Exit.
 9    CONTINUE

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
