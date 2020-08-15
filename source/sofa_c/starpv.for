      SUBROUTINE iau_STARPV ( RA, DEC, PMR, PMD, PX, RV, PV, J )
*+
*  - - - - - - - - - - -
*   i a u _ S T A R P V
*  - - - - - - - - - - -
*
*  Convert star catalog coordinates to position+velocity vector.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given (Note 1):
*     RA       d         right ascension (radians)
*     DEC      d         declination (radians)
*     PMR      d         RA proper motion (radians/year)
*     PMD      d         Dec proper motion (radians/year)
*     PX       d         parallax (arcseconds)
*     RV       d         radial velocity (km/s, positive = receding)
*
*  Returned (Note 2):
*     PV       d(3,2)    pv-vector (au, au/day)
*     J        i         status:
*                           0 = no warnings
*                           1 = distance overridden (Note 6)
*                           2 = excessive velocity (Note 7)
*                           4 = solution didn't converge (Note 8)
*                        else = binary logical OR of the above
*
*  Notes:
*
*  1) The star data accepted by this routine are "observables" for an
*     imaginary observer at the solar-system barycenter.  Proper motion
*     and radial velocity are, strictly, in terms of barycentric
*     coordinate time, TCB.  For most practical applications, it is
*     permissible to neglect the distinction between TCB and ordinary
*     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
*     limited by the intrinsic accuracy of the proper-motion and radial-
*     velocity data;  moreover, the pv-vector is likely to be merely an
*     intermediate result, so that a change of time unit would cancel
*     out overall.
*
*     In accordance with normal star-catalog conventions, the object's
*     right ascension and declination are freed from the effects of
*     secular aberration.  The frame, which is aligned to the catalog
*     equator and equinox, is Lorentzian and centered on the SSB.
*
*  2) The resulting position and velocity pv-vector is with respect to
*     the same frame and, like the catalog coordinates, is freed from
*     the effects of secular aberration.  Should the "coordinate
*     direction", where the object was located at the catalog epoch, be
*     required, it may be obtained by calculating the magnitude of the
*     position vector PV(1-3,1) dividing by the speed of light in au/day
*     to give the light-time, and then multiplying the space velocity
*     PV(1-3,2) by this light-time and adding the result to PV(1-3,1).
*
*     Summarizing, the pv-vector returned is for most stars almost
*     identical to the result of applying the standard geometrical
*     "space motion" transformation.  The differences, which are the
*     subject of the Stumpff paper referenced below, are:
*
*     (i) In stars with significant radial velocity and proper motion,
*     the constantly changing light-time distorts the apparent proper
*     motion.  Note that this is a classical, not a relativistic,
*     effect.
*
*     (ii) The transformation complies with special relativity.
*
*  3) Care is needed with units.  The star coordinates are in radians
*     and the proper motions in radians per Julian year, but the
*     parallax is in arcseconds; the radial velocity is in km/s, but
*     the pv-vector result is in au and au/day.
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
*     of which is an arbitrary (large) value (see the constant PXMIN).
*     When the distance is overridden in this way, the status, initially
*     zero, has 1 added to it.
*
*  7) If the space velocity is a significant fraction of c (see the
*     constant VMAX), it is arbitrarily set to zero.  When this action
*     occurs, 2 is added to the status.
*
*  8) The relativistic adjustment involves an iterative calculation.
*     If the process fails to converge within a set number (IMAX) of
*     iterations, 4 is added to the status.
*
*  9) The inverse transformation is performed by the routine iau_PVSTAR.
*
*  Called:
*     iau_S2PV     spherical coordinates to pv-vector
*     iau_PM       modulus of p-vector
*     iau_ZP       zero p-vector
*     iau_PN       decompose p-vector into modulus and direction
*     iau_PDP      scalar product of two p-vectors
*     iau_SXP      multiply p-vector by scalar
*     iau_PMP      p-vector minus p-vector
*     iau_PPP      p-vector plus p-vector
*
*  Reference:
*
*     Stumpff, P., Astron.Astrophys. 144, 232-240 (1985).
*
*  This revision:  2017 March 16
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RA, DEC, PMR, PMD, PX, RV, PV(3,2)
      INTEGER J

*  Smallest allowed parallax
      DOUBLE PRECISION PXMIN
      PARAMETER ( PXMIN = 1D-7 )

*  Largest allowed speed (fraction of c)
      DOUBLE PRECISION VMAX
      PARAMETER ( VMAX = 0.5D0 )

*  Julian years to days
      DOUBLE PRECISION Y2D
      PARAMETER ( Y2D = 365.25D0 )

*  Radians to arcseconds
      DOUBLE PRECISION DR2AS
      PARAMETER ( DR2AS = 206264.8062470963551564734D0 )

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

*  Maximum number of iterations for relativistic solution
      INTEGER I,IMAX
      PARAMETER ( IMAX = 100 )

      INTEGER IWARN
      DOUBLE PRECISION W, R, RD, RAD, DECD, V, X(3), USR(3), UST(3),
     :                 VSR, VST, BETST, BETSR, BETT, BETR, OD, ODEL,
     :                 DD, DDEL, ODD, ODDEL, D, DEL, UR(3), UT(3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Distance (au).
      IF ( PX.GE.PXMIN ) THEN
         W = PX
         IWARN = 0
      ELSE
         W = PXMIN
         IWARN = 1
      END IF
      R = DR2AS / W

*  Radial velocity (au/day).
      RD = D2S * RV * 1D3 / AUM

*  Proper motion (radian/day).
      RAD = PMR / Y2D
      DECD = PMD / Y2D

*  To pv-vector (au,au/day).
      CALL iau_S2PV ( RA, DEC, R, RAD, DECD, RD, PV )

*  If excessive velocity, arbitrarily set it to zero.
      CALL iau_PM ( PV(1,2), V )
      IF ( V/C .GT. VMAX ) THEN
         CALL iau_ZP ( PV(1,2) )
         IWARN = IWARN + 2
      END IF

*  Isolate the radial component of the velocity (au/day).
      CALL iau_PN ( PV(1,1), W, X )
      CALL iau_PDP ( X, PV(1,2), VSR )
      CALL iau_SXP ( VSR, X, USR )

*  Isolate the transverse component of the velocity (au/day).
      CALL iau_PMP ( PV(1,2), USR, UST )
      CALL iau_PM ( UST, VST )

*  Special-relativity dimensionless parameters.
      BETSR = VSR / C
      BETST = VST / C

*  Determine the inertial-to-observed relativistic correction terms.
      OD = 0D0
      ODEL = 0D0
      ODD = 0D0
      ODDEL = 0D0
      BETT = BETST
      BETR = BETSR
      DO 1 I=1,IMAX
         D = 1D0 + BETR
         W = BETR*BETR + BETT*BETT
         DEL = - W / ( SQRT(1D0-W) + 1D0 )
         BETR = D*BETSR + DEL
         BETT = D*BETST
         IF ( I .GT. 1 ) THEN
            DD = ABS(D-OD)
            DDEL = ABS(DEL-ODEL)
            IF ( I.GT.2 .AND.
     :           DD.GE.ODD .AND.
     :           DDEL.GE.ODDEL ) GO TO 2
            IF ( I .GE. IMAX ) IWARN = IWARN + 4
            ODD = DD
            ODDEL = DDEL
         END IF
         OD = D
         ODEL = DEL
 1    CONTINUE
 2    CONTINUE

*  Replace observed radial velocity with inertial value.
      IF ( BETSR .NE. 0D0 ) THEN
         W = D + DEL/BETSR
      ELSE
         W = 1D0
      END IF
      CALL iau_SXP ( W, USR, UR )

*  Replace observed tangential velocity with inertial value.
      CALL iau_SXP ( D, UST, UT )

*  Combine the two to obtain the inertial space velocity.
      CALL iau_PPP ( UR, UT, PV(1,2) )

*  Return the status.
      J = IWARN

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
