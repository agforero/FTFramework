      SUBROUTINE iau_GC2GDE ( A, F, XYZ, ELONG, PHI, HEIGHT, J )
*+
*  - - - - - - - - - - -
*   i a u _ G C 2 G D E
*  - - - - - - - - - - -
*
*  Transform geocentric coordinates to geodetic for a reference
*  ellipsoid of specified form.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     A         d      equatorial radius (Notes 2,4)
*     F         d      flattening (Note 3)
*     XYZ       d(3)   geocentric vector (Note 4)
*
*  Returned:
*     ELONG     d      longitude (radians, east +ve)
*     PHI       d      latitude (geodetic, radians)
*     HEIGHT    d      height above ellipsoid (geodetic, Note 4)
*     J         i      status:  0 = OK
*                              -1 = illegal F
*                              -2 = illegal A
*
*  Notes:
*
*  1) This routine is closely based on the GCONV2H subroutine by
*     Toshio Fukushima (see reference).
*
*  2) The equatorial radius, A, can be in any units, but meters is
*     the conventional choice.
*
*  3) The flattening, F, is (for the Earth) a value around 0.00335,
*     i.e. around 1/298.
*
*  4) The equatorial radius, A, and the geocentric vector, XYZ,
*     must be given in the same units, and determine the units of
*     the returned height, HEIGHT.
*
*  5) If an error occurs (J<0), ELONG, PHI and HEIGHT are unchanged.
*
*  6) The inverse transformation is performed in the routine iau_GD2GCE.
*
*  7) The transformation for a standard ellipsoid (such as WGS84) can
*     more conveniently be performed by calling iau_GC2GD, which uses a
*     numerical code (1 for WGS84) to identify the required A and F
*     values.
*
*  Reference:
*
*     Fukushima, T., "Transformation from Cartesian to geodetic
*     coordinates accelerated by Halley's method", J.Geodesy (2006)
*     79: 689-693
*
*  This revision:  2014 November 7
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION A, F, XYZ(3), ELONG, PHI, HEIGHT
      INTEGER J

*  Pi
      DOUBLE PRECISION DPI
      PARAMETER ( DPI = 3.141592653589793238462643D0 )

      DOUBLE PRECISION AEPS2, E2, E4T, EC2, EC, B, X, Y, Z, P2, ABSZ, P,
     :                 S0, PN, ZC, C0, C02, C03, S02, S03, A02, A0, A03,
     :                 D0, F0, B0, S1, CC, S12, CC2


*  -------------
*  Preliminaries
*  -------------

*  Validate ellipsoid parameters.
      IF ( F.LT.0D0 .OR. F.GE.1D0 ) THEN
         J = -1
         GO TO 9999
      ELSE IF ( A .LE. 0D0 ) THEN
         J = -2
         GO TO 9999
      END IF

*  Functions of ellipsoid parameters (with further validation of F).
      AEPS2 = A*A*1D-32
      E2 = (2D0-F)*F
      E4T = E2*E2*1.5D0
      EC2 = 1D0-E2
      IF ( EC2 .LE. 0D0 ) THEN
         J = -1
         GO TO 9999
      END IF
      EC = SQRT(EC2)
      B = A*EC

*  Cartesian components.
      X = XYZ(1)
      Y = XYZ(2)
      Z = XYZ(3)

*  Distance from polar axis squared.
      P2 = X*X + Y*Y

*  Longitude.
      IF ( P2.GT.0D0 ) THEN
         ELONG = ATAN2(Y,X)
      ELSE
         ELONG = 0D0
      END IF

*  Unsigned z-coordinate.
      ABSZ = ABS(Z)

*  Proceed unless polar case.
      IF ( P2.GT.AEPS2 ) THEN

*     Distance from polar axis.
         P = SQRT(P2)

*     Normalization.
         S0 = ABSZ/A
         PN = P/A
         ZC = EC*S0

*     Prepare Newton correction factors.
         C0 = EC*PN
         C02 = C0*C0
         C03 = C02*C0
         S02 = S0*S0
         S03 = S02*S0
         A02 = C02+S02
         A0 = SQRT(A02)
         A03 = A02*A0
         D0 = ZC*A03 + E2*S03
         F0 = PN*A03 - E2*C03

*     Prepare Halley correction factor.
         B0 = E4T*S02*C02*PN*(A0-EC)
         S1 = D0*F0 - B0*S0
         CC = EC*(F0*F0-B0*C0)

*     Evaluate latitude and height.
         PHI = ATAN(S1/CC)
         S12 = S1*S1
         CC2 = CC*CC
         HEIGHT = (P*CC+ABSZ*S1-A*SQRT(EC2*S12+CC2))/SQRT(S12+CC2)
      ELSE

*     Exception: pole.
         PHI = DPI/2D0
         HEIGHT = ABSZ-B
      END IF

*  Restore sign of latitude.
      IF ( Z.LT.0D0 ) PHI = -PHI

*  OK status.
      J = 0

*  Finished.
 9999 CONTINUE

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
