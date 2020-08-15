      SUBROUTINE iau_FK5HZ ( R5, D5, DATE1, DATE2, RH, DH )
*+
*  - - - - - - - - - -
*   i a u _ F K 5 H Z
*  - - - - - - - - - -
*
*  Transform an FK5 (J2000.0) star position into the system of the
*  Hipparcos catalogue, assuming zero Hipparcos proper motion.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     R5              d      FK5 RA (radians), equinox J2000.0, at date
*     D5              d      FK5 Dec (radians), equinox J2000.0, at date
*     DATE1,DATE2     d      TDB date (Notes 1,2)
*
*  Returned:
*     RH              d      Hipparcos RA (radians)
*     DH              d      Hipparcos Dec (radians)
*
*  Notes:
*
*  1) This routine converts a star position from the FK5 system to
*     the Hipparcos system, in such a way that the Hipparcos proper
*     motion is zero.  Because such a star has, in general, a non-zero
*     proper motion in the FK5 system, the routine requires the date
*     at which the position in the FK5 system was determined.
*
*  2) The TDB date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TDB)=2450123.7 could be expressed in any of these ways, among
*     others:
*
*            DATE1          DATE2
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
*  3) The FK5 to Hipparcos transformation is modeled as a pure
*     rotation and spin;  zonal errors in the FK5 catalogue are
*     not taken into account.
*
*  4) The position returned by this routine is in the Hipparcos
*     reference system but at date DATE1+DATE2.
*
*  5) See also iau_FK52H, iau_H2FK5, iau_HFK5Z.
*
*  Called:
*     iau_S2C      spherical coordinates to unit vector
*     iau_FK5HIP   FK5 to Hipparcos rotation and spin
*     iau_SXP      multiply p-vector by scalar
*     iau_RV2M     r-vector to r-matrix
*     iau_TRXP     product of transpose of r-matrix and p-vector
*     iau_PXP      vector product of two p-vectors
*     iau_C2S      p-vector to spherical
*     iau_ANP      normalize angle into range 0 to 2pi
*
*  Reference:
*
*     F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).
*
*  This revision:  2012 September 5
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R5, D5, DATE1, DATE2, RH, DH

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian year
      DOUBLE PRECISION DJY
      PARAMETER ( DJY = 365.25D0 )

      DOUBLE PRECISION T, P5E(3), R5H(3,3), S5H(3), VST(3), RST(3,3),
     :                 P5(3), PH(3), W

      DOUBLE PRECISION iau_ANP

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval from given date to fundamental epoch J2000.0 (JY).
      T = - ( ( DATE1-DJ00 ) + DATE2 ) / DJY

*  FK5 barycentric position vector.
      CALL iau_S2C ( R5, D5, P5E )

*  FK5 to Hipparcos orientation matrix and spin vector.
      CALL iau_FK5HIP ( R5H, S5H )

*  Accumulated Hipparcos wrt FK5 spin over that interval.
      CALL iau_SXP ( T, S5H, VST )

*  Express the accumulated spin as a rotation matrix.
      CALL iau_RV2M ( VST, RST )

*  Derotate the vector's FK5 axes back to date.
      CALL iau_TRXP ( RST, P5E, P5 )

*  Rotate the vector into the Hipparcos system.
      CALL iau_RXP ( R5H, P5, PH )

*  Hipparcos vector to spherical.
      CALL iau_C2S ( PH, W, DH )
      RH = iau_ANP ( W )

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
