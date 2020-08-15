      SUBROUTINE iau_H2FK5 ( RH, DH, DRH, DDH, PXH, RVH,
     :                       R5, D5, DR5, DD5, PX5, RV5 )
*+
*  - - - - - - - - - -
*   i a u _ H 2 F K 5
*  - - - - - - - - - -
*
*  Transform Hipparcos star data into the FK5 (J2000.0) system.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given (all Hipparcos, epoch J2000.0):
*     RH        d      RA (radians)
*     DH        d      Dec (radians)
*     DRH       d      proper motion in RA (dRA/dt, rad/Jyear)
*     DDH       d      proper motion in Dec (dDec/dt, rad/Jyear)
*     PXH       d      parallax (arcsec)
*     RVH       d      radial velocity (km/s, positive = receding)
*
*  Returned (all FK5, equinox J2000.0, epoch J2000.0):
*     R5        d      RA (radians)
*     D5        d      Dec (radians)
*     DR5       d      proper motion in RA (dRA/dt, rad/Jyear)
*     DD5       d      proper motion in Dec (dDec/dt, rad/Jyear)
*     PX5       d      parallax (arcsec)
*     RV5       d      radial velocity (km/s, positive = receding)
*
*  Notes:
*
*  1) This routine transforms Hipparcos star positions and proper
*     motions into FK5 J2000.0.
*
*  2) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
*     and are per year rather than per century.
*
*  3) The FK5 to Hipparcos transformation is modeled as a pure rotation
*     and spin;  zonal errors in the FK5 catalog are not taken into
*     account.
*
*  4) See also iau_FK52H, iau_FK5HZ, iau_HFK5Z.
*
*  Called:
*     iau_STARPV   star catalog data to space motion pv-vector
*     iau_FK5HIP   FK5 to Hipparcos rotation and spin
*     iau_RV2M     r-vector to r-matrix
*     iau_RXP      product of r-matrix and p-vector
*     iau_TRXP     product of transpose of r-matrix and p-vector
*     iau_PXP      vector product of two p-vectors
*     iau_PMP      p-vector minus p-vector
*     iau_PVSTAR   space motion pv-vector to star catalog data
*
*  Reference:
*
*     F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
*
*  This revision:  2017 October 12
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RH, DH, DRH, DDH, PXH, RVH,
     :                 R5, D5, DR5, DD5, PX5, RV5

      DOUBLE PRECISION PVH(3,2), R5H(3,3), S5H(3), SH(3), WXP(3),
     :                 VV(3), PV5(3,2)
      INTEGER J, I

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Hipparcos barycentric position/velocity pv-vector (normalized).
      CALL iau_STARPV ( RH, DH, DRH, DDH, PXH, RVH, PVH, J )

*  FK5 to Hipparcos orientation matrix and spin vector.
      CALL iau_FK5HIP ( R5H, S5H )

*  Make spin units per day instead of per year.
      DO 1 I=1,3
         S5H(I) = S5H(I) / 365.25D0
 1    CONTINUE

*  Orient the spin into the Hipparcos system.
      CALL iau_RXP ( R5H, S5H, SH )

*  De-orient the Hipparcos position into the FK5 system.
      CALL iau_TRXP ( R5H, PVH(1,1), PV5(1,1) )

*  Apply spin to the position giving an extra space motion component.
      CALL iau_PXP ( PVH(1,1), SH, WXP )

*  Subtract this component from the Hipparcos space motion.
      CALL iau_PMP ( PVH(1,2), WXP, VV )

*  De-orient the Hipparcos space motion into the FK5 system.
      CALL iau_TRXP ( R5H, VV, PV5(1,2) )

*  FK5 pv-vector to spherical.
      CALL iau_PVSTAR ( PV5, R5, D5, DR5, DD5, PX5, RV5, J )

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
