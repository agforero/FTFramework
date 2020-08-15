      SUBROUTINE iau_APER13 ( UT11, UT12, ASTROM )
*+
*  - - - - - - - - - - -
*   i a u _ A P E R 1 3
*  - - - - - - - - - - -
*
*  In the star-independent astrometry parameters, update only the
*  Earth rotation angle.  The caller provides UT1 (n.b. not UTC).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     UT11     d      UT1 as a 2-part...
*     UT12     d      ...Julian Date (Note 1)
*     ASTROM   d(30)  star-independent astrometry parameters:
*               (1)      not used
*               (2-4)    not used
*               (5-7)    not used
*               (8)      not used
*               (9-11)   not used
*               (12)     not used
*               (13-21)  not used
*               (22)     longitude + s' (radians)
*               (23)     not used
*               (24)     not used
*               (25)     not used
*               (26)     not used
*               (27)     not used
*               (28)     not used
*               (29)     not used
*               (30)     not used
*
*  Returned:
*     ASTROM   d(30)  star-independent astrometry parameters:
*               (1)      unchanged
*               (2-4)    unchanged
*               (5-7)    unchanged
*               (8)      unchanged
*               (9-11)   unchanged
*               (12)     unchanged
*               (13-21)  unchanged
*               (22)     unchanged
*               (23)     unchanged
*               (24)     unchanged
*               (25)     unchanged
*               (26)     unchanged
*               (27)     unchanged
*               (28)     "local" Earth rotation angle (radians)
*               (29)     unchanged
*               (30)     unchanged
*
*  Notes:
*
*  1) The UT1 date (n.b. not UTC) UT11+UT12 is a Julian Date,
*     apportioned in any convenient way between the arguments UT11 and
*     UT12.  For example, JD(UT1)=2450123.7 could be expressed in any
*     of these ways, among others:
*
*            UT11           UT12
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in cases
*     where the loss of several decimal digits of resolution is
*     acceptable.  The J2000 and MJD methods are good compromises
*     between resolution and convenience.  The date & time method is
*     best matched to the algorithm used:  maximum precision is
*     delivered when the UT11 argument is for 0hrs UT1 on the day in
*     question and the UT12 argument lies in the range 0 to 1, or vice
*     versa.
*
*  2) If the caller wishes to provide the Earth rotation angle itself,
*     the routine iau_APER can be used instead.  One use of this
*     technique is to substitute Greenwich apparent sidereal time and
*     thereby to support equinox based transformations directly.
*
*  3) This is one of several routines that inserts into the ASTROM
*     array star-independent parameters needed for the chain of
*     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
*
*     The various routines support different classes of observer and
*     portions of the transformation chain:
*
*            routines           observer        transformation
*
*        iau_APCG iau_APCG13    geocentric      ICRS <-> GCRS
*        iau_APCI iau_APCI13    terrestrial     ICRS <-> CIRS
*        iau_APCO iau_APCO13    terrestrial     ICRS <-> observed
*        iau_APCS iau_APCS13    space           ICRS <-> GCRS
*        iau_APER iau_APER13    terrestrial     update Earth rotation
*        iau_APIO iau_APIO13    terrestrial     CIRS <-> observed
*
*     Those with names ending in "13" use contemporary SOFA models to
*     compute the various ephemerides.  The others accept ephemerides
*     supplied by the caller.
*
*     The transformation from ICRS to GCRS covers space motion,
*     parallax, light deflection, and aberration.  From GCRS to CIRS
*     comprises frame bias and precession-nutation.  From CIRS to
*     observed takes account of Earth rotation, polar motion, diurnal
*     aberration and parallax (unless subsumed into the ICRS <-> GCRS
*     transformation), and atmospheric refraction.
*
*  Called:
*     iau_APER     astrometry parameters: update ERA
*     iau_ERA00    Earth rotation angle, IAU 2000
*
*  This revision:   2017 March 12
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION UT11, UT12, ASTROM(30)

      DOUBLE PRECISION iau_ERA00

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      CALL iau_APER ( iau_ERA00 ( UT11, UT12 ), ASTROM )

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
