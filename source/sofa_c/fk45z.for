      SUBROUTINE iau_FK45Z ( R1950, D1950, BEPOCH, R2000, D2000 )
*+
*  - - - - - - - - - -
*   i a u _ F K 4 5 Z
*  - - - - - - - - - -
*
*  Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
*  proper motion in the FK5 system.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  This routine converts a star's catalog data from the old FK4
*  (Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system,
*  in such a way that the FK5 proper motion is zero.  Because such a
*  star has, in general, a non-zero proper motion in the FK4 system,
*  the routine requires the epoch at which the position in the FK4
*  system was determined.
*
*  Given:
*     R1950,D1950      d     B1950.0 FK4 RA,Dec at epoch (rad)
*     BEPOCH           d     Besselian epoch (e.g. 1979.3D0)
*
*  Returned:
*     R2000,D2000      d     J2000.0 FK5 RA,Dec (rad)
*
*  Notes:
*
*  1) The epoch BEPOCH is strictly speaking Besselian, but if a Julian
*     epoch is supplied the result will be affected only to a negligible
*     extent.
*
*  2) The method is from Appendix 2 of Aoki et al. (1983), but using the
*     constants of Seidelmann (1992).  See the routine iau_FK425 for a
*     general introduction to the FK4 to FK5 conversion.
*
*  3) Conversion from equinox B1950.0 FK4 to equinox J2000.0 FK5 only is
*     provided for.  Conversions for different starting and/or ending
*     epochs would require additional treatment for precession, proper
*     motion and E-terms.
*
*  4) In the FK4 catalog the proper motions of stars within 10 degrees
*     of the poles do not embody differential E-terms effects and
*     should, strictly speaking, be handled in a different manner from
*     stars outside these regions.  However, given the general lack of
*     homogeneity of the star data available for routine astrometry, the
*     difficulties of handling positions that may have been determined
*     from astrometric fields spanning the polar and non-polar regions,
*     the likelihood that the differential E-terms effect was not taken
*     into account when allowing for proper motion in past astrometry,
*     and the undesirability of a discontinuity in the algorithm, the
*     decision has been made in this SOFA algorithm to include the
*     effects of differential E-terms on the proper motions for all
*     stars, whether polar or not.  At epoch 2000.0, and measuring "on
*     the sky" rather than in terms of RA change, the errors resulting
*     from this simplification are less than 1 milliarcsecond in
*     position and 1 milliarcsecond per century in proper motion.
*
*  References:
*
*     Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
*     FK4-based positions of stars to epoch J2000.0 positions in
*     accordance with the new IAU resolutions".  Astron.Astrophys.
*     128, 263-267.
*
*     Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
*     Astronomical Almanac", ISBN 0-935702-68-7.
*
*  Called:
*     iau_ANP      normalize angle into range 0 to 2pi
*     iau_C2S      p-vector to spherical
*     iau_EPB2JD   Besselian epoch to Julian date
*     iau_EPJ      Julian date to Julian epoch
*     iau_PDP      scalar product of two p-vectors
*     iau_PMP      p-vector minus p-vector
*     iau_PPSP     p-vector plus scaled p-vector
*     iau_PVU      update a pv-vector
*     iau_S2C      spherical to p-vector
*
*  This revision:   2018 January 11
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION R1950, D1950, BEPOCH, R2000, D2000

*  2pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Radians per year to arcsec per century
      DOUBLE PRECISION PMF
      PARAMETER ( PMF = 100D0*60D0*60D0*360D0/D2PI )

*  Position and position+velocity vectors
      DOUBLE PRECISION R0(3), A1(3), P1(3), P2(3), PV1(3,2), PV2(3,2)

*  Miscellaneous
      DOUBLE PRECISION W, DJM0, DJM
      INTEGER K, J, I

*  Functions
      DOUBLE PRECISION iau_EPJ, iau_ANP

*
*  CANONICAL CONSTANTS
*

*  Vectors A and Adot (Seidelmann 3.591-2)
      DOUBLE PRECISION A(3), AD(3)
      DATA A,AD/ -1.62557D-6,  -0.31919D-6, -0.13843D-6,
     :           +1.245D-3,    -1.580D-3,   -0.659D-3 /

*  3x2 matrix of p-vectors (cf. Seidelmann 3.591-4, matrix M)
      DOUBLE PRECISION EM(3,3,2)
      DATA EM /
     :   +0.9999256782D0,     -0.0111820611D0,     -0.0048579477D0,
     :   +0.0111820610D0,     +0.9999374784D0,     -0.0000271765D0,
     :   +0.0048579479D0,     -0.0000271474D0,     +0.9999881997D0,
     :   -0.000551D0,         -0.238565D0,         +0.435739D0,
     :   +0.238514D0,         -0.002667D0,         -0.008541D0,
     :   -0.435623D0,         +0.012254D0,         +0.002117D0
     :        /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Spherical coordinates to p-vector.
      CALL iau_S2C ( R1950, D1950, R0 )

*  Adjust p-vector A to give zero proper motion in FK5.
      W = ( BEPOCH - 1950D0 ) / PMF
      CALL iau_PPSP ( A, W, AD, A1 )

*  Remove E-terms.
      CALL iau_PDP ( R0, A1, W )
      CALL iau_PPSP ( A1, -W, R0, P1 )
      CALL iau_PMP ( R0, P1, P2 )

*  Convert to Fricke system pv-vector (cf. Seidelmann 3.591-3).
      DO 3 K = 1,2
         DO 2 J=1,3
            W = 0D0
            DO 1 I=1,3
               W = W + EM(I,J,K)*P2(I)
 1          CONTINUE
            PV1(J,K) = W
 2       CONTINUE
 3    CONTINUE

*  Allow for fictitious proper motion.
      CALL iau_EPB2JD ( BEPOCH, DJM0, DJM )
      W = ( iau_EPJ(DJM0,DJM) - 2000D0 ) / PMF
      CALL iau_PVU ( W, PV1, PV2 )

*  Revert to spherical coordinates.
      CALL iau_C2S ( PV2, W, D2000 )
      R2000 = iau_ANP ( W )

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
