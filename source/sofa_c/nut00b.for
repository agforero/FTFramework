      SUBROUTINE iau_NUT00B ( DATE1, DATE2, DPSI, DEPS )
*+
*  - - - - - - - - - - -
*   i a u _ N U T 0 0 B
*  - - - - - - - - - - -
*
*  Nutation, IAU 2000B model.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2   d   TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     DPSI,DEPS     d   nutation, luni-solar + planetary (Note 2)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in cases
*     where the loss of several decimal digits of resolution is
*     acceptable.  The J2000 method is best matched to the way the
*     argument is handled internally and will deliver the optimum
*     resolution.  The MJD method and the date & time methods are both
*     good compromises between resolution and convenience.
*
*  2) The nutation components in longitude and obliquity are in radians
*     and with respect to the equinox and ecliptic of date.  The
*     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
*     value of 84381.448 arcsec.  (The errors that result from using
*     this routine with the IAU 2006 value of 84381.406 arcsec can be
*     neglected.)
*
*     The nutation model consists only of luni-solar terms, but includes
*     also a fixed offset which compensates for certain long-period
*     planetary terms (Note 7).
*
*  3) This routine is an implementation of the IAU 2000B abridged
*     nutation model formally adopted by the IAU General Assembly in
*     2000.  The routine computes the MHB_2000_SHORT luni-solar nutation
*     series (Luzum 2001), but without the associated corrections for
*     the precession rate adjustments and the offset between the GCRS
*     and J2000.0 mean poles.
*
*  4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
*     terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
*     77 terms, plus additional simplifications, yet still delivers
*     results of 1 mas accuracy at present epochs.  This combination of
*     accuracy and size makes the IAU 2000B abridged nutation model
*     suitable for most practical applications.
*
*     The routine delivers a pole accurate to 1 mas from 1900 to 2100
*     (usually better than 1 mas, very occasionally just outside 1 mas).
*     The full IAU 2000A model, which is implemented in the routine
*     iau_NUT00A (q.v.), delivers considerably greater accuracy at
*     current epochs;  however, to realize this improved accuracy,
*     corrections for the essentially unpredictable free-core-nutation
*     (FCN) must also be included.
*
*  5) The present routine provides classical nutation.  The
*     MHB_2000_SHORT algorithm, from which it is adapted, deals also
*     with (i) the offsets between the GCRS and mean poles and (ii) the
*     adjustments in longitude and obliquity due to the changed
*     precession rates.  These additional functions, namely frame bias
*     and precession adjustments, are supported by the SOFA routines
*     iau_BI00 and iau_PR00.
*
*  6) The MHB_2000_SHORT algorithm also provides "total" nutations,
*     comprising the arithmetic sum of the frame bias, precession
*     adjustments, and nutation (luni-solar + planetary).  These total
*     nutations can be used in combination with an existing IAU 1976
*     precession implementation, such as iau_PMAT76, to deliver GCRS-to-
*     true predictions of mas accuracy at current epochs.  However, for
*     symmetry with the iau_NUT00A routine (q.v. for the reasons), the
*     SOFA routines do not generate the "total nutations" directly.
*     Should they be required, they could of course easily be generated
*     by calling iau_BI00, iau_PR00 and the present routine and adding
*     the results.
*
*  7) The IAU 2000B model includes "planetary bias" terms that are fixed
*     in size but compensate for long-period nutations.  The amplitudes
*     quoted in McCarthy & Luzum (2003), namely Dpsi = -1.5835 mas and
*     Depsilon = +1.6339 mas, are optimized for the "total nutations"
*     method described in Note 6.  The Luzum (2001) values used in this
*     SOFA implementation, namely -0.135 mas and +0.388 mas, are
*     optimized for the "rigorous" method, where frame bias, precession
*     and nutation are applied separately and in that order.  During the
*     interval 1995-2050, the SOFA implementation delivers a maximum
*     error of 1.001 mas (not including FCN).
*
*  References:
*
*     Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions
*     for the precession quantities based upon the IAU /1976/ system of
*     astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)
*
*     Luzum, B., private communication, 2001 (Fortran code
*     MHB_2000_SHORT)
*
*     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
*     precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.
*     85, 37-49 (2003)
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
*
*  This revision:  2009 December 15
*
*  SOFA release 2020-07-21
*
*  Copyright (C) 2020 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Milliarcseconds to radians
      DOUBLE PRECISION DMAS2R
      PARAMETER ( DMAS2R = DAS2R / 1D3 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Miscellaneous
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 DPSILS, DEPSLS, DPSIPL, DEPSPL
      INTEGER I, J

*  -------------------------
*  Luni-Solar nutation model
*  -------------------------

*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 77 )

*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)

*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)

*  ---------------------------------------
*  Fixed offset in lieu of planetary terms (radians)
*  ---------------------------------------

      DOUBLE PRECISION DPPLAN, DEPLAN
      PARAMETER ( DPPLAN = - 0.135D0 * DMAS2R,
     :            DEPLAN = + 0.388D0 * DMAS2R )

*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------

*
*  Luni-Solar argument multipliers:
*
*               L     L'    F     D     Om

      DATA ( ( NALS(I,J), I=1,5 ), J= 1,10 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=11,20 ) /
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=21,30 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=31,40 ) /
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=41,50 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=51,60 ) /
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=61,70 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=71,77 ) /
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2 /

*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*

      DATA ( ( CLS(I,J), I=1,6 ), J= 1,10 ) /
     : -172064161D0, -174666D0,  33386D0, 92052331D0,  9086D0, 15377D0,
     :  -13170906D0,   -1675D0, -13696D0,  5730336D0, -3015D0, -4587D0,
     :   -2276413D0,    -234D0,   2796D0,   978459D0,  -485D0,  1374D0,
     :    2074554D0,     207D0,   -698D0,  -897492D0,   470D0,  -291D0,
     :    1475877D0,   -3633D0,  11817D0,    73871D0,  -184D0, -1924D0,
     :    -516821D0,    1226D0,   -524D0,   224386D0,  -677D0,  -174D0,
     :     711159D0,      73D0,   -872D0,    -6750D0,     0D0,   358D0,
     :    -387298D0,    -367D0,    380D0,   200728D0,    18D0,   318D0,
     :    -301461D0,     -36D0,    816D0,   129025D0,   -63D0,   367D0,
     :     215829D0,    -494D0,    111D0,   -95929D0,   299D0,   132D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=11,20 ) /
     :     128227D0,     137D0,    181D0,   -68982D0,    -9D0,    39D0,
     :     123457D0,      11D0,     19D0,   -53311D0,    32D0,    -4D0,
     :     156994D0,      10D0,   -168D0,    -1235D0,     0D0,    82D0,
     :      63110D0,      63D0,     27D0,   -33228D0,     0D0,    -9D0,
     :     -57976D0,     -63D0,   -189D0,    31429D0,     0D0,   -75D0,
     :     -59641D0,     -11D0,    149D0,    25543D0,   -11D0,    66D0,
     :     -51613D0,     -42D0,    129D0,    26366D0,     0D0,    78D0,
     :      45893D0,      50D0,     31D0,   -24236D0,   -10D0,    20D0,
     :      63384D0,      11D0,   -150D0,    -1220D0,     0D0,    29D0,
     :     -38571D0,      -1D0,    158D0,    16452D0,   -11D0,    68D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=21,30 ) /
     :      32481D0,       0D0,      0D0,   -13870D0,     0D0,     0D0,
     :     -47722D0,       0D0,    -18D0,      477D0,     0D0,   -25D0,
     :     -31046D0,      -1D0,    131D0,    13238D0,   -11D0,    59D0,
     :      28593D0,       0D0,     -1D0,   -12338D0,    10D0,    -3D0,
     :      20441D0,      21D0,     10D0,   -10758D0,     0D0,    -3D0,
     :      29243D0,       0D0,    -74D0,     -609D0,     0D0,    13D0,
     :      25887D0,       0D0,    -66D0,     -550D0,     0D0,    11D0,
     :     -14053D0,     -25D0,     79D0,     8551D0,    -2D0,   -45D0,
     :      15164D0,      10D0,     11D0,    -8001D0,     0D0,    -1D0,
     :     -15794D0,      72D0,    -16D0,     6850D0,   -42D0,    -5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=31,40 ) /
     :      21783D0,       0D0,     13D0,     -167D0,     0D0,    13D0,
     :     -12873D0,     -10D0,    -37D0,     6953D0,     0D0,   -14D0,
     :     -12654D0,      11D0,     63D0,     6415D0,     0D0,    26D0,
     :     -10204D0,       0D0,     25D0,     5222D0,     0D0,    15D0,
     :      16707D0,     -85D0,    -10D0,      168D0,    -1D0,    10D0,
     :      -7691D0,       0D0,     44D0,     3268D0,     0D0,    19D0,
     :     -11024D0,       0D0,    -14D0,      104D0,     0D0,     2D0,
     :       7566D0,     -21D0,    -11D0,    -3250D0,     0D0,    -5D0,
     :      -6637D0,     -11D0,     25D0,     3353D0,     0D0,    14D0,
     :      -7141D0,      21D0,      8D0,     3070D0,     0D0,     4D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=41,50 ) /
     :      -6302D0,     -11D0,      2D0,     3272D0,     0D0,     4D0,
     :       5800D0,      10D0,      2D0,    -3045D0,     0D0,    -1D0,
     :       6443D0,       0D0,     -7D0,    -2768D0,     0D0,    -4D0,
     :      -5774D0,     -11D0,    -15D0,     3041D0,     0D0,    -5D0,
     :      -5350D0,       0D0,     21D0,     2695D0,     0D0,    12D0,
     :      -4752D0,     -11D0,     -3D0,     2719D0,     0D0,    -3D0,
     :      -4940D0,     -11D0,    -21D0,     2720D0,     0D0,    -9D0,
     :       7350D0,       0D0,     -8D0,      -51D0,     0D0,     4D0,
     :       4065D0,       0D0,      6D0,    -2206D0,     0D0,     1D0,
     :       6579D0,       0D0,    -24D0,     -199D0,     0D0,     2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=51,60 ) /
     :       3579D0,       0D0,      5D0,    -1900D0,     0D0,     1D0,
     :       4725D0,       0D0,     -6D0,      -41D0,     0D0,     3D0,
     :      -3075D0,       0D0,     -2D0,     1313D0,     0D0,    -1D0,
     :      -2904D0,       0D0,     15D0,     1233D0,     0D0,     7D0,
     :       4348D0,       0D0,    -10D0,      -81D0,     0D0,     2D0,
     :      -2878D0,       0D0,      8D0,     1232D0,     0D0,     4D0,
     :      -4230D0,       0D0,      5D0,      -20D0,     0D0,    -2D0,
     :      -2819D0,       0D0,      7D0,     1207D0,     0D0,     3D0,
     :      -4056D0,       0D0,      5D0,       40D0,     0D0,    -2D0,
     :      -2647D0,       0D0,     11D0,     1129D0,     0D0,     5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=61,70 ) /
     :      -2294D0,       0D0,    -10D0,     1266D0,     0D0,    -4D0,
     :       2481D0,       0D0,     -7D0,    -1062D0,     0D0,    -3D0,
     :       2179D0,       0D0,     -2D0,    -1129D0,     0D0,    -2D0,
     :       3276D0,       0D0,      1D0,       -9D0,     0D0,     0D0,
     :      -3389D0,       0D0,      5D0,       35D0,     0D0,    -2D0,
     :       3339D0,       0D0,    -13D0,     -107D0,     0D0,     1D0,
     :      -1987D0,       0D0,     -6D0,     1073D0,     0D0,    -2D0,
     :      -1981D0,       0D0,      0D0,      854D0,     0D0,     0D0,
     :       4026D0,       0D0,   -353D0,     -553D0,     0D0,  -139D0,
     :       1660D0,       0D0,     -5D0,     -710D0,     0D0,    -2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=71,77 ) /
     :      -1521D0,       0D0,      9D0,      647D0,     0D0,     4D0,
     :       1314D0,       0D0,      0D0,     -700D0,     0D0,     0D0,
     :      -1283D0,       0D0,      0D0,      672D0,     0D0,     0D0,
     :      -1331D0,       0D0,      8D0,      663D0,     0D0,     4D0,
     :       1383D0,       0D0,     -2D0,     -594D0,     0D0,    -2D0,
     :       1405D0,       0D0,      4D0,     -610D0,     0D0,     2D0,
     :       1290D0,       0D0,      0D0,     -556D0,     0D0,     0D0 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------

*
*  Fundamental (Delaunay) arguments from Simon et al. (1994)
*

*  Mean anomaly of the Moon.
      EL  = MOD (         485868.249036D0 +
     :            ( + 1717915923.2178D0 ) * T, TURNAS ) * DAS2R

*  Mean anomaly of the Sun.
      ELP = MOD (        1287104.79305D0 +
     :            (  + 129596581.0481D0 ) * T, TURNAS ) * DAS2R

*  Mean argument of the latitude of the Moon.
      F   = MOD (         335779.526232D0 +
     :            ( + 1739527262.8478D0 ) * T, TURNAS ) * DAS2R

*  Mean elongation of the Moon from the Sun.
      D   = MOD (        1072260.70369D0 +
     :            ( + 1602961601.2090D0 ) * T, TURNAS ) * DAS2R

*  Mean longitude of the ascending node of the Moon.
      OM  = MOD (         450160.398036D0 +
     :            (    - 6962890.5431D0 ) * T, TURNAS ) * DAS2R

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG

 100  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R

*  -----------------------------
*  IN LIEU OF PLANETARY NUTATION
*  -----------------------------

*  Fixed offset to correct for missing terms in truncated series.
      DPSIPL = DPPLAN
      DEPSPL = DEPLAN

*  -------
*  RESULTS
*  -------

*  Add luni-solar and planetary components.
      DPSI = DPSILS + DPSIPL
      DEPS = DEPSLS + DEPSPL

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
