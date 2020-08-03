C$Procedure ZZSGP4 ( SGP4 wrapper )

      SUBROUTINE ZZSGP4 ( GEOPHS, ELEMS, OPMODE, T, STATE )

C$ Abstract
C
C     Umbrella for the SGP4 initializer and evaluator routines.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     None.
C
C$ Declarations

      IMPLICIT NONE
      INCLUDE 'zzsgp4.inc'

      DOUBLE PRECISION         GEOPHS    ( 8 )
      DOUBLE PRECISION         ELEMS     ( 10 )
      INTEGER                  OPMODE
      DOUBLE PRECISION         T
      DOUBLE PRECISION         STATE     ( 6 )

C$ Brief_I/O
C
C     Variable  I/O  Entry Point
C     --------  ---  --------------------------------------------------
C     GEOPHS     I   XXSGP4I
C     ELEMS      I   XXSGP4I
C     OPMODE     I   XXSGP4I
C     T          I   XXSGP4E
C     STATE      O   XXSGP4E
C
C$ Detailed_Input
C
C     See Individual Entry points.
C
C$ Detailed_Output
C
C     See Individual Entry points.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) 'SPICE(BOGUSENTRY)' will signal if ZZSGP4 is called.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine wraps XXSGP4E and XXSGP4I. As entry points to
C     this routine, they share the local memory space.
C
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine
C     specific arithmetic implementation.
C
C     Use a set of TLEs to calculate a collection of states for a
C     time interval centered at the TLE set epoch.
C
C           PROGRAM ZZSGP4_T
C
C     C
C     C     Read a data file containing sets of TLEs, then calculate
C     C     states at -1440 to 1440 minutes from the epoch of each
C     C     TLE set in steps of 10 minutes.
C     C
C     C     Example cases listed in sgp4-ver.tle.
C     C
C     C     1 00005U 58002B   00179.78495062  .00000023
C     C       00000-0  28098-4 0  4753
C     C     2 00005  34.2682 348.7242 1859667 331.7664
C     C       19.3264 10.82419157413667
C     C
C     C     1 04632U 70093B   04031.91070959 -.00000084
C     C       00000-0  10000-3 0  9955
C     C     2 04632  11.4628 273.1101 1450506 207.6000
C     C       143.9350  1.20231981 44145
C     C
C     C     1 06251U 62025E   06176.82412014  .00008885
C     C       00000-0  12808-3 0  3985
C     C     2 06251  58.0579  54.0425 0030035 139.1568
C     C       221.1854 15.56387291  6774
C     C
C     C     1 08195U 75081A   06176.33215444  .00000099
C     C       00000-0  11873-3 0   813
C     C     2 08195  64.1586 279.0717 6877146 264.7651
C     C       20.2257  2.00491383225656
C     C
C     C     1 09880U 77021A   06176.56157475  .00000421
C     C       00000-0  10000-3 0  9814
C     C     2 09880  64.5968 349.3786 7069051 270.0229
C     C       16.3320  2.00813614112380
C     C
C     C     1 09998U 74033F   05148.79417928 -.00000112
C     C       00000-0  00000+0 0  4480
C     C     2 09998   9.4958 313.1750 0270971 327.5225
C     C       30.8097  1.16186785 45878
C     C
C     C     1 11801U          80230.29629788  .01431103
C     C       00000-0  14311-1      13
C     C     2 11801  46.7916 230.4354 7318036  47.4722
C     C       10.4117  2.28537848    13
C     C
C     C     1 14128U 83058A   06176.02844893 -.00000158
C     C       00000-0  10000-3 0  9627
C     C     2 14128  11.4384  35.2134 0011562  26.4582
C     C       333.5652  0.98870114 46093
C     C
C     C     1 16925U 86065D   06151.67415771  .02550794
C     C       -30915-6  18784-3 0  4486
C     C     2 16925  62.0906 295.0239 5596327 245.1593
C     C       47.9690  4.88511875148616
C     C
C     C     1 20413U 83020D   05363.79166667  .00000000
C     C       00000-0  00000+0 0  7041
C     C     2 20413  12.3514 187.4253 7864447 196.3027
C     C       356.5478  0.24690082  7978
C     C
C     C     1 21897U 92011A   06176.02341244 -.00001273
C     C       00000-0 -13525-3 0  3044
C     C     2 21897  62.1749 198.0096 7421690 253.0462
C     C       20.1561  2.01269994104880
C     C
C     C     1 22312U 93002D   06094.46235912  .99999999
C     C       81888-5  49949-3 0  3953
C     C     2 22312  62.1486  77.4698 0308723 267.9229
C     C       88.7392 15.95744531 98783
C     C
C     C     1 22674U 93035D   06176.55909107  .00002121
C     C       00000-0  29868-3 0  6569
C     C     2 22674  63.5035 354.4452 7541712 253.3264
C     C       18.7754  1.96679808 93877
C     C
C     C     1 23177U 94040C   06175.45752052  .00000386
C     C       00000-0  76590-3 0    95
C     C     2 23177   7.0496 179.8238 7258491 296.0482
C     C       8.3061  2.25906668 97438
C     C
C     C     1 23333U 94071A   94305.49999999 -.00172956
C     C       26967-3  10000-3 0    15
C     C     2 23333  28.7490   2.3720 9728298  30.4360
C     C       1.3500  0.07309491    70
C     C
C     C     1 23599U 95029B   06171.76535463  .00085586
C     C       12891-6  12956-2 0  2905
C     C     2 23599   6.9327   0.2849 5782022 274.4436
C     C       25.2425  4.47796565123555
C     C
C     C     1 24208U 96044A   06177.04061740 -.00000094
C     C       00000-0  10000-3 0  1600
C     C     2 24208   3.8536  80.0121 0026640 311.0977
C     C       48.3000  1.00778054 36119
C     C
C     C     1 25954U 99060A   04039.68057285 -.00000108
C     C       00000-0  00000-0 0  6847
C     C     2 25954   0.0004 243.8136 0001765  15.5294
C     C       22.7134  1.00271289 15615
C     C
C     C     1 26900U 01039A   06106.74503247  .00000045
C     C       00000-0  10000-3 0  8290
C     C     2 26900   0.0164 266.5378 0003319  86.1794
C     C       182.2590  1.00273847 16981
C     C
C     C     1 26975U 78066F   06174.85818871  .00000620
C     C       00000-0  10000-3 0  6809
C     C     2 26975  68.4714 236.1303 5602877 123.7484
C     C       302.5767  2.05657553 67521
C     C
C     C     1 28057U 03049A   06177.78615833  .00000060
C     C       00000-0  35940-4 0  1836
C     C     2 28057  98.4283 247.6961 0000884  88.1964
C     C       271.9322 14.35478080140550
C     C
C     C     1 28129U 03058A   06175.57071136 -.00000104
C     C       00000-0  10000-3 0   459
C     C     2 28129  54.7298 324.8098 0048506 266.2640
C     C       93.1663  2.00562768 18443
C     C
C     C     1 28350U 04020A   06167.21788666  .16154492
C     C       76267-5  18678-3 0  8894
C     C     2 28350  64.9977 345.6130 0024870 260.7578
C     C       99.9590 16.47856722116490
C     C
C     C     1 28623U 05006B   06177.81079184  .00637644
C     C       69054-6  96390-3 0  6000
C     C     2 28623  28.5200 114.9834 6249053 170.2550
C     C       212.8965  3.79477162 12753
C     C
C     C     1 28626U 05008A   06176.46683397 -.00000205
C     C       00000-0  10000-3 0  2190
C     C     2 28626   0.0019 286.9433 0000335  13.7918
C     C       55.6504  1.00270176  4891
C     C
C     C     1 28872U 05037B   05333.02012661  .25992681
C     C       00000-0  24476-3 0  1534
C     C     2 28872  96.4736 157.9986 0303955 244.0492
C     C       110.6523 16.46015938 10708
C     C
C     C     1 29141U 85108AA  06170.26783845  .99999999
C     C       00000-0  13519-0 0   718
C     C     2 29141  82.4288 273.4882 0015848 277.2124
C     C       83.9133 15.93343074  6828
C     C
C     C     1 29238U 06022G   06177.28732010  .00766286
C     C       10823-4  13334-2 0   101
C     C     2 29238  51.5595 213.7903 0202579  95.2503
C     C       267.9010 15.73823839  1061
C     C
C     C     1 88888U          80275.98708465  .00073094
C     C       13844-3  66816-4 0    87
C     C     2 88888  72.8435 115.9689 0086731  52.6988
C     C       110.5714 16.05824518  1058
C     C
C
C           IMPLICIT NONE
C
C           INCLUDE 'zzsgp4.inc'
C
C           CHARACTER*(72)           LINES  ( 2 )
C           CHARACTER*(72)           TLEDAT
C
C           INTEGER                  FRSTYR
C           INTEGER                  I
C           INTEGER                  OPMODE
C
C           DOUBLE PRECISION         DELT
C           DOUBLE PRECISION         ELEMS  ( 10 )
C           DOUBLE PRECISION         EPOCH
C           DOUBLE PRECISION         GEOPHS ( 8 )
C           DOUBLE PRECISION         STATE  ( 6 )
C           DOUBLE PRECISION         TF
C           DOUBLE PRECISION         TIME
C           DOUBLE PRECISION         TS
C
C           LOGICAL                  EOF
C
C     C
C     C     SPICELIB routines.
C     C
C           LOGICAL                  FAILED
C
C     C
C     C     Load a leapseconds kernel for time conversion. Required
C     C     by the SPK 10 evaluator.
C     C
C           CALL FURNSH ( '/kernels/gen/lsk/naif0011.tls' )
C
C     C
C     C     Define the geophysical quantities using the values
C     C     from geophysical.ker.
C     C
C           GEOPHS( K_J2 ) =    1.082616D-3
C           GEOPHS( K_J3 ) =   -2.53881D-6
C           GEOPHS( K_J4 ) =   -1.65597D-6
C           GEOPHS( K_KE ) =    7.43669161D-2
C           GEOPHS( K_QO ) =  120.0D0
C           GEOPHS( K_SO ) =   78.0D0
C           GEOPHS( K_ER ) = 6378.135D0
C           GEOPHS( K_AE ) =    1.0D0
C
C           TLEDAT = 'sgp4-ver1.tle'
C
C     C
C     C     Error subsystem to report to ensure execution continues
C     C     if an error signals.
C     C
C           CALL ERRACT( 'SET', 'REPORT')
C
C     C
C     C     Use Spacetrack #3 algorithm to calculate sidereal time.
C     C
C           OPMODE = AFSPC
C
C     C
C     C     Identify the earliest year for the elements.
C     C
C           FRSTYR = 1958
C
C     C
C     C     Start and final offsets from TLE epochs. [-1440, 1400]
C     C     minutes.
C     C
C           TS     =  -1440.0D0
C           TF     =   1440.0D0
C
C     C
C     C     Step size for elements output 10 minutes.
C     C
C           DELT   = 10.D0
C
C     C
C     C     Read the TLE data file.
C     C
C           CALL RDTEXT ( TLEDAT, LINES(1), EOF )
C           CALL RDTEXT ( TLEDAT, LINES(2), EOF )
C
C     C
C     C     Loop over data file until end-of-file.
C     C
C           DO WHILE ( .NOT. EOF )
C
C     C
C     C        Parse the elements to something SPICE can use.
C     C
C              CALL GETELM ( FRSTYR, LINES, EPOCH, ELEMS )
C
C              WRITE(*, FMT='(A72)') LINES(1)
C              WRITE(*, FMT='(A72)') LINES(2)
C              WRITE(*,*) ' '
C
C     C
C     C        Initialize SGP4 calculations based on values in
C     C        GEOPHS, ELEMS, and AFSPC.
C     C
C              CALL XXSGP4I ( GEOPHS, ELEMS, OPMODE )
C
C     C
C     C        Start time keyed in minutes from TLE epoch.
C     C
C              TIME   = TS
C
C              DO WHILE ( TIME .LE. DABS(TF) .AND. (.NOT. FAILED()) )
C
C     C
C     C           Calculate the STATE at TIME.
C     C
C                 CALL XXSGP4E ( TIME, STATE )
C
C     C
C     C           If the propagation succeeded, output the STATE.
C     C
C                 IF ( .NOT. FAILED() ) THEN
C
C                    WRITE(*, FMT ='(7F17.8)' ) TIME,
C          .                                    (STATE(I),I=1,6)
C
C                 END IF
C
C     C
C     C           Increment the evaluation time by one step.
C     C
C                 TIME = TIME + DELT
C
C              END DO
C
C              WRITE(*,*) ' '
C
C     C
C     C        reset the error subsystem for the next loop.
C     C
C              CALL RESET()
C
C     C
C     C        Read the next two lines (if any) from the TLE
C     C        data file.
C     C
C              CALL RDTEXT ( TLEDAT, LINES(1), EOF )
C              CALL RDTEXT ( TLEDAT, LINES(2), EOF )
C
C           END DO
C
C           END
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C   [1] Hoots, F. R., and Roehrich, R. L. 1980. "Models for
C       Propagation of the NORAD Element Sets." Spacetrack Report #3.
C       U.S. Air Force: Aerospace Defense Command.
C
C   [2] Hoots, Felix R. "Spacetrack Report #6: Models for Propagation
C       of Space Command Element Sets." Space Command,
C       U. S. Air Force, CO.
C
C   [3] Hoots, Felix R., P. W. Schumacher, and R. A. Glover. 2004.
C       History of Analytical Orbit Modeling in the U. S. Space
C       Surveillance System. Journal of Guidance, Control, and
C       Dynamics. 27(2):174-185.
C
C   [4] Vallado, David, Crawford, Paul, Hujsak, Richard,
C       and Kelso, T.S. 2006. Revisiting Spacetrack Report #3. Paper
C       AIAA 2006-6753 presented at the AIAA/AAS Astrodynamics
C       Specialist Conference, August 21-24, 2006. Keystone, CO.
C
C$ Author_and_Institution
C
C     David Vallado   (AGI)
C     E. D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, SEP-15-2014 (EDW)
C
C        Based on routine SGP4, 28-JUN-2005, Vallado 2006 [4].
C
C-&

C$ Index_Entries
C
C  SGP4
C
C-&

C
C     Local Variables
C

      DOUBLE PRECISION         A
      DOUBLE PRECISION         ALTA
      DOUBLE PRECISION         ALTP
      DOUBLE PRECISION         ARGPDOT
      DOUBLE PRECISION         AYCOF
      DOUBLE PRECISION         CC1
      DOUBLE PRECISION         CC4
      DOUBLE PRECISION         CC5
      DOUBLE PRECISION         CON41
      DOUBLE PRECISION         ETA
      DOUBLE PRECISION         INCLO
      DOUBLE PRECISION         MDOT
      DOUBLE PRECISION         NODEDOT
      DOUBLE PRECISION         OMGCOF
      DOUBLE PRECISION         SINMAO
      DOUBLE PRECISION         T2COF
      DOUBLE PRECISION         T3COF
      DOUBLE PRECISION         T4COF
      DOUBLE PRECISION         T5COF
      DOUBLE PRECISION         X1MTH2
      DOUBLE PRECISION         X7THM1
      DOUBLE PRECISION         XLCOF
      DOUBLE PRECISION         XMCOF
      DOUBLE PRECISION         XNODCF

C
C     DS values
C

      DOUBLE PRECISION         E3
      DOUBLE PRECISION         EE2
      DOUBLE PRECISION         PEO
      DOUBLE PRECISION         PGHO
      DOUBLE PRECISION         PHO
      DOUBLE PRECISION         PINCO
      DOUBLE PRECISION         PLO
      DOUBLE PRECISION         GSTO
      DOUBLE PRECISION         XFACT
      DOUBLE PRECISION         XGH2
      DOUBLE PRECISION         XGH3
      DOUBLE PRECISION         XGH4
      DOUBLE PRECISION         XH2
      DOUBLE PRECISION         XH3
      DOUBLE PRECISION         XI2
      DOUBLE PRECISION         XI3
      DOUBLE PRECISION         XL2
      DOUBLE PRECISION         XL3
      DOUBLE PRECISION         XL4
      DOUBLE PRECISION         XLAMO
      DOUBLE PRECISION         ZMOL
      DOUBLE PRECISION         ZMOS
      DOUBLE PRECISION         ATIME
      DOUBLE PRECISION         XLI
      DOUBLE PRECISION         XNI

      DOUBLE PRECISION         AINV
      DOUBLE PRECISION         AM
      DOUBLE PRECISION         AO
      DOUBLE PRECISION         ARGPM
      DOUBLE PRECISION         ARGPO
      DOUBLE PRECISION         ARGPP
      DOUBLE PRECISION         AXNL
      DOUBLE PRECISION         AYNL
      DOUBLE PRECISION         BETAL
      DOUBLE PRECISION         BSTAR
      DOUBLE PRECISION         CC1SQ
      DOUBLE PRECISION         CC2
      DOUBLE PRECISION         CC3
      DOUBLE PRECISION         CNOD
      DOUBLE PRECISION         CNODM
      DOUBLE PRECISION         COEF
      DOUBLE PRECISION         COEF1
      DOUBLE PRECISION         CON42
      DOUBLE PRECISION         COS2U
      DOUBLE PRECISION         COSEO1
      DOUBLE PRECISION         COSI
      DOUBLE PRECISION         COSIM
      DOUBLE PRECISION         COSIO
      DOUBLE PRECISION         COSIO2
      DOUBLE PRECISION         COSIO4
      DOUBLE PRECISION         COSIP
      DOUBLE PRECISION         COSISQ
      DOUBLE PRECISION         COSOMM
      DOUBLE PRECISION         COSSU
      DOUBLE PRECISION         COSU
      DOUBLE PRECISION         D2
      DOUBLE PRECISION         D2201
      DOUBLE PRECISION         D2211
      DOUBLE PRECISION         D3
      DOUBLE PRECISION         D3210
      DOUBLE PRECISION         D3222
      DOUBLE PRECISION         D4
      DOUBLE PRECISION         D4410
      DOUBLE PRECISION         D4422
      DOUBLE PRECISION         D5220
      DOUBLE PRECISION         D5232
      DOUBLE PRECISION         D5421
      DOUBLE PRECISION         D5433
      DOUBLE PRECISION         DAY
      DOUBLE PRECISION         DEDT
      DOUBLE PRECISION         DEL1
      DOUBLE PRECISION         DEL2
      DOUBLE PRECISION         DEL3
      DOUBLE PRECISION         DELM
      DOUBLE PRECISION         DELMO
      DOUBLE PRECISION         DELOMG
      DOUBLE PRECISION         DIDT
      DOUBLE PRECISION         DMDT
      DOUBLE PRECISION         DNDT
      DOUBLE PRECISION         DNODT
      DOUBLE PRECISION         DOMDT
      DOUBLE PRECISION         ECCM
      DOUBLE PRECISION         ECCO
      DOUBLE PRECISION         ECCP
      DOUBLE PRECISION         ECCSQ
      DOUBLE PRECISION         ECOSE
      DOUBLE PRECISION         EETA
      DOUBLE PRECISION         EL2
      DOUBLE PRECISION         EMSQ
      DOUBLE PRECISION         EO1
      DOUBLE PRECISION         EPOCH
      DOUBLE PRECISION         ER
      DOUBLE PRECISION         ESINE
      DOUBLE PRECISION         ETASQ
      DOUBLE PRECISION         GAM
      DOUBLE PRECISION         INCLM
      DOUBLE PRECISION         J2
      DOUBLE PRECISION         J3
      DOUBLE PRECISION         J3OJ2
      DOUBLE PRECISION         J4
      DOUBLE PRECISION         KPS
      DOUBLE PRECISION         MM
      DOUBLE PRECISION         MO
      DOUBLE PRECISION         MP
      DOUBLE PRECISION         MR
      DOUBLE PRECISION         MV
      DOUBLE PRECISION         NO
      DOUBLE PRECISION         NODEM
      DOUBLE PRECISION         NODEO
      DOUBLE PRECISION         NODEP
      DOUBLE PRECISION         OMEOSQ
      DOUBLE PRECISION         OMGADF
      DOUBLE PRECISION         PERIGE
      DOUBLE PRECISION         PINVSQ
      DOUBLE PRECISION         PL
      DOUBLE PRECISION         POSQ
      DOUBLE PRECISION         PSISQ
      DOUBLE PRECISION         QZMS24
      DOUBLE PRECISION         QZMS2T
      DOUBLE PRECISION         RDOTL
      DOUBLE PRECISION         RL
      DOUBLE PRECISION         RP
      DOUBLE PRECISION         RTEMSQ
      DOUBLE PRECISION         RTEOSQ
      DOUBLE PRECISION         RVDOT
      DOUBLE PRECISION         RVDOTL
      DOUBLE PRECISION         S1
      DOUBLE PRECISION         S2
      DOUBLE PRECISION         S3
      DOUBLE PRECISION         S4
      DOUBLE PRECISION         S5
      DOUBLE PRECISION         S6
      DOUBLE PRECISION         S7
      DOUBLE PRECISION         SE2
      DOUBLE PRECISION         SE3
      DOUBLE PRECISION         SFOUR
      DOUBLE PRECISION         SGH2
      DOUBLE PRECISION         SGH3
      DOUBLE PRECISION         SGH4
      DOUBLE PRECISION         SH2
      DOUBLE PRECISION         SH3
      DOUBLE PRECISION         SI2
      DOUBLE PRECISION         SI3
      DOUBLE PRECISION         SIN2U
      DOUBLE PRECISION         SINEO1
      DOUBLE PRECISION         SINI
      DOUBLE PRECISION         SINIM
      DOUBLE PRECISION         SINIO
      DOUBLE PRECISION         SINIP
      DOUBLE PRECISION         SINOMM
      DOUBLE PRECISION         SINSU
      DOUBLE PRECISION         SINU
      DOUBLE PRECISION         SL2
      DOUBLE PRECISION         SL3
      DOUBLE PRECISION         SL4
      DOUBLE PRECISION         SNOD
      DOUBLE PRECISION         SNODM
      DOUBLE PRECISION         SS
      DOUBLE PRECISION         SS1
      DOUBLE PRECISION         SS2
      DOUBLE PRECISION         SS3
      DOUBLE PRECISION         SS4
      DOUBLE PRECISION         SS5
      DOUBLE PRECISION         SS6
      DOUBLE PRECISION         SS7
      DOUBLE PRECISION         SU
      DOUBLE PRECISION         SZ1
      DOUBLE PRECISION         SZ11
      DOUBLE PRECISION         SZ12
      DOUBLE PRECISION         SZ13
      DOUBLE PRECISION         SZ2
      DOUBLE PRECISION         SZ21
      DOUBLE PRECISION         SZ22
      DOUBLE PRECISION         SZ23
      DOUBLE PRECISION         SZ3
      DOUBLE PRECISION         SZ31
      DOUBLE PRECISION         SZ32
      DOUBLE PRECISION         SZ33
      DOUBLE PRECISION         T2
      DOUBLE PRECISION         T3
      DOUBLE PRECISION         T4
      DOUBLE PRECISION         TC
      DOUBLE PRECISION         TEM5
      DOUBLE PRECISION         TEMP
      DOUBLE PRECISION         TEMP1
      DOUBLE PRECISION         TEMP2
      DOUBLE PRECISION         TEMP3
      DOUBLE PRECISION         TEMP4
      DOUBLE PRECISION         TEMPA
      DOUBLE PRECISION         TEMPE
      DOUBLE PRECISION         TEMPL
      DOUBLE PRECISION         TSI
      DOUBLE PRECISION         TUMIN
      DOUBLE PRECISION         TVEC   ( 8 )
      DOUBLE PRECISION         TZERO
      DOUBLE PRECISION         U
      DOUBLE PRECISION         UX
      DOUBLE PRECISION         UY
      DOUBLE PRECISION         UZ
      DOUBLE PRECISION         VX
      DOUBLE PRECISION         VY
      DOUBLE PRECISION         VZ
      DOUBLE PRECISION         X2O3
      DOUBLE PRECISION         XHDOT1
      DOUBLE PRECISION         XINC
      DOUBLE PRECISION         XINCP
      DOUBLE PRECISION         XKE
      DOUBLE PRECISION         XL
      DOUBLE PRECISION         XLM
      DOUBLE PRECISION         XMDF
      DOUBLE PRECISION         XMX
      DOUBLE PRECISION         XMY
      DOUBLE PRECISION         XN
      DOUBLE PRECISION         XNODDF
      DOUBLE PRECISION         XNODE
      DOUBLE PRECISION         XPIDOT
      DOUBLE PRECISION         Z1
      DOUBLE PRECISION         Z11
      DOUBLE PRECISION         Z12
      DOUBLE PRECISION         Z13
      DOUBLE PRECISION         Z2
      DOUBLE PRECISION         Z21
      DOUBLE PRECISION         Z22
      DOUBLE PRECISION         Z23
      DOUBLE PRECISION         Z3
      DOUBLE PRECISION         Z31
      DOUBLE PRECISION         Z32
      DOUBLE PRECISION         Z33

      INTEGER                  ITER
      INTEGER                  SVMODE
      INTEGER                  IREZ

      LOGICAL                  DOINIT
      LOGICAL                  DOSIMP
      LOGICAL                  DODEEP

C
C     SPICELIB routines.
C

      DOUBLE PRECISION         PI
      DOUBLE PRECISION         TWOPI

      LOGICAL                  FAILED
      LOGICAL                  RETURN

      SAVE                     A
      SAVE                     ALTA
      SAVE                     ALTP
      SAVE                     ARGPDOT
      SAVE                     ARGPO
      SAVE                     ATIME
      SAVE                     AYCOF
      SAVE                     BSTAR
      SAVE                     CC1
      SAVE                     CC4
      SAVE                     CC5
      SAVE                     CON41
      SAVE                     D2
      SAVE                     D2201
      SAVE                     D2211
      SAVE                     D3
      SAVE                     D3210
      SAVE                     D3222
      SAVE                     D4
      SAVE                     D4410
      SAVE                     D4422
      SAVE                     D5220
      SAVE                     D5232
      SAVE                     D5421
      SAVE                     D5433
      SAVE                     DEDT
      SAVE                     DEL1
      SAVE                     DEL2
      SAVE                     DEL3
      SAVE                     DELMO
      SAVE                     DIDT
      SAVE                     DMDT
      SAVE                     DNODT
      SAVE                     DODEEP
      SAVE                     DOMDT
      SAVE                     DOSIMP
      SAVE                     E3
      SAVE                     ECCO
      SAVE                     EE2
      SAVE                     ER
      SAVE                     ETA
      SAVE                     GSTO
      SAVE                     INCLO
      SAVE                     IREZ
      SAVE                     J2
      SAVE                     J3OJ2
      SAVE                     MDOT
      SAVE                     MO
      SAVE                     NO
      SAVE                     NODEDOT
      SAVE                     NODEO
      SAVE                     OMGCOF
      SAVE                     PEO
      SAVE                     PGHO
      SAVE                     PHO
      SAVE                     PINCO
      SAVE                     PLO
      SAVE                     SE2
      SAVE                     SE3
      SAVE                     SGH2
      SAVE                     SGH3
      SAVE                     SGH4
      SAVE                     SH2
      SAVE                     SH3
      SAVE                     SI2
      SAVE                     SI3
      SAVE                     SINMAO
      SAVE                     SL2
      SAVE                     SL3
      SAVE                     SL4
      SAVE                     SVMODE
      SAVE                     T2COF
      SAVE                     T3COF
      SAVE                     T4COF
      SAVE                     T5COF
      SAVE                     X1MTH2
      SAVE                     X7THM1
      SAVE                     XFACT
      SAVE                     XGH2
      SAVE                     XGH3
      SAVE                     XGH4
      SAVE                     XH2
      SAVE                     XH3
      SAVE                     XI2
      SAVE                     XI3
      SAVE                     XKE
      SAVE                     XL2
      SAVE                     XL3
      SAVE                     XL4
      SAVE                     XLAMO
      SAVE                     XLCOF
      SAVE                     XLI
      SAVE                     XMCOF
      SAVE                     XNI
      SAVE                     XNODCF
      SAVE                     ZMOL
      SAVE                     ZMOS

      CALL CHKIN  ( 'ZZSGP4' )
      CALL SETMSG ( 'The routine ZZSGP4 is an umbrella for '
     .//            'the SGP4 initializer and propagator entry '
     .//            'points. Do not call ZZSGP4. It is '
     .//            'likely that a programming error has '
     .//            'been made.' )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)'  )
      CALL CHKOUT ( 'ZZSGP4' )

      RETURN



C$Procedure XXSGP4I ( SGP4 initializer )

      ENTRY XXSGP4I ( GEOPHS, ELEMS, OPMODE )

C$ Abstract
C
C     This subroutine initializes variables for SGP4.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     None.
C
C$ Declarations
C
C     Refer to Declarations section in ZZSGP4.
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     GEOPHS     I   Geophysical constants
C     ELEMS      I   Two-line element data
C     OPMODE     I   Flag to indicate operation mode for GMST.
C
C$ Detailed_Input
C
C     GEOPHS      is a collection of 8 geophysical constants needed
C                 for computing a state.  The order of these
C                 constants must be:
C
C                 GEOPHS(1) = J2 gravitational harmonic for earth
C                 GEOPHS(2) = J3 gravitational harmonic for earth
C                 GEOPHS(3) = J4 gravitational harmonic for earth
C
C                 These first three constants are dimensionless.
C
C                 GEOPHS(4) = KE: Square root of the GM for earth where
C                             GM is expressed in earth radii cubed per
C                             minutes squared.
C
C                 GEOPHS(5) = QO: Low altitude bound for atmospheric
C                             model in km.
C
C                 GEOPHS(6) = SO: High altitude bound for atmospheric
C                             model in km.
C
C                 GEOPHS(7) = RE: Equatorial radius of the earth in km.
C
C
C                 GEOPHS(8) = AE: Distance units/earth radius
C                             (normally 1)
C
C                 Below are currently recommended values for these
C                 items:
C
C                   J2 =    1.082616D-3
C                   J3 =   -2.53881D-6
C                   J4 =   -1.65597D-6
C
C                 The next item is the square root of GM for the
C                 earth given in units of earth-radii**1.5/Minute
C
C                   KE =    7.43669161D-2
C
C                 The next two items give the top and
C                 bottom of the atmospheric drag model
C                 used by the type 10 ephemeris type.
C                 Don't adjust these unless you understand
C                 the full implications of such changes.
C
C                   QO =  120.0D0
C                   SO =   78.0D0
C
C                 The following is the equatorial radius
C                 of the earth as used by NORAD in km.
C
C                   ER = 6378.135D0
C
C                 The value of AE is the number of
C                 distance units per earth radii used by
C                 the NORAD state propagation software.
C                 The value should be 1 unless you've got
C                 a very good understanding of the NORAD
C                 routine SGP4 and the affect of changing
C                 this value..
C
C                   AE =    1.0D0
C
C     ELEMS       is an array containing two-line element data
C                 as prescribed below. The elements XNDD6O and BSTAR
C                 must already be scaled by the proper exponent stored
C                 in the two line elements set.  Moreover, the
C                 various items must be converted to the units shown
C                 here.
C
C                    ELEMS (  1 ) = XNDT2O in radians/minute**2
C                    ELEMS (  2 ) = XNDD6O in radians/minute**3
C                    ELEMS (  3 ) = BSTAR
C                    ELEMS (  4 ) = XINCL  in radians
C                    ELEMS (  5 ) = XNODEO in radians
C                    ELEMS (  6 ) = EO
C                    ELEMS (  7 ) = OMEGAO in radians
C                    ELEMS (  8 ) = XMO    in radians
C                    ELEMS (  9 ) = XNO    in radians/minute
C                    ELEMS ( 10 ) = EPOCH of the elements in seconds
C                                   past ephemeris epoch J2000.
C
C     OPMODE         Flag indicating which technique
C                    to use to calculate sidereal time.
C
C$ Detailed_Output
C
C     None.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) SPICE(SUBORBITAL) signals if radius of perigee has
C        value less-than 1.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is based on the SGP4INIT code by David Vallado
C     corresponding to "Revisiting Spacetrack Report #3".
C     The intent is to maintain the original Vallado algorithm,
C     changing code only to meet NAIF format standards and to
C     integrate with SPICELIB.
C
C        1) Implemented error checks using SPICE error subsystem.
C           On detecting an error, control returns to the calling
C           routine. This behavior differs from the original
C           version.
C
C        2) Comments prefixed with "SGP4FIX" indicate a comment
C           from the code by Vallado et. al concerning a correction
C           to the STR#3 code.
C
C        3) Eliminated the use of COMMON blocks.
C
C        Removed getgravconst call, replaced with GEOPHS array.
C
C        xBStar,
C        xEcco,
C        Epoch,
C        xArgpo,
C        xInclo,
C        xMo,
C        xNo,
C        xnodeo replaced with ELEMS array.
C        radiusearthkm replaced with ER
C
C$ Examples
C
C     Refer to Examples section in ZZSGP4.
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C   [1] Hoots, F. R., and Roehrich, R. L. 1980. "Models for
C       Propagation of the NORAD Element Sets." Spacetrack Report #3.
C       U.S. Air Force: Aerospace Defense Command.
C
C   [2] Hoots, Felix R. "Spacetrack Report #6: Models for Propagation
C       of Space Command Element Sets." Space Command,
C       U. S. Air Force, CO.
C
C   [3] Hoots, Felix R., P. W. Schumacher, and R. A. Glover. 2004.
C       History of Analytical Orbit Modeling in the U. S. Space
C       Surveillance System. Journal of Guidance, Control, and
C       Dynamics. 27(2):174-185.
C
C   [4] Vallado, David, Crawford, Paul, Hujsak, Richard,
C       and Kelso, T.S. 2006. Revisiting Spacetrack Report #3. Paper
C       AIAA 2006-6753 presented at the AIAA/AAS Astrodynamics
C       Specialist Conference, August 21-24, 2006. Keystone, CO.
C
C$ Author_and_Institution
C
C     David Vallado   (AGI)
C     E. D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, DEC-11-2014 (EDW)
C
C        Based on routine SGP4INIT, 28-JUN-2005, Vallado 2006 [4].
C
C-&

C$ Index_Entries
C
C  SGP4
C
C-&

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'XXSGP4I' )

C
C     Initialize.
C
      DODEEP = .FALSE.
      DOSIMP = .FALSE.
      SVMODE = OPMODE

C
C     This code block replaces the call:
C
C     sgp4fix - note the following variables are also passed directly
C     via sgp4 common. It is possible to streamline the XXSGP4I call
C     by deleting the "x" variables, but the user would need to set
C     the common values first. we include the additional assignment
C     in case twoline2rv is not used.
C
C        bstar  = xbstar
C        ecco   = xecco
C        argpo  = xargpo
C        inclo  = xinclo
C        mo     = xmo
C        no     = xno
C        nodeo  = xnodeo

      BSTAR = ELEMS ( KBSTAR )
      INCLO = ELEMS ( KINCL  )
      NODEO = ELEMS ( KNODE0 )
      ECCO  = ELEMS ( KECC   )
      ARGPO = ELEMS ( KOMEGA )
      MO    = ELEMS ( KMO    )
      NO    = ELEMS ( KNO    )

C
C       Remember that sgp4 uses units of days from 0 jan 1950
C       (sgp4epoch) and minutes from the epoch (time)
C
C       2433281.5 JD TDB = 1949-12-31 00:00:00.000000 TDB
C       2400000.5 JD TDB = 1858-11-17 00:00:00.000000 TDB
C
C       2433281.5 - 2400000.5 = 33281.0
C

C
C     Convert the J2000 TDB representation of the epoch to
C     JD UTC then calculate the offset from the JD 2433281.5 UTC
C     reference.
C

      TVEC(1) = ELEMS ( KEPOCH )
      CALL TTRANS ( 'TDB', 'JDUTC', TVEC )
      EPOCH = TVEC(1) - 2433281.5D0

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'XXSGP4I' )
         RETURN
      END IF

C
C     This code block replaces the call:
C
C     CALL getgravconst( whichconst, tumin,
C     .                  mu, radiusearthkm, xke,
C     .                  j2, j3, j4, j3oj2 )
C
      J2    = GEOPHS(K_J2)
      J3    = GEOPHS(K_J3)
      J4    = GEOPHS(K_J4)
      ER    = GEOPHS(K_ER)
      XKE   = GEOPHS(K_KE)
      TUMIN = 1.D0/GEOPHS(K_KE)
      J3OJ2 = J3/J2

C
C     The following assignment and IF block is taken
C     from TWOLINE2RVSGP4.
C
      A = (NO*TUMIN)**(-2.0D0/3.0D0)

      IF (DABS(ECCO-1.0D0) .GT. 0.000001D0) THEN
         ALTP= (A*(1.0D0-ECCO))-1.0D0
         ALTA= (A*(1.0D0+ECCO))-1.0D0
      ELSE
         ALTA= 999999.9D0
         ALTP= 2.0D0* (4.0D0/(NO*NO)**(1.0D0/3.0D0))
      END IF

      SS     = 78.0D0/ER + 1.0D0
      QZMS2T = ((120.0D0-78.0D0)/ER) ** 4
      X2O3   =  2.0D0 / 3.0D0

C
C     sgp4fix divisor for divide by zero check on inclination
C     the old check used 1.0D0 + cos(pi-1.0D-9), but then compared
C     it to 1.5D-12, so the threshold was changed to 1.5D-12 for
C     consistency.
C
        TEMP4 =   1.5D-12
        TZERO = 0.0D0
        DOINIT= .TRUE.

        CALL ZZINIL( GEOPHS, OPMODE, ECCO,  EPOCH, INCLO,  NO,
     .               AINV,   AO,     CON41, CON42, COSIO,  COSIO2,
     .               ECCSQ,  OMEOSQ, POSQ,  RP,    RTEOSQ, SINIO,
     .               GSTO )

        IF ( FAILED() ) THEN
           CALL CHKOUT ( 'XXSGP4I' )
           RETURN
        END IF


C
C       Check RP for a reasonable value. The propagator may not
C       calculate correct state values for RP < 1.
C
        IF (RP .LT. 1.0D0) THEN

           CALL SETMSG ( 'TLE elements suborbital.'     )
           CALL SIGERR ( 'SPICE(SUBORBITAL)'            )
           CALL CHKOUT ( 'XXSGP4I'                      )
           RETURN

        END IF

C
C       If nodeo and No are gtr 0
C
        IF( (OMEOSQ .GE. 0.0D0) .OR. (NO .GE. 0.0D0) ) THEN

            DOSIMP = .FALSE.

            IF (RP .LT. (220.0D0/ER+1.0D0)) THEN
                DOSIMP = .TRUE.
            END IF

            SFOUR  = SS
            QZMS24 = QZMS2T
            PERIGE = (RP-1.0D0)*ER

C
C           For perigees below 156 km, S and Qoms2t are altered.
C
            IF(PERIGE .LT. 156.0D0) THEN
                SFOUR = PERIGE-78.0D0

                IF(PERIGE .LE. 98.0D0) THEN
                    SFOUR = 20.0D0
                END IF

                QZMS24 = ( (120.0D0-SFOUR)/ER )**4
                SFOUR  = SFOUR/ER + 1.0D0

            END IF

            PINVSQ = 1.0D0/POSQ

            TSI    = 1.0D0/(AO-SFOUR)
            ETA    = AO*ECCO*TSI
            ETASQ  = ETA*ETA
            EETA   = ECCO*ETA
            PSISQ  = DABS(1.0D0-ETASQ)
            COEF   = QZMS24*TSI**4
            COEF1  = COEF/PSISQ**3.5D0
            CC2    = COEF1*NO* (AO* (1.0D0+1.5D0*ETASQ+EETA*
     .               (4.0D0+ETASQ) )+0.375D0*
     .         J2*TSI/PSISQ*CON41*(8.0D0+3.0D0*ETASQ*(8.0D0+ETASQ)))
            CC1    = BSTAR*CC2
            CC3    = 0.0D0

            IF(ECCO .GT. 1.0D-4) THEN
                CC3 = -2.0D0*COEF*TSI*J3OJ2*NO*SINIO/ECCO
            END IF

            X1MTH2 = 1.0D0-COSIO2
            CC4    = 2.0D0*NO*COEF1*AO*OMEOSQ*(ETA*(2.0D0+0.5D0*ETASQ)
     .              +ECCO*(0.5D0 + 2.0D0*ETASQ) - J2*TSI / (AO*PSISQ)*
     .              (-3.0D0*CON41*(1.0D0-2.0D0*
     .       EETA+ETASQ*(1.5D0-0.5D0*EETA))+0.75D0*X1MTH2*(2.0D0*ETASQ
     .       -EETA*(1.0D0+ETASQ))*DCOS(2.0D0*ARGPO)))
            CC5    = 2.0D0*COEF1*AO*OMEOSQ* (1.0D0 + 2.75D0*
     .               (ETASQ + EETA) + EETA*ETASQ )
            COSIO4 = COSIO2*COSIO2
            TEMP1  = 1.5D0*J2*PINVSQ*NO
            TEMP2  = 0.5D0*TEMP1*J2*PINVSQ
            TEMP3  = -0.46875D0*J4*PINVSQ*PINVSQ*NO
            MDOT   = NO + 0.5D0*TEMP1*RTEOSQ*CON41 + 0.0625D0*TEMP2*
     .               RTEOSQ*(13.0D0 - 78.0D0*COSIO2 + 137.0D0*COSIO4)
            ARGPDOT= -0.5D0*TEMP1*CON42 + 0.0625D0*TEMP2*
     .               (7.0D0 - 114.0D0*COSIO2 +
     .        395.0D0*COSIO4)+TEMP3*(3.0D0-36.0D0*COSIO2+49.0D0*COSIO4)
            XHDOT1 = -TEMP1*COSIO
            NODEDOT = XHDOT1+(0.5D0*TEMP2*(4.0D0-19.0D0*COSIO2)+
     .                 2.0D0*TEMP3*(3.0D0 - 7.0D0*COSIO2))*COSIO
            XPIDOT = ARGPDOT+NODEDOT
            OMGCOF = BSTAR*CC3*DCOS(ARGPO)
            XMCOF  = 0.0D0

            IF(ECCO .GT. 1.0D-4) THEN
                XMCOF = -X2O3*COEF*BSTAR/EETA
            END IF

            XNODCF = 3.5D0*OMEOSQ*XHDOT1*CC1
            T2COF  = 1.5D0*CC1

C
C           sgp4fix for divide by zero with xinco = 180 deg.
C
            IF (DABS(COSIO+1.D0).GT. 1.5D-12) THEN
                XLCOF  = -0.25D0*J3OJ2*SINIO*
     .                   (3.0D0+5.0D0*COSIO)/(1.0D0+COSIO)
            ELSE
                XLCOF  = -0.25D0*J3OJ2*SINIO*
     .                   (3.0D0+5.0D0*COSIO)/TEMP4
            END IF

            AYCOF  = -0.5D0*J3OJ2*SINIO
            DELMO  = (1.0D0+ETA*DCOS(MO))**3
            SINMAO = DSIN(MO)
            X7THM1 = 7.0D0*COSIO2-1.0D0

C
C           Deep Space Initialization
C
            IF ((TWOPI()/NO) .GE. 225.0D0) THEN

                DODEEP = .TRUE.
                DOSIMP = .TRUE.
                TC     = 0.0D0
                INCLM  = INCLO

C
C               Common.
C
                CALL ZZDSCM ( EPOCH, ECCO, ARGPO, TC, INCLO,
     .                  NODEO,  NO,
     .                  SNODM,  CNODM, SINIM, COSIM, SINOMM, COSOMM,
     .                  DAY,    E3,    EE2,   ECCM,  EMSQ,   GAM,
     .                  PEO,    PGHO,  PHO,   PINCO, PLO,
     .                  RTEMSQ, SE2,   SE3,   SGH2,  SGH3,   SGH4,
     .                  SH2,    SH3,   SI2,   SI3,   SL2,    SL3,
     .                  SL4,    S1,    S2,    S3,    S4,     S5,
     .                  S6,     S7,    SS1,   SS2,   SS3,    SS4,
     .                  SS5,    SS6,   SS7,   SZ1,   SZ2,    SZ3,
     .                  SZ11,   SZ12,  SZ13,  SZ21,  SZ22,   SZ23,
     .                  SZ31,   SZ32,  SZ33,  XGH2,  XGH3,   XGH4,
     .                  XH2,    XH3,   XI2,   XI3,   XL2,    XL3,
     .                  XL4,    XN,    Z1,    Z2,    Z3,     Z11,
     .                  Z12,    Z13,   Z21,   Z22,   Z23,    Z31,
     .                  Z32,    Z33,   ZMOL,  ZMOS )

C
C               Long period perturbations.
C
                CALL ZZDSPR( OPMODE, E3, EE2, PEO, PGHO, PHO,
     .                  PINCO,
     .                  PLO,  SE2,   SE3,   SGH2, SGH3, SGH4,
     .                  SH2,  SH3,   SI2,   SI3,  SL2,  SL3,
     .                  SL4,  TZERO, XGH2,  XGH3, XGH4, XH2,
     .                  XH3,  XI2,   XI3,   XL2,  XL3,  XL4,
     .                  ZMOL, ZMOS,  INCLM, DOINIT,
     .                  ECCO, INCLO, NODEO, ARGPO, MO )

                ARGPM  = 0.0D0
                NODEM  = 0.0D0
                MM     = 0.0D0

C
C               Initialization
C
                CALL ZZDSIN ( GEOPHS,
     .                   COSIM,  EMSQ,  ARGPO, S1,    S2,    S3,
     .                   S4,     S5,    SINIM, SS1,   SS2,   SS3,
     .                   SS4,    SS5,   SZ1,   SZ3,   SZ11,  SZ13,
     .                   SZ21,   SZ23,  SZ31,  SZ33,  TZERO, TC,
     .                   GSTO,   MO,    MDOT,  NO,    NODEO, NODEDOT,
     .                   XPIDOT, Z1,    Z3,    Z11,   Z13,   Z21,
     .                   Z23,    Z31,   Z33,   ECCO,  ECCSQ,
     .                   ECCM,   ARGPM, INCLM, MM,    XN,    NODEM,
     .                   IREZ,   ATIME, D2201, D2211, D3210, D3222,
     .                   D4410,  D4422, D5220, D5232, D5421, D5433,
     .                   DEDT,   DIDT,  DMDT,  DNDT,  DNODT, DOMDT,
     .                   DEL1,   DEL2,  DEL3,  XFACT, XLAMO, XLI,
     .                   XNI )

            END IF

C
C           Set variables if not deep space or rp < 220
C
            IF ( .NOT. DOSIMP ) THEN
               CC1SQ = CC1*CC1
               D2    = 4.0D0*AO*TSI*CC1SQ
               TEMP  = D2*TSI*CC1 / 3.0D0
               D3    = (17.0D0*AO + SFour) * TEMP
               D4    = 0.5D0*TEMP*AO*TSI*
     .                  (221.0D0*AO + 31.0D0*SFour)*CC1
               T3COF = D2 + 2.0D0*CC1SQ
               T4COF = 0.25D0* (3.0D0*D3+CC1*(12.0D0*D2+10.0D0*CC1SQ) )
               T5COF = 0.2D0* (3.0D0*D4 + 12.0D0*CC1*D3 + 6.0D0*D2*D2 +
     .                  15.0D0*CC1SQ* (2.0D0*D2 + CC1SQ) )

            END IF

         END IF

      DOINIT = .FALSE.

      CALL CHKOUT ( 'XXSGP4I' )
      RETURN



C$Procedure XXSGP4E ( SGP4 evaluator )

      ENTRY XXSGP4E ( T, STATE )

C$ Abstract
C
C     This procedure is the SGP4 prediction model from Space Command.
C     This is an updated and combined version of SGP4 and SDP4
C     originally published separately in Spacetrack report #3 [1].
C     This version follows the methodology from the 2006 AIAA
C     Vallado paper [4] describing the history and development of
C     \the code.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     TLE
C     Two line elements
C
C$ Declarations
C
C     Refer to Declarations section in ZZSGP4.
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     T          I   Time to evaluate state in minutes from epoch
C     STATE      O   Evaluated state
C
C$ Detailed_Input
C
C
C     T          The time in minutes from the elements epoch
C                at which to calculate a state from the
C                TLE set.
C
C$ Detailed_Output
C
C     STATE      the state produced by evaluating the input elements
C                at the input epoch T. Units are km and km/sec.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) SPICE(BADMEANMOTION) Mean motion less than 0.0
C
C     2) SPICE(BADMECCENTRICITY) signals when the mean eccentricity is
C        not bounded by [-0.001, 1].
C
C     3) SPICE(BADMSEMIMAJOR) signals when the mean semimajor axis
C        has value less-than .95.
C
C     4) SPICE(BADPECCENTRICITY) signals when the perturbed eccentricity
C       is not bounded by [0, 1].
C
C     5) SPICE(BADSEMILATUS) signals when the semi-latus rectum
C        has value less-than 0.
C
C     6) SPICE(ORBITDECAY) signals if the scaled orbit radial
C        distance has value less-than 1.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is based on the SGP4 code by David Vallado
C     corresponding to "Revisiting Spacetrack Report #3."
C     The intent is to maintain the original Vallado algorithm,
C     changing code only to meet NAIF format standards and to
C     integrate with SPICELIB.
C
C        1) Implemented error checks using SPICE error subsystem.
C           On detecting an error, control returns to the calling
C           routine. This behavior differs from the original
C           version.
C
C        2) Comments prefixed with "SGP4FIX" indicate a comment
C           from the code by Vallado et. al concerning a correction
C           to the STR#3 code.
C
C        3) Eliminated the use of COMMON blocks.
C
C        Removed getgravconst call, replaced with GEOPHS array.
C
C        Capitalize all variables.
C
C        radiusearthkm replaced with ER
C        VKmPerSec     replaced with KPS
C        r, v          replaced with STATE
C        method        replaced with DODEEP
C        isimp         replaced with DOSIMP
C        Error         eliminated
C        whichconst    eliminated, function provided by GEOPHS
C
C$ Examples
C
C     Refer to Examples section in ZZSGP4.
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C   [1] Hoots, F. R., and Roehrich, R. L. 1980. "Models for
C       Propagation of the NORAD Element Sets." Spacetrack Report #3.
C       U.S. Air Force: Aerospace Defense Command.
C
C   [2] Hoots, Felix R. "Spacetrack Report #6: Models for Propagation
C       of Space Command Element Sets." Space Command,
C       U. S. Air Force, CO.
C
C   [3] Hoots, Felix R., P. W. Schumacher, and R. A. Glover. 2004.
C       History of Analytical Orbit Modeling in the U. S. Space
C       Surveillance System. Journal of Guidance, Control, and
C       Dynamics. 27(2):174-185.
C
C   [4] Vallado, David, Crawford, Paul, Hujsak, Richard,
C       and Kelso, T.S. 2006. Revisiting Spacetrack Report #3. Paper
C       AIAA 2006-6753 presented at the AIAA/AAS Astrodynamics
C       Specialist Conference, August 21-24, 2006. Keystone, CO.
C
C$ Author_and_Institution
C
C     David Vallado   (AGI)
C     E. D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, SEP-15-2014 (EDW)
C
C        Based on routine SGP4, 28-JUN-2005 [4].
C
C-&

C$ Index_Entries
C
C  SGP4
C
C-&

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'XXSGP4E' )

C
C     Local constants. Keep compiler ok for warnings on
C     uninitialized variables
C
      X2O3   = 2.0D0/3.0D0
      MR     = 0.0D0
      COSEO1 = 1.0D0
      SINEO1 = 0.0D0

C
C     Set mathematical constants.
C
C     This code block replaces the call:
C
C     sgp4fix identify constants and allow alternate values.
C
C     CALL getgravconst( whichconst, tumin,
C     .                  mu, radiusearthkm, xke,
C     .                  j2, j3, j4, j3oj2 )
C

C
C     sgp4fix divisor for divide by zero check on inclination
C     the old check used 1.0D0 + cos(pi-1.0D-9), but then compared it to
C     1.5D-12, so the threshold was changed to 1.5D-12 for consistency.
C
      TEMP4 =  1.5D-12
      KPS   =  ER * XKE/60.0D0

C
C     UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG
C
      XMDF   = MO + MDOT*T
      OMGADF = ARGPO + ARGPDOT*T
      XNODDF = NODEO + NODEDOT*T
      ARGPM  = OMGADF
      MM     = XMDF
      T2     = T*T
      NODEM  = XNODDF + XNODCF*T2
      TEMPA  = 1.0D0 - CC1*T
      TEMPE  = BSTAR*CC4*T
      TEMPL  = T2COF*T2


      IF ( .NOT. DOSIMP ) THEN

         DELOMG = OMGCOF*T
         DELM   = XMCOF*(( 1.0D0+ETA*DCOS(XMDF) )**3-DELMO)
         TEMP   = DELOMG + DELM
         MM     = XMDF + TEMP
         ARGPM  = OMGADF - TEMP
         T3     = T2*T
         T4     = T3*T
         TEMPA  = TEMPA - D2*T2 - D3*T3 - D4*T4
         TEMPE  = TEMPE + BSTAR*CC5*(DSIN(MM) - SINMAO)
         TEMPL  = TEMPL + T3COF*T3 + T4*(T4COF + T*T5COF)

      END IF

      XN    = NO
      ECCM  = ECCO
      INCLM = INCLO

      IF (DODEEP) THEN
         TC     = T

         CALL ZZDSPC (  IREZ,  D2201, D2211,   D3210, D3222, D4410,
     .                  D4422, D5220, D5232,   D5421, D5433, DEDT,
     .                  DEL1,  DEL2,  DEL3,    DIDT,  DMDT,  DNODT,
     .                  DOMDT, ARGPO, ARGPDOT, T,     TC,    GSTO,
     .                  XFACT, XLAMO, NO,      ATIME, ECCM,
     .                  ARGPM, INCLM, XLI,     MM,    XNI,
     .                  NODEM, DNDT,  XN  )

      END IF

C
C     Mean motion less than 0.0.
C
      IF(XN .LE. 0.0D0) THEN

         CALL SETMSG ( 'Mean motion less-than zero. '
     .    //           'This error may indicate a bad TLE set.'     )
         CALL SIGERR ( 'SPICE(BADMEANMOTION)'                       )
         CALL CHKOUT ( 'XXSGP4E'                                    )
         RETURN

      END IF

      AM   = (XKE/XN)**X2O3*TEMPA**2
      XN   = XKE/AM**1.5D0
      ECCM = ECCM-TEMPE

C
C     Fix tolerance for error recognition. Vallado code used
C     a lower limit of -0.001. This value apparently prevents
C     an error signal due to roundoff error.
C
      IF (ECCM .GE. 1.0D0 .OR. ECCM .LT. -0.001D0) THEN

         CALL SETMSG ( 'Mean eccentricity value, #, beyond '
     .    //           'allowed bounds [-0.001,1.0). This '
     .    //           'error may indicate a bad TLE set.'          )
         CALL ERRDP  ( '#', Eccm                                    )
         CALL SIGERR ( 'SPICE(BADMECCENTRICITY)'                     )
         CALL CHKOUT ( 'XXSGP4E'                                    )
         RETURN

      END IF

      IF ( AM .LT. 0.95 ) THEN

         CALL SETMSG ( 'Mean semi-major axis value, #, below allowed '
     .    //           'minimum of 0.95. This error may indicate '
     .    //           'a bad TLE set or a decayed orbit.'          )
         CALL ERRDP  ( '#', Eccm                                    )
         CALL SIGERR ( 'SPICE(BADMSEMIMAJOR)'                        )
         CALL CHKOUT ( 'XXSGP4E'                                    )
         RETURN

      END IF


C
C     sgp4fix change test condition for eccentricity
C
      IF (ECCM .LT. 1.0D-6) ECCM = 1.0D-6

      MM     = MM+NO*TEMPL
      XLM    = MM+ARGPM+NODEM
      EMSQ   = ECCM*ECCM
      TEMP   = 1.0D0 - EMSQ
      NODEM  = DMOD( NODEM, TWOPI())
      ARGPM  = DMOD( ARGPM, TWOPI())
      XLM    = DMOD( XLM,   TWOPI())
      MM     = DMOD( XLM - ARGPM - NODEM, TWOPI())

C
C     Compute extra mean quantities
C
      SINIM  = DSIN(INCLM)
      COSIM  = DCOS(INCLM)

C
C     Add lunar-solar periodics
C
      ECCP   = ECCM
      XINCP  = INCLM
      ARGPP  = ARGPM
      NODEP  = NODEM
      MP     = MM
      SINIP  = SINIM
      COSIP  = COSIM

C
C     Use deep space perturbation if indicated.
C

      IF (DODEEP) THEN

         CALL ZZDSPR ( SVMODE, E3,  EE2,   PEO,    PGHO, PHO, PINCO,
     .                 PLO,  SE2,   SE3,   SGH2,   SGH3, SGH4,
     .                 SH2,  SH3,   SI2,   SI3,    SL2,  SL3,
     .                 SL4,  T,     XGH2,  XGH3,   XGH4, XH2,
     .                 XH3,  XI2,   XI3,   XL2,    XL3,  XL4,
     .                 ZMOL, ZMOS,  INCLO, .FALSE.,
     .                 ECCP, XINCP, NODEP, ARGPP, MP )

         IF(XINCP .LT. 0.0D0) THEN
            XINCP  = -XINCP
            NODEP  = NODEP + PI()
            ARGPP  = ARGPP - PI()
         END IF

         IF ( (ECCP .LT. 0.0D0) .OR.
     .        (ECCP .GT. 1.0D0) ) THEN

            CALL SETMSG ( 'Perturbed eccentricity value, #, beyond '
     .       //           'allowed bounds [0,1]. This error may '
     .       //           'indicate a bad TLE set.'                )
            CALL ERRDP  ( '#', ECCP                                )
            CALL SIGERR ( 'SPICE(BADPECCENTRICITY)'                )
            CALL CHKOUT ( 'XXSGP4E'                                )
            RETURN

         END IF

      END IF

C
C     Update for long period periodics if a deep space trajectory.
C

      IF (DODEEP) THEN

          SINIP =  DSIN(XINCP)
          COSIP =  DCOS(XINCP)
          AYCOF = -0.5D0*J3OJ2*SINIP

C
C         sgp4fix for divide by zero with xincp = 180 deg
C

          IF (DABS(COSIP+1.D0).GT. 1.5D-12) THEN
             XLCOF  = -0.25D0*J3OJ2*SINIP*
     .                 (3.0D0+5.0D0*COSIP)/(1.0D0+COSIP)
          ELSE
             XLCOF  = -0.25D0*J3OJ2*SINIP*
     .                 (3.0D0+5.0D0*COSIP)/TEMP4
          END IF

      END IF

      AXNL = ECCP*DCOS(ARGPP)
      TEMP = 1.0D0 / (AM*(1.0D0-ECCP*ECCP))
      AYNL = ECCP*DSIN(ARGPP) + TEMP*AYCOF
      XL   = MP + ARGPP + NODEP + TEMP*XLCOF*AXNL

C
C     Solve Kepler's equation.
C
      U    = DMOD(XL-NODEP,TWOPI())
      EO1  = U
      ITER = 0

C
C     sgp4fix for Kepler iteration the following iteration needs
C     better limits on corrections
C

      TEMP = 9999.9D0

      DO WHILE ((TEMP.GE.1.0D-12).AND.(ITER.LT.10))

         ITER  = ITER+1
         SINEO1= DSIN(EO1)
         COSEO1= DCOS(EO1)
         TEM5  = 1.0D0 - COSEO1*AXNL - SINEO1*AYNL
         TEM5  = (U - AYNL*COSEO1 + AXNL*SINEO1 - EO1) / TEM5
         TEMP  = DABS(TEM5)

C
C        Stop excessive correction.
C
         IF(TEMP.GT.1.0D0) TEM5=TEM5/TEMP

         EO1   = EO1+TEM5

      END DO

C
C     Short period preliminary quantities.
C
      ECOSE = AXNL*COSEO1+AYNL*SINEO1
      ESINE = AXNL*SINEO1-AYNL*COSEO1
      EL2   = AXNL*AXNL+AYNL*AYNL
      PL    = AM*(1.0D0-EL2)

C
C     Error check for semi-latus rectum < 0.0
C
      IF ( PL .LT. 0.0D0 ) THEN

         CALL SETMSG ( 'Semi-latus rectum less-than zero.'          )
         CALL SIGERR ( 'SPICE(BADSEMILATUS)'                        )
         CALL CHKOUT ( 'XXSGP4E'                                    )
         RETURN

      END IF

      RL    = AM*(1.0D0-ECOSE)
      RDOTL = DSQRT(AM)*ESINE/RL
      RVDOTL= DSQRT(PL)/RL
      BETAL = DSQRT(1.0D0-EL2)
      TEMP  = ESINE/(1.0D0+BETAL)
      SINU  = AM/RL*(SINEO1-AYNL-AXNL*TEMP)
      COSU  = AM/RL*(COSEO1-AXNL+AYNL*TEMP)
      SU    = DATAN2(SINU,COSU)
      SIN2U = (COSU+COSU)*SINU
      COS2U = 1.0D0-2.0D0*SINU*SINU
      TEMP  = 1.0D0/PL
      TEMP1 = 0.5D0*J2*TEMP
      TEMP2 = TEMP1*TEMP

C
C     Update for short period periodics if a deep space trajectory.
C
      IF (DODEEP) THEN
          COSISQ = COSIP*COSIP
          CON41  = 3.0D0*COSISQ - 1.0D0
          X1MTH2 = 1.0D0 - COSISQ
          X7THM1 = 7.0D0*COSISQ - 1.0D0
      END IF

      MR    = RL*(1.0D0 - 1.5D0*TEMP2*BETAL*CON41) +
     .           0.5D0*TEMP1*X1MTH2*COS2U
      SU    = SU - 0.25D0*TEMP2*X7THM1*SIN2U
      XNODE = NODEP + 1.5D0*TEMP2*COSIP*SIN2U
      XINC  = XINCP + 1.5D0*TEMP2*COSIP*SINIP*COS2U
      MV    = RDOTL - XN*TEMP1*X1MTH2*SIN2U / XKE
      RVDOT = RVDOTL + XN*TEMP1* (X1MTH2*COS2U+1.5D0*CON41) / XKE

C
C     Orientation vectors.
C

      SINSU=  DSIN(SU)
      COSSU=  DCOS(SU)
      SNOD =  DSIN(XNODE)
      CNOD =  DCOS(XNODE)
      SINI =  DSIN(XINC)
      COSI =  DCOS(XINC)
      XMX  = -SNOD*COSI
      XMY  =  CNOD*COSI
      UX   =  XMX*SINSU + CNOD*COSSU
      UY   =  XMY*SINSU + SNOD*COSSU
      UZ   =  SINI*SINSU
      VX   =  XMX*COSSU - CNOD*SINSU
      VY   =  XMY*COSSU - SNOD*SINSU
      VZ   =  SINI*COSSU

C
C     Position and velocity.
C

      STATE(1) = MR*UX * ER
      STATE(2) = MR*UY * ER
      STATE(3) = MR*UZ * ER
      STATE(4) = (MV*UX + RVDOT*VX) * KPS
      STATE(5) = (MV*UY + RVDOT*VY) * KPS
      STATE(6) = (MV*UZ + RVDOT*VZ) * KPS

C
C     sgp4fix for decaying satellites
C
C     Place this test here to ensure evaluation of STATE.
C     The result may be physically invalid.
C
      IF (MR .LT. 1.0D0) THEN

         CALL SETMSG ( 'Satellite has decayed.' )
         CALL SIGERR ( 'SPICE(ORBITDECAY)'      )
         CALL CHKOUT ( 'XXSGP4E'                 )
         RETURN

      END IF

      CALL CHKOUT ( 'XXSGP4E' )

      RETURN

      END

