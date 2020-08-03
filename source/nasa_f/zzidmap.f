C$Procedure ZZIDMAP ( Private --- SPICE body ID/name assignments )

      SUBROUTINE  ZZIDMAP( BLTCOD, BLTNAM )
      IMPLICIT NONE

C$ Abstract
C
C     The default SPICE body/ID mapping assignments available
C     to the SPICE library.
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
C     naif_ids.req
C
C$ Keywords
C
C     Body mappings.
C
C$ Declarations

      INCLUDE              'zzbodtrn.inc'

      INTEGER              BLTCOD(NPERM)
      CHARACTER*(MAXL)     BLTNAM(NPERM)

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     BLTCOD     O  List of default integer ID codes
C     BLTNAM     O  List of default names
C     NPERM      P  Number of name/ID mappings
C
C$ Detailed_Input
C
C     None.
C
C$ Detailed_Output
C
C     BLTCOD     The array of NPERM elements listing the body ID codes.
C
C     BLTNAM     The array of NPERM elements listing the body names 
C                corresponding to the ID entry in BLTCOD
C
C$ Parameters
C
C     NPERM      The length of both BLTCOD, BLTNAM
C                (read from zzbodtrn.inc).
C
C$ Exceptions
C
C     None.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     Each ith entry of BLTCOD maps to the ith entry of BLTNAM.
C
C$ Examples
C
C     Simple to use, a call the ZZIDMAP returns the arrays defining the
C     name/ID mappings.
C
C
C        INCLUDE            'zzbodtrn.inc'
C
C        INTEGER             ID  ( NPERM )
C        CHARACTER*(MAXL)    NAME( NPERM )
C
C        CALL ZZIDMAP( ID, NAME )
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     E.D. Wright, 04-APR-2017 (JPL)
C
C$ Version
C
C-    SPICELIB 1.0.9 04-APR-2017 (EDW)
C
C        Added information stating the frames subsystem performs
C        frame ID-name mappings and the DSK subsystem performs
C        surface ID-name mappings.
C
C        Edited body/ID assignment format to indicate whitespace
C        between 'NAME' and Comments.
C
C     Added:
C
C             -302   HELIOS 2
C             -301   HELIOS 1
C             -198   NASA-ISRO SAR MISSION
C             -198   NISAR
C             -159   EURC
C             -159   EUROPA CLIPPER
C             -152   CH2
C             -152   CHANDRAYAAN-2
C             -143   TRACE GAS ORBITER
C             -143   TGO
C             -143   EXOMARS 2016 TGO
C             -117   EDL DEMONSTRATOR MODULE
C             -117   EDM
C             -117   EXOMARS 2016 EDM
C              -76   CURIOSITY
C              -69   PSYC
C              -66   MCOB
C              -66   MARCO-B
C              -65   MCOA
C              -65   MARCO-A
C              -62   EMM
C              -62   EMIRATES MARS MISSION
C              -49   LUCY
C              -28   JUPITER ICY MOONS EXPLORER
C              -28   JUICE
C              553   DIA
C          2000016   PSYCHE
C          2101955   BENNU
C
C     Removed assignments:
C
C             -159   EUROPA ORBITER
C              -69   MPO
C              -69   MERCURY PLANETARY ORBITER
C
C     Modified assignments
C
C             -121   MERCURY PLANETARY ORBITER
C             -121   MPO
C             -121   BEPICOLOMBO MPO
C              -68   MERCURY MAGNETOSPHERIC ORBITER
C              -68   MMO
C              -68   BEPICOLOMBO MMO
C
C-    SPICELIB 1.0.8 06-MAY-2014 (EDW)
C
C         Edited text comments in Asteroids section and Comets section.
C
C         Eliminated "PI" IAU Number from "CHARON" description.
C
C         HYROKKIN (644) spelling corrected to HYRROKKIN.
C
C     Added:
C
C             -750   SPRINT-AS
C             -189   NSYT
C             -189   INSIGHT
C             -170   JWST
C             -170   JAMES WEBB SPACE TELESCOPE
C             -144   SOLO
C             -144   SOLAR ORBITER
C              -96   SPP
C              -96   SOLAR PROBE PLUS
C              -64   ORX
C              -64   OSIRIS-REX
C              -54   ARM
C              -54   ASTEROID RETRIEVAL MISSION
C              -12   LADEE
C               -3   MOM
C               -3   MARS ORBITER MISSION
C                0   SOLAR_SYSTEM_BARYCENTER
C                1   MERCURY_BARYCENTER
C                2   VENUS_BARYCENTER
C                3   EARTH_BARYCENTER
C                4   MARS_BARYCENTER
C                5   JUPITER_BARYCENTER
C                6   SATURN_BARYCENTER
C                7   URANUS_BARYCENTER
C                8   NEPTUNE_BARYCENTER
C                9   PLUTO_BARYCENTER
C              644   HYRROKKIN
C              904   KERBEROS
C              905   STYX
C          1003228   C/2013 A1
C          1003228   SIDING SPRING
C          2000002   PALLAS
C          2000511   DAVIDA
C
C     Removed assignments:
C
C             -486   HERSCHEL
C             -489   PLANCK
C             -187   SOLAR PROBE
C
C-    SPICELIB 1.0.7 20-MAY-2010 (EDW)
C
C        Edit to vehicle ID list to correct -76 not in proper
C        numerical (descending) order.
C
C     Added:
C
C               -5   AKATSUKI
C               -5   VCO
C             -121   BEPICOLOMBO
C             -177   GRAIL-A
C             -181   GRAIL-B
C             -202   MAVEN
C             -205   SOIL MOISTURE ACTIVE AND PASSIVE
C             -205   SMAP
C             -362   RADIATION BELT STORM PROBE A
C             -362   RBSP_A
C             -363   RADIATION BELT STORM PROBE B
C             -363   RBSP_B
C              550   HERSE
C              653   AEGAEON
C          1000093   TEMPEL_1
C          2000021   LUTETIA
C          2004179   TOUTATIS
C
C-    SPICELIB 1.0.6 08-APR-2009 (EDW)
C
C     Added:
C
C               -5   PLC
C               -5   PLANET-C
C              -68   MMO
C              -68   MERCURY MAGNETOSPHERIC ORBITER
C              -69   MPO
C              -69   MERCURY PLANETARY ORBITER
C          2002867   STEINS
C             -140   EPOCH
C             -140   DIXI
C
C-    SPICELIB 1.0.5 09-JAN-2008 (EDW)
C
C     Added:
C
C              -18   LCROSS
C              -29   NEXT
C              -86   CH1
C              -86   CHANDRAYAAN-1
C             -131   KAGUYA
C             -140   EPOXI
C             -151   CHANDRA
C             -187   SOLAR PROBE
C              636   AEGIR
C              637   BEBHIONN
C              638   BERGELMIR
C              639   BESTLA
C              640   FARBAUTI
C              641   FENRIR
C              642   FORNJOT
C              643   HATI
C              644   HYROKKIN
C              645   KARI
C              646   LOGE
C              647   SKOLL
C              648   SURTUR
C              649   ANTHE
C              650   JARNSAXA
C              651   GREIP
C              652   TARQEQ
C              809   HALIMEDE
C              810   PSAMATHE
C              811   SAO
C              812   LAOMEDEIA
C              813   NESO
C
C     NAIF modified the Jovian system listing to conform to the
C     current (as of this date) name/body mapping.
C
C              540   MNEME
C              541   AOEDE
C              542   THELXINOE
C              543   ARCHE
C              544   KALLICHORE
C              545   HELIKE
C              546   CARPO
C              547   EUKELADE
C              548   CYLLENE
C              549   KORE
C
C     Removed assignments:
C
C             -172   SPACETECH-3 COMBINER
C             -174   PLUTO-KUIPER EXPRESS
C             -175   PLUTO-KUIPER EXPRESS SIMULATION
C             -205   SPACETECH-3 COLLECTOR
C              514   1979J2
C              515   1979J1
C              516   1979J3
C              610   1980S1
C              611   1980S3
C              612   1980S6
C              613   1980S13
C              614   1980S25
C              615   1980S28
C              616   1980S27
C              617   1980S26
C              706   1986U7
C              707   1986U8
C              708   1986U9
C              709   1986U4
C              710   1986U6
C              711   1986U3
C              712   1986U1
C              713   1986U2
C              714   1986U5
C              715   1985U1
C              718   1986U10
C              901   1978P1
C
C     Spelling correction:
C
C        MAGACLITE to MEGACLITE
C
C     Rename:
C
C        ERRIAPO to ERRIAPUS
C        STV-1 to STV51
C        STV-2 to STV52
C        STV-3 to STV53
C
C
C-    SPICELIB 1.0.4 01-NOV-2006 (EDW)
C
C     NAIF removed several provisional name/ID mappings from
C     the Jovian system listing:
C
C     539         'HEGEMONE'              JXXXIX
C     540         'MNEME'                 JXL
C     541         'AOEDE'                 JXLI
C     542         'THELXINOE'             JXLII
C     543         'ARCHE'                 JXLIII
C     544         'KALLICHORE'            JXLIV
C     545         'HELIKE'                JXLV
C     546         'CARPO'                 JXLVI
C     547         'EUKELADE'              JXLVII
C     548         'CYLLENE'               JXLVIII
C
C     The current mapping set for the range 539-561:
C
C              540   ARCHE
C              541   EUKELADE
C              546   HELIKE
C              547   AOEDE
C              548   HEGEMONE
C              551   KALLICHORE
C              553   CYLLENE
C              560   CARPO
C              561   MNEME
C
C     The new mapping leaves the IDs 539, 542-545, 549, 550, 552,
C     554-559 unassigned.
C
C     Added:
C
C              635   DAPHNIS
C              722   FRANCISCO
C              723   MARGARET
C              724   FERDINAND
C              725   PERDITA
C              726   MAB
C              727   CUPID
C              -61   JUNO
C              -76   MSL
C              -76   MARS SCIENCE LABORATORY
C             -212   STV-1
C             -213   STV-2
C             -214   STV-3
C              902   NIX
C              903   HYDRA
C             -85    LRO
C             -85    LUNAR RECON ORBITER
C             -85    LUNAR RECONNAISSANCE ORBITER
C
C     Spelling correction
C
C              632   METHODE to METHONE
C
C-    SPICELIB 1.0.3 14-NOV-2005 (EDW)
C
C     Added:
C
C              539   HEGEMONE
C              540   MNEME
C              541   AOEDE
C              542   THELXINOE
C              543   ARCHE
C              544   KALLICHORE
C              545   HELIKE
C              546   CARPO
C              547   EUKELADE
C              548   CYLLENE
C              631   NARVI
C              632   METHODE
C              633   PALLENE
C              634   POLYDEUCES
C          2025143   ITOKAWA
C              -98   NEW HORIZONS
C             -248   VENUS EXPRESS, VEX
C             -500   RSAT, SELENE Relay Satellite, SELENE Rstar, Rstar
C             -502   VSAT, SELENE VLBI Radio Satellite,
C                    SELENE VRAD Satellite, SELENE Vstar
C           399064   DSS-64
C
C      Change in spelling:
C
C              623   SUTTUNG to SUTTUNGR
C              627   SKADI   to SKATHI
C              630   THRYM   to THRYMR
C
C-    SPICELIB 1.0.2 20-DEC-2004 (EDW)
C
C     Added:
C
C           Due to the previous definition of Parkes with DSS-05,
C           the Parkes ID remains 399005.
C
C             -486   HERSCHEL
C             -489   PLANCK
C           399049   DSS-49
C           399055   DSS-55
C             -203   DAWN
C          1000012   67P/CHURYUMOV-GERASIMENKO (1969 R1)
C          1000012   CHURYUMOV-GERASIMENKO
C          398989    NOTO
C             -84    PHOENIX
C            -131    SELENE
C            -238    SMART-1, S1, SM1, SMART1
C            -130    HAYABUSA
C
C-    SPICELIB 1.0.1 19-DEC-2003 (EDW)
C
C     Added:
C              -79   SPITZER
C          2000216   KLEOPATRA
C
C-    SPICELIB 1.0.0 27-JUL-2003 (EDW)
C
C     Added:
C              -47   GNS
C              -74   MRO
C              -74   MARS RECON ORBITER
C             -130   MUSES-C
C             -142   TERRA
C             -154   AQUA
C             -159   EUROPA ORBITER
C             -190   SIM
C             -198   INTEGRAL
C             -227   KEPLER
C             -234   STEREO AHEAD
C             -235   STEREO BEHIND
C             -253   OPPORTUNITY
C             -254   SPIRIT
C              528   AUTONOE
C              529   THYONE
C              530   HERMIPPE
C              531   AITNE
C              532   EURYDOME
C              533   EUANTHE
C              534   EUPORIE
C              535   ORTHOSIE
C              536   SPONDE
C              537   KALE
C              538   PASITHEE
C              619   YMIR
C              620   PAALIAQ
C              621   TARVOS
C              622   IJIRAQ
C              623   SUTTUNG
C              624   KIVIUQ
C              625   MUNDILFARI
C              626   ALBIORIX
C              627   SKADI
C              628   ERRIAPO
C              629   SIARNAQ
C              630   THRYM
C              718   PROSPERO
C              719   SETEBOS
C              720   STEPHANO
C              721   TRINCULO
C           398990   NEW NORCIA
C          2431011   DACTYL
C          2000001   CERES
C          2000004   VESTA
C
C     Renamed:
C
C              -25   LPM to
C              -25   LP
C
C             -180   MUSES-C to
C             -130   MUSES-B
C
C             -172   STARLIGHT COMBINER to
C             -172   SPACETECH-3 COMBINER
C
C             -205   STARLIGHT COLLECTOR to
C             -205   SPACETECH-3 COLLECTOR
C
C      Removed:
C             -172   SLCOMB
C
C
C-&

C$ Index_Entries
C
C     body ID mapping
C
C-&

C
C     A script generates this file. Do not edit by hand.
C     Edit the creation script to modify the contents of
C     ZZIDMAP.
C

      BLTCOD(1) =   0
      BLTNAM(1) =  'SOLAR_SYSTEM_BARYCENTER'

      BLTCOD(2) =   0
      BLTNAM(2) =  'SSB'

      BLTCOD(3) =   0
      BLTNAM(3) =  'SOLAR SYSTEM BARYCENTER'

      BLTCOD(4) =   1
      BLTNAM(4) =  'MERCURY_BARYCENTER'

      BLTCOD(5) =   1
      BLTNAM(5) =  'MERCURY BARYCENTER'

      BLTCOD(6) =   2
      BLTNAM(6) =  'VENUS_BARYCENTER'

      BLTCOD(7) =   2
      BLTNAM(7) =  'VENUS BARYCENTER'

      BLTCOD(8) =   3
      BLTNAM(8) =  'EARTH_BARYCENTER'

      BLTCOD(9) =   3
      BLTNAM(9) =  'EMB'

      BLTCOD(10) =   3
      BLTNAM(10) =  'EARTH MOON BARYCENTER'

      BLTCOD(11) =   3
      BLTNAM(11) =  'EARTH-MOON BARYCENTER'

      BLTCOD(12) =   3
      BLTNAM(12) =  'EARTH BARYCENTER'

      BLTCOD(13) =   4
      BLTNAM(13) =  'MARS_BARYCENTER'

      BLTCOD(14) =   4
      BLTNAM(14) =  'MARS BARYCENTER'

      BLTCOD(15) =   5
      BLTNAM(15) =  'JUPITER_BARYCENTER'

      BLTCOD(16) =   5
      BLTNAM(16) =  'JUPITER BARYCENTER'

      BLTCOD(17) =   6
      BLTNAM(17) =  'SATURN_BARYCENTER'

      BLTCOD(18) =   6
      BLTNAM(18) =  'SATURN BARYCENTER'

      BLTCOD(19) =   7
      BLTNAM(19) =  'URANUS_BARYCENTER'

      BLTCOD(20) =   7
      BLTNAM(20) =  'URANUS BARYCENTER'

      BLTCOD(21) =   8
      BLTNAM(21) =  'NEPTUNE_BARYCENTER'

      BLTCOD(22) =   8
      BLTNAM(22) =  'NEPTUNE BARYCENTER'

      BLTCOD(23) =   9
      BLTNAM(23) =  'PLUTO_BARYCENTER'

      BLTCOD(24) =   9
      BLTNAM(24) =  'PLUTO BARYCENTER'

      BLTCOD(25) =   10
      BLTNAM(25) =  'SUN'

      BLTCOD(26) =   199
      BLTNAM(26) =  'MERCURY'

      BLTCOD(27) =   299
      BLTNAM(27) =  'VENUS'

      BLTCOD(28) =   399
      BLTNAM(28) =  'EARTH'

      BLTCOD(29) =   301
      BLTNAM(29) =  'MOON'

      BLTCOD(30) =   499
      BLTNAM(30) =  'MARS'

      BLTCOD(31) =   401
      BLTNAM(31) =  'PHOBOS'

      BLTCOD(32) =   402
      BLTNAM(32) =  'DEIMOS'

      BLTCOD(33) =   599
      BLTNAM(33) =  'JUPITER'

      BLTCOD(34) =   501
      BLTNAM(34) =  'IO'

      BLTCOD(35) =   502
      BLTNAM(35) =  'EUROPA'

      BLTCOD(36) =   503
      BLTNAM(36) =  'GANYMEDE'

      BLTCOD(37) =   504
      BLTNAM(37) =  'CALLISTO'

      BLTCOD(38) =   505
      BLTNAM(38) =  'AMALTHEA'

      BLTCOD(39) =   506
      BLTNAM(39) =  'HIMALIA'

      BLTCOD(40) =   507
      BLTNAM(40) =  'ELARA'

      BLTCOD(41) =   508
      BLTNAM(41) =  'PASIPHAE'

      BLTCOD(42) =   509
      BLTNAM(42) =  'SINOPE'

      BLTCOD(43) =   510
      BLTNAM(43) =  'LYSITHEA'

      BLTCOD(44) =   511
      BLTNAM(44) =  'CARME'

      BLTCOD(45) =   512
      BLTNAM(45) =  'ANANKE'

      BLTCOD(46) =   513
      BLTNAM(46) =  'LEDA'

      BLTCOD(47) =   514
      BLTNAM(47) =  'THEBE'

      BLTCOD(48) =   515
      BLTNAM(48) =  'ADRASTEA'

      BLTCOD(49) =   516
      BLTNAM(49) =  'METIS'

      BLTCOD(50) =   517
      BLTNAM(50) =  'CALLIRRHOE'

      BLTCOD(51) =   518
      BLTNAM(51) =  'THEMISTO'

      BLTCOD(52) =   519
      BLTNAM(52) =  'MAGACLITE'

      BLTCOD(53) =   520
      BLTNAM(53) =  'TAYGETE'

      BLTCOD(54) =   521
      BLTNAM(54) =  'CHALDENE'

      BLTCOD(55) =   522
      BLTNAM(55) =  'HARPALYKE'

      BLTCOD(56) =   523
      BLTNAM(56) =  'KALYKE'

      BLTCOD(57) =   524
      BLTNAM(57) =  'IOCASTE'

      BLTCOD(58) =   525
      BLTNAM(58) =  'ERINOME'

      BLTCOD(59) =   526
      BLTNAM(59) =  'ISONOE'

      BLTCOD(60) =   527
      BLTNAM(60) =  'PRAXIDIKE'

      BLTCOD(61) =   528
      BLTNAM(61) =  'AUTONOE'

      BLTCOD(62) =   529
      BLTNAM(62) =  'THYONE'

      BLTCOD(63) =   530
      BLTNAM(63) =  'HERMIPPE'

      BLTCOD(64) =   531
      BLTNAM(64) =  'AITNE'

      BLTCOD(65) =   532
      BLTNAM(65) =  'EURYDOME'

      BLTCOD(66) =   533
      BLTNAM(66) =  'EUANTHE'

      BLTCOD(67) =   534
      BLTNAM(67) =  'EUPORIE'

      BLTCOD(68) =   535
      BLTNAM(68) =  'ORTHOSIE'

      BLTCOD(69) =   536
      BLTNAM(69) =  'SPONDE'

      BLTCOD(70) =   537
      BLTNAM(70) =  'KALE'

      BLTCOD(71) =   538
      BLTNAM(71) =  'PASITHEE'

      BLTCOD(72) =   539
      BLTNAM(72) =  'HEGEMONE'

      BLTCOD(73) =   540
      BLTNAM(73) =  'MNEME'

      BLTCOD(74) =   541
      BLTNAM(74) =  'AOEDE'

      BLTCOD(75) =   542
      BLTNAM(75) =  'THELXINOE'

      BLTCOD(76) =   543
      BLTNAM(76) =  'ARCHE'

      BLTCOD(77) =   544
      BLTNAM(77) =  'KALLICHORE'

      BLTCOD(78) =   545
      BLTNAM(78) =  'HELIKE'

      BLTCOD(79) =   546
      BLTNAM(79) =  'CARPO'

      BLTCOD(80) =   547
      BLTNAM(80) =  'EUKELADE'

      BLTCOD(81) =   548
      BLTNAM(81) =  'CYLLENE'

      BLTCOD(82) =   549
      BLTNAM(82) =  'KORE'

      BLTCOD(83) =   550
      BLTNAM(83) =  'HERSE'

      BLTCOD(84) =   553
      BLTNAM(84) =  'DIA'

      BLTCOD(85) =   699
      BLTNAM(85) =  'SATURN'

      BLTCOD(86) =   601
      BLTNAM(86) =  'MIMAS'

      BLTCOD(87) =   602
      BLTNAM(87) =  'ENCELADUS'

      BLTCOD(88) =   603
      BLTNAM(88) =  'TETHYS'

      BLTCOD(89) =   604
      BLTNAM(89) =  'DIONE'

      BLTCOD(90) =   605
      BLTNAM(90) =  'RHEA'

      BLTCOD(91) =   606
      BLTNAM(91) =  'TITAN'

      BLTCOD(92) =   607
      BLTNAM(92) =  'HYPERION'

      BLTCOD(93) =   608
      BLTNAM(93) =  'IAPETUS'

      BLTCOD(94) =   609
      BLTNAM(94) =  'PHOEBE'

      BLTCOD(95) =   610
      BLTNAM(95) =  'JANUS'

      BLTCOD(96) =   611
      BLTNAM(96) =  'EPIMETHEUS'

      BLTCOD(97) =   612
      BLTNAM(97) =  'HELENE'

      BLTCOD(98) =   613
      BLTNAM(98) =  'TELESTO'

      BLTCOD(99) =   614
      BLTNAM(99) =  'CALYPSO'

      BLTCOD(100) =   615
      BLTNAM(100) =  'ATLAS'

      BLTCOD(101) =   616
      BLTNAM(101) =  'PROMETHEUS'

      BLTCOD(102) =   617
      BLTNAM(102) =  'PANDORA'

      BLTCOD(103) =   618
      BLTNAM(103) =  'PAN'

      BLTCOD(104) =   619
      BLTNAM(104) =  'YMIR'

      BLTCOD(105) =   620
      BLTNAM(105) =  'PAALIAQ'

      BLTCOD(106) =   621
      BLTNAM(106) =  'TARVOS'

      BLTCOD(107) =   622
      BLTNAM(107) =  'IJIRAQ'

      BLTCOD(108) =   623
      BLTNAM(108) =  'SUTTUNGR'

      BLTCOD(109) =   624
      BLTNAM(109) =  'KIVIUQ'

      BLTCOD(110) =   625
      BLTNAM(110) =  'MUNDILFARI'

      BLTCOD(111) =   626
      BLTNAM(111) =  'ALBIORIX'

      BLTCOD(112) =   627
      BLTNAM(112) =  'SKATHI'

      BLTCOD(113) =   628
      BLTNAM(113) =  'ERRIAPUS'

      BLTCOD(114) =   629
      BLTNAM(114) =  'SIARNAQ'

      BLTCOD(115) =   630
      BLTNAM(115) =  'THRYMR'

      BLTCOD(116) =   631
      BLTNAM(116) =  'NARVI'

      BLTCOD(117) =   632
      BLTNAM(117) =  'METHONE'

      BLTCOD(118) =   633
      BLTNAM(118) =  'PALLENE'

      BLTCOD(119) =   634
      BLTNAM(119) =  'POLYDEUCES'

      BLTCOD(120) =   635
      BLTNAM(120) =  'DAPHNIS'

      BLTCOD(121) =   636
      BLTNAM(121) =  'AEGIR'

      BLTCOD(122) =   637
      BLTNAM(122) =  'BEBHIONN'

      BLTCOD(123) =   638
      BLTNAM(123) =  'BERGELMIR'

      BLTCOD(124) =   639
      BLTNAM(124) =  'BESTLA'

      BLTCOD(125) =   640
      BLTNAM(125) =  'FARBAUTI'

      BLTCOD(126) =   641
      BLTNAM(126) =  'FENRIR'

      BLTCOD(127) =   642
      BLTNAM(127) =  'FORNJOT'

      BLTCOD(128) =   643
      BLTNAM(128) =  'HATI'

      BLTCOD(129) =   644
      BLTNAM(129) =  'HYRROKKIN'

      BLTCOD(130) =   645
      BLTNAM(130) =  'KARI'

      BLTCOD(131) =   646
      BLTNAM(131) =  'LOGE'

      BLTCOD(132) =   647
      BLTNAM(132) =  'SKOLL'

      BLTCOD(133) =   648
      BLTNAM(133) =  'SURTUR'

      BLTCOD(134) =   649
      BLTNAM(134) =  'ANTHE'

      BLTCOD(135) =   650
      BLTNAM(135) =  'JARNSAXA'

      BLTCOD(136) =   651
      BLTNAM(136) =  'GREIP'

      BLTCOD(137) =   652
      BLTNAM(137) =  'TARQEQ'

      BLTCOD(138) =   653
      BLTNAM(138) =  'AEGAEON'

      BLTCOD(139) =   799
      BLTNAM(139) =  'URANUS'

      BLTCOD(140) =   701
      BLTNAM(140) =  'ARIEL'

      BLTCOD(141) =   702
      BLTNAM(141) =  'UMBRIEL'

      BLTCOD(142) =   703
      BLTNAM(142) =  'TITANIA'

      BLTCOD(143) =   704
      BLTNAM(143) =  'OBERON'

      BLTCOD(144) =   705
      BLTNAM(144) =  'MIRANDA'

      BLTCOD(145) =   706
      BLTNAM(145) =  'CORDELIA'

      BLTCOD(146) =   707
      BLTNAM(146) =  'OPHELIA'

      BLTCOD(147) =   708
      BLTNAM(147) =  'BIANCA'

      BLTCOD(148) =   709
      BLTNAM(148) =  'CRESSIDA'

      BLTCOD(149) =   710
      BLTNAM(149) =  'DESDEMONA'

      BLTCOD(150) =   711
      BLTNAM(150) =  'JULIET'

      BLTCOD(151) =   712
      BLTNAM(151) =  'PORTIA'

      BLTCOD(152) =   713
      BLTNAM(152) =  'ROSALIND'

      BLTCOD(153) =   714
      BLTNAM(153) =  'BELINDA'

      BLTCOD(154) =   715
      BLTNAM(154) =  'PUCK'

      BLTCOD(155) =   716
      BLTNAM(155) =  'CALIBAN'

      BLTCOD(156) =   717
      BLTNAM(156) =  'SYCORAX'

      BLTCOD(157) =   718
      BLTNAM(157) =  'PROSPERO'

      BLTCOD(158) =   719
      BLTNAM(158) =  'SETEBOS'

      BLTCOD(159) =   720
      BLTNAM(159) =  'STEPHANO'

      BLTCOD(160) =   721
      BLTNAM(160) =  'TRINCULO'

      BLTCOD(161) =   722
      BLTNAM(161) =  'FRANCISCO'

      BLTCOD(162) =   723
      BLTNAM(162) =  'MARGARET'

      BLTCOD(163) =   724
      BLTNAM(163) =  'FERDINAND'

      BLTCOD(164) =   725
      BLTNAM(164) =  'PERDITA'

      BLTCOD(165) =   726
      BLTNAM(165) =  'MAB'

      BLTCOD(166) =   727
      BLTNAM(166) =  'CUPID'

      BLTCOD(167) =   899
      BLTNAM(167) =  'NEPTUNE'

      BLTCOD(168) =   801
      BLTNAM(168) =  'TRITON'

      BLTCOD(169) =   802
      BLTNAM(169) =  'NEREID'

      BLTCOD(170) =   803
      BLTNAM(170) =  'NAIAD'

      BLTCOD(171) =   804
      BLTNAM(171) =  'THALASSA'

      BLTCOD(172) =   805
      BLTNAM(172) =  'DESPINA'

      BLTCOD(173) =   806
      BLTNAM(173) =  'GALATEA'

      BLTCOD(174) =   807
      BLTNAM(174) =  'LARISSA'

      BLTCOD(175) =   808
      BLTNAM(175) =  'PROTEUS'

      BLTCOD(176) =   809
      BLTNAM(176) =  'HALIMEDE'

      BLTCOD(177) =   810
      BLTNAM(177) =  'PSAMATHE'

      BLTCOD(178) =   811
      BLTNAM(178) =  'SAO'

      BLTCOD(179) =   812
      BLTNAM(179) =  'LAOMEDEIA'

      BLTCOD(180) =   813
      BLTNAM(180) =  'NESO'

      BLTCOD(181) =   999
      BLTNAM(181) =  'PLUTO'

      BLTCOD(182) =   901
      BLTNAM(182) =  'CHARON'

      BLTCOD(183) =   902
      BLTNAM(183) =  'NIX'

      BLTCOD(184) =   903
      BLTNAM(184) =  'HYDRA'

      BLTCOD(185) =   904
      BLTNAM(185) =  'KERBEROS'

      BLTCOD(186) =   905
      BLTNAM(186) =  'STYX'

      BLTCOD(187) =   -1
      BLTNAM(187) =  'GEOTAIL'

      BLTCOD(188) =   -3
      BLTNAM(188) =  'MOM'

      BLTCOD(189) =   -3
      BLTNAM(189) =  'MARS ORBITER MISSION'

      BLTCOD(190) =   -5
      BLTNAM(190) =  'AKATSUKI'

      BLTCOD(191) =   -5
      BLTNAM(191) =  'VCO'

      BLTCOD(192) =   -5
      BLTNAM(192) =  'PLC'

      BLTCOD(193) =   -5
      BLTNAM(193) =  'PLANET-C'

      BLTCOD(194) =   -6
      BLTNAM(194) =  'P6'

      BLTCOD(195) =   -6
      BLTNAM(195) =  'PIONEER-6'

      BLTCOD(196) =   -7
      BLTNAM(196) =  'P7'

      BLTCOD(197) =   -7
      BLTNAM(197) =  'PIONEER-7'

      BLTCOD(198) =   -8
      BLTNAM(198) =  'WIND'

      BLTCOD(199) =   -12
      BLTNAM(199) =  'VENUS ORBITER'

      BLTCOD(200) =   -12
      BLTNAM(200) =  'P12'

      BLTCOD(201) =   -12
      BLTNAM(201) =  'PIONEER 12'

      BLTCOD(202) =   -12
      BLTNAM(202) =  'LADEE'

      BLTCOD(203) =   -13
      BLTNAM(203) =  'POLAR'

      BLTCOD(204) =   -18
      BLTNAM(204) =  'MGN'

      BLTCOD(205) =   -18
      BLTNAM(205) =  'MAGELLAN'

      BLTCOD(206) =   -18
      BLTNAM(206) =  'LCROSS'

      BLTCOD(207) =   -20
      BLTNAM(207) =  'P8'

      BLTCOD(208) =   -20
      BLTNAM(208) =  'PIONEER-8'

      BLTCOD(209) =   -21
      BLTNAM(209) =  'SOHO'

      BLTCOD(210) =   -23
      BLTNAM(210) =  'P10'

      BLTCOD(211) =   -23
      BLTNAM(211) =  'PIONEER-10'

      BLTCOD(212) =   -24
      BLTNAM(212) =  'P11'

      BLTCOD(213) =   -24
      BLTNAM(213) =  'PIONEER-11'

      BLTCOD(214) =   -25
      BLTNAM(214) =  'LP'

      BLTCOD(215) =   -25
      BLTNAM(215) =  'LUNAR PROSPECTOR'

      BLTCOD(216) =   -27
      BLTNAM(216) =  'VK1'

      BLTCOD(217) =   -27
      BLTNAM(217) =  'VIKING 1 ORBITER'

      BLTCOD(218) =   -28
      BLTNAM(218) =  'JUPITER ICY MOONS EXPLORER'

      BLTCOD(219) =   -28
      BLTNAM(219) =  'JUICE'

      BLTCOD(220) =   -29
      BLTNAM(220) =  'STARDUST'

      BLTCOD(221) =   -29
      BLTNAM(221) =  'SDU'

      BLTCOD(222) =   -29
      BLTNAM(222) =  'NEXT'

      BLTCOD(223) =   -30
      BLTNAM(223) =  'VK2'

      BLTCOD(224) =   -30
      BLTNAM(224) =  'VIKING 2 ORBITER'

      BLTCOD(225) =   -30
      BLTNAM(225) =  'DS-1'

      BLTCOD(226) =   -31
      BLTNAM(226) =  'VG1'

      BLTCOD(227) =   -31
      BLTNAM(227) =  'VOYAGER 1'

      BLTCOD(228) =   -32
      BLTNAM(228) =  'VG2'

      BLTCOD(229) =   -32
      BLTNAM(229) =  'VOYAGER 2'

      BLTCOD(230) =   -40
      BLTNAM(230) =  'CLEMENTINE'

      BLTCOD(231) =   -41
      BLTNAM(231) =  'MEX'

      BLTCOD(232) =   -41
      BLTNAM(232) =  'MARS EXPRESS'

      BLTCOD(233) =   -44
      BLTNAM(233) =  'BEAGLE2'

      BLTCOD(234) =   -44
      BLTNAM(234) =  'BEAGLE 2'

      BLTCOD(235) =   -46
      BLTNAM(235) =  'MS-T5'

      BLTCOD(236) =   -46
      BLTNAM(236) =  'SAKIGAKE'

      BLTCOD(237) =   -47
      BLTNAM(237) =  'PLANET-A'

      BLTCOD(238) =   -47
      BLTNAM(238) =  'SUISEI'

      BLTCOD(239) =   -47
      BLTNAM(239) =  'GNS'

      BLTCOD(240) =   -47
      BLTNAM(240) =  'GENESIS'

      BLTCOD(241) =   -48
      BLTNAM(241) =  'HUBBLE SPACE TELESCOPE'

      BLTCOD(242) =   -48
      BLTNAM(242) =  'HST'

      BLTCOD(243) =   -49
      BLTNAM(243) =  'LUCY'

      BLTCOD(244) =   -53
      BLTNAM(244) =  'MARS PATHFINDER'

      BLTCOD(245) =   -53
      BLTNAM(245) =  'MPF'

      BLTCOD(246) =   -53
      BLTNAM(246) =  'MARS ODYSSEY'

      BLTCOD(247) =   -53
      BLTNAM(247) =  'MARS SURVEYOR 01 ORBITER'

      BLTCOD(248) =   -54
      BLTNAM(248) =  'ARM'

      BLTCOD(249) =   -54
      BLTNAM(249) =  'ASTEROID RETRIEVAL MISSION'

      BLTCOD(250) =   -55
      BLTNAM(250) =  'ULYSSES'

      BLTCOD(251) =   -58
      BLTNAM(251) =  'VSOP'

      BLTCOD(252) =   -58
      BLTNAM(252) =  'HALCA'

      BLTCOD(253) =   -59
      BLTNAM(253) =  'RADIOASTRON'

      BLTCOD(254) =   -61
      BLTNAM(254) =  'JUNO'

      BLTCOD(255) =   -62
      BLTNAM(255) =  'EMM'

      BLTCOD(256) =   -62
      BLTNAM(256) =  'EMIRATES MARS MISSION'

      BLTCOD(257) =   -64
      BLTNAM(257) =  'ORX'

      BLTCOD(258) =   -64
      BLTNAM(258) =  'OSIRIS-REX'

      BLTCOD(259) =   -65
      BLTNAM(259) =  'MCOA'

      BLTCOD(260) =   -65
      BLTNAM(260) =  'MARCO-A'

      BLTCOD(261) =   -66
      BLTNAM(261) =  'VEGA 1'

      BLTCOD(262) =   -66
      BLTNAM(262) =  'MCOB'

      BLTCOD(263) =   -66
      BLTNAM(263) =  'MARCO-B'

      BLTCOD(264) =   -67
      BLTNAM(264) =  'VEGA 2'

      BLTCOD(265) =   -68
      BLTNAM(265) =  'MERCURY MAGNETOSPHERIC ORBITER'

      BLTCOD(266) =   -68
      BLTNAM(266) =  'MMO'

      BLTCOD(267) =   -68
      BLTNAM(267) =  'BEPICOLOMBO MMO'

      BLTCOD(268) =   -69
      BLTNAM(268) =  'PSYC'

      BLTCOD(269) =   -70
      BLTNAM(269) =  'DEEP IMPACT IMPACTOR SPACECRAFT'

      BLTCOD(270) =   -74
      BLTNAM(270) =  'MRO'

      BLTCOD(271) =   -74
      BLTNAM(271) =  'MARS RECON ORBITER'

      BLTCOD(272) =   -76
      BLTNAM(272) =  'CURIOSITY'

      BLTCOD(273) =   -76
      BLTNAM(273) =  'MSL'

      BLTCOD(274) =   -76
      BLTNAM(274) =  'MARS SCIENCE LABORATORY'

      BLTCOD(275) =   -77
      BLTNAM(275) =  'GLL'

      BLTCOD(276) =   -77
      BLTNAM(276) =  'GALILEO ORBITER'

      BLTCOD(277) =   -78
      BLTNAM(277) =  'GIOTTO'

      BLTCOD(278) =   -79
      BLTNAM(278) =  'SPITZER'

      BLTCOD(279) =   -79
      BLTNAM(279) =  'SPACE INFRARED TELESCOPE FACILITY'

      BLTCOD(280) =   -79
      BLTNAM(280) =  'SIRTF'

      BLTCOD(281) =   -81
      BLTNAM(281) =  'CASSINI ITL'

      BLTCOD(282) =   -82
      BLTNAM(282) =  'CAS'

      BLTCOD(283) =   -82
      BLTNAM(283) =  'CASSINI'

      BLTCOD(284) =   -84
      BLTNAM(284) =  'PHOENIX'

      BLTCOD(285) =   -85
      BLTNAM(285) =  'LRO'

      BLTCOD(286) =   -85
      BLTNAM(286) =  'LUNAR RECON ORBITER'

      BLTCOD(287) =   -85
      BLTNAM(287) =  'LUNAR RECONNAISSANCE ORBITER'

      BLTCOD(288) =   -86
      BLTNAM(288) =  'CH1'

      BLTCOD(289) =   -86
      BLTNAM(289) =  'CHANDRAYAAN-1'

      BLTCOD(290) =   -90
      BLTNAM(290) =  'CASSINI SIMULATION'

      BLTCOD(291) =   -93
      BLTNAM(291) =  'NEAR EARTH ASTEROID RENDEZVOUS'

      BLTCOD(292) =   -93
      BLTNAM(292) =  'NEAR'

      BLTCOD(293) =   -94
      BLTNAM(293) =  'MO'

      BLTCOD(294) =   -94
      BLTNAM(294) =  'MARS OBSERVER'

      BLTCOD(295) =   -94
      BLTNAM(295) =  'MGS'

      BLTCOD(296) =   -94
      BLTNAM(296) =  'MARS GLOBAL SURVEYOR'

      BLTCOD(297) =   -95
      BLTNAM(297) =  'MGS SIMULATION'

      BLTCOD(298) =   -96
      BLTNAM(298) =  'SPP'

      BLTCOD(299) =   -96
      BLTNAM(299) =  'SOLAR PROBE PLUS'

      BLTCOD(300) =   -97
      BLTNAM(300) =  'TOPEX/POSEIDON'

      BLTCOD(301) =   -98
      BLTNAM(301) =  'NEW HORIZONS'

      BLTCOD(302) =   -107
      BLTNAM(302) =  'TROPICAL RAINFALL MEASURING MISSION'

      BLTCOD(303) =   -107
      BLTNAM(303) =  'TRMM'

      BLTCOD(304) =   -112
      BLTNAM(304) =  'ICE'

      BLTCOD(305) =   -116
      BLTNAM(305) =  'MARS POLAR LANDER'

      BLTCOD(306) =   -116
      BLTNAM(306) =  'MPL'

      BLTCOD(307) =   -117
      BLTNAM(307) =  'EDL DEMONSTRATOR MODULE'

      BLTCOD(308) =   -117
      BLTNAM(308) =  'EDM'

      BLTCOD(309) =   -117
      BLTNAM(309) =  'EXOMARS 2016 EDM'

      BLTCOD(310) =   -121
      BLTNAM(310) =  'MERCURY PLANETARY ORBITER'

      BLTCOD(311) =   -121
      BLTNAM(311) =  'MPO'

      BLTCOD(312) =   -121
      BLTNAM(312) =  'BEPICOLOMBO MPO'

      BLTCOD(313) =   -127
      BLTNAM(313) =  'MARS CLIMATE ORBITER'

      BLTCOD(314) =   -127
      BLTNAM(314) =  'MCO'

      BLTCOD(315) =   -130
      BLTNAM(315) =  'MUSES-C'

      BLTCOD(316) =   -130
      BLTNAM(316) =  'HAYABUSA'

      BLTCOD(317) =   -131
      BLTNAM(317) =  'SELENE'

      BLTCOD(318) =   -131
      BLTNAM(318) =  'KAGUYA'

      BLTCOD(319) =   -135
      BLTNAM(319) =  'DRTS-W'

      BLTCOD(320) =   -140
      BLTNAM(320) =  'EPOCH'

      BLTCOD(321) =   -140
      BLTNAM(321) =  'DIXI'

      BLTCOD(322) =   -140
      BLTNAM(322) =  'EPOXI'

      BLTCOD(323) =   -140
      BLTNAM(323) =  'DEEP IMPACT FLYBY SPACECRAFT'

      BLTCOD(324) =   -142
      BLTNAM(324) =  'TERRA'

      BLTCOD(325) =   -142
      BLTNAM(325) =  'EOS-AM1'

      BLTCOD(326) =   -143
      BLTNAM(326) =  'TRACE GAS ORBITER'

      BLTCOD(327) =   -143
      BLTNAM(327) =  'TGO'

      BLTCOD(328) =   -143
      BLTNAM(328) =  'EXOMARS 2016 TGO'

      BLTCOD(329) =   -144
      BLTNAM(329) =  'SOLO'

      BLTCOD(330) =   -144
      BLTNAM(330) =  'SOLAR ORBITER'

      BLTCOD(331) =   -146
      BLTNAM(331) =  'LUNAR-A'

      BLTCOD(332) =   -150
      BLTNAM(332) =  'CASSINI PROBE'

      BLTCOD(333) =   -150
      BLTNAM(333) =  'HUYGENS PROBE'

      BLTCOD(334) =   -150
      BLTNAM(334) =  'CASP'

      BLTCOD(335) =   -151
      BLTNAM(335) =  'AXAF'

      BLTCOD(336) =   -151
      BLTNAM(336) =  'CHANDRA'

      BLTCOD(337) =   -152
      BLTNAM(337) =  'CH2'

      BLTCOD(338) =   -152
      BLTNAM(338) =  'CHANDRAYAAN-2'

      BLTCOD(339) =   -154
      BLTNAM(339) =  'AQUA'

      BLTCOD(340) =   -159
      BLTNAM(340) =  'EURC'

      BLTCOD(341) =   -159
      BLTNAM(341) =  'EUROPA CLIPPER'

      BLTCOD(342) =   -164
      BLTNAM(342) =  'YOHKOH'

      BLTCOD(343) =   -164
      BLTNAM(343) =  'SOLAR-A'

      BLTCOD(344) =   -165
      BLTNAM(344) =  'MAP'

      BLTCOD(345) =   -166
      BLTNAM(345) =  'IMAGE'

      BLTCOD(346) =   -170
      BLTNAM(346) =  'JWST'

      BLTCOD(347) =   -170
      BLTNAM(347) =  'JAMES WEBB SPACE TELESCOPE'

      BLTCOD(348) =   -177
      BLTNAM(348) =  'GRAIL-A'

      BLTCOD(349) =   -178
      BLTNAM(349) =  'PLANET-B'

      BLTCOD(350) =   -178
      BLTNAM(350) =  'NOZOMI'

      BLTCOD(351) =   -181
      BLTNAM(351) =  'GRAIL-B'

      BLTCOD(352) =   -183
      BLTNAM(352) =  'CLUSTER 1'

      BLTCOD(353) =   -185
      BLTNAM(353) =  'CLUSTER 2'

      BLTCOD(354) =   -188
      BLTNAM(354) =  'MUSES-B'

      BLTCOD(355) =   -189
      BLTNAM(355) =  'NSYT'

      BLTCOD(356) =   -189
      BLTNAM(356) =  'INSIGHT'

      BLTCOD(357) =   -190
      BLTNAM(357) =  'SIM'

      BLTCOD(358) =   -194
      BLTNAM(358) =  'CLUSTER 3'

      BLTCOD(359) =   -196
      BLTNAM(359) =  'CLUSTER 4'

      BLTCOD(360) =   -198
      BLTNAM(360) =  'INTEGRAL'

      BLTCOD(361) =   -198
      BLTNAM(361) =  'NASA-ISRO SAR MISSION'

      BLTCOD(362) =   -198
      BLTNAM(362) =  'NISAR'

      BLTCOD(363) =   -200
      BLTNAM(363) =  'CONTOUR'

      BLTCOD(364) =   -202
      BLTNAM(364) =  'MAVEN'

      BLTCOD(365) =   -203
      BLTNAM(365) =  'DAWN'

      BLTCOD(366) =   -205
      BLTNAM(366) =  'SOIL MOISTURE ACTIVE AND PASSIVE'

      BLTCOD(367) =   -205
      BLTNAM(367) =  'SMAP'

      BLTCOD(368) =   -212
      BLTNAM(368) =  'STV51'

      BLTCOD(369) =   -213
      BLTNAM(369) =  'STV52'

      BLTCOD(370) =   -214
      BLTNAM(370) =  'STV53'

      BLTCOD(371) =   -226
      BLTNAM(371) =  'ROSETTA'

      BLTCOD(372) =   -227
      BLTNAM(372) =  'KEPLER'

      BLTCOD(373) =   -228
      BLTNAM(373) =  'GLL PROBE'

      BLTCOD(374) =   -228
      BLTNAM(374) =  'GALILEO PROBE'

      BLTCOD(375) =   -234
      BLTNAM(375) =  'STEREO AHEAD'

      BLTCOD(376) =   -235
      BLTNAM(376) =  'STEREO BEHIND'

      BLTCOD(377) =   -236
      BLTNAM(377) =  'MESSENGER'

      BLTCOD(378) =   -238
      BLTNAM(378) =  'SMART1'

      BLTCOD(379) =   -238
      BLTNAM(379) =  'SM1'

      BLTCOD(380) =   -238
      BLTNAM(380) =  'S1'

      BLTCOD(381) =   -238
      BLTNAM(381) =  'SMART-1'

      BLTCOD(382) =   -248
      BLTNAM(382) =  'VEX'

      BLTCOD(383) =   -248
      BLTNAM(383) =  'VENUS EXPRESS'

      BLTCOD(384) =   -253
      BLTNAM(384) =  'OPPORTUNITY'

      BLTCOD(385) =   -253
      BLTNAM(385) =  'MER-1'

      BLTCOD(386) =   -254
      BLTNAM(386) =  'SPIRIT'

      BLTCOD(387) =   -254
      BLTNAM(387) =  'MER-2'

      BLTCOD(388) =   -301
      BLTNAM(388) =  'HELIOS 1'

      BLTCOD(389) =   -302
      BLTNAM(389) =  'HELIOS 2'

      BLTCOD(390) =   -362
      BLTNAM(390) =  'RADIATION BELT STORM PROBE A'

      BLTCOD(391) =   -362
      BLTNAM(391) =  'RBSP_A'

      BLTCOD(392) =   -363
      BLTNAM(392) =  'RADIATION BELT STORM PROBE B'

      BLTCOD(393) =   -363
      BLTNAM(393) =  'RBSP_B'

      BLTCOD(394) =   -500
      BLTNAM(394) =  'RSAT'

      BLTCOD(395) =   -500
      BLTNAM(395) =  'SELENE Relay Satellite'

      BLTCOD(396) =   -500
      BLTNAM(396) =  'SELENE Rstar'

      BLTCOD(397) =   -500
      BLTNAM(397) =  'Rstar'

      BLTCOD(398) =   -502
      BLTNAM(398) =  'VSAT'

      BLTCOD(399) =   -502
      BLTNAM(399) =  'SELENE VLBI Radio Satellite'

      BLTCOD(400) =   -502
      BLTNAM(400) =  'SELENE VRAD Satellite'

      BLTCOD(401) =   -502
      BLTNAM(401) =  'SELENE Vstar'

      BLTCOD(402) =   -502
      BLTNAM(402) =  'Vstar'

      BLTCOD(403) =   -550
      BLTNAM(403) =  'MARS-96'

      BLTCOD(404) =   -550
      BLTNAM(404) =  'M96'

      BLTCOD(405) =   -550
      BLTNAM(405) =  'MARS 96'

      BLTCOD(406) =   -550
      BLTNAM(406) =  'MARS96'

      BLTCOD(407) =   -750
      BLTNAM(407) =  'SPRINT-A'

      BLTCOD(408) =   50000001
      BLTNAM(408) =  'SHOEMAKER-LEVY 9-W'

      BLTCOD(409) =   50000002
      BLTNAM(409) =  'SHOEMAKER-LEVY 9-V'

      BLTCOD(410) =   50000003
      BLTNAM(410) =  'SHOEMAKER-LEVY 9-U'

      BLTCOD(411) =   50000004
      BLTNAM(411) =  'SHOEMAKER-LEVY 9-T'

      BLTCOD(412) =   50000005
      BLTNAM(412) =  'SHOEMAKER-LEVY 9-S'

      BLTCOD(413) =   50000006
      BLTNAM(413) =  'SHOEMAKER-LEVY 9-R'

      BLTCOD(414) =   50000007
      BLTNAM(414) =  'SHOEMAKER-LEVY 9-Q'

      BLTCOD(415) =   50000008
      BLTNAM(415) =  'SHOEMAKER-LEVY 9-P'

      BLTCOD(416) =   50000009
      BLTNAM(416) =  'SHOEMAKER-LEVY 9-N'

      BLTCOD(417) =   50000010
      BLTNAM(417) =  'SHOEMAKER-LEVY 9-M'

      BLTCOD(418) =   50000011
      BLTNAM(418) =  'SHOEMAKER-LEVY 9-L'

      BLTCOD(419) =   50000012
      BLTNAM(419) =  'SHOEMAKER-LEVY 9-K'

      BLTCOD(420) =   50000013
      BLTNAM(420) =  'SHOEMAKER-LEVY 9-J'

      BLTCOD(421) =   50000014
      BLTNAM(421) =  'SHOEMAKER-LEVY 9-H'

      BLTCOD(422) =   50000015
      BLTNAM(422) =  'SHOEMAKER-LEVY 9-G'

      BLTCOD(423) =   50000016
      BLTNAM(423) =  'SHOEMAKER-LEVY 9-F'

      BLTCOD(424) =   50000017
      BLTNAM(424) =  'SHOEMAKER-LEVY 9-E'

      BLTCOD(425) =   50000018
      BLTNAM(425) =  'SHOEMAKER-LEVY 9-D'

      BLTCOD(426) =   50000019
      BLTNAM(426) =  'SHOEMAKER-LEVY 9-C'

      BLTCOD(427) =   50000020
      BLTNAM(427) =  'SHOEMAKER-LEVY 9-B'

      BLTCOD(428) =   50000021
      BLTNAM(428) =  'SHOEMAKER-LEVY 9-A'

      BLTCOD(429) =   50000022
      BLTNAM(429) =  'SHOEMAKER-LEVY 9-Q1'

      BLTCOD(430) =   50000023
      BLTNAM(430) =  'SHOEMAKER-LEVY 9-P2'

      BLTCOD(431) =   1000001
      BLTNAM(431) =  'AREND'

      BLTCOD(432) =   1000002
      BLTNAM(432) =  'AREND-RIGAUX'

      BLTCOD(433) =   1000003
      BLTNAM(433) =  'ASHBROOK-JACKSON'

      BLTCOD(434) =   1000004
      BLTNAM(434) =  'BOETHIN'

      BLTCOD(435) =   1000005
      BLTNAM(435) =  'BORRELLY'

      BLTCOD(436) =   1000006
      BLTNAM(436) =  'BOWELL-SKIFF'

      BLTCOD(437) =   1000007
      BLTNAM(437) =  'BRADFIELD'

      BLTCOD(438) =   1000008
      BLTNAM(438) =  'BROOKS 2'

      BLTCOD(439) =   1000009
      BLTNAM(439) =  'BRORSEN-METCALF'

      BLTCOD(440) =   1000010
      BLTNAM(440) =  'BUS'

      BLTCOD(441) =   1000011
      BLTNAM(441) =  'CHERNYKH'

      BLTCOD(442) =   1000012
      BLTNAM(442) =  '67P/CHURYUMOV-GERASIMENKO (1969 R1)'

      BLTCOD(443) =   1000012
      BLTNAM(443) =  'CHURYUMOV-GERASIMENKO'

      BLTCOD(444) =   1000013
      BLTNAM(444) =  'CIFFREO'

      BLTCOD(445) =   1000014
      BLTNAM(445) =  'CLARK'

      BLTCOD(446) =   1000015
      BLTNAM(446) =  'COMAS SOLA'

      BLTCOD(447) =   1000016
      BLTNAM(447) =  'CROMMELIN'

      BLTCOD(448) =   1000017
      BLTNAM(448) =  'D''ARREST'

      BLTCOD(449) =   1000018
      BLTNAM(449) =  'DANIEL'

      BLTCOD(450) =   1000019
      BLTNAM(450) =  'DE VICO-SWIFT'

      BLTCOD(451) =   1000020
      BLTNAM(451) =  'DENNING-FUJIKAWA'

      BLTCOD(452) =   1000021
      BLTNAM(452) =  'DU TOIT 1'

      BLTCOD(453) =   1000022
      BLTNAM(453) =  'DU TOIT-HARTLEY'

      BLTCOD(454) =   1000023
      BLTNAM(454) =  'DUTOIT-NEUJMIN-DELPORTE'

      BLTCOD(455) =   1000024
      BLTNAM(455) =  'DUBIAGO'

      BLTCOD(456) =   1000025
      BLTNAM(456) =  'ENCKE'

      BLTCOD(457) =   1000026
      BLTNAM(457) =  'FAYE'

      BLTCOD(458) =   1000027
      BLTNAM(458) =  'FINLAY'

      BLTCOD(459) =   1000028
      BLTNAM(459) =  'FORBES'

      BLTCOD(460) =   1000029
      BLTNAM(460) =  'GEHRELS 1'

      BLTCOD(461) =   1000030
      BLTNAM(461) =  'GEHRELS 2'

      BLTCOD(462) =   1000031
      BLTNAM(462) =  'GEHRELS 3'

      BLTCOD(463) =   1000032
      BLTNAM(463) =  'GIACOBINI-ZINNER'

      BLTCOD(464) =   1000033
      BLTNAM(464) =  'GICLAS'

      BLTCOD(465) =   1000034
      BLTNAM(465) =  'GRIGG-SKJELLERUP'

      BLTCOD(466) =   1000035
      BLTNAM(466) =  'GUNN'

      BLTCOD(467) =   1000036
      BLTNAM(467) =  'HALLEY'

      BLTCOD(468) =   1000037
      BLTNAM(468) =  'HANEDA-CAMPOS'

      BLTCOD(469) =   1000038
      BLTNAM(469) =  'HARRINGTON'

      BLTCOD(470) =   1000039
      BLTNAM(470) =  'HARRINGTON-ABELL'

      BLTCOD(471) =   1000040
      BLTNAM(471) =  'HARTLEY 1'

      BLTCOD(472) =   1000041
      BLTNAM(472) =  'HARTLEY 2'

      BLTCOD(473) =   1000042
      BLTNAM(473) =  'HARTLEY-IRAS'

      BLTCOD(474) =   1000043
      BLTNAM(474) =  'HERSCHEL-RIGOLLET'

      BLTCOD(475) =   1000044
      BLTNAM(475) =  'HOLMES'

      BLTCOD(476) =   1000045
      BLTNAM(476) =  'HONDA-MRKOS-PAJDUSAKOVA'

      BLTCOD(477) =   1000046
      BLTNAM(477) =  'HOWELL'

      BLTCOD(478) =   1000047
      BLTNAM(478) =  'IRAS'

      BLTCOD(479) =   1000048
      BLTNAM(479) =  'JACKSON-NEUJMIN'

      BLTCOD(480) =   1000049
      BLTNAM(480) =  'JOHNSON'

      BLTCOD(481) =   1000050
      BLTNAM(481) =  'KEARNS-KWEE'

      BLTCOD(482) =   1000051
      BLTNAM(482) =  'KLEMOLA'

      BLTCOD(483) =   1000052
      BLTNAM(483) =  'KOHOUTEK'

      BLTCOD(484) =   1000053
      BLTNAM(484) =  'KOJIMA'

      BLTCOD(485) =   1000054
      BLTNAM(485) =  'KOPFF'

      BLTCOD(486) =   1000055
      BLTNAM(486) =  'KOWAL 1'

      BLTCOD(487) =   1000056
      BLTNAM(487) =  'KOWAL 2'

      BLTCOD(488) =   1000057
      BLTNAM(488) =  'KOWAL-MRKOS'

      BLTCOD(489) =   1000058
      BLTNAM(489) =  'KOWAL-VAVROVA'

      BLTCOD(490) =   1000059
      BLTNAM(490) =  'LONGMORE'

      BLTCOD(491) =   1000060
      BLTNAM(491) =  'LOVAS 1'

      BLTCOD(492) =   1000061
      BLTNAM(492) =  'MACHHOLZ'

      BLTCOD(493) =   1000062
      BLTNAM(493) =  'MAURY'

      BLTCOD(494) =   1000063
      BLTNAM(494) =  'NEUJMIN 1'

      BLTCOD(495) =   1000064
      BLTNAM(495) =  'NEUJMIN 2'

      BLTCOD(496) =   1000065
      BLTNAM(496) =  'NEUJMIN 3'

      BLTCOD(497) =   1000066
      BLTNAM(497) =  'OLBERS'

      BLTCOD(498) =   1000067
      BLTNAM(498) =  'PETERS-HARTLEY'

      BLTCOD(499) =   1000068
      BLTNAM(499) =  'PONS-BROOKS'

      BLTCOD(500) =   1000069
      BLTNAM(500) =  'PONS-WINNECKE'

      BLTCOD(501) =   1000070
      BLTNAM(501) =  'REINMUTH 1'

      BLTCOD(502) =   1000071
      BLTNAM(502) =  'REINMUTH 2'

      BLTCOD(503) =   1000072
      BLTNAM(503) =  'RUSSELL 1'

      BLTCOD(504) =   1000073
      BLTNAM(504) =  'RUSSELL 2'

      BLTCOD(505) =   1000074
      BLTNAM(505) =  'RUSSELL 3'

      BLTCOD(506) =   1000075
      BLTNAM(506) =  'RUSSELL 4'

      BLTCOD(507) =   1000076
      BLTNAM(507) =  'SANGUIN'

      BLTCOD(508) =   1000077
      BLTNAM(508) =  'SCHAUMASSE'

      BLTCOD(509) =   1000078
      BLTNAM(509) =  'SCHUSTER'

      BLTCOD(510) =   1000079
      BLTNAM(510) =  'SCHWASSMANN-WACHMANN 1'

      BLTCOD(511) =   1000080
      BLTNAM(511) =  'SCHWASSMANN-WACHMANN 2'

      BLTCOD(512) =   1000081
      BLTNAM(512) =  'SCHWASSMANN-WACHMANN 3'

      BLTCOD(513) =   1000082
      BLTNAM(513) =  'SHAJN-SCHALDACH'

      BLTCOD(514) =   1000083
      BLTNAM(514) =  'SHOEMAKER 1'

      BLTCOD(515) =   1000084
      BLTNAM(515) =  'SHOEMAKER 2'

      BLTCOD(516) =   1000085
      BLTNAM(516) =  'SHOEMAKER 3'

      BLTCOD(517) =   1000086
      BLTNAM(517) =  'SINGER-BREWSTER'

      BLTCOD(518) =   1000087
      BLTNAM(518) =  'SLAUGHTER-BURNHAM'

      BLTCOD(519) =   1000088
      BLTNAM(519) =  'SMIRNOVA-CHERNYKH'

      BLTCOD(520) =   1000089
      BLTNAM(520) =  'STEPHAN-OTERMA'

      BLTCOD(521) =   1000090
      BLTNAM(521) =  'SWIFT-GEHRELS'

      BLTCOD(522) =   1000091
      BLTNAM(522) =  'TAKAMIZAWA'

      BLTCOD(523) =   1000092
      BLTNAM(523) =  'TAYLOR'

      BLTCOD(524) =   1000093
      BLTNAM(524) =  'TEMPEL_1'

      BLTCOD(525) =   1000093
      BLTNAM(525) =  'TEMPEL 1'

      BLTCOD(526) =   1000094
      BLTNAM(526) =  'TEMPEL 2'

      BLTCOD(527) =   1000095
      BLTNAM(527) =  'TEMPEL-TUTTLE'

      BLTCOD(528) =   1000096
      BLTNAM(528) =  'TRITTON'

      BLTCOD(529) =   1000097
      BLTNAM(529) =  'TSUCHINSHAN 1'

      BLTCOD(530) =   1000098
      BLTNAM(530) =  'TSUCHINSHAN 2'

      BLTCOD(531) =   1000099
      BLTNAM(531) =  'TUTTLE'

      BLTCOD(532) =   1000100
      BLTNAM(532) =  'TUTTLE-GIACOBINI-KRESAK'

      BLTCOD(533) =   1000101
      BLTNAM(533) =  'VAISALA 1'

      BLTCOD(534) =   1000102
      BLTNAM(534) =  'VAN BIESBROECK'

      BLTCOD(535) =   1000103
      BLTNAM(535) =  'VAN HOUTEN'

      BLTCOD(536) =   1000104
      BLTNAM(536) =  'WEST-KOHOUTEK-IKEMURA'

      BLTCOD(537) =   1000105
      BLTNAM(537) =  'WHIPPLE'

      BLTCOD(538) =   1000106
      BLTNAM(538) =  'WILD 1'

      BLTCOD(539) =   1000107
      BLTNAM(539) =  'WILD 2'

      BLTCOD(540) =   1000108
      BLTNAM(540) =  'WILD 3'

      BLTCOD(541) =   1000109
      BLTNAM(541) =  'WIRTANEN'

      BLTCOD(542) =   1000110
      BLTNAM(542) =  'WOLF'

      BLTCOD(543) =   1000111
      BLTNAM(543) =  'WOLF-HARRINGTON'

      BLTCOD(544) =   1000112
      BLTNAM(544) =  'LOVAS 2'

      BLTCOD(545) =   1000113
      BLTNAM(545) =  'URATA-NIIJIMA'

      BLTCOD(546) =   1000114
      BLTNAM(546) =  'WISEMAN-SKIFF'

      BLTCOD(547) =   1000115
      BLTNAM(547) =  'HELIN'

      BLTCOD(548) =   1000116
      BLTNAM(548) =  'MUELLER'

      BLTCOD(549) =   1000117
      BLTNAM(549) =  'SHOEMAKER-HOLT 1'

      BLTCOD(550) =   1000118
      BLTNAM(550) =  'HELIN-ROMAN-CROCKETT'

      BLTCOD(551) =   1000119
      BLTNAM(551) =  'HARTLEY 3'

      BLTCOD(552) =   1000120
      BLTNAM(552) =  'PARKER-HARTLEY'

      BLTCOD(553) =   1000121
      BLTNAM(553) =  'HELIN-ROMAN-ALU 1'

      BLTCOD(554) =   1000122
      BLTNAM(554) =  'WILD 4'

      BLTCOD(555) =   1000123
      BLTNAM(555) =  'MUELLER 2'

      BLTCOD(556) =   1000124
      BLTNAM(556) =  'MUELLER 3'

      BLTCOD(557) =   1000125
      BLTNAM(557) =  'SHOEMAKER-LEVY 1'

      BLTCOD(558) =   1000126
      BLTNAM(558) =  'SHOEMAKER-LEVY 2'

      BLTCOD(559) =   1000127
      BLTNAM(559) =  'HOLT-OLMSTEAD'

      BLTCOD(560) =   1000128
      BLTNAM(560) =  'METCALF-BREWINGTON'

      BLTCOD(561) =   1000129
      BLTNAM(561) =  'LEVY'

      BLTCOD(562) =   1000130
      BLTNAM(562) =  'SHOEMAKER-LEVY 9'

      BLTCOD(563) =   1000131
      BLTNAM(563) =  'HYAKUTAKE'

      BLTCOD(564) =   1000132
      BLTNAM(564) =  'HALE-BOPP'

      BLTCOD(565) =   1003228
      BLTNAM(565) =  'C/2013 A1'

      BLTCOD(566) =   1003228
      BLTNAM(566) =  'SIDING SPRING'

      BLTCOD(567) =   9511010
      BLTNAM(567) =  'GASPRA'

      BLTCOD(568) =   2431010
      BLTNAM(568) =  'IDA'

      BLTCOD(569) =   2431011
      BLTNAM(569) =  'DACTYL'

      BLTCOD(570) =   2000001
      BLTNAM(570) =  'CERES'

      BLTCOD(571) =   2000002
      BLTNAM(571) =  'PALLAS'

      BLTCOD(572) =   2000004
      BLTNAM(572) =  'VESTA'

      BLTCOD(573) =   2000016
      BLTNAM(573) =  'PSYCHE'

      BLTCOD(574) =   2000021
      BLTNAM(574) =  'LUTETIA'

      BLTCOD(575) =   2000216
      BLTNAM(575) =  'KLEOPATRA'

      BLTCOD(576) =   2000433
      BLTNAM(576) =  'EROS'

      BLTCOD(577) =   2000511
      BLTNAM(577) =  'DAVIDA'

      BLTCOD(578) =   2000253
      BLTNAM(578) =  'MATHILDE'

      BLTCOD(579) =   2002867
      BLTNAM(579) =  'STEINS'

      BLTCOD(580) =   2009969
      BLTNAM(580) =  '1992KD'

      BLTCOD(581) =   2009969
      BLTNAM(581) =  'BRAILLE'

      BLTCOD(582) =   2004015
      BLTNAM(582) =  'WILSON-HARRINGTON'

      BLTCOD(583) =   2004179
      BLTNAM(583) =  'TOUTATIS'

      BLTCOD(584) =   2025143
      BLTNAM(584) =  'ITOKAWA'

      BLTCOD(585) =   2101955
      BLTNAM(585) =  'BENNU'

      BLTCOD(586) =   398989
      BLTNAM(586) =  'NOTO'

      BLTCOD(587) =   398990
      BLTNAM(587) =  'NEW NORCIA'

      BLTCOD(588) =   399001
      BLTNAM(588) =  'GOLDSTONE'

      BLTCOD(589) =   399002
      BLTNAM(589) =  'CANBERRA'

      BLTCOD(590) =   399003
      BLTNAM(590) =  'MADRID'

      BLTCOD(591) =   399004
      BLTNAM(591) =  'USUDA'

      BLTCOD(592) =   399005
      BLTNAM(592) =  'DSS-05'

      BLTCOD(593) =   399005
      BLTNAM(593) =  'PARKES'

      BLTCOD(594) =   399012
      BLTNAM(594) =  'DSS-12'

      BLTCOD(595) =   399013
      BLTNAM(595) =  'DSS-13'

      BLTCOD(596) =   399014
      BLTNAM(596) =  'DSS-14'

      BLTCOD(597) =   399015
      BLTNAM(597) =  'DSS-15'

      BLTCOD(598) =   399016
      BLTNAM(598) =  'DSS-16'

      BLTCOD(599) =   399017
      BLTNAM(599) =  'DSS-17'

      BLTCOD(600) =   399023
      BLTNAM(600) =  'DSS-23'

      BLTCOD(601) =   399024
      BLTNAM(601) =  'DSS-24'

      BLTCOD(602) =   399025
      BLTNAM(602) =  'DSS-25'

      BLTCOD(603) =   399026
      BLTNAM(603) =  'DSS-26'

      BLTCOD(604) =   399027
      BLTNAM(604) =  'DSS-27'

      BLTCOD(605) =   399028
      BLTNAM(605) =  'DSS-28'

      BLTCOD(606) =   399033
      BLTNAM(606) =  'DSS-33'

      BLTCOD(607) =   399034
      BLTNAM(607) =  'DSS-34'

      BLTCOD(608) =   399042
      BLTNAM(608) =  'DSS-42'

      BLTCOD(609) =   399043
      BLTNAM(609) =  'DSS-43'

      BLTCOD(610) =   399045
      BLTNAM(610) =  'DSS-45'

      BLTCOD(611) =   399046
      BLTNAM(611) =  'DSS-46'

      BLTCOD(612) =   399049
      BLTNAM(612) =  'DSS-49'

      BLTCOD(613) =   399053
      BLTNAM(613) =  'DSS-53'

      BLTCOD(614) =   399054
      BLTNAM(614) =  'DSS-54'

      BLTCOD(615) =   399055
      BLTNAM(615) =  'DSS-55'

      BLTCOD(616) =   399061
      BLTNAM(616) =  'DSS-61'

      BLTCOD(617) =   399063
      BLTNAM(617) =  'DSS-63'

      BLTCOD(618) =   399064
      BLTNAM(618) =  'DSS-64'

      BLTCOD(619) =   399065
      BLTNAM(619) =  'DSS-65'

      BLTCOD(620) =   399066
      BLTNAM(620) =  'DSS-66'



      RETURN
      END

