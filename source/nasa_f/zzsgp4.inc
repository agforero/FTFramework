C$Procedure ZZSGP4 ( SGP4 parameters )
C
C$ Abstract
C
C      Parameter assignments for SGP4 algorithm as expressed
C      by Vallado [2].
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
C     None.
C
C$ Brief_I/O
C
C     None.
C
C$ Detailed_Input
C
C     None.
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
C     None.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     None.
C
C$ Examples
C
C     J2    = GEOPHS(K_J2)
C     J3    = GEOPHS(K_J3)
C     J4    = GEOPHS(K_J4)
C     ER    = GEOPHS(K_ER)
C     XKE   = GEOPHS(K_KE)
C
C     TUMIN = 1.D0/XKE
C     J3OJ2 = J3/J2
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
C   [2] Vallado, David, Crawford, Paul, Hujsak, Richard, and Kelso, T.S.
C       2006. Revisiting Spacetrack Report #3. Paper AIAA 2006-6753
C       presented at the AIAA/AAS Astrodynamics Specialist Conference,
C       August 21-24, 2006. Keystone, CO.
C
C$ Author_and_Institution
C
C     E. D. Wright    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0 22-JUL-2014 (EDW)
C
C-&

C$ Index_Entries
C
C  SGP4
C
C-&


C
C      WGS gravitational constants IDs.
C

      INTEGER                       NGRAVS
      PARAMETER                   ( NGRAVS = 3 )

      INTEGER                       NGRAVC
      PARAMETER                   ( NGRAVC = 8 )


      INTEGER                       WGS721
      PARAMETER                   ( WGS721 = 1 )

      INTEGER                       WGS72
      PARAMETER                   ( WGS72  = 2 )

      INTEGER                       WGS84
      PARAMETER                   ( WGS84  = 3 )


C
C      Gravitational constant indices.
C

      INTEGER                       P_RAD
      PARAMETER                   ( P_RAD =  1 )

      INTEGER                       P_XKE
      PARAMETER                   ( P_XKE =  2 )

      INTEGER                       P_MU
      PARAMETER                   ( P_MU =   3 )

      INTEGER                       P_TUMN
      PARAMETER                   ( P_TUMN = 4 )

      INTEGER                       P_J2
      PARAMETER                   ( P_J2 =   5 )

      INTEGER                       P_J3
      PARAMETER                   ( P_J3 =   6 )

      INTEGER                       P_J4
      PARAMETER                   ( P_J4 =   7 )

      INTEGER                       P_J3J2
      PARAMETER                   ( P_J3J2 = 8 )


C
C     The following parameters give the location in the GEOPHS
C     array of the various geophysical parameters needed for
C     the two line element sets.
C
C     K_J2  --- location of J2
C     K_J3  --- location of J3
C     K_J4  --- location if J4
C     K_KE  --- location of KE = sqrt(GM) in earth-radii**1.5/MIN
C     K_QO  --- upper bound of atmospheric model in KM
C     K_SO  --- lower bound of atmospheric model in KM
C     K_ER  --- earth equatorial radius in KM.
C     K_AE  --- distance units/earth radius
C

      INTEGER                       K_J2
      PARAMETER                   ( K_J2 = 1 )

      INTEGER                       K_J3
      PARAMETER                   ( K_J3 = 2 )

      INTEGER                       K_J4
      PARAMETER                   ( K_J4 = 3 )

      INTEGER                       K_KE
      PARAMETER                   ( K_KE = 4 )

      INTEGER                       K_QO
      PARAMETER                   ( K_QO = 5 )

      INTEGER                       K_SO
      PARAMETER                   ( K_SO = 6 )

      INTEGER                       K_ER
      PARAMETER                   ( K_ER = 7 )

      INTEGER                       K_AE
      PARAMETER                   ( K_AE = 8 )

      INTEGER                       NGEO
      PARAMETER                   ( NGEO = K_AE )

C
C     Operation mode values, OPMODE.
C

      INTEGER                       AFSPC
      PARAMETER                   ( AFSPC  = 1 )

      INTEGER                       IMPRVD
      PARAMETER                   ( IMPRVD = 2 )

C
C     An enumeration of the various components of the
C     elements array---ELEMS
C
C     KNDT20  --- location of NDT20
C     KNDD60  --- location of NDD60
C     KBSTAR  --- location of BSTAR
C     KINCL   --- location of INCL
C     KNODE0  --- location of NODE0
C     KECC    --- location of ECC
C     KOMEGA  --- location of OMEGA
C     KMO     --- location of MO
C     KNO     --- location of NO
C

      INTEGER                       KNDT20
      PARAMETER                   ( KNDT20 = 1 )

      INTEGER                       KNDD60
      PARAMETER                   ( KNDD60 = 2 )

      INTEGER                       KBSTAR
      PARAMETER                   ( KBSTAR = 3 )

      INTEGER                       KINCL
      PARAMETER                   ( KINCL  = 4 )

      INTEGER                       KNODE0
      PARAMETER                   ( KNODE0 = 5 )

      INTEGER                       KECC
      PARAMETER                   ( KECC   = 6 )

      INTEGER                       KOMEGA
      PARAMETER                   ( KOMEGA = 7 )

      INTEGER                       KMO
      PARAMETER                   ( KMO    = 8 )

      INTEGER                       KNO
      PARAMETER                   ( KNO    = 9 )

      INTEGER                       KEPOCH
      PARAMETER                   ( KEPOCH = 10 )

      INTEGER                       NELEMS
      PARAMETER                   ( NELEMS = KEPOCH )


