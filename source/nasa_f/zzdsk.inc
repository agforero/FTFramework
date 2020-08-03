
C
C     File: zzdsk.inc
C
C
C     Version 4.0.0 13-NOV-2015 (NJB)
C
C        Changed parameter LBTLEN to CVTLEN.
C        Added parameter LMBCRV.
C        
C     Version 3.0.0 05-NOV-2015 (NJB)
C
C        Added parameters
C
C           CTRCOR
C           ELLCOR
C           GUIDED
C           LBTLEN
C           PNMBRL
C           TANGNT
C           TMTLEN
C           UMBRAL
C
C     Version 2.0.0 04-MAR-2015 (NJB)
C
C        Removed declaration of parameter SHPLEN.
C        This name is already in use in the include
C        file gf.inc.
C
C     Version 1.0.0 26-JAN-2015 (NJB)
C

C
C     Parameters supporting METHOD string parsing:
C
      INTEGER               DSKSHP
      PARAMETER           ( DSKSHP = 2 )

      INTEGER               ELLSHP
      PARAMETER           ( ELLSHP = 1 )

C
C     Local method length.
C
      INTEGER               MTHLEN
      PARAMETER           ( MTHLEN = 500 )

C
C     Length of sub-point type string.
C
      INTEGER               SUBLEN
      PARAMETER           ( SUBLEN = 20 )

C
C     Length of curve type string.
C
      INTEGER               CVTLEN
      PARAMETER           ( CVTLEN = 20 )

C
C     Limb type parameter codes.
C
      INTEGER               TANGNT
      PARAMETER           ( TANGNT = 1 )

      INTEGER               GUIDED
      PARAMETER           ( GUIDED = 2 )

C
C     Length of terminator type string.
C
      INTEGER               TMTLEN
      PARAMETER           ( TMTLEN = 20 )

C
C     Terminator type and limb parameter codes.
C
      INTEGER               LMBCRV
      PARAMETER           ( LMBCRV = 0 )

      INTEGER               UMBRAL
      PARAMETER           ( UMBRAL = 1 )

      INTEGER               PNMBRL
      PARAMETER           ( PNMBRL = 2 )

C
C     Length of aberration correction locus string.
C
      INTEGER               ACLLEN
      PARAMETER           ( ACLLEN = 25 )

C
C     Aberration correction locus codes.
C
      INTEGER               CTRCOR
      PARAMETER           ( CTRCOR = 1 )

      INTEGER               ELLCOR
      PARAMETER           ( ELLCOR = 2 )

C
C     End of include file zzdsk.inc
C
