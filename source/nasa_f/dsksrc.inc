
C
C     File: dsksrc.inc
C
C
C     Version 1.0.0 19-FEB-2016 (NJB)

C
C     This file declares sizes of certain output array arguments
C     returned by DSKXSI. These arrays contain information on the
C     source of surface data on which an intercept was found.
C
C     See the include files
C
C        dla.inc
C        dskdsc.inc
C
C     for the declarations of DLA and DSK descriptor sizes,
C     respectively.
C
C     Caution: the sizes declared here may be increased
C     in a future SPICE Toolkit version. These sizes
C     are correct for the N0066 Toolkit.
C     
C

C
C     Size of double precision source information component:
C
      INTEGER               DCSIZE
      PARAMETER           ( DCSIZE = 1 )

C     Size of integer source information component:
C
      INTEGER               ICSIZE
      PARAMETER           ( ICSIZE = 1 )

C
C     The contents of the components are DSK data type-dependent.
C
C         Type 2
C         ======
C
C                    Integer component        Double Precision component
C
C         Element    Contents
C         -------    ---------------------------------------------------
C            1       Plate ID                 Unused
C
C
C
C     End of include file dsksrc.inc
C
