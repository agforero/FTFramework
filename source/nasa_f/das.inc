
C
C     Include file das.inc
C
C     This include file declares public parameters for the DAS 
C     subsystem.
C
C        Version 1.0.0 10-FEB-2017 (NJB)
C
C     Parameter declarations follow.
C

C
C     DAS file table size:
C
C        The parameter name is FTSIZE. The value of the parameter is
C        defined in the include file
C
C           zzddhman.inc
C
C        That value is duplicated here, since zzddhman.inc contains
C        other declarations that conflict with some of those in DAS
C        routines.
C
      INTEGER               FTSIZE
      PARAMETER           ( FTSIZE = 5000 )

C 
C     Capacity of DAS data records:
C
C        -- NWD double precision numbers.
C        -- NWI integers.
C        -- NWC characters.
 
C     These parameters are named to enhance ease of maintenance of
C     the code; the values should not be changed.
 
      INTEGER               NWD
      PARAMETER           ( NWD     =  128 )
 
      INTEGER               NWI
      PARAMETER           ( NWI     =  256 )
 
      INTEGER               NWC
      PARAMETER           ( NWC     = 1024 )

C
C     End of include file das.inc
C

