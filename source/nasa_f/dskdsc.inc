C
C     Include file dskdsc.inc
C
C     This include file declares parameters for DSK segment descriptors.
C
C-       SPICELIB Version 1.0.0 08-FEB-2017 (NJB)
C
C           Updated version info.
C
C           22-JAN-2016 (NJB)
C
C              Added parameter for data class 2. Changed name of data
C              class 1 parameter. Corrected data class descriptions.
C
C           13-MAY-2010 (NJB)
C
C              Descriptor now contains two ID codes, one for the
C              surface, one for the associated ephemeris object. This
C              supports association of multiple surfaces with one
C              ephemeris object without creating file management
C              issues.
C
C              Room was added for coordinate system definition
C              parameters.
C
C               Flag arrays and model ID/component entries were deleted.
C
C            11-SEP-2008 (NJB)
C
C
C     DSK segment descriptors are implemented as an array of d.p.
C     numbers.  Note that each integer descriptor datum occupies one
C     d.p. value. 
C
C     
C     
C
C     Segment descriptor parameters
C
C     Each segment descriptor occupies a contiguous
C     range of DAS d.p. addresses.
C
C        The DSK segment descriptor layout is:
C
C           +---------------------+
C           | Surface ID code     |      
C           +---------------------+
C           | Center ID code      |      
C           +---------------------+
C           | Data class code     |  
C           +---------------------+
C           | Data type           |  
C           +---------------------+
C           | Ref frame code      |  
C           +---------------------+
C           | Coord sys code      |  
C           +---------------------+
C           | Coord sys parameters|  {10 elements}
C           +---------------------+
C           | Min coord 1         |  
C           +---------------------+
C           | Max coord 1         |  
C           +---------------------+
C           | Min coord 2         |  
C           +---------------------+
C           | Max coord 2         |  
C           +---------------------+
C           | Min coord 3         |  
C           +---------------------+
C           | Max coord 3         |  
C           +---------------------+
C           | Start time          |  
C           +---------------------+
C           | Stop time           |  
C           +---------------------+
C
C     Parameters defining offsets for segment descriptor elements
C     follow.
C
C
C     Surface ID code:
C
      INTEGER               SRFIDX
      PARAMETER           ( SRFIDX = 1 )

C
C     Central ephemeris object NAIF ID:
C
      INTEGER               CTRIDX
      PARAMETER           ( CTRIDX = SRFIDX + 1 )

C
C     Data class:
C
C     The "data class" is a code indicating the category of
C     data contained in the segment.
C     
      INTEGER               CLSIDX
      PARAMETER           ( CLSIDX = CTRIDX + 1 )

C
C     Data type:
C
      INTEGER               TYPIDX
      PARAMETER           ( TYPIDX = CLSIDX + 1 )

C
C     Frame ID:
C
      INTEGER               FRMIDX
      PARAMETER           ( FRMIDX = TYPIDX + 1 )

C
C     Coordinate system code:
C
      INTEGER               SYSIDX
      PARAMETER           ( SYSIDX = FRMIDX + 1 )
      
C
C     Coordinate system parameter start index:
C
      INTEGER               PARIDX
      PARAMETER           ( PARIDX = SYSIDX + 1 )

C
C     Number of coordinate system parameters:
C
      INTEGER               NSYPAR
      PARAMETER           ( NSYPAR = 10 )

C
C     Ranges for coordinate bounds:
C
      INTEGER               MN1IDX
      PARAMETER           ( MN1IDX = PARIDX + NSYPAR )

      INTEGER               MX1IDX
      PARAMETER           ( MX1IDX = MN1IDX + 1 )

      INTEGER               MN2IDX
      PARAMETER           ( MN2IDX = MX1IDX + 1 )

      INTEGER               MX2IDX
      PARAMETER           ( MX2IDX = MN2IDX + 1 )

      INTEGER               MN3IDX
      PARAMETER           ( MN3IDX = MX2IDX + 1 )

      INTEGER               MX3IDX
      PARAMETER           ( MX3IDX = MN3IDX + 1 )
      
C
C     Coverage time bounds:
C
      INTEGER               BTMIDX
      PARAMETER           ( BTMIDX = MX3IDX + 1 )

      INTEGER               ETMIDX
      PARAMETER           ( ETMIDX = BTMIDX + 1 )

C
C     Descriptor size (24):
C
      INTEGER               DSKDSZ
      PARAMETER           ( DSKDSZ = ETMIDX )



C
C     Data class values:
C
C        Class 1 indicates a surface that can be represented as a
C                single-valued function of its domain coordinates.
C
C                An example is a surface defined by a function that
C                maps each planetodetic longitude and latitude pair to
C                a unique altitude.
C       
      INTEGER               SVFCLS
      PARAMETER           ( SVFCLS = 1 )
C
C        Class 2 indicates a general surface. Surfaces that
C                have multiple points for a given pair of domain
C                coordinates---for example, multiple radii for a given
C                latitude and longitude---belong to class 2.
C     
      INTEGER               GENCLS
      PARAMETER           ( GENCLS = 2 )

C
C
C     Coordinate system values:
C
C        The coordinate system code indicates the system to which the
C        tangential coordinate bounds belong. 
C
C        Code 1 refers to the planetocentric latitudinal system. 
C
C        In this system, the first tangential coordinate is longitude
C        and the second tangential coordinate is latitude. The third
C        coordinate is radius.
C         
C
      INTEGER               LATSYS
      PARAMETER           ( LATSYS = 1 )


C
C        Code 2 refers to the cylindrical system. 
C
C        In this system, the first tangential coordinate is radius and
C        the second tangential coordinate is longitude. The third,
C        orthogonal coordinate is Z.
C         
C
      INTEGER               CYLSYS
      PARAMETER           ( CYLSYS = 2 )

 
C
C        Code 3 refers to the rectangular system. 
C
C        In this system, the first tangential coordinate is X and
C        the second tangential coordinate is Y. The third,
C        orthogonal coordinate is Z.
C         
C
      INTEGER               RECSYS
      PARAMETER           ( RECSYS = 3 )


C
C        Code 4 refers to the planetodetic/geodetic system. 
C
C        In this system, the first tangential coordinate is longitude
C        and the second tangential coordinate is planetodetic
C        latitude. The third, orthogonal coordinate is altitude.
C         
C
      INTEGER               PDTSYS
      PARAMETER           ( PDTSYS = 4 )



C
C     End of include file dskdsc.inc
C

