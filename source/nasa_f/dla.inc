
C
C     Include file dla.inc
C
C     This include file declares parameters for DLA format 
C     version zero.
C
C        Version 3.0.1 17-OCT-2016 (NJB)
C
C           Corrected comment: VERIDX is now described as a DAS
C           integer address rather than a d.p. address.
C
C        Version 3.0.0 20-JUN-2006 (NJB)
C
C           Changed name of parameter DSCSIZ to DLADSZ.
C
C        Version 2.0.0 09-FEB-2005 (NJB)
C
C           Changed descriptor layout to make backward pointer 
C           first element.  Updated DLA format version code to 1.
C
C           Added parameters for format version and number of bytes per
C           DAS comment record.
C
C        Version 1.0.0 28-JAN-2004 (NJB)
C
      

C
C     DAS integer address of DLA version code.
C
      INTEGER               VERIDX
      PARAMETER           ( VERIDX = 1 )

C
C     Linked list parameters
C
C     Logical arrays (aka "segments") in a DAS linked array (DLA) file
C     are organized as a doubly linked list.  Each logical array may
C     actually consist of character, double precision, and integer
C     components.  A component of a given data type occupies a
C     contiguous range of DAS addresses of that type.  Any or all
C     array components may be empty.
C
C     The segment descriptors in a SPICE DLA (DAS linked array) file
C     are connected by a doubly linked list.  Each node of the list is
C     represented by a pair of integers acting as forward and backward
C     pointers.  Each pointer pair occupies the first two integers of a
C     segment descriptor in DAS integer address space.  The DLA file
C     contains pointers to the first integers of both the first and
C     last segment descriptors. 
C
C     At the DLA level of a file format implementation, there is
C     no knowledge of the data contents.  Hence segment descriptors
C     provide information only about file layout (in contrast with
C     the DAF system).  Metadata giving specifics of segment contents
C     are stored within the segments themselves in DLA-based file
C     formats.
C
C
C     Parameter declarations follow.
C
C     DAS integer addresses of first and last segment linked list
C     pointer pairs.  The contents of these pointers
C     are the DAS addresses of the first integers belonging
C     to the first and last link pairs, respectively.
C
C     The acronyms "LLB" and "LLE" denote "linked list begin"
C     and "linked list end" respectively.
C     
      INTEGER               LLBIDX
      PARAMETER           ( LLBIDX = VERIDX + 1 )

      INTEGER               LLEIDX
      PARAMETER           ( LLEIDX = LLBIDX + 1 )

C
C     Null pointer parameter. 
C
      INTEGER               NULPTR
      PARAMETER           ( NULPTR = -1 )

C
C     Segment descriptor parameters
C
C     Each segment descriptor occupies a contiguous
C     range of DAS integer addresses.
C
C        The segment descriptor layout is:
C
C           +---------------+
C           | BACKWARD PTR  | Linked list backward pointer
C           +---------------+
C           | FORWARD PTR   | Linked list forward pointer
C           +---------------+
C           | BASE INT ADDR | Base DAS integer address
C           +---------------+
C           | INT COMP SIZE | Size of integer segment component
C           +---------------+
C           | BASE DP ADDR  | Base DAS d.p. address
C           +---------------+
C           | DP COMP SIZE  | Size of d.p. segment component
C           +---------------+
C           | BASE CHR ADDR | Base DAS character address
C           +---------------+
C           | CHR COMP SIZE | Size of character segment component
C           +---------------+
C
C     Parameters defining offsets for segment descriptor elements
C     follow.
C
      INTEGER               BWDIDX
      PARAMETER           ( BWDIDX = 1 )

      INTEGER               FWDIDX
      PARAMETER           ( FWDIDX = BWDIDX + 1 )

      INTEGER               IBSIDX
      PARAMETER           ( IBSIDX = FWDIDX + 1 )

      INTEGER               ISZIDX
      PARAMETER           ( ISZIDX = IBSIDX + 1 )

      INTEGER               DBSIDX
      PARAMETER           ( DBSIDX = ISZIDX + 1 )

      INTEGER               DSZIDX
      PARAMETER           ( DSZIDX = DBSIDX + 1 )

      INTEGER               CBSIDX
      PARAMETER           ( CBSIDX = DSZIDX + 1 )

      INTEGER               CSZIDX
      PARAMETER           ( CSZIDX = CBSIDX + 1 )

C
C     Descriptor size:
C
      INTEGER               DLADSZ
      PARAMETER           ( DLADSZ = CSZIDX )


C
C     Other DLA parameters:
C
C
C     DLA format version.  (This number is expected to occur very
C     rarely at integer address VERIDX in uninitialized DLA files.)
C
      INTEGER               FMTVER
      PARAMETER           ( FMTVER = 1000000 )

C
C     Characters per DAS comment record.
C
      INTEGER               NCHREC
      PARAMETER           ( NCHREC = 1024 )


C
C     End of include file dla.inc
C

