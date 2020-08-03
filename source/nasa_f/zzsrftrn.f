C$Procedure ZZSRFTRN ( Surface name/ID mapping umbrella )

      SUBROUTINE ZZSRFTRN ( BODYID, SRFNAM, SURFID,  
     .                      USRCTR, FOUND,  UPDATE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Umbrella routine for surface name/ID mapping entry points
C     and the data structures they use.
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
C     NAIF_IDS
C
C$ Keywords
C
C     CONVERSION
C     DSK
C     ID
C     NAME
C     STRING
C     SURFACE
C
C$ Declarations

      IMPLICIT NONE
      
      INCLUDE 'srftrn.inc'
      INCLUDE 'zzctr.inc'
      INCLUDE 'zzsrftrn.inc'

      INTEGER               BODYID
      CHARACTER*(*)         SRFNAM
      INTEGER               SURFID
      INTEGER               USRCTR ( 2 )
      LOGICAL               FOUND
      LOGICAL               UPDATE
  
C$ Brief_I/O
C
C     Variable  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     BODYID     I   ZZSRFC2N, ZZSRFN2C
C     SRFNAM    I-O  ZZSRFC2N, ZZSRFN2C
C     SURFID    I-O  ZZSRFC2N, ZZSRFN2C
C     USRCTR    I-O  ZZSRFTRK
C     FOUND      O   ZZSRFC2N, ZZSRFN2C
C     UPDATE     O   ZZSRFTRK
C
C$ Detailed_Input
C
C     See the entry points for descriptions of their inputs.
C
C$ Detailed_Output
C 
C     See the entry points for descriptions of their outputs.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If this routine is called directly, the error 
C        SPICE(BOGUSENTRY) is signaled.
C
C$ Files
C    
C     Surface-to-name mappings may be defined at run time by loading
C     text kernels containing kernel variable assignments of the form
C
C        NAIF_SURFACE_NAME += ( <surface name 1>, ... )
C        NAIF_SURFACE_CODE += ( <surface code 1>, ... )
C        NAIF_SURFACE_BODY += ( <body code 1>,    ... )
C
C     Here the set of the three Ith list items on the right hand side
C     of the three assignments define the Ith surface name/ID mapping.
C
C     The same effect can be achieved using assignments formatted as
C     follows:
C
C        NAIF_SURFACE_NAME += <surface name 1>
C        NAIF_SURFACE_CODE += <surface code 1>
C        NAIF_SURFACE_BODY += <body code 1>
C
C        NAIF_SURFACE_NAME += <surface name 2>
C        NAIF_SURFACE_CODE += <surface code 2>
C        NAIF_SURFACE_BODY += <body code 2>
C
C           ...
C
C     Note the use of the
C
C        +=
C
C     operator; this operator appends to rather than overwrites the
C     kernel variable named on the left hand side of the assignment.
C
C$ Particulars
C
C     This umbrella routine contains declarations of data structures
C     that support surface name/ID mapping.
C
C     This umbrella routine contains the following entry points:
C
C        ZZSRFC2N   {Surface code to name}
C        ZZSRFN2C   {Surface name to code}
C        ZZSRFTRK   {Track surface map updates}
C
C
C$ Examples
C
C     See the routines
C
C        SRFC2S
C        SRFCSS
C        SRFS2C 
C        SRFSCC
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
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 01-APR-2016 (NJB) (EDW) (BVS)
C
C-&

C$ Index_Entries
C
C     surface name id mapping umbrella
C
C-&


C
C     SPICELIB functions
C
      INTEGER               ZZHASH2
      INTEGER               ZZHASHI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
C
C     Hash control area items.
C
      INTEGER               SIZIDX
      PARAMETER           ( SIZIDX =   0 )
 
      INTEGER               FREIDX
      PARAMETER           ( FREIDX =  -1 )

C
C     Local variables
C

C
C     Data structures in this package
C     ===============================
C
C     The kernel variable table
C     -------------------------
C
C     This table contains file scope arrays populated with the values
C     of the kernel variables that define the surface name/ID mapping.
C     These arrays contain
C
C        - surface names
C        - surface ID codes
C        - bodies associated with surface name/ID pairs
C        - an array of normalized names. These names are
C          upper case, left-justified, and compressed so that
C          the names contain no consecutive, embedded blanks
C
C
C     The surface ID table
C     --------------------
C
C     This table enables pairs of surface IDs and body IDs to be mapped
C     to surface names. The table consists of
C
C        - a singly linked list pool
C        - a list head array
C        - a pointer array that maps pool nodes to 
C          indices in the kernel variable table
C
C     The pointer array maps each node belonging to a collision list in
C     the pool to the index in the kernel table of the associated
C     values. The kernel table values used by this mapping are 
C
C        - the  surface ID code
C        - the body ID code 
C        - the original surface name
C
C     An integer hash function is used to map each surface ID to the
C     index in the list head array where the index of the head node for
C     that surface ID is located.
C
C     The layout of the structure is:
C
C                                                 Kernel variable table
C                                                 (only the portions 
C                                                 used here are shown)
C
C                                                         body IDs
C                                                            |
C   +-- integer_hash( surface ID )                           |  original
C   |                                            surface IDs |  surface
C   |                                                   |    |   names
C   |                                                   |    |    |
C   |  list heads     list pool   pointer array         |    |    |
C   |  +---------+    +--------+  +-----------+        +--+ +--+ +--+
C   |  |         |    |        |  |           |   +--->|  | |  | |  |
C   |  +---------+    +--------+  +-----------+   |    +--+ +--+ +--+
C   +->|         |-+  | ^   *  |->|           |---+    |  | |  | |  |
C      +---------+ |  +-|---|--+  +-----------+        +--+ +--+ +--+
C                  |    |   |
C          ...     |    |...|          ...                   ...
C                  |    |   |
C      +---------+ |  +-|---|--+  +-----------+        +--+ +--+ +--+
C      |         | +->| *   |  |->|           |---+ +->|  | |  | |  |
C      +---------+    +-----|--+  +-----------+   | |  +--+ +--+ +--+
C          ...           ...|          ...      +-|-+        ...
C      +---------+    +-----|--+  +-----------+ | |    +--+ +--+ +--+
C      |         |    |     v  |->|           |-+ +--->|  | |  | |  |
C      +---------+    +--------+  +-----------+        +--+ +--+ +--+ 
C
C      ----------------------------------------        --------------
C                       NROOM                              MXNSRF
C
C
C      The diagram above is not to scale: the arrays on the left have
C      available length NROOM, while the arrays on the right have
C      length MXNSRF. 
C
C      Note that the pool array is dimensioned (LBSNGL:NROOM). Elements
C      at indices -1 and 0 contain the size and location of the 
C      first free element, respectively.
C
C
C     The surface name table
C     ----------------------
C
C     This table enables pairs of surface names and body IDs to be
C     mapped to surface IDs. The structure is parallel to that of
C     the surface ID table; it contains
C
C        - a singly linked list pool
C        - a list head array
C        - a pointer array that maps pool nodes to 
C          indices in the kernel variable table
C
C     The pointer array maps each node belonging to a collision list in
C     the pool to the index in the kernel table of the associated
C     values. The kernel table values used by this mapping are 
C
C        - the normalized surface name
C        - the surface ID code
C        - the body ID code 
C
C     An string hash function is used to map each surface name to the
C     index in the list head array where the index of the head node for
C     that surface ID is located. 
C
C     The hash function is applied to the input string after it has
C     been normalized and then had all embedded blanks compressed out.
C     This allows the hash function terminate when it encounters the
C     first blank in the input string, while taking into account all
C     non-blank characters in the string. This makes it efficient while
C     enabling it to discriminate well between strings that may have
C     initial words in common. These compressed strings are not used
C     for any other purpose than hashing. For detection of the correct
C     matching elements in the kernel table, the normalized version
C     of the input string (which may contain blanks) is used.
C
C     The layout of the structure is:
C
C
C                                                 Kernel variable table
C                                                 (only the portions 
C                                                 used here are shown)
C
C                                                         body IDs
C                                                            |  
C   +-- string_hash(surface name)                            |  
C   |                                         normalized     |  surface
C   |                                         surface names  |   IDs
C   |                                                   |    |    |
C   |                                                   |    |    |
C   |  list heads     list pool   pointer array         |    |    |
C   |  +---------+    +--------+  +-----------+        +--+ +--+ +--+
C   |  |         |    |        |  |           |   +--->|  | |  | |  |
C   |  +---------+    +--------+  +-----------+   |    +--+ +--+ +--+
C   +->|         |-+  | ^   *  |->|           |---+    |  | |  | |  |
C      +---------+ |  +-|---|--+  +-----------+        +--+ +--+ +--+
C                  |    |   |
C          ...     |    |...|          ...                   ...
C                  |    |   |
C      +---------+ |  +-|---|--+  +-----------+        +--+ +--+ +--+
C      |         | +->| *   |  |->|           |---+ +->|  | |  | |  |
C      +---------+    +-----|--+  +-----------+   | |  +--+ +--+ +--+
C          ...           ...|          ...      +-|-+        ...
C      +---------+    +-----|--+  +-----------+ | |    +--+ +--+ +--+
C      |         |    |     v  |->|           |-+ +--->|  | |  | |  |
C      +---------+    +--------+  +-----------+        +--+ +--+ +--+
C
C
C      ----------------------------------------        --------------
C                       NROOM                              MXNSRF
C
C
C      The diagram above is not to scale: the arrays on the left have
C      available length NROOM, while the arrays on the right have
C      length MXNSRF. 
C
C      Note that the pool array is dimensioned (LBSNGL:NROOM). Elements
C      at indices -1 and 0 contain the size and location of the 
C      first free element, respectively.
C
C
C
C
C     Declarations of data structures
C     ===============================
C
C        Kernel variable table
C        =====================
C
C        Input names:                 KERNAM        
C        Input surface IDs:           KERSID
C        Input body IDs:              KERBID
C
      INTEGER               KERSID ( MXNSRF )
      INTEGER               KERBID ( MXNSRF )
      CHARACTER*(SFNMLN)    KERNAM ( MXNSRF )
      INTEGER               NKVAR
C
C        Normalized names:            NRMNAM
C
C     Each of these surface names is prefixed with an 11-character
C     string containing the associated body ID.
C

      CHARACTER*(SFNMLN)    NORNAM ( MXNSRF )

C
C
C        Surface ID table
C        ================
C
C        Surface ID list heads:       SIDHLS
C        Surface ID pool:             SIDPOL
C        Surface ID name pointers:    SIDIDX
C
C
      INTEGER               SIDHLS ( NROOM )
      INTEGER               SIDPOL ( LBSNGL : NROOM )
      INTEGER               SIDIDX ( NROOM )


C        Surface Name table
C        ==================
C
C        Surface name list heads:     SNMHLS
C        Surface name pool:           SNMPOL
C        Surface name ID pointers:    SNMIDP
C
C
      INTEGER               SNMHLS ( NROOM )
      INTEGER               SNMPOL ( LBSNGL : NROOM )
      INTEGER               SNMIDX ( NROOM )


C
C     Other local declarations:
C
      CHARACTER*(SFNMLN)    NSRFNM
      CHARACTER*(SFNMLN)    SQSHNM

      INTEGER               ITEMAT
      INTEGER               LOOKAT
      INTEGER               NODE

C
C     POLCTR tracks the state of the kernel pool.
C     SRFCTR tracks the state of the surface mapping
C     kernel variables.
C
      INTEGER               POLCTR ( CTRSIZ )
      INTEGER               SRFCTR ( CTRSIZ )

      LOGICAL               EXTKER
      LOGICAL               PASS1
      LOGICAL               LUPDTE

C
C
C     Saved variables
C     
      SAVE


C
C     Initial values
C
      DATA                  EXTKER / .FALSE.    /
      DATA                  POLCTR / CTRSIZ * 0 /
      DATA                  SRFCTR / CTRSIZ * 0 /
      DATA                  PASS1  / .TRUE.     /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZSRFTRN'                                    )
      CALL SETMSG ( 'ZZSRFTRN is an umbrella routine. It should '
     .//            'never be called directly.'                   )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)'                           )
      CALL CHKOUT ( 'ZZSRFTRN'                                    )
      RETURN




C$Procedure ZZSRFN2C ( Surface name to ID code mapping )

      ENTRY ZZSRFN2C ( SRFNAM, BODYID, SURFID, FOUND )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Map a surface name and body ID code to a surface ID code.
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
C     NAIF_IDS
C
C$ Keywords
C
C     CONVERSION
C     DSK
C     ID
C     NAME
C     STRING
C     SURFACE
C
C$ Declarations
C
C     CHARACTER*(*)         SRFNAM
C     INTEGER               BODYID
C     INTEGER               SURFID
C     LOGICAL               FOUND
C  
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     SRFNAM     I   Surface name.
C     BODYID     I   Body ID code.
C     SURFID     O   Surface ID code.
C     FOUND      O   Found flag.
C
C$ Detailed_Input
C
C     SRFNAM     is the name of surface to be translated.
C
C     BODYID     is the ID code of a body with which the surface
C                designated by SRFNAM is associated.
C
C$ Detailed_Output
C 
C     SURFID     is the surface ID code associated with the input 
C                surface name and body ID code. 
C
C                If multiple assignments for the input surface
C                name and body ID are present in the kernel pool,
C                the latest one is used to determine SURFID.
C
C     FOUND      is a logical flag that is .TRUE. if the inputs
C                were mapped to a surface ID code and .FALSE.
C                otherwise.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs when the data structures of this package
C         are initialized, it will be signaled by a routine in the
C         call tree of this routine.
C
C$ Files
C    
C     Surface-to-name mappings may be defined at run time by loading
C     text kernels containing kernel variable assignments of the form
C
C        NAIF_SURFACE_NAME += ( <surface name 1>, ... )
C        NAIF_SURFACE_CODE += ( <surface code 1>, ... )
C        NAIF_SURFACE_BODY += ( <body code 1>,    ... )
C
C     Here the set of the three Ith list items on the right hand side
C     of the three assignments define the Ith surface name/ID mapping.
C
C     The same effect can be achieved using assignments formatted as
C     follows:
C
C        NAIF_SURFACE_NAME += <surface name 1>
C        NAIF_SURFACE_CODE += <surface code 1>
C        NAIF_SURFACE_BODY += <body code 1>
C
C        NAIF_SURFACE_NAME += <surface name 2>
C        NAIF_SURFACE_CODE += <surface code 2>
C        NAIF_SURFACE_BODY += <body code 2>
C
C           ...
C
C     Note the use of the
C
C        +=
C
C     operator; this operator appends to rather than overwrites the
C     kernel variable named on the left hand side of the assignment.
C
C$ Particulars
C
C     This routine maps pairs of surface names and body ID codes to
C     surface ID codes. It relies on the mapping variables in the
C     kernel pool.
C
C     On the first pass through this routine, this routine
C     initializes the shared data structures of this package,
C     if the initialization has not already been done.
C
C$ Examples
C
C     See the routines
C
C        SRFS2C 
C        SRFSCC
C
C$ Restrictions
C
C     This routine must not be called by user applications.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 01-APR-2016 (NJB) (EDW) (BVS)
C
C-&

C$ Index_Entries
C
C     map surface name and body id to surface id
C
C-&

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSRFN2C' )

C
C     No result has been found.
C
      FOUND = .FALSE.


      IF ( PASS1 ) THEN
C
C        Initialize the surface kernel variable update counter
C        and the local pool counter. Note that this routine
C        is a "subsystem" as seen by its callers and a "user"
C        with respect to the kernel pool. Hence the different
C        initializations.
C
         CALL ZZCTRSIN ( SRFCTR )
         CALL ZZCTRUIN ( POLCTR ) 

C
C        Initialize local data structures. The first instance of this
C        call also sets a watch on the surface mapping kernel
C        variables.
C
         CALL ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                   EXTKER, NKVAR,  SNMHLS, SNMPOL, 
     .                   SNMIDX, SIDHLS, SIDPOL, SIDIDX ) 
C
C        Sync POLCTR with the kernel pool counter.
C        
         CALL ZZCVPOOL ( 'ZZSRFTRN', POLCTR, LUPDTE )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSRFN2C' )
            RETURN
         END IF

         PASS1  = .FALSE.

      END IF

C
C     Determine whether the data structures need to be updated
C     due to a change in the kernel pool contents. 
C
      CALL ZZCVPOOL ( 'ZZSRFTRN', POLCTR, LUPDTE )

      IF ( LUPDTE ) THEN         
C
C        Conservatively increment the ZZSRFTRN state counter in
C        expectation of successful update.
C
         CALL ZZCTRINC ( SRFCTR ) 
C
C        Initialize local data structures.
C
         CALL ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                   EXTKER, NKVAR,  SNMHLS, SNMPOL, 
     .                   SNMIDX, SIDHLS, SIDPOL, SIDIDX ) 

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSRFN2C' )
            RETURN
         END IF

      END IF

C
C     No translation can be done if the surface mapping variables
C     are not in the pool.
C
      IF ( .NOT. EXTKER ) THEN
         CALL CHKOUT ( 'ZZSRFN2C' )
         RETURN
      END IF

C
C     Get a "normalized" copy of the input name: left-justified,
C     compressed, upper case.
C
      CALL LJUCRS ( 1, SRFNAM, NSRFNM )

C
C     Get a "squished" version of the above name: a version
C     containing no blanks.
C
      CALL CMPRSS ( ' ', 0, NSRFNM, SQSHNM )

C
C     Find the hash value of the squished input name.
C
      LOOKAT = ZZHASH2( SQSHNM, SNMPOL(SIZIDX) )
      NODE   = SNMHLS ( LOOKAT )

      FOUND = .FALSE.
      
      IF ( NODE .GT. 0 ) THEN
C
C        Start at the head node and check each normalized name saved
C           for this hash value until we find a name and body ID that
C           match or run out of items in the collision list.
C
         DO WHILE (  ( NODE .GT. 0 ) .AND. ( .NOT. FOUND )  )

            FOUND  =       ( NSRFNM  .EQ. NORNAM( SNMIDX(NODE) ) )
     .               .AND. ( BODYID  .EQ. KERBID( SNMIDX(NODE) ) )
            
            ITEMAT = NODE
            NODE   = SNMPOL ( NODE )

         END DO
C
C        ITEMAT is the value of the last node checked, or
C        0 if the list is empty.
C       
      END IF

      IF ( FOUND ) THEN

         SURFID = KERSID( SNMIDX(ITEMAT) )

      END IF

      CALL CHKOUT ( 'ZZSRFN2C' )
      RETURN

 



C$Procedure ZZSRFC2N ( Surface ID code to name mapping )

      ENTRY ZZSRFC2N ( SURFID, BODYID, SRFNAM, FOUND )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Map a surface ID code and body ID code to a surface name.
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
C     NAIF_IDS
C
C$ Keywords
C
C     CONVERSION
C     DSK
C     ID
C     NAME
C     STRING
C     SURFACE
C
C$ Declarations
C
C     INTEGER               SURFID
C     INTEGER               BODYID
C     CHARACTER*(*)         SRFNAM
C     LOGICAL               FOUND
C  
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     SRFNAM     I   Surface name.
C     BODYID     I   Body ID code.
C     SURFID     O   Surface ID code.
C     FOUND      O   Found flag.
C
C$ Detailed_Input
C
C     SURFID     is a surface ID code.
C
C     BODYID     is the ID code of a body with which the surface
C                designated by SURFID is associated.
C
C$ Detailed_Output
C 
C     SRFNAM     is the name of the surface corresponding to the 
C                input surface ID code and body ID code.
C
C                If multiple assignments for the input surface
C                name and body ID are present in the kernel pool,
C                the latest one is used to determine SRFNAM.
C
C     FOUND      is a logical flag that is .TRUE. if the inputs
C                were mapped to a surface name and .FALSE.
C                otherwise.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs when the data structures of this package
C         are initialized, it will be signaled by a routine in the
C         call tree of this routine.
C
C$ Files
C    
C     Surface-to-name mappings may be defined at run time by loading
C     text kernels containing kernel variable assignments of the form
C
C        NAIF_SURFACE_NAME += ( <surface name 1>, ... )
C        NAIF_SURFACE_CODE += ( <surface code 1>, ... )
C        NAIF_SURFACE_BODY += ( <body code 1>,    ... )
C
C     Here the set of the three Ith list items on the right hand side
C     of the three assignments define the Ith surface name/ID mapping.
C
C     The same effect can be achieved using assignments formatted as
C     follows:
C
C        NAIF_SURFACE_NAME += <surface name 1>
C        NAIF_SURFACE_CODE += <surface code 1>
C        NAIF_SURFACE_BODY += <body code 1>
C
C        NAIF_SURFACE_NAME += <surface name 2>
C        NAIF_SURFACE_CODE += <surface code 2>
C        NAIF_SURFACE_BODY += <body code 2>
C
C           ...
C
C     Note the use of the
C
C        +=
C
C     operator; this operator appends to rather than overwrites the
C     kernel variable named on the left hand side of the assignment.
C
C$ Particulars
C
C     This routine maps pairs of surface ID codes and body ID codes to
C     surface names. It relies on the mapping variables in the kernel
C     pool.
C
C     On the first pass through this routine, this routine
C     initializes the shared data structures of this package,
C     if the initialization has not already been done.
C
C$ Examples
C
C     See the routines
C
C        SRFC2S
C        SRFCSS
C
C$ Restrictions
C
C     This routine must not be called by user applications.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 01-APR-2016 (NJB) (EDW) (BVS)
C
C-&

C$ Index_Entries
C
C     map surface id and body id to surface name
C
C-&


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSRFC2N' )

C
C     No result has been found.
C
      FOUND = .FALSE.


      IF ( PASS1 ) THEN
C
C        Initialize the surface kernel variable update counter
C        and the local pool counter. Note that this routine
C        is a "subsystem" as seen by its callers and a "user"
C        with respect to the kernel pool. Hence the different
C        initializations.
C
         CALL ZZCTRSIN ( SRFCTR )
         CALL ZZCTRUIN ( POLCTR ) 

C
C        Initialize local data structures. The first instance of this
C        call also sets a watch on the surface mapping kernel
C        variables.
C
         CALL ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                   EXTKER, NKVAR,  SNMHLS, SNMPOL, 
     .                   SNMIDX, SIDHLS, SIDPOL, SIDIDX ) 
C
C        Sync SRFCTR with the kernel pool counter.
C        
         CALL ZZCVPOOL ( 'ZZSRFTRN', POLCTR, LUPDTE )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSRFC2N' )
            RETURN
         END IF


         PASS1  = .FALSE.

      END IF

C
C     Determine whether the data structures need to be updated
C     due to a change in the kernel pool contents.
C
      CALL ZZCVPOOL ( 'ZZSRFTRN', POLCTR, LUPDTE )

      IF ( LUPDTE ) THEN
C
C        Conservatively increment the ZZSRFTRN state counter in
C        expectation of successful update.
C
         CALL ZZCTRINC ( SRFCTR ) 

C
C        Initialize local data structures.
C
         CALL ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                   EXTKER, NKVAR,  SNMHLS, SNMPOL, 
     .                   SNMIDX, SIDHLS, SIDPOL, SIDIDX ) 

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSRFC2N' )
            RETURN
         END IF

      END IF

C
C     No translation can be done if the surface mapping variables
C     are not in the pool.
C
      IF ( .NOT. EXTKER ) THEN
         CALL CHKOUT ( 'ZZSRFC2N' )
         RETURN
      END IF

C
C     Find the hash value of the squished input name.
C
      LOOKAT = ZZHASHI( SURFID, SIDPOL(SIZIDX) )
      NODE   = SIDHLS ( LOOKAT )
 
      FOUND  = .FALSE.
      
      IF ( NODE .GT. 0 ) THEN
C
C        Start at the head node and check each normalized name saved
C           for this hash value until we find a name and body ID that
C           match or run out of items in the collision list.
C
         DO WHILE (  ( NODE .GT. 0 ) .AND. ( .NOT. FOUND )  )

            FOUND  =       ( SURFID  .EQ. KERSID( SIDIDX(NODE) ) )
     .               .AND. ( BODYID  .EQ. KERBID( SIDIDX(NODE) ) )
            
            ITEMAT = NODE
            NODE   = SIDPOL ( NODE )

         END DO
C
C        ITEMAT is the value of the last node checked, or
C        0 if the list is empty.
C       
      END IF

      IF ( FOUND ) THEN

         SRFNAM = KERNAM( SIDIDX(ITEMAT) )

      END IF

      CALL CHKOUT ( 'ZZSRFC2N' )
      RETURN



 


C$Procedure ZZSRFTRK ( Surface mapping tracker )

      ENTRY ZZSRFTRK ( USRCTR, UPDATE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Check a user counter against the surface mapping counter
C     maintained by this package. Indicate whether the caller 
C     needs to update variables associated with the user counter.
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
C     NAIF_IDS
C
C$ Keywords
C
C     CONVERSION
C     DSK
C     ID
C     NAME
C     STRING
C     SURFACE
C
C$ Declarations
C
C     INTEGER               USRCTR ( * )
C     LOGICAL               UPDATE
C  
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     USRCTR    I-O  User counter.
C     UPDATE     O   Update flag.
C
C$ Detailed_Input
C
C     USRCTR     is a counter passed in by a calling routine. Normally
C                USRCTR would be used to enable the application calling
C                this routine to determine whether local associations
C                of surface names and IDs are up to date.
C 
C
C$ Detailed_Output
C
C     USRCTR     is the user counter, updated if necessary to match the
C                counter maintained by this package.
C
C
C     UPDATE     is a logical flag indicating whether USRCTR was
C                updated.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs when the data structures of this package
C         are initialized, it will be signaled by a routine in the
C         call tree of this routine.
C
C$ Files
C    
C     Surface-to-name mappings may be defined at run time by loading
C     text kernels containing kernel variable assignments of the form
C
C        NAIF_SURFACE_NAME += ( <surface name 1>, ... )
C        NAIF_SURFACE_CODE += ( <surface code 1>, ... )
C        NAIF_SURFACE_BODY += ( <body code 1>,    ... )
C
C     Here the set of the three Ith list items on the right hand side
C     of the three assignments define the Ith surface name/ID mapping.
C
C     The same effect can be achieved using assignments formatted as
C     follows:
C
C        NAIF_SURFACE_NAME += <surface name 1>
C        NAIF_SURFACE_CODE += <surface code 1>
C        NAIF_SURFACE_BODY += <body code 1>
C
C        NAIF_SURFACE_NAME += <surface name 2>
C        NAIF_SURFACE_CODE += <surface code 2>
C        NAIF_SURFACE_BODY += <body code 2>
C
C           ...
C
C     Note the use of the
C
C        +=
C
C     operator; this operator appends to rather than overwrites the
C     kernel variable named on the left hand side of the assignment.
C
C$ Particulars
C
C     This routine allows SPICELIB routines to determine whether the
C     the surface mapping has been updated.
C
C     On every pass through this routine, this routine tests whether
C     the local mapping data structures need to be updated. If they do,
C     the data structures are re-initialized using the surface mapping
C     variables, if these are present in the kernel pool.
C
C$ Examples
C
C     See use of this routine in SPICELIB high-level geometry 
C     routines such as
C
C        ILLUMF
C        SINCPT
C        SUBPNT
C
C$ Restrictions
C
C     This routine must not be called by user applications.     
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C     E.D. Wright     (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 01-APR-2016 (NJB) (EDW) (BVS)
C
C-&

C$ Index_Entries
C
C     track surface mapping variable updates
C
C-&

C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF


      IF ( PASS1 ) THEN
C
C        Check in because ZZSRFKER can fail.
C
         CALL CHKIN ( 'ZZSRFTRK' )
C
C        Initialize the surface kernel variable update counter
C        and the local pool counter. Note that this routine
C        is a "subsystem" as seen by its callers and a "user"
C        with respect to the kernel pool. Hence the different
C        initializations.
C
         CALL ZZCTRSIN ( SRFCTR )
         CALL ZZCTRUIN ( POLCTR )         

C
C        Initialize local data structures. The first instance of this
C        call also sets a watch on the surface mapping kernel
C        variables.
C
         CALL ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                   EXTKER, NKVAR,  SNMHLS, SNMPOL, 
     .                   SNMIDX, SIDHLS, SIDPOL, SIDIDX ) 
C
C        Sync SRFCTR with the kernel pool counter.
C        
         CALL ZZCVPOOL ( 'ZZSRFTRN', POLCTR, LUPDTE )
C
C        Check out here since this routine doesn't check out
C        before its normal exit.
C
         CALL CHKOUT ( 'ZZSRFTRK' )

         IF ( FAILED() ) THEN
            RETURN
         END IF

         PASS1  = .FALSE.

      END IF
 
C
C     Check for updates to the kernel pool variables.
C
      CALL ZZCVPOOL ( 'ZZSRFTRN', POLCTR, LUPDTE )      

      IF ( LUPDTE ) THEN
C
C        Check in because ZZSRFKER can fail.
C
         CALL CHKIN ( 'ZZSRFTRK' )
         
C
C        Conservatively increment the ZZSRFTRN state counter in
C        expectation of successful update.
C
         CALL ZZCTRINC ( SRFCTR ) 

C
C        Update kernel pool mapping lists and hashes.
C  
         CALL ZZSRFKER ( KERNAM, NORNAM, KERSID, KERBID, 
     .                   EXTKER, NKVAR,  SNMHLS, SNMPOL,
     .                   SNMIDX, SIDHLS, SIDPOL, SIDIDX ) 

         CALL CHKOUT ( 'ZZSRFTRK' )

         IF ( FAILED() ) THEN
            RETURN
         END IF

      END IF

C
C     Check the input counter against the ZZSRFTRN counter;
C     sync the user counter.
C
      CALL ZZCTRCHK ( SRFCTR, USRCTR, UPDATE )
      
      RETURN
      END 

