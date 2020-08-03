C$Procedure ZZSRFINI ( Private --- Surface-Code Hash Initialization )
 
      SUBROUTINE ZZSRFINI ( NORNAM, CODES,  BODIES,  
     .                      NVALS,  MAXVAL, SNMHLS, SNMPOL,
     .                      SNMIDX, SIDHLS, SIDPOL, SIDIDX )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines.  Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Initialize the name-based and ID-based hashes used for efficient
C     access to surface-name mapping arrays. This routine should be
C     called by ZZSRFTRN and ZZSRFKER only.
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
C     UTILITY
C
C$ Declarations
 
      IMPLICIT NONE
 
      INCLUDE              'srftrn.inc'
      INCLUDE              'zzsrftrn.inc'
 
      CHARACTER*(*)         NORNAM (          * )
      INTEGER               CODES  (          * )
      INTEGER               BODIES (          * )
      INTEGER               NVALS
      INTEGER               MAXVAL
      INTEGER               SNMHLS (          * )
      INTEGER               SNMPOL ( LBSNGL : * )
      INTEGER               SNMIDX (          * )
      INTEGER               SIDHLS (          * )
      INTEGER               SIDPOL ( LBSNGL : * )
      INTEGER               SIDIDX (          * )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     NORNAM     I   Array of normalized surface names
C     CODES      I   Array of surface ID codes for NAMES/NORNAM
C     BODIES     I   Array of body ID codes
C     NVALS      I   Length of NAMES, NORNAM, and CODES arrays
C     MAXVAL     I   Size of the hash arrays
C     SNMHLS     O   Surface name-based hash head node pointer list
C     SNMPOL     O   Surface name-based hash node collision list
C     SNMIDX     O   Surface name-based hash index storage array 
C     SIDHLS     O   Surface ID-based hash head node pointer list
C     SIDPOL     O   Surface ID-based hash node collision list
C     SIDIDX     O   Surface ID-based hash index storage array
C     LBSNGL     P   Lower bound of hash pool arrays
C     SFNMLN     P   Maximum length of surface name strings
C
C$ Detailed_Input
C
C     NORNAM    is the array of normalized surface names, made from
C               elements of NAMES by upper-casing, left-justifying, and
C               compressing groups of spaces to a single space. This
C               represents the canonical member of the equivalence
C               class to which each parallel entry in NAMES belongs.
C
C               This array is parallel to CODES and BODIES.
C
C     CODES     is the array of surface codes extracted. This array is
C               parallel to NAMES and NORNAM.
C
C     BODIES    is the array of body ID codes associated with the input
C               surface names.  
C
C     NVALS     is the number of items contained in NAMES, NORNAM,
C               CODES.
C
C     MAXVAL    is the output hash size.
C
C$ Detailed_Output
C
C     All output arrays must be declared with the dimension MAXVAL.
C     MAXVAL must be greater than or equal to NVALS.
C
C     SNMHLS
C     SNMPOL    are the surface name-based hash head node pointer and
C               collision lists. Together with the arrays SNMIDX,
C               NORNAM and BODIES, they enable mapping pairs of
C               normalized surface names and body ID codes to surface
C               ID codes.
C
C     SNMIDX    is the surface name-based hash index storage array. It
C               maps nodes in the name collision list to entries in the 
C               parallel NORNAM and BODIES arrays.
C
C     SIDHLS
C     SIDPOL    are the surface ID-based hash head node pointer and  
C               collision lists. Together with the arrays SIDIDX,
C               CODES and BODIES, they enable mapping pairs of
C               surface ID codes and body ID codes to surface
C               names.
C
C     SIDIDX    is the surface ID-based hash index storage array. It
C               maps nodes in the ID collision list to entries in the 
C               parallel CODES and BODIES arrays.
C 
C$ Parameters
C
C     LBSNGL    is the lower bound of the hashes' collision list array.
C
C     SFNMLN    is the maximum length of a surface name. Defined in the
C               include file 'srftrn.inc'.
C
C$ Exceptions
C
C     1) If the input number of bodies NVALS is not less than or equal
C        to the size of the output hash, the error 'SPICE(BUG1)' will be
C        signaled.
C     
C     2) If registering an ID in the output ID-based hash fails, the
C        error 'SPICE(BUG2)' will be signaled.
C
C     3) If registering a name in the output name-based hash fails,
C        the error 'SPICE(BUG3)' will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This is a utility routine used for initializing the hashes
C     facilitating efficient surface name-ID translation in ZZSRFTRN.
C
C     The order of mappings in the input arrays determines their
C     priority, with the mapping having the lowest priority being first
C     and the mapping with the highest priority being last.

C     If more than one entry with a particular normalized name and body
C     ID is present in the input arrays, only the latest entry is
C     registered in the name-based hash.
C
C     If more than one entry with a particular surface ID and body ID
C     is present in the input arrays, only the latest entry that maps
C     to a not-yet-registered normalized name is registered in the
C     ID-based hash. Registering IDs only for not-yet-registered names
C     achieves masking all IDs with the lower priority in cases when a
C     single normalized name and body ID map to more than one surface
C     ID.
C     
C$ Examples
C
C     See the routine ZZSRFTRN.
C
C$ Restrictions
C
C     1)  This routine is intended only for use by ZZSRFTRN and
C         ZZSRFKER.
C
C     2)  All output hash arrays must be declared with the same
C         dimension which is greater than or equal to MAXVAL.
C
C     3)  The order of mappings in the input arrays determines the
C         priority, with the mapping with the lowest priority being the
C         first and the mapping with the highest priority being the
C         last.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     B.V. Semenov       (JPL)
C     M.J. Spencer       (JPL)
C     W.L. Taber         (JPL)
C     F.S. Turner        (JPL)
C     E.D. Wright        (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 03-DEC-2015 (NJB) (BVS) (EDW)
C
C-&
 

C
C     SPICELIB functions.
C
      INTEGER               ZZHASH2
      INTEGER               ZZHASHI

C
C     Hash control area items.
C
      INTEGER               SIZIDX
      PARAMETER           ( SIZIDX =   0 )
 
      INTEGER               FREIDX
      PARAMETER           ( FREIDX =  -1 )

C
C     Local Variables
C 
      CHARACTER*(SFNMLN)    SQSHNM

      INTEGER               HEAD
      INTEGER               I
      INTEGER               ITEMAT
      INTEGER               LOOKAT
      INTEGER               NODE

      LOGICAL               FULL
      LOGICAL               IDNEW
      LOGICAL               LFOUND
      LOGICAL               NAMNEW
C
C     Consistency check.
C
      IF ( MAXVAL .LT. NVALS ) THEN

         CALL CHKIN  ( 'ZZSRFINI'                                )
         CALL SETMSG ( 'There is an inconsistency between the '  
     .   //            'number of input bodies and the size of ' 
     .   //            'the output hashes. The number of input ' 
     .   //             'bodies was #. The size of the output '  
     .   //            'hashes was #.'                           )
         CALL ERRINT ( '#', NVALS                                )
         CALL ERRINT ( '#', MAXVAL                               )
         CALL SIGERR ( 'SPICE(BUG1)'                             )
         CALL CHKOUT ( 'ZZSRFINI'                                )
         RETURN

      END IF

C
C     Initialize output hashes. Set all collision list pointers
C     to 0, which is the null value.
C
      CALL ZZHSIINI ( MAXVAL, SIDHLS, SIDPOL )
      CALL ZZHSCINI ( MAXVAL, SNMHLS, SNMPOL )

      CALL CLEARI ( SIDPOL(0), SIDPOL(1))
      CALL CLEARI ( SNMPOL(0), SNMPOL(1))
       
C
C     Loop through the input arrays to populate hashes. We do it
C     backwards to pick and register only the highest priority (latest)
C     values for each pair of normalized surface name and body ID code.
C
      DO I = NVALS, 1, -1 
C
C        Register this normalized surface name and body ID, but only if
C        this pair is not already in the hash.
C
C        We must traverse the collision list for the normalized surface
C        name "manually," since we have to check the body ID for each
C        matching name.
C
C        Use hash function to get index of the head node.
C
         CALL CMPRSS ( ' ', 0, NORNAM(I), SQSHNM )

         LOOKAT = ZZHASH2( SQSHNM, SNMPOL(SIZIDX) )

         HEAD   = SNMHLS ( LOOKAT )

C
C        Indicate name and body were not found to begin with.
C
         LFOUND  = .FALSE.
         ITEMAT  = 0
         NAMNEW  = .TRUE.

C
C        See if this normalized name and corresponding body ID are,
C        respectively, in the normalized name list and body ID list.
C        Note that the body ID list is not a parallel array to the
C        normalized name array: we use the name pool pointer array
C        SNMIDX to indicate the location of the body ID corresponding
C        to a name.
C
         NODE = HEAD

         IF ( NODE .GT. 0 ) THEN
C
C           Start at the head node and check each normalized name saved
C           for this hash value until we find a name and body ID that
C           match or run out of items in the collision list.
C
            DO WHILE (  ( NODE .GT. 0 ) .AND. ( .NOT. LFOUND )  )

               LFOUND =       ( NORNAM( SNMIDX(NODE) ) .EQ. NORNAM(I) )
     .                  .AND. ( BODIES( SNMIDX(NODE) ) .EQ. BODIES(I) )

               ITEMAT = NODE
               NODE   = SNMPOL ( NODE )

            END DO
C
C           ITEMAT is the value of the last node in the list, or
C           0 if the list is empty.
C
            NAMNEW = .NOT. LFOUND

         END IF


         IF ( NAMNEW ) THEN
C
C           We need to add the current normalized name and BODY ID
C           to the hash. Make sure there's room.
C
            FULL = ( SNMPOL(FREIDX) .GT. SNMPOL(SIZIDX) )
           
            IF ( FULL ) THEN

               CALL CHKIN  ( 'ZZSRFINI'                          )
               CALL SETMSG ( 'Could not add name # body ID # '
     .         //            'to the hash.'                      )
               CALL ERRCH  ( '#', NORNAM(I)                      )
               CALL ERRINT ( '#', BODIES(I)                      )
               CALL SIGERR ( 'SPICE(BUG2)'                       )
               CALL CHKOUT ( 'ZZSRFINI'                          ) 
               RETURN

            ELSE
C
C              Store the item at the first free location in
C              the collision pool.
C               
               NODE           = SNMPOL(FREIDX)
               SNMPOL(FREIDX) = SNMPOL(FREIDX) + 1

               IF ( HEAD .GT. 0 ) THEN
C
C                 Link the new entry at the tail of the applicable
C                 collision list. The index of the tail node is ITEMAT.
C
                  SNMPOL(ITEMAT) = NODE

               ELSE
C
C                 Insert the new head node into the head list.
C
                  SNMHLS(LOOKAT) = NODE

               END IF
C
C              Set the index in the data arrays for the new pool
C              entry.
C
               SNMIDX(NODE) = I

            END IF            

C
C           NAMNEW indicates that the Ith normalized name and body ID
C           pair was not in the hash prior to the above block of code.
C
C           We may have a situation when a single normalized surface
C           name and body ID pair maps to more than one surface ID. In
C           such cases we want to completely mask all surface IDs with
C           the lower priority. This is easy to do by simply not
C           attempting to register any more surface IDs if the name is
C           already registered.
C
C           Register this surface ID and body ID pair, but only if it
C           is not already in the hash.
C
C           We must traverse the collision list for the normalized
C           surface name "manually," since we have to check the body ID
C           for each matching surface ID.
C
C           Use hash function to get index of the head node.
C  
            LOOKAT = ZZHASHI( CODES(I), SIDPOL(SIZIDX) )

            HEAD   = SIDHLS ( LOOKAT )

C
C           Indicate surface ID and body were not found to begin with.
C
            LFOUND  = .FALSE.
            ITEMAT  = 0
            IDNEW   = .TRUE.

C
C           See if this surface ID and corresponding body ID are,
C           respectively, in the surface ID list and body ID list.
C
            NODE = HEAD

            IF ( NODE .GT. 0 ) THEN
C
C              Start at the head node and check each surface ID saved
C              for this hash value until we find a surface ID and body
C              ID that match or run out of items in this collision
C              list.
C
               DO WHILE (  ( NODE .GT. 0 ) .AND. ( .NOT. LFOUND )  )

                  LFOUND =      ( CODES (SIDIDX(NODE)) .EQ. CODES (I) )
     .                    .AND. ( BODIES(SIDIDX(NODE)) .EQ. BODIES(I) )

                  ITEMAT = NODE
                  NODE   = SIDPOL ( NODE )

               END DO
C
C              ITEMAT is the value of the last node in the list, or
C              0 if the list is empty.
C
               IDNEW = .NOT. LFOUND

            END IF


            IF ( IDNEW ) THEN
C
C              We need to add the current surface ID and BODY ID
C              to the hash. Make sure there's room.
C
               FULL = ( SIDPOL(FREIDX) .GT. SIDPOL(SIZIDX) )
           
               IF ( FULL ) THEN

                  CALL CHKIN  ( 'ZZSRFINI'                             )
                  CALL SETMSG ( 'Could not add surface ID # body ID # '
     .            //            'to the hash.'                         )
                  CALL ERRINT ( '#', CODES(I)                          )
                  CALL ERRINT ( '#', BODIES(I)                         )
                  CALL SIGERR ( 'SPICE(BUG3)'                          )
                  CALL CHKOUT ( 'ZZSRFINI'                             )
                  RETURN

               ELSE
C
C                 Store the item at the first free location in the
C                 collision pool.
C               
                  NODE           = SIDPOL(FREIDX)
                  SIDPOL(FREIDX) = SIDPOL(FREIDX) + 1

                  IF ( HEAD .GT. 0 ) THEN
C
C                    Link the new entry at the tail of the applicable
C                    collision list. The index of the tail node is
C                    ITEMAT.
C
                     SIDPOL(ITEMAT) = NODE

                  ELSE
C
C                    Insert the new head node into the head list.
C
                     SIDHLS(LOOKAT) = NODE

                  END IF
C
C                 Set the index in the data arrays for the new pool
C                 entry.
C
                  SIDIDX(NODE) = I

               END IF            

            END IF
C
C           We've processed the new (surface ID, body ID) pair.
C
         END IF
C
C        We've processed the Ith mapping between (surface name, body
C        ID) and (surface ID, body ID).
C 
      END DO
 
      RETURN
      END
