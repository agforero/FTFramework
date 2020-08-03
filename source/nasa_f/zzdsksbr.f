C$Procedure ZZDSKSBR ( DSK, remove entries from API segment buffer )
 
      SUBROUTINE ZZDSKSBR ( NEEDED, 
     .                      MAXBOD, STSIZE, BTBODY, BTNBOD, BTSEGP,  
     .                      BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                      STOFF,  STCTR,  STRAD                   )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Remove entries from DSK API segment buffer data structures
C     to make room for a new entry.
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
C     DSK
C
C$ Keywords
C
C     DLA
C     DSK
C     TOPOGRAPHY
C     UTILITY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'

      INTEGER               NEEDED
      INTEGER               MAXBOD
      INTEGER               STSIZE 
      INTEGER               BTBODY ( * )
      INTEGER               BTNBOD
      INTEGER               BTSEGP ( * )
      INTEGER               BTSTSZ ( * )
      INTEGER               STHAN  ( * )
      DOUBLE PRECISION      STDSCR ( DSKDSZ, * )
      INTEGER               STDLAD ( DLADSZ, * )
      INTEGER               STFREE
      DOUBLE PRECISION      STOFF  ( 3, * )
      DOUBLE PRECISION      STCTR  ( 3, * )
      DOUBLE PRECISION      STRAD  ( * )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     NEEDED     I   Number of segment table entries needed.
C     MAXBOD     I   Maximum size of body table.
C     STSIZE     I   Maximum size of segment table.
C     BTBODY    I-O  Body table's body ID array.
C     BTNBOD    I-O  Number of entries in body table.
C     BTSEGP    I-O  Array of pointers from bodies to segment table.
C     BTSTSZ    I-O  Array of sizes of body segment lists.
C     STHAN     I-O  Array of handles of segments.
C     STDSCR    I-O  Array of DSK descriptors of segments.
C     STDLAD    I-O  Array of DLA descriptors of segments.
C     STFREE    I-O  Index of first free entry in segment table.
C     STOFF     I-O  Offsets of segment frame centers from body.
C     STCTR     I-O  Centers of bounding spheres for segments.
C     STRAD     I-O  Radii of bounding spheres for segments.
C
C$ Detailed_Input
C
C
C     NEEDED     is the number of segment table entries needed to store
C                the segment list for a new body.
C
C     MAXBOD     is the maximum size of the body table. BTBODY, BTSEGP,
C                and BTSTSZ must each be declared by the caller with
C                size at least MAXBOD.
C
C     STSIZE     is the size of the segment table. STHAN, STDSCR,
C                STDLAD, STOFF, STCTR, and STRAD must be declared large
C                enough to hold STSIZE entries.
C
C     BTBODY     is a table of body IDs. 
C
C     BTNBOD     is the number of body IDs currently in BTBODY.
C
C     BTSEGP     is an array of start indices in the segment table
C                of existing entries.
C
C     BTSTSZ     is an array of segment list sizes. The Ith entry
C                of BTSTSZ is the length of the segment list for the
C                Ith body.
C
C     STHAN      is an array of DAS handles. Each entry of STHAN that
C                is in use contains the handle of the DAS file 
C                containing the DSK segment to which that entry
C                corresponds. The Ith entries of STHAN, STDSCR, STDLAD,
C                STOFF, STCTR, and STRAD correspond to the same DSK
C                segment.
C
C     STDSCR     is an array of DSK descriptors.
C
C     STDLAD     is an array of DLA descriptors.
C
C     STFREE     is the index of the first free entry in the segment
C                table. 
C
C     STOFF      is an array of offsets of segment frame centers
C                from the central bodies of the segment. These offsets
C                are expressed in the reference frame of the segment.
C                They are constant vectors. Units are km.
C
C     STCTR      is an array of centers of outer bounding spheres for
C                segments. Each segment has an outer bounding sphere
C                that completely encloses that segment. Each center
C                is a vector expressed in the reference frame of the 
C                the segment, and it is an offset from the frame's
C                center. Units are km.
C
C     STRAD      is an array of radii of outer bounding spheres of
C                segments. Units are km.
C
C
C$ Detailed_Output
C
C
C     BTBODY     is the input body ID table, modified by this routine.
C                If it was necessary to delete a body from the table to
C                make room, that has been done.
C
C
C     BTNBOD     is the number of bodies in the body table.
C
C
C     BTSEGP     is an array of start indices in the segment table
C                of existing entries, updated to reflect deletion of
C                segments.
C
C
C     BTSTSZ     is the array of sizes of segment lists of bodies. If
C                it was necessary to make room, entries of BTSTSZ
C                will have been deleted.
C
C
C     STHAN      is the segment handle array. If it was necessary to
C                make room, entries of STHAN will have been deleted.
C
C
C     STDSCR     is the segment DSK descriptor array. If it was
C                necessary to make room, entries of STDESR will have
C                been deleted.
C
C
C     STDLAD     is the segment DLA descriptor array. If it was
C                necessary to make room, entries of STDLAD will have
C                been deleted.
C
C
C     STFREE     is the index of the first free element in each 
C                segment table array.
C
C
C     STOFF      is the segment frame center offset array.
C                If it was necessary to make room, entries of
C                STOFF will have been deleted.
C
C
C     STCTR      is the segment bounding sphere center array.
C                If it was necessary to make room, entries of
C                STCTR will have been deleted.
C
C
C     STRAD      is the segment bounding sphere radius array.
C                If it was necessary to make room, entries of
C                STRAD will have been deleted. 
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If MAXBOD is less than 1, the error SPICE(INVALIDTABLESIZE)
C         is signaled.
C
C     2)  If the segment table doesn't have at least NEEDED entries,
C         the error SPICE(SEGTABLETOOSMALL) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine makes room for new entries in the data structures
C     used by SPICELIB geometry APIs to perform efficient computations
C     using DSK data. It does so by deleting existing entries, if
C     necessary. The umbrella routine in which these structures
C     are declared is ZZDSKSBF.
C
C     In a sense, ZZDSKSBF sits "above" the ZZDSKBSR subsystem.
C     High-level geometry routines access the data stored in the
C     ZZDSKSBF arrays directly; they don't call entry points of
C     ZZDSKBSR. 
C
C     The ZZDSKSBF subsystem maintains two logical tables: a body
C     table and a segment table. Both tables are designed to store
C     a fixed maximum number of items; there is no "virtual storage"
C     concept that applies. 
C
C     Items in both tables have priority assigned on a first-in,
C     first-out basis. When it becomes necessary to remove an
C     item to make room for another, the lowest-indexed item is
C     removed first, and the tables are compressed so that the
C     remaining items are contiguous.
C
C     When the state of the underlying ZZDSKBSR system changes due to
C     loading or unloading of DSK files, the ZZDSKSBF body and segment
C     tables must be re-initialized. Entry points of ZZDSKSBF are
C     expected to perform this action by calling ZZDSKSBI.
C
C$ Examples
C
C     See usage in ZZDSKSBF.
C
C$ Restrictions
C
C     1) This is a private routine. 
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 11-JUL-2016 (NJB) 
C
C        30-JAN-2016 (NJB) 
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     remove entries from DSK API segment buffer
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C
      INTEGER               AVAIL
      INTEGER               I
      INTEGER               J
      INTEGER               NBTDEL
      INTEGER               NSTDEL


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKSBR' )
 

      IF ( NEEDED .GT. STSIZE ) THEN
C
C        There's not enough room in the entire segment table.
C
         CALL SETMSG ( 'Size of segment table is #; number of '
     .   //            'entries requested is #. The segment '    
     .   //            'table is supposed to be declared '
     .   //            'with sufficient size to accommodate all '
     .   //            'loaded DSK segments.'                    )
         CALL ERRINT ( '#', STSIZE                               )
         CALL ERRINT ( '#', NEEDED                               )
         CALL SIGERR ( 'SPICE(SEGTABLETOOSMALL)'                 )
         CALL CHKOUT ( 'ZZDSKSBR'                                )
         RETURN

      END IF

C
C     We can't make room in a body table of zero size.
C
      IF ( MAXBOD .LT. 1 ) THEN

         CALL SETMSG ( 'Body table size must be at least 1 but is #.' )
         CALL ERRINT ( '#', MAXBOD                                    )
         CALL SIGERR ( 'SPICE(INVALIDTABLESIZE)'                      )
         CALL CHKOUT ( 'ZZDSKSBR'                                     )
         RETURN

      END IF

C
C     AVAIL is the number of entries currently available.
C
      AVAIL = STSIZE - STFREE + 1

      IF ( AVAIL .LT. NEEDED ) THEN
C
C        We need to make room in the segment table.
C
C        The entries at the end of the body table have the highest
C        priority. We scan forward through this table, summing the
C        entry counts for each body, until we have enough entries. Let
C        NE represent the number of available entries. NE is initially
C        the number of unused entries.
C
         NSTDEL = 0
         I      = 1

         DO WHILE ( ( I .LE. BTNBOD ) .AND. ( AVAIL .LT. NEEDED ) )
C
C           Add the segment count for the Ith body to the total.
C
            NSTDEL = NSTDEL + BTSTSZ(I) 
            AVAIL  = AVAIL  + NSTDEL
            I      = I      + 1

         END DO

C
C        Backstop: we should always have enough room in the segment
C        table at this point.
C
         IF ( AVAIL .LT. NEEDED ) THEN

            CALL SETMSG ( 'The requested number of segment entries is '
     .      //            '#; the size STSIZE of the input segment  '
     .      //            'table is #. This error should have been '
     .      //            'trapped before this point.'                )
            CALL ERRINT ( '#', NEEDED                                 )
            CALL ERRINT ( '#', STSIZE                                 )
            CALL SIGERR ( 'SPICE(BUG)'                                )
            CALL CHKOUT ( 'ZZDSKSBR'                                  )
            RETURN

         END IF
         
C
C        At this point, I is the index of the first retained body,
C        unless all were deleted, in which case I is BTNBOD+1. We need
C        to delete the segment table entries of the bodies indexed from
C        1 to I-1.
C
         NBTDEL = I - 1

         IF ( NSTDEL .GT. 0 ) THEN
C
C           Adjust the tables to be consistent with the deletions.
C
C           Shift the body table and update the body table pointers.
C           
            DO I = NBTDEL+1, BTNBOD

               J         = I - NBTDEL
               BTBODY(J) = BTBODY(I)
               BTSTSZ(J) = BTSTSZ(I)
               BTSEGP(J) = BTSEGP(I) - NSTDEL
               
            END DO

C
C           Update the body table count.
C
            BTNBOD = BTNBOD - NBTDEL

C
C           Shift the segment table entries forward by NSTDEL to make
C           room at the rear of the table.
C
            DO I = NSTDEL+1, STFREE-1

               J        = I - NSTDEL
               STHAN(J) = STHAN(I)

               CALL MOVED( STDSCR(1,I), DSKDSZ, STDSCR(1,J) )
               CALL MOVEI( STDLAD(1,I), DLADSZ, STDLAD(1,J) )
               CALL MOVED( STOFF (1,I), 3,      STOFF (1,J) )
               CALL MOVED( STCTR (1,I), 3,      STCTR (1,J) )

               STRAD(J) = STRAD(I)

            END DO

            STFREE = STFREE - NSTDEL
   
         END IF

      END IF

      CALL CHKOUT ( 'ZZDSKSBR' )
      RETURN
      END
