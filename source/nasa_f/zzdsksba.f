C$Procedure ZZDSKSBA ( DSK, add entry to API segment buffer )
 
      SUBROUTINE ZZDSKSBA ( BODYID, 
     .                      MAXBOD, STSIZE, BTBODY, BTNBOD, BTSEGP, 
     .                      BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                      STOFF,  STCTR,  STRAD                  )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Add an entry for a specified body to the DSK API segment buffer
C     data structure.
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

      INTEGER               BODYID
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
C     BODYID     I   Body ID code.
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
C     BODYID     is the integer ID code of a body for which a new 
C                entry is to be made. 
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
C
C$ Detailed_Output
C
C     BTBODY     is the input body ID table, modified by this routine.
C                If BODYID was not in the table on input, it has been
C                appended to the table. 
C
C                If it was necessary to delete a body from the table
C                to make room for the body designated by BODYID, that
C                has been done.
C
C
C     BTNBOD     is the number of bodies in the body table. Depending
C                on what deletions may have been necessary, this 
C                number may be greater than, equal to, or less than
C                its value on input.
C
C
C     BTSEGP     is an array of start indices in the segment table
C                of existing entries, updated to reflect appending of
C                segments corresponding to BODYID.
C
C
C     BTSTSZ     is the array of sizes of segment lists of bodies.
C                BTSTSZ contains an entry for the body designated by
C                BODYID.
C
C                If it was necessary to make room, other entries of
C                BTSTSZ may have been deleted.
C
C
C     STHAN      is the segment handle array, updated to include
C                an entry for each loaded DSK segment associated 
C                with BODYID.
C
C                If it was necessary to make room, other entries of
C                STHAN may have been deleted.
C
C
C     STDSCR     is the segment DSK descriptor array, updated to include
C                an entry for each loaded DSK segment associated 
C                with BODYID. 
C
C                Segment entries are created in the order segments are
C                found by ZZDSKSNS, so the highest-priority segment's
C                entry is at the lowest index in STDSCR.
C
C                If it was necessary to make room, other entries of
C                STDESR may have been deleted.
C
C
C     STDLAD     is the segment DLA descriptor array, updated to include
C                an entry for each loaded DSK segment associated 
C                with BODYID.
C
C                If it was necessary to make room, other entries of
C                STDLAD may have been deleted.
C
C
C     STFREE     is the index of the first free element in each 
C                segment table array.
C
C
C     STOFF      is the segment frame center offset array, updated to
C                include an entry for each loaded DSK segment
C                associated with BODYID.
C
C                If it was necessary to make room, other entries of
C                STOFF may have been deleted.
C
C
C     STCTR      is the segment bounding sphere center array, updated
C                to include an entry for each loaded DSK segment
C                associated with BODYID.
C
C                If it was necessary to make room, other entries of
C                STCTR may have been deleted.
C
C
C     STRAD      is the segment bounding sphere radius array, updated
C                to include an entry for each loaded DSK segment
C                associated with BODYID.
C
C                If it was necessary to make room, other entries of
C                STRAD may have been deleted.
C
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an attempt is made to add an entry for a body already in
C         the body table, the error SPICE(INVALIDADD) is signaled.
C
C     2)  If a DSK segment has an unrecognized coordinate system, an
C         error will be signaled by a routine in the call tree
C         of this routine.
C
C     3)  If the center of a frame of a DSK segment cannot be 
C         obtained from the frame's ID code, the error
C         SPICE(NOFRAMEINFO) is signaled.
C
C     4)  If the name of a frame of a DSK segment cannot be 
C         obtained from the frame's ID code, the error
C         SPICE(NOFRAMENAME) is signaled.
C
C     5)  If an error occurs while this routine attempts to
C         obtain segment information from ZZDSKBSR, an error
C         will be signaled by a routine in the call tree
C         of this routine.
C
C     6)  If the offset of a frame's center from the frame's
C         body is not computable at the TDB epoch at the
C         midpoint of the segment's time coverage interval,
C         an error will be signaled by a routine in the 
C         SPK subsystem.
C
C     7)  If the segment table cannot accommodate all segments for
C         the specified body, the error SPICE(SEGMENTTABLEFULL) will
C         be signaled.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C     If any loaded DSK segment has a reference frame that is not
C     centered at the segment's central (target) body, SPK data are
C     required to compute the offset between the frame's center and
C     the segment's center.
C
C     Frame kernels may be required in order to look up a segment's
C     frame center offset. In some cases, additional kernels such
C     as CK kernels and SCLK kernels could be required to support
C     the offset vector lookup.
C  
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine maintains data structures used by SPICELIB 
C     geometry APIs to perform efficient computations using 
C     DSK data. The umbrella routine in which these structures
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
C     New entries are added to the tables by appending. 
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
C-    SPICELIB Version 1.0.0, 13-JUL-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     add entry to DSK API segment buffer
C
C-&

C
C     SPICELIB functions
C
      INTEGER               ISRCHI

      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               ZZDSKSBD

C
C     External routines
C
      EXTERNAL              ZZDSKBDC 

C
C     Local parameters
C
      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

C
C     Local variables
C
      CHARACTER*(FRNMLN)    FRNAME

      DOUBLE PRECISION      DSKDSC ( DSKDSZ )
      DOUBLE PRECISION      ET
      DOUBLE PRECISION      LT

      INTEGER               AVAIL
      INTEGER               DLADSC ( DLADSZ )
      INTEGER               FRMCTR
      INTEGER               HANDLE
      INTEGER               I
      INTEGER               J
      INTEGER               NSEG
      INTEGER               SEGCLD
      INTEGER               SEGCLS
      INTEGER               SEGCTR
      INTEGER               SEGFID

      LOGICAL               FRMFND
      LOGICAL               SEGFND
      LOGICAL               STATUS


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKSBA' )

C
C     Check the body table for presence of the new body. It's
C     an error to call this routine for a body that's already
C     present. (Such a call likely indicates the tables were
C     not re-initialized after a BSR state change.)
C
      I = ISRCHI( BODYID, BTNBOD, BTBODY )

      IF ( I .GT. 0 ) THEN

         CALL SETMSG ( 'Body # is already present in the DSK segment '
     .   //            'buffer body table. The table must be '
     .   //            're-initialized before this body can be added.' )
         CALL ERRINT ( '#', BODYID                                     )
         CALL SIGERR ( 'SPICE(INVALIDADD)'                             )
         CALL CHKOUT ( 'ZZDSKSBA'                                      )
         RETURN

      END IF
 
C
C     Make sure the BSR segment list for the body is up to 
C     date.
C     
      CALL ZZDSKBBL ( BODYID )

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZDSKSBA' )
         RETURN
      END IF

C
C     Count the segments in the BSR system for the body.
C
      NSEG   = 0
      STATUS = ZZDSKSBD( BODYID )
     
      CALL ZZDSKBSS( BODYID )
      CALL ZZDSKSNS( ZZDSKBDC, HANDLE, DLADSC, DSKDSC, SEGFND )

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZDSKSBA' )
         RETURN
      END IF

      DO WHILE ( SEGFND )

         NSEG = NSEG + 1

         CALL ZZDSKSNS( ZZDSKBDC, HANDLE, DLADSC, DSKDSC, SEGFND )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZDSKSBA' )
            RETURN
         END IF

      END DO

C
C     Check the number of segments for BODY against the size of the
C     segment table. If the table isn't big enough, we can't make 
C     room by deleting existing entries. This is a backstop check;
C     this situation should not occur if STSIZE is consistent with
C     the value in ZZDSKBSR.
C
      IF ( NSEG .GT. STSIZE ) THEN
         
         CALL SETMSG ( 'The number of segments for body # is #; '
     .   //            'the size STSIZE of the input segment table '
     .   //            'is #.'                                       )
         CALL ERRINT ( '#', BODYID                                   )
         CALL ERRINT ( '#', NSEG                                     )
         CALL ERRINT ( '#', STSIZE                                   )
         CALL SIGERR ( 'SPICE(SEGMENTTABLEFULL)'                     )
         CALL CHKOUT ( 'ZZDSKSBA'                                    )
         RETURN

      END IF

C
C     If we don't have enough room to store new entries in the body
C     table or in the segment table, make room.
C
      AVAIL = STSIZE - STFREE + 1

      IF (  ( BTNBOD .EQ. MAXBOD ) .OR. ( AVAIL .LT. NSEG )  ) THEN

         CALL ZZDSKSBR ( NSEG,
     .                   MAXBOD, STSIZE, BTBODY, BTNBOD, BTSEGP, 
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZDSKSBA' )
            RETURN
         END IF

      END IF

C
C     Append the new body ID to the body table. We've ensured there's
C     room in the table.
C
      BTNBOD           = BTNBOD + 1
      BTBODY( BTNBOD ) = BODYID
      BTSEGP( BTNBOD ) = STFREE
      BTSTSZ( BTNBOD ) = NSEG

C
C     Make a second pass through the BSR segment list, this time
C     accumulating segments in the input segment table as we go.
C
      STATUS = ZZDSKSBD( BODYID )

      CALL ZZDSKBSS( BODYID )
      CALL ZZDSKSNS( ZZDSKBDC, HANDLE, DLADSC, DSKDSC, SEGFND )

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZDSKSBA' )
         RETURN
      END IF

      DO WHILE ( SEGFND )
C
C        Insert handle and descriptor data for the current segment at
C        index STFREE in the segment table.
C
         STHAN( STFREE ) = HANDLE

         CALL MOVEI ( DLADSC, DLADSZ, STDLAD(1,STFREE) )
         CALL MOVED ( DSKDSC, DSKDSZ, STDSCR(1,STFREE) )

         STFREE          = STFREE + 1

         CALL ZZDSKSNS( ZZDSKBDC, HANDLE, DLADSC, DSKDSC, SEGFND )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZDSKSBA' )
            RETURN
         END IF

      END DO

C
C     Compute bounding spheres and frame center offsets for each
C     segment.
C
      DO I = 1, NSEG
C
C        J is the index in the segment table of the Ith segment
C        for BODYID.
C
         J  = BTSEGP( BTNBOD ) + I - 1

         CALL ZZSEGBOX ( STDSCR(1,J), STCTR(1,J), STRAD(J) )
         
         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZDSKSBA' )
            RETURN
         END IF

C
C        Obtain the center of the frame for the Ith segment.
C
         SEGFID = NINT( STDSCR(FRMIDX,J) )

         CALL FRINFO ( SEGFID, FRMCTR, SEGCLS, SEGCLD, FRMFND )
         
         IF ( .NOT. FRMFND ) THEN

            CALL SETMSG ( 'Could not look up frame info for '
     .      //            'segment frame having ID #.'       )
            CALL ERRINT ( '#', SEGFID                        )
            CALL SIGERR ( 'SPICE(NOFRAMEINFO)'               )
            CALL CHKOUT ( 'ZZDSKSBA'                         )
            RETURN

         END IF
C
C        If the frame center is not the same as the central
C        body, compute the offset between the two. Otherwise
C        set the offset to zero.
C        
         SEGCTR = NINT( STDSCR(CTRIDX,J) )

         IF ( SEGCTR .EQ. FRMCTR ) THEN

            CALL CLEARD ( 3, STOFF(1,J) )

         ELSE

            CALL FRMNAM ( SEGFID, FRNAME )

            IF ( FRNAME .EQ. ' ' ) THEN

               CALL SETMSG ( 'Could not look up frame info for '
     .         //            'segment frame having ID #.'       )
               CALL ERRINT ( '#', SEGFID                        )
               CALL SIGERR ( 'SPICE(NOFRAMENAME)'               )
               CALL CHKOUT ( 'ZZDSKSBA'                         )
               RETURN

            END IF

C
C           Note that SPK data must be available at the midpoint of the
C           DSK coverage epoch in order for the following call to work.
C
            ET = ( STDSCR(BTMIDX,J) + STDSCR(ETMIDX,J) ) / 2

            CALL SPKGPS ( FRMCTR, ET, FRNAME, SEGCTR, STOFF(1,J), LT )
           
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZDSKSBA' )
               RETURN
            END IF 

         END IF

      END DO


      CALL CHKOUT ( 'ZZDSKSBA' )
      RETURN
      END
