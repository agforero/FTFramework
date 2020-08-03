C$Procedure      ZZDSKBSR ( DSK, buffer segments for readers )
 
      SUBROUTINE ZZDSKBSR ( FNAME,
     .                      BODYID,
     .                      HANDLE,
     .                      CMPFUN,
     .                      USRCTR,
     .                      UPDATE,
     .                      DLADSC,
     .                      DSKDSC,
     .                      FOUND   )

      IMPLICIT NONE
 
C$ Abstract
C
C     Load and unload DSK files for use by the readers. Buffer segments
C     for readers.
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
C     DAS
C
C$ Keywords
C
C     TOPOGRAPHY
C
C$ Declarations

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'zzctr.inc'
 
      CHARACTER*(*)         FNAME
      INTEGER               BODYID
      INTEGER               HANDLE
      LOGICAL               CMPFUN
      EXTERNAL              CMPFUN
      INTEGER               USRCTR ( * )
      LOGICAL               UPDATE
      INTEGER               DLADSC ( * )
      DOUBLE PRECISION      DSKDSC ( * )
      LOGICAL               FOUND


      INTEGER               FTSIZE
      PARAMETER           ( FTSIZE =  5000 )
 
      INTEGER               BTSIZE
      PARAMETER           ( BTSIZE =  100 )
  
      INTEGER               LBPOOL
      PARAMETER           ( LBPOOL =   -5 )

      INTEGER               STSIZE
      PARAMETER           ( STSIZE = 10000 )
 
C$ Brief_I/O
C
C     Variable  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     FNAME      I   ZZDSKLSF
C     BODYID     I   ZZDSKBSS
C     HANDLE    I,O  ZZDSKLSF, ZZDSKUPF, ZZDSKSNS
C     USRCTR    I,O  ZZDSKCHK
C     UPDATE     O   ZZDSKCHK
C     DESCR      O   ZZDSKSNS
C     FOUND      O   ZZDSKSNS
C
C$ Detailed_Input
C
C     FNAME      is the name of a binary DSK file to be loaded.
C
C     HANDLE     on input is the handle of a binary DSK file to be
C                unloaded.
C
C
C     The purpose of entry points ZZDSKBSS and ZZDSKSNS is to search for
C     segments in DSK files matching certain criteria. 
C
C     USRCTR     on input is the value of a DSK loaded kernel counter
C                maintained by a DSK routine or subsystem.
C
C
C$ Detailed_Output
C
C     HANDLE     on output is the handle of the DSK file
C                containing a located segment.
C
C     USRCTR     on output is the value of a DSK loaded kernel counter
C                maintained by the ZZDSKBSR subsystem.
C
C     UPDATE     is a logical flag indicating whether the state of 
C                the loaded DSK set has changed since the ZZDSKBSR
C                state counter was equal to a given user counter.
C
C     DESCR      is the packed descriptor of a located segment.
C
C     FOUND      indicates whether a requested segment was found or not.
C
C$ Parameters
C
C     FTSIZE     is the maximum number of shape files that can
C                be loaded by ZZDSKLSF at any given time for use by the
C                readers.
C
C     BTSIZE     is the maximum number of bodies whose segments
C                are buffered by ZZDSKSNS.
C
C     STSIZE     is the maximum number of segments that can be buffered
C                at any given time by ZZDSKSNS.
C
C$ Exceptions
C
C     1) If ZZDSKBSR is called directly, the error SPICE(DSKBOGUSENTRY)
C        is signaled.
C
C     2) See entry points ZZDSKLSF, ZZDSKUPF, ZZDSKBSS, and ZZDSKSNS
C        for exceptions specific to them.
C
C$ Files
C
C     DSK shape files are indicated by filename before loading
C     (see ZZDSKLSF) and handle after loading (all other places).
C
C$ Particulars
C
C     ZZDSKBSR serves as an umbrella, allowing data to be shared by its
C     entry points:
C
C        ZZDSKLSF       Load shape file.
C        ZZDSKUPF       Unload shape file.
C        ZZDSKBSS       Begin search for segment.
C        ZZDSKSNS       Select next segment.
C        ZZDSKCHK       Check for change in loaded kernel set.
C
C     Before a file can be read by the DSK readers, it must be
C     loaded by ZZDSKLSF, which among other things load the file into
C     the DAS subsystem.
C
C     Up to FTSIZE files may be loaded for use simultaneously, and a 
C     file only has to be loaded once to become a potential search 
C     target for any number of subsequent reads.
C
C     Once a DSK has been loaded, it is assigned a file
C     handle, which is used to keep track of the file internally, and
C     which is used by the calling program to refer to the file in all
C     subsequent calls to DSK routines.
C
C     A file may be removed from the list of files for potential
C     searching by unloading it via a call to ZZDSKUPF.
C 
C     ZZDSKBSS and ZZDSKSNS are used together to search through loaded
C     files for segments.
C
C     ZZDSKBSS sets up the search.
C
C     ZZDSKSNS finds segments matching the search criteria set up by
C     ZZDSKBSS. Last-loaded files get searched first, and individual
C     files are searched backwards.
C
C     When an applicable segment is found, ZZDSKSNS returns that
C     segment's descriptor and identifier, along with the handle of the
C     file containing the segment.
C
C     Subsequent calls to ZZDSKSNS continue the search, picking up
C     where the previous call to this routine left off.
C
C     ZZDSKSNS uses information on loaded files to manage a buffer of
C     saved segment descriptors and identifiers. The buffer is used to
C     speed up access time by minimizing file reads.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     1) If Fortran I/O errors occur while searching a loaded DSK
C        file, the internal state of this suite of routines may
C        be corrupted.  It may be possible to correct the state
C        by unloading the pertinent DSK files and then re-loading
C        them.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C     R.E. Thurman   (JPL)
C     I.M. Underwood (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) 
C
C        Updated version info.
C
C        31-MAR-2016 (NJB) 
C
C           Deleted unused variable references. 
C           Added header for ZZDSKCHK.
C
C           Cleaned up in-line comments.
C
C        19-SEP-2014 (NJB) 
C
C           Specification change: now treats the target body as the
C           primary search key. Segments lists are associated with
C           bodies rather than surfaces.
C
C           The input argument SURFID has been replaced with BODYID.
C
C           A counter system has been implemented to enable SPICE
C           routines to detect changes in the state of the loaded
C           kernel set. The user interface entry point for checking for
C           state updates in ZZDSKCHK.
C
C        21-MAY-2013 (NJB) 
C
C           Removed debugging WRITE statement from entry point
C           ZZDSKSNS. Edited headers to remove long lines.
C
C        02-APR-2010 (NJB) (RET) (IMU)
C
C           Original version.
C
C-&
 
C$ Index_Entries
C
C     buffer dsk segments for readers
C
C-&
 
 
 
 
C$ Revisions
C
C     None.
C
C-&
 
C
C     SPICELIB functions
C
      INTEGER               INTMAX
      INTEGER               ISRCHI
      INTEGER               LNKNFN
      INTEGER               LNKNXT
      INTEGER               LNKPRV
 
      LOGICAL               FAILED
      LOGICAL               RETURN
  
C
C     Local parameters


      INTEGER               SLEN
      PARAMETER           ( SLEN   = 40 )
 
C
C     Constants used in the doubly linked list structure:
C     
      INTEGER               FORWRD
      PARAMETER           ( FORWRD  = 1 )
 
C
C     Local variables
C
 
C
C     DSK loaded kernel set state change counter:
C
      INTEGER               DSKCTR ( CTRSIZ )

C
C     The file table contains the handle and file number of each file
C     that has been loaded for use with the DSK readers. File
C     numbers begin at one, and are incremented until they reach a 
C     value of INTMAX() - 1, at which point they are mapped to the
C     range 1:NFT, where NFT is the number of loaded DSK files.
C
C     A file number is similar to a file handle, but it is assigned
C     and used exclusively by this module. The purpose of file numbers
C     is to keep track of the order in which files are loaded and the
C     order in which they are searched.
C
C     All names begin with FT.
C
C        HAN      Handle
C        NUM      File number
C
C     NFT is the number of currently loaded DSK files. NEXT is
C     incremented whenever a new file is loaded to give the file
C     number for that file. FINDEX is the index of whatever file is
C     of current interest.
C
C     New files are added at the end of the table. As files are
C     removed, succeeding files are moved forward to take up the
C     slack. This keeps the table ordered by file number.
C
      INTEGER               NFT
      INTEGER               FTHAN    ( FTSIZE )
      INTEGER               FTNUM    ( FTSIZE )
      INTEGER               NEXT
      INTEGER               FINDEX
 
C
C     The body table contains the beginning of the list of the
C     stored segments for each body and the
C     expense at which that list was constructed. (The expense of an
C     body list is the number of segment descriptors examined
C     during the construction of the list.) It also contains the
C     highest and lowest file numbers searched during the construction
C     of the list.
C
C     All names begin with BT.
C
C        INS      Body number
C        EXP      Expense
C        HFS      Highest file (number) searched
C        LFS      Lowest  file (number) searched
C        BEG      Beginning of segment list
C
C     NIT is the number of bodies for which segments are currently
C     being stored in the table. BINDEX is the index of whatever
C     body is of current interest at any given time.
C
C     New bodies are added at the end of the table. As bodies
C     are removed, the last body is moved forward to take up the
C     slack. This keeps the entries in the table contiguous.
C
      INTEGER               BINDEX
      INTEGER               BTBEG    ( BTSIZE )
      INTEGER               BTEXP    ( BTSIZE )
      INTEGER               BTHFS    ( BTSIZE )
      INTEGER               BTBOD    ( BTSIZE )
      INTEGER               BTLFS    ( BTSIZE )
      INTEGER               NBT
 
C
C     The segment table contains the handle, descriptor, and identifier
C     for each segment that has been found so far.
C
C     The segment table is implemented as a set of arrays indexed by
C     a SPICE doubly linked list structure.  For each body
C     in the body table, there is a segment table list; each 
C     node of a list points to data associated with a segment.  In 
C     each list, the head node corresponds to the highest-priority 
C     segment in that list, and segment priority decreases in the 
C     forward direction.
C
C     All names begin with ST.
C
C        DLAD     DLA segment descriptor
C        DSKD     DSK segment descriptor
C        HAN      Handle
C        POOL     Doubly linked list pool.
C
C     New segments are added to the front or end of an body list
C     as appropriate, according to the rules spelled out under
C     entry point ZZDSKSNS.
C

      DOUBLE PRECISION      STDSKD  ( DSKDSZ,     STSIZE )

      INTEGER               STHAN   (             STSIZE )
      INTEGER               STDLAD  ( DLADSZ,     STSIZE )
      INTEGER               STPOOL ( 2,  LBPOOL : STSIZE )
 
C
C     Other local variables
C 
      CHARACTER*(SLEN)      DOING
      CHARACTER*(SLEN)      STACK    ( 2 )
      CHARACTER*(SLEN)      STATUS
      CHARACTER*(SLEN)      URGENT

      DOUBLE PRECISION      DSKLDS ( DSKDSZ )

  
      INTEGER               CHEAP
      INTEGER               COST
      INTEGER               DLALDS ( DLADSZ )
      INTEGER               DLANXT ( DLADSZ )
      INTEGER               DLAPRV ( DLADSZ )
      INTEGER               HEAD
      INTEGER               I
      INTEGER               J
      INTEGER               MINEXP
      INTEGER               NEW
      INTEGER               NODE
      INTEGER               NXTSEG
      INTEGER               P
      INTEGER               PRVNOD
      INTEGER               SAVEP
      INTEGER               SVBODY
      INTEGER               TAIL
      INTEGER               TOP
      LOGICAL               BEGSCH
      LOGICAL               FND
      LOGICAL               PASS1

C
C     Saved variables
C
      SAVE                  BEGSCH
      SAVE                  FINDEX
      SAVE                  FND
      SAVE                  FTHAN
      SAVE                  FTNUM
      SAVE                  BINDEX
      SAVE                  BTBEG
      SAVE                  BTEXP
      SAVE                  BTHFS
      SAVE                  BTBOD
      SAVE                  BTLFS
      SAVE                  DSKCTR
      SAVE                  NEXT
      SAVE                  NFT
      SAVE                  NBT
      SAVE                  PASS1
      SAVE                  SAVEP
      SAVE                  SVBODY
      SAVE                  STATUS
      SAVE                  STDLAD
      SAVE                  STDSKD
      SAVE                  STHAN
      SAVE                  STPOOL
      SAVE                  TOP
 
C
C     Initial values
C
      DATA                  NFT       / 0             /
      DATA                  NBT       / 0             /
      DATA                  NEXT      / 0             /
      DATA                  PASS1     / .TRUE.        /
      DATA                  SAVEP     / 0             /
      DATA                  STATUS    / 'BOGUS ENTRY' /
 
 
C
C     Nobody has any business calling ZZDSKBSR directly.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZDSKBSR' )
      CALL SIGERR ( 'SPICE(DSKBOGUSENTRY)' )
      CALL CHKOUT ( 'ZZDSKBSR' )
 
      RETURN
 
 
C$Procedure ZZDSKLSF ( DSK, load shape file )
 
      ENTRY ZZDSKLSF ( FNAME, HANDLE )
 
C$ Abstract
C
C     Load a DSK shape file for use by the DSK readers.  Return that
C     file's handle, to be used by other DSK routines to refer to the
C     file.
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
C     DAS
C
C$ Keywords
C
C     TOPOGRAPHY
C
C$ Declarations
C
C     CHARACTER*(*)         FNAME
C     INTEGER               HANDLE
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     FNAME      I   Name of the DSK file to be loaded.
C     HANDLE     O   Loaded file's handle.
C     FTSIZE     P   Maximum number of loaded DSK files.
C
C$ Detailed_Input
C
C     FNAME      is the name of a DSK file to be loaded.
C
C$ Detailed_Output
C
C     HANDLE     is an integer handle assigned to the file upon loading.
C                Almost every other DSK routine will subsequently use
C                this number to refer to the file.
C
C$ Parameters
C
C     FTSIZE     is the maximum number of DSK files that may 
C                be loaded simultaneously under any circumstances.
C                FTSIZE is currently set to match the maximum number
C                of DAS files that may be loaded simultaneously.
C
C$ Exceptions
C
C     1) If an attempt is made to open more DAS files than is specified
C        by the parameter FTSIZE in DASAH, an error is signaled by a
C        routine in the call tree of this routine.
C
C     2) If an attempt is made to load more files than is specified
C        by the local parameter FTSIZE, and if the DAS system has 
C        room to load another file, the error SPICE(DSKTOOMANYFILES)
C        signaled. The current setting of FTSIZE does not allow this
C        situation to arise:  the DAS system will trap the error 
C        before this routine has the chance.
C
C     3) If the file specified by FNAME can not be opened, an error
C        is signaled by a routine that this routine calls.
C
C     4) If the file specified by FNAME has already been loaded,
C        it will become the "last-loaded" file.  The readers
C        search the last-loaded file first.
C
C$ Files
C
C     The DSK file specified by FNAME is loaded. The file is assigned
C     an integer handle by ZZDSKLSF. Other DSK routines will refer to
C     this file by its handle.
C
C$ Particulars
C
C     See Particulars above, in ZZDSKBSR.
C
C     If there is room for a new file, ZZDSKLSF opens the file for
C     reading.
C
C     DSK readers search files loaded with ZZDSKLSF in the reverse order
C     in which they were loaded.  That is, last-loaded files are
C     searched first.
C
C$ Examples
C
C     See the Example above, in ZZDSKBSR.
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
C     N.J. Bachman   (JPL)
C     R.E. Thurman   (JPL)
C     I.M. Underwood (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) 
C
C        Updated version info.
C
C        19-SEP-2014 (NJB) 
C
C           Specification change: now treats the target body as the
C           primary search key. Segments lists are associated with
C           bodies rather than surfaces.
C
C           Now initializes the DSK state change counter on the first
C           pass and updates it on every subsequent call.

C       21-MAY-2013 (NJB) 
C
C           Edited headers to remove long lines.
C
C       02-APR-2010 (NJB) (RET) (IMU)
C
C-&
 
C$ Index_Entries
C
C     load dsk shape file
C
C-&
 
 
C$ Revisions
C
C     None.
C-&
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKLSF' )

      IF ( PASS1 ) THEN
C
C        Initialize the BSR counter.
C
         CALL ZZCTRSIN ( DSKCTR )

         PASS1 = .FALSE.

      END IF

C
C     Increment the BSR counter regardless of whether 
C     the load operation is successful.
C
      CALL ZZCTRINC ( DSKCTR )
    
C
C     Don't allow a search to continue after loading a file; a new
C     search should be re-started.
C
      STATUS = 'BOGUS ENTRY'
 
C
C     Nothing works unless at least one file has been loaded, so
C     this is as good a place as any to initialize the free list
C     whenever the body table is empty.
C
      IF ( NBT .EQ. 0 ) THEN
         CALL LNKINI ( STSIZE, STPOOL )
      END IF
  
C
C     To load a new file, first try to open it for reading.
C
      CALL DASOPR ( FNAME, HANDLE )
 
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZDSKLSF' )
         RETURN
      END IF
 
C
C     Determine if the file is already in the table.
C
      FINDEX = ISRCHI ( HANDLE, NFT, FTHAN )
  
      IF ( FINDEX .GT. 0 ) THEN
C
C        The last call we made to DASOPR added another DAS link to
C        the DSK file.  Remove this link.
C
         CALL DASCLS ( HANDLE )

C
C        Handle is already in the table.  Remove it.
C
         NFT = NFT - 1
 
         DO I = FINDEX, NFT
            FTHAN( I ) = FTHAN( I + 1 )
            FTNUM( I ) = FTNUM( I + 1 )
         END DO
 
C
C        Unlink any segments that came from this file.
C
         I = 1
 
         DO WHILE ( I .LE. NBT )
 
            P = BTBEG( I )
 
            DO WHILE ( P .GT. 0 )
C
C              Find the successor of P, if any.
C 
               NXTSEG = LNKNXT ( P, STPOOL )
 
               IF ( STHAN( P ) .EQ. HANDLE ) THEN
C
C                 The segment corresponding to node P came from
C                 the file we're unloading.  Delete the node for
C                 P from the segment list for body I; if P happens
C                 to be the head node for body I's segment list,
C                 make the successor of P the head of the list.
C
                  CALL LNKFSL ( P, P, STPOOL )

                  IF ( P .EQ. BTBEG(I) ) THEN
                     BTBEG( I ) = NXTSEG 
                  END IF
 
               END IF 
C
C              Update P.
C 
               P = NXTSEG
 
            END DO
 
C
C           If the list for this body is now empty, shorten the
C           current table by one: put all the entries for the last
C           body in the table into the space occupied by the
C           one we've deleted.
C
            IF ( BTBEG( I ) .LE. 0 ) THEN

               BTBOD( I ) = BTBOD( NBT )
               BTEXP( I ) = BTEXP( NBT )
               BTHFS( I ) = BTHFS( NBT )
               BTLFS( I ) = BTLFS( NBT )
               BTBEG( I ) = BTBEG( NBT )
 
               NBT = NBT - 1
 
            ELSE
 
               I = I + 1
 
            END IF
 
         END DO
 
      ELSE
C
C        This is a new file.  Make sure that there are unused slots
C        in the file table.
C
         IF ( NFT .EQ. FTSIZE ) THEN
 
            CALL DASCLS ( HANDLE )
 
            CALL SETMSG ( 'Number of files loaded is at a maximum, ' //
     .                    'as specified by the parameter FTSIZE, '   //
     .                    'the value of which is #. You will need '  //
     .                    'to load fewer files. Consider unloading ' //
     .                    'any files that are not needed.'            )
            CALL ERRINT ( '#', FTSIZE                                 )
            CALL SIGERR ( 'SPICE(DSKTOOMANYFILES)'                    )
            CALL CHKOUT ( 'ZZDSKLSF'                                  )
            RETURN
 
         END IF
 
      END IF
 
C
C     Determine the next file number.
C
C     Programmer's note: this section is normally not reached.
C     It should be tested by temporarily setting the comparison
C     value to a smaller number, for example 2*FTSIZE.
C     
      IF ( NEXT .LT. INTMAX()-1 ) THEN

         NEXT = NEXT + 1

      ELSE
C
C        The user is to be congratulated:  we've run out of file 
C        numbers.
C
C        Re-set the valid file numbers so they lie in the range 1:NFT, 
C        with the Ith file in the file table having file number I.
C        First update the LFS and HFS components of the body table
C        according to this mapping.  
C  
C        Set any body table entries that are lower than FTNUM(1) 
C        to zero.  
C
         DO I = 1, NBT 
C
C           Re-map the HFS table for the Ith body.
C
            J = ISRCHI ( BTHFS(I), NFT, FTNUM )

C
C           Either the highest file searched for body I is the Jth 
C           file in the file table, or the file is not in the table.
C           In both cases, J is the correct value to assign to 
C           BTHFS(I).
C
            BTHFS(I) = J

C
C           When the highest file searched for body I is not in the
C           file table, the highest file searched has been unloaded.
C           Note that this assignment makes all files appear to be
C           "new" when a lookup for body I is performed.
C
C           Re-map the LFS table for the Ith body.
C
            J = ISRCHI ( BTLFS(I), NFT, FTNUM )

            IF ( J .GT. 0 ) THEN
C
C              The lowest file searched for body I is the Jth file
C              in the file table.
C
               BTLFS(I) = J

            ELSE
C
C              The lowest file searched for body I is not in the 
C              file table.  This occurs when the lowest file searched 
C              has been unloaded.  Zero out both the lowest and
C              highest file searched to force reconstruction of the 
C              list.
C
               BTLFS(I) = 0
               BTHFS(I) = 0

            END IF

         END DO

C
C        Re-map the file number table itself.
C
         DO I = 1, NFT
 
            FTNUM(I) = I

         END DO

C
C        Assign a new file number.
C
         NEXT = NFT + 1

      END IF

C
C     Now add this file to file table.
C
      NFT        = NFT  + 1
      FTHAN(NFT) = HANDLE
      FTNUM(NFT) = NEXT

      CALL CHKOUT ( 'ZZDSKLSF' )
      RETURN
 
 


C$Procedure ZZDSKUSF ( DSK, Unload shape file )
 
      ENTRY ZZDSKUSF ( HANDLE )
 
C$ Abstract
C
C     Unload a DSK shape file so that it will no longer be searched
C     by the readers.
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
C     DAS
C
C$ Keywords
C
C     TOPOGRAPHY
C
C$ Declarations
C
C     INTEGER               HANDLE
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK file to be unloaded
C
C$ Detailed_Input
C
C     HANDLE     Integer handle assigned to the file upon loading.
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
C     1) Unloading a file that has not been loaded is a no-op.
C        No error is signaled.
C
C$ Files
C
C     The file referred to by HANDLE is unloaded.
C
C$ Particulars
C
C     See Particulars section above, in ZZDSKBSR.
C
C     Unloading a file with ZZDSKUSF removes that file from
C     consideration by the DSK readers. In doing so, it frees up space
C     for another file to be loaded.
C
C$ Examples
C
C     See the Example above, in ZZDSKBSR.
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
C     N.J. Bachman   (JPL)
C     R.E. Thurman   (JPL)
C     I.M. Underwood (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) 
C
C        Updated version info.
C
C        19-SEP-2014 (NJB) 
C
C           Specification change: now treats the target body as the
C           primary search key. Segments lists are associated with
C           bodies rather than surfaces.
C
C           Now initializes the DSK state change counter on the first
C           pass and updates it on every subsequent call.

C        21-MAY-2013 (NJB) 
C
C           Edited headers to remove long lines.
C           Removed debugging output.
C
C        02-APR-2010 (NJB) (RET) (IMU)
C
C-&
 
C$ Index_Entries
C
C     unload ck shape file
C
C-&
 
 
 
C$ Revisions
C
C     None.
C
C-&
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKUSF' )

      IF ( PASS1 ) THEN
C
C        Initialize the BSR counter.
C
         CALL ZZCTRSIN ( DSKCTR )

         PASS1 = .FALSE.

      END IF

C
C     Increment the BSR counter regardless of whether 
C     the load operation is successful.
C
      CALL ZZCTRINC ( DSKCTR )


C
C     Don't allow a search to continue after unloading a file; a new
C     search should be re-started.
C
      STATUS = 'BOGUS ENTRY'

C
C     All of the stored segments from the file must be removed 
C     from the segment table (by returning the corresponding nodes 
C     to the segment table pool.)
C
C     Don't do anything if the given handle is not in the file table.
C
      FINDEX = ISRCHI ( HANDLE, NFT, FTHAN )

      IF ( FINDEX .EQ. 0 ) THEN
         CALL CHKOUT ( 'ZZDSKUSF' )
         RETURN
      END IF
C
C
C     First get rid of the entry in the file table. Close the file
C     before wiping out the handle.
C
      CALL DASCLS ( FTHAN(FINDEX) )

 
      NFT = NFT - 1

      DO I = FINDEX, NFT
         FTHAN( I ) = FTHAN( I + 1 )
         FTNUM( I ) = FTNUM( I + 1 )
      END DO
 
C
C     Check each body list individually. Note that the first
C     node on each list, having no predecessor, must be handled
C     specially.
C
      I = 1
 
      DO WHILE ( I .LE. NBT )
 
         P = BTBEG( I )
 
         DO WHILE ( P .GT. 0 )
 
            NXTSEG = LNKNXT ( P, STPOOL )

            IF ( STHAN(P) .EQ. HANDLE ) THEN
       
               IF ( P .EQ. BTBEG(I) ) THEN
                  BTBEG( I ) = NXTSEG
               END IF
C
C              Free this segment table entry.
C            
               CALL LNKFSL ( P, P, STPOOL )

            END IF

            P = NXTSEG
 
         END DO
 
C
C        If the list for this body is now empty, shorten the
C        current table by one: put all the entries for the last
C        body in the table into the space occupied by the
C        one we've deleted.
C
         IF ( BTBEG(I) .LE. 0 ) THEN

            IF ( I .NE. NBT ) THEN

               BTBOD (I) = BTBOD (NBT)
               BTEXP (I) = BTEXP (NBT)
               BTHFS (I) = BTHFS (NBT)
               BTLFS (I) = BTLFS (NBT)
               BTBEG (I) = BTBEG (NBT)

            END IF

            NBT = NBT - 1
                        
         ELSE

            I = I + 1

         END IF

      END DO


      CALL CHKOUT ( 'ZZDSKUSF' )
      RETURN

 
 
 
C$Procedure ZZDSKBSS ( DSK, begin search for segment )
 
      ENTRY ZZDSKBSS ( BODYID )
 
C$ Abstract
C
C     Initiate search through loaded files to find segments
C     satisfying search criteria.
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
C     DAS
C
C$ Keywords
C
C     TOPOGRAPHY
C
C$ Declarations
C
C     INTEGER               BODYID
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C
C$ Detailed_Input
C
C     ZZDSKBSS sets up a search for segments. The four quantities below
C     establish the search criteria.
C
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
C     1)  If no files have been loaded, the error
C         SPICE(NOLOADEDDSKFILES) is signaled.
C
C$ Files
C
C     All files loaded by ZZDSKLSF are potential search targets for
C     ZZDSKSNS.
C
C$ Particulars
C
C     ZZDSKBSS sets up a search for segments by ZZDSKSNS. It records the
C     body and time to be searched for, and whether to require
C     segments containing angular velocity data. If angular velocity
C     data are required, only segments containing angular velocity
C     data will be returned by ZZDSKSNS. If angular velocity data are
C     not required, segments returned by ZZDSKSNS may or may not contain
C     angular velocity data.
C
C     ZZDSKBSS determines the first task that ZZDSKSNS will have to
C     perform if it is called to get an applicable segment.
C
C$ Examples
C
C     See Examples in ZZDSKBSR.
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
C     N.J. Bachman   (JPL)
C     R.E. Thurman   (JPL)
C     I.M. Underwood (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) 
C
C        Updated version info.
C
C        19-SEP-2014 (NJB) 
C
C           Specification change: now treats the target body as the
C           primary search key. Segments lists are associated with
C           bodies rather than surfaces.
C
C           Now initializes the DSK state change counter on the first
C           pass.
C
C        21-MAY-2013 (NJB) 
C
C           Edited headers to remove long lines.
C
C        02-APR-2010 (NJB) (RET) (IMU)
C
C
C-&
 
C$ Index_Entries
C
C     begin search for dsk segment
C
C-&
 
 
 
C$ Revisions
C
C     None.
C
C-&
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKBSS' )


      IF ( PASS1 ) THEN
C
C        Initialize the BSR counter.
C
         CALL ZZCTRSIN ( DSKCTR )

         PASS1 = .FALSE.

      END IF


C
C     Make a saved copy of the body ID code.  
C
      SVBODY = BODYID
 
C
C     There must be at least one file loaded.
C
      IF ( NFT .EQ. 0 ) THEN
 
         CALL SETMSG ( 'At least one DSK file needs must be loaded ' //
     .                 'by ZZDSKLSF before beginning a search.'       )
         CALL SIGERR ( 'SPICE(NOLOADEDDSKFILES)'                      )
         CALL CHKOUT ( 'ZZDSKBSS'                                     )
         RETURN
 
      END IF
 
C
C     The stack of suspended tasks is empty.
C
      TOP = 0

C
C     Is the body already in the body table?  The answer
C     determines what the first task for ZZDSKSNS will be.
C
      BINDEX = ISRCHI ( SVBODY, NBT, BTBOD )

 
      IF ( BINDEX .EQ. 0 ) THEN

         STATUS = 'NEW BODY'

      ELSE
C
C        Set the status so that ZZDSKSNS will determine whether to check
C        the segment list or search new files.
C
         STATUS = '?'

      END IF

C
C     The saved segment list pointer is no longer valid.
C
      SAVEP  = -1
 
      CALL CHKOUT ( 'ZZDSKBSS' )
      RETURN
 
 
 
C$Procedure ZZDSKSNS ( DSK, Select next segment )
 
      ENTRY ZZDSKSNS ( CMPFUN, HANDLE, DLADSC, DSKDSC, FOUND )
 
C$ Abstract
C
C     Search through loaded files to find a segment matching the
C     requested body, time, and need for angular velocity.
C     Buffer segment descriptors, identifiers, and handles in the
C     process to minimize file reads.
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
C     DAS
C
C$ Keywords
C
C     TOPOGRAPHY
C
C$ Declarations
C
C     LOGICAL               CMPFUN
C     EXTERNAL              CMPFUN
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C     LOGICAL               FOUND
C
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     CMPFUN     I   Comparison function.
C     HANDLE     O   Handle of file containing the applicable segment.
C     DLSDSC     O   DLA descriptor of the applicable segment.
C     DSKDSC     O   DSK descriptor of the applicable segment.
C     FOUND      O   True if a segment was found.
C
C$ Detailed_Input
C
C     CMPFUN
C
C$ Detailed_Output
C
C     DLADSC
C
C     DSKDS
C
C
C     HANDLE     is an integer handle of the file containing the
C                segment matching the body and time
C                specifications made in the last call to ZZDSKBSS.
C
C     FOUND      is true if an applicable segment was found.  False
C                otherwise.  If FOUND is false, the values of the
C                other arguments are meaningless.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If ZZDSKSNS is called without ZZDSKBSS ever having been
C        called, the error 'SPICE(CALLZZDSKBSSFIRST)' is signaled.
C
C     2) If no segment is found that matches the search criteria,
C        FOUND is set to false, but the values of HANDLE, DESCR,
C        and SEGID will be meaningless.
C
C$ Files
C
C     All files loaded by ZZDSKLSF are potential search targets for
C     ZZDSKSNS. The files are all referred to by their integer handles.
C
C$ Particulars
C
C     ZZDSKSNS is used to locate segments based on the search criteria
C     established by the most recent call to ZZDSKBSS.  When a segment
C     is found it will have the following characteristics:
C
C        1) Its body will match the body specified in the
C           call to ZZDSKBSS.
C
C        2) Its time interval will intersect the time interval
C
C              [SCLKDP - TOL, SCLKDP + TOL],
C
C           where SCLKDP and TOL were specified in the call to ZZDSKBSS.
C
C        3) If there is a need for angular velocity data, as specified
C           by NEEDAV in the call to ZZDSKBSS, a returned segment
C           will contain angular velocity data. If there is no need
C           for such data, the returned segment may or may not contain
C           angular velocity data.
C
C     The first call to ZZDSKSNS following a call to ZZDSKBSS starts a
C     search through loaded files and either returns the first
C     applicable segment, or indicates that no segment was found.
C
C     ZZDSKSNS searches through last-loaded files first. Individual
C     files are searched backwards, so that segments that were inserted
C     last into the file get checked first.
C
C     Subsequent calls to ZZDSKSNS pick up the search exactly where the
C     previous calls left off. If a segment is not found, future calls
C     will also indicate that no segment could be found, until a new
C     search is begun.
C
C     ZZDSKSNS also buffers segment descriptors and identifiers, to
C     attempt to minimize file reads.
C
C$ Examples
C
C     See Examples in ZZDSKBSR.
C
C$ Restrictions
C
C     1) This subroutine assumes that a search has been initiated by
C        a call to ZZDSKBSS.
C
C     2) When a DSK file is loaded or unloaded, a new search must 
C        be started via a call to ZZDSKBSS before this routine may
C        be called.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C     R.E. Thurman   (JPL)
C     I.M. Underwood (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) 
C
C        Updated version info.
C
C        31-MAR-2016 (NJB) 
C
C           References to unneeded variables were removed.
C
C        19-SEP-2014 (NJB) 
C
C           Specification change: now treats the target body as the
C           primary search key. Segment lists are associated with
C           bodies rather than surfaces.
C
C           Now initializes the DSK state change counter on the first
C           pass.
C
C        21-MAY-2013 (NJB) 
C
C            Removed debugging WRITE statement. Edited headers to
C            remove long lines.
C
C        02-APR-2010 (NJB) (RET) (IMU)
C
C-&
 
C$ Index_Entries
C
C     select next dsk segment
C
C-&
 
 
C$ Revisions
C
C     None.
C
C-&
 
 
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKSNS' )



      IF ( PASS1 ) THEN
C
C        Initialize the BSR counter.
C
         CALL ZZCTRSIN ( DSKCTR )

         PASS1 = .FALSE.

      END IF

C
C     Nothing's been found yet.
C
      FOUND = .FALSE.

C
C     Initialize the segment list pointer to the saved value from
C     the previous pass through this routine, if any.
C
      P = SAVEP

C
C     ZZDSKSNS buffers segment descriptors and identifiers, to
C     attempt to minimize file reads. Buffering segments involves
C     maintaining three tables:  the file table, the body table,
C     and the segment table. ZZDSKSNS is broken down into various tasks,
C     described in the code below, which perform these manipulations.
C
C     A description of the components of each table is provided in
C     the declarations section of ZZDSKBSR.
C
C     Basically, the buffering is performed as follows: once a request
C     for a segment for a particular body is made, if there are
C     no adequate entries in the buffer already, a search is made
C     through loaded files for applicable segments.  Every segment
C     pertaining to that body in a searched file is buffered,
C     before a check of the current buffer is made.  If the search
C     doesn't turn up a segment matching the specified search criteria
C     the next file is searched and new segments are added to the list,
C     and so on.
C
C     The information in the segment table (ST) is stored in a
C     doubly-linked list. Each node in the list contains several
C     individual pieces of data, which are stored in parallel
C     arrays.
C
C     In the following loop, we will try to simplify things by
C     doing exactly one thing on each pass through the loop.
C     After each pass, the status of the loop (STATUS) will be
C     adjusted to reflect the next thing that needs to be done.
C     The first task is set by ZZDSKBSS.
C
C     Occasionally, the current task will have to be interrupted
C     until another task can be carried out. (For example, when
C     collecting new segments, an interrupt might place a segment
C     at the front or end of the current body list; when placing
C     the segment on the list, a second interrupt might free
C     room in the segment table in order to allow the addition
C     to proceed.) In this case, the current task will be saved and
C     restored after the more urgent task has been completed.
C
C     The loop can terminate in only one of two ways (unless an error
C     occurs). First, if an applicable segment is found in the segment 
C     table, the handle, descriptor, and identifier for the segment 
C     are returned immediately.  Second, if the table does not contain 
C     an applicable segment, and if no files remain to be searched, 
C     the loop terminates normally, and no data are returned.
C
C     The status is saved on exit, however, so that subsequent calls
C     will resume a search exactly where previous calls left off.
C
C     Each status is described below.
C
C     'NEW BODY'
C
C        This indicates that the specified body has
C        no segments stored for it at all. It must be added to the
C        body table.  (This is followed immediately by an
C        OLD FILES search, in which every file loaded is considered an
C        old file.)
C
C     'NEW FILES'
C
C        This indicates that at least one new file has been added
C        since the last time the segment list for the specified
C        body was searched. Find the oldest of these new files,
C        and begin a NEW SEGMENTS search in forward order for
C        segments to add to the front of the list.
C
C     'NEW SEGMENTS'
C
C        Continue a NEW FILES search, adding segments for the specified
C        body to the front of the list.
C
C     'OLD FILES'
C
C        This indicates that although the list has been searched
C        and found to contain no applicable segment, some of the
C        older files remain to be searched. Find the newest of these
C        old files, and begin an OLD SEGMENTS search in backward order.
C
C     'OLD SEGMENTS'
C
C        Continue an OLD FILES search, adding segments for the specified
C        body to the end of the list.
C
C     'CHECK LIST'
C
C        This indicates that the list is ready to be searched,
C        either because no new files have been added, or because
C        segments from a new file or an old file have recently
C        been added.
C
C        The list is never checked until all new files have been
C        searched.
C
C        If an applicable segment is found, it is returned.
C
C     'MAKE ROOM' (Interrupt)
C
C        This indicates that one of the bodies must be removed,
C        along with its stored segments, to make room for another
C        body or segment.  The body (other than the
C        specified body) with the smallest expense is selected
C        for this honor.
C
C     'ADD TO FRONT' (Interrupt)
C
C        This indicates that a segment has been found (during the
C        course of a NEW FILES search) and must be added to the front
C        of the list.
C
C     'ADD TO END' (Interrupt)
C
C        This indicates that a segment has been found (during the
C        course of an OLD FILES search) and must be added to the end
C        of the list.
C
C     'SUSPEND'
C
C        This indicates that the current task (DOING) should be
C        interrupted until a more urgent task (URGENT) can be
C        carried out. The current task is placed on a stack for
C        safekeeping.
C
C     'RESUME'
C
C        This indicates that the most recently interrupted task
C        should be resumed immediately.
C
C     '?'
C
C        This indicates that the next task is not immediately
C        apparent: if new files exist, they should be searched;
C        otherwise the list should be checked.
C
C     'HOPELESS'
C
C        This indicates that the table does not contain an applicable
C        segment, and no files remain to be searched.
C
C      'BOGUS ENTRY'
C
C        This is the initial value of STATUS and indicates that no
C        call to ZZDSKBSS was ever made. If this is the case then an
C        error will be signaled.
C
      
      IF ( STATUS .EQ. 'BOGUS ENTRY' )  THEN
 
         CALL SETMSG ( 'Must begin a search by calling ZZDSKBSS '
     .   //            'first.'                                   )
         CALL SIGERR ( 'SPICE(CALLZZDSKBSSFIRST)'                 )
         CALL CHKOUT ( 'ZZDSKSNS'                                 )
         RETURN
 
      END IF

      
      DO WHILE ( STATUS .NE. 'HOPELESS' )
C
C        If new files have been added, they have to be searched.
C        Otherwise, go right to the list of stored segments.
C
         IF ( STATUS .EQ. '?' ) THEN
C
C           There are two ways to get to this point.
C
C           1)  Status may have been set to '?' by ZZDSKBSS.
C
C           2)  Status was set to '?' by the NEW SEGMENTS block
C               of code as the result of finishing the read of
C               a new file.
C
 
            IF ( BTHFS( BINDEX ) .LT. FTNUM( NFT ) ) THEN
 
               STATUS = 'NEW FILES'

            ELSE

C              If the segment list for this body is empty, make 
C              sure the expense is set to 0.
C 
               IF ( BTBEG(BINDEX) .LE. 0 ) THEN
                  BTEXP(BINDEX) = 0
               END IF

C
C              Prepare to look at the first segment in the list for
C              this body.
C 
               P      =  BTBEG( BINDEX )
               STATUS = 'CHECK LIST'
 
            END IF
 

         ELSE IF ( STATUS .EQ. 'NEW BODY' ) THEN
C
C           New bodies are added to the end of the body 
C           table. If the table is full, one of the current occupants 
C           must be removed to make room for the new one.
C
C           Setting LFS to one more than the highest current file 
C           number means the 'OLD FILES' search that follows will
C           begin with the last-loaded file.
C
C           There is one way to get here:
C
C           1)  The variable STATUS was set to NEW BODY prior 
C               in ZZDSKBSS.
C
C           Find the cheapest slot in the body table to store
C           the initial information about this body.
C
C           NOTE:  This used to be handled by the MAKE ROOM section.
C           However, trying to handle this special case there was
C           just more trouble than it was worth.
C
            IF ( NBT .LT. BTSIZE ) THEN
C
C              If the body table isn't full, the cheapest place is
C              just the next unused row of the table.
C
               NBT   = NBT + 1
               CHEAP = NBT
 
            ELSE
C
C              The body table is full.  Find the least
C              expensive body in the table and remove it.
C
               CHEAP  = 1
               MINEXP = BTEXP(1)
 
               DO I = 2, NBT
 
                  IF ( BTEXP(I) .LT. MINEXP ) THEN
                     CHEAP  = I
                     MINEXP = BTEXP(I)
                  END IF
 
               END DO

C
C              If there are any segments associated with the
C              least expensive body, we put them back on the free
C              list.
C
               HEAD = BTBEG(CHEAP)

               IF ( HEAD .GT. 0 ) THEN

                  TAIL =  - LNKPRV ( HEAD, STPOOL )
                  CALL LNKFSL ( HEAD, TAIL, STPOOL )

               END IF

            END IF

C
C           Set up a table entry for the new body. 
C 
            BTBOD (CHEAP) = SVBODY
            BTEXP (CHEAP) = 0
            BTHFS (CHEAP) = FTNUM(NFT)
            BTLFS (CHEAP) = FTNUM(NFT) + 1
            BTBEG (CHEAP) = 0
            BINDEX        = CHEAP

C
C           Now search all of the files for segments relating to
C           this body.
C
            STATUS = 'OLD FILES' 


         ELSE IF ( STATUS .EQ. 'NEW FILES' ) THEN
C
C           When new files exist, they should be searched in forward
C           order, beginning with the oldest new file not yet searched.
C           All new files must be searched before the list can be
C           checked, to ensure that the best (newest) segments are
C           being used.
C
C           Begin a forward search, and prepare to look for individual
C           segments from the file.
C
C           The only way to get here is to have STATUS set to
C           the value NEW FILES in the STATUS .EQ. '?' block
C           of the IF structure.
C
C           Find the next file to search; set FINDEX to the
C           corresponding file table entry.
 
            FINDEX = 1
 
            DO WHILE ( BTHFS(BINDEX) .GE. FTNUM(FINDEX) )
 
               FINDEX = FINDEX + 1
 
            END DO
 
            BTHFS( BINDEX ) = FTNUM( FINDEX )
 
C
C           Start a forward search through the current file.
C
            BEGSCH = .TRUE.
            CALL DLABFS (  FTHAN( FINDEX ),  DLALDS,  FND  )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZDSKSNS' )
               RETURN
            END IF
            
            STATUS = 'NEW SEGMENTS'
 
C
C           The cost of the list contributed by the new file is
C           zero so far.
C
            COST   = 0


         ELSE IF ( STATUS .EQ. 'NEW SEGMENTS' ) THEN
C
C           New files are searched in forward order. Segments, when
C           found, are inserted at the front of the list. 
C
C           Each segment examined, whether applicable or not, adds to
C           the expense of the list.
C
C           The only ways to get here are:
C
C               1) Enter from the NEW FILES block of the IF structure.
C               2) Re-enter from this block if there are more segments
C                  to examine in the current file and the last segment
C                  seen wasn't for the body of interest.
C               3) Enter from the RESUME state after adding a segment
C                  to the front of the list for the current body.
C
            IF ( BEGSCH ) THEN
C
C              We already have a FND value, and if FND is true, a
C              DLA descriptor.
C
               BEGSCH = .FALSE.

            ELSE
C
C              Use the current DLA descriptor to look up the next one.
C
               CALL DLAFNS ( FTHAN(FINDEX), DLALDS, DLANXT, FND )

               IF ( FND ) THEN
                  CALL MOVEI  ( DLANXT, DLADSZ, DLALDS )
               END IF

            END IF

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZDSKSNS' )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN
C
C              We're out of segments in the current file.  Decide
C              whether we need to examine another new file, or 
C              whether we're ready to check the list.
C
               STATUS = '?'
               BTEXP( BINDEX ) = BTEXP( BINDEX ) + COST
 
            ELSE
C
C              Get the DSK segment descriptor for the current
C              segment.
C
               CALL DSKGD( FTHAN(FINDEX), DLALDS, DSKLDS )

               IF (  NINT( DSKLDS(CTRIDX) ) .EQ. SVBODY ) THEN
C
C                 The segment is for the body of interest. Add this
C                 segment to the front of the list.
C
                  DOING  = 'NEW SEGMENTS'
                  URGENT = 'ADD TO FRONT'
                  STATUS = 'SUSPEND'
 
               END IF
 
               COST = COST + 1
 
            END IF
C
C           If we haven't reset the status, we'll return for another
C           'NEW SEGMENTS' pass.
C           

         ELSE IF ( STATUS .EQ. 'OLD FILES' ) THEN
C
C           When old files must be searched (because the segments in 
C           the list are inadequate), they should be searched in
C           backward order, beginning with the newest old file not 
C           yet searched.  The segment list will be re-checked 
C           after each file is searched.  If a match is found, 
C           the search terminates, so some old files may not be 
C           searched.
C
C           Begin a backwards search, and prepare to look for 
C           individual segments from the file.
C
C           You can get to this block in two ways.
C
C           1) We can have a NEW BODY.
C
C           2) We have checked the current list (CHECK LIST) for
C              this body, didn't find an applicable segment and
C              have some files left that have not been searched.
 
            FINDEX = NFT
 
            DO WHILE ( BTLFS( BINDEX ) .LE. FTNUM( FINDEX ) )
               FINDEX = FINDEX - 1
            END DO
 
            BEGSCH = .TRUE.
            CALL DLABBS (  FTHAN( FINDEX ),  DLALDS,  FND )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZDSKSNS' )
               RETURN
            END IF

            STATUS = 'OLD SEGMENTS'
 
C
C           The next thing we'll do is search through all the segments
C           of this file for those that applicable to this body.
C           The cost of the list contributed by the current file is
C           zero so far.
C
            COST = 0
  
C
C        Old files are searched in backward order. Segments, when
C        found, are inserted at the end of the list.  
C
C        Each segment examined, whether applicable or not, adds to
C        the expense of the list.
C
         ELSE IF ( STATUS .EQ. 'OLD SEGMENTS' ) THEN
C
C           There is only one way to get here---from the
C           block 'OLD FILES'.  Note we do not add to the
C           expense of the list for this body until we've
C           completely searched this file.
C 
            IF ( BEGSCH ) THEN
C
C              We already have a value of FND, and if FND is true,
C              a DLA segment from the current file.
C
               BEGSCH = .FALSE.

            ELSE
C
C              Look up the previous segment from this file.
C
               CALL DLAFPS ( FTHAN(FINDEX), DLALDS, DLAPRV, FND )

               IF ( FND ) THEN

                  CALL MOVEI ( DLAPRV, DLADSZ, DLALDS )

               END IF

            END IF

 
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZDSKSNS' )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN
C
C              All of the segments in this file have been exhausted.
C              Change the lowest file searched indicator for this
C              body to be the current file, and go check the
C              current list.
C
               BTLFS( BINDEX ) =  FTNUM( FINDEX )
               BTEXP( BINDEX ) =  BTEXP( BINDEX ) + COST
               P               =  BTBEG( BINDEX )
               STATUS          = 'CHECK LIST'
 
            ELSE
C
C              Get the DSK descriptor for this segment.
C
               CALL DSKGD ( FTHAN(FINDEX), DLALDS, DSKLDS )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZDSKSNS' )
                  RETURN
               END IF

               IF (  NINT( DSKLDS(CTRIDX) )  .EQ.  SVBODY  ) THEN
C
C                 This is a segment for the body of interest.
C
                  DOING  = 'OLD SEGMENTS'
                  URGENT = 'ADD TO END'
                  STATUS = 'SUSPEND'
 
               END IF
 
               COST = COST + 1
 
            END IF
 
         ELSE IF ( STATUS .EQ. 'CHECK LIST' ) THEN
C
C           Okay, all the new files (and maybe an old file or two)
C           have been searched. Time to look at the list of segments
C           stored for the body, to see if there is one applicable
C           to the specified epoch and need for angular velocity data.
C
C           If so, return it.  If not, try another old file.  If there
C           are no more old files, give up the ghost.
C
C           There are two ways to get to this point.
C
C           1) From the '?' block.
C           2) From the 'OLD SEGMENTS' block.
C
C           For every segment examined, adjust the re-use interval
C           associated with the current body.
C
C           P always points to the current segment in the list. Reject
C           a segment if there is a need for angular velocity data and
C           the segment doesn't have it.
C
C
            DO WHILE ( P .GT. 0 )
              
               IF ( CMPFUN( STHAN(P), STDLAD(1,P), STDSKD(1,P) )  ) THEN

                  CALL MOVEI ( STDLAD(1,P), DLADSZ, DLADSC )

                  CALL MOVED ( STDSKD(1,P), DSKDSZ, DSKDSC )

                  HANDLE =  STHAN ( P )
                  FOUND  = .TRUE.
C
C                 Go ahead and move the pointer up before returning
C                 so that the search for the next applicable segment
C                 will start at the right place.
C
                  SAVEP = STPOOL ( FORWRD, P )

                  CALL CHKOUT ( 'ZZDSKSNS' )
                  RETURN
 
               END IF


C
C              Get the next node.  We avoid LNKNXT here in order
C              to speed up the operation.
C 
               P = STPOOL ( FORWRD, P )
 
            END DO

C
C           If we're still here we didn't have information for this
C           body in the segment list.
C
C           If there are more files, search them.
C           Otherwise, things are hopeless, set the status that way.
C 
            IF ( BTLFS( BINDEX ) .GT. FTNUM( 1 ) ) THEN
               STATUS = 'OLD FILES'
            ELSE
               STATUS = 'HOPELESS'
            END IF


         ELSE IF ( STATUS .EQ. 'MAKE ROOM' ) THEN
C
C           When adding a new segment to a full table, one of the 
C           current bodies must be dropped.  The ideal 
C           candidate is the one whose list was constructed at the 
C           lowest expense.  The candidate should be removed from 
C           the body table, and its list transferred to the 
C           segment table pool. 
C
C           There is ``room'' if the segment table pool contains at 
C           least one free node.
C
C           It is possible that a single body requires more than the
C           entire segment table for its own segments. Two things might
C           happen in such a case:
C
C              1) If the list under consideration was being added to at
C                 the end, then a search is continued without buffering
C                 any segments.
C
C              2) If the list was being added to at the beginning, then
C                 that means there was a NEW FILES search going on, and
C                 so a brand new list is constructed for the body,
C                 much as in a 'NEW BODY' task.
C
C           There are two different ways to get to this point.
C
C              1) From 'ADD TO FRONT' if the segment table pool is full.
C              2) From 'ADD TO END' if the segment table pool is full.
C
C           Try to make room by deleting a segment list.  CHEAP will 
C           be the index of the "cheapest" segment list in the 
C           body table.
C
            MINEXP = INTMAX()
            CHEAP  = 0
 

      
            DO I = 1, NBT
 
               IF ( I .NE. BINDEX ) THEN

                  IF (      ( BTEXP(I) .LT. MINEXP   ) 
     .                 .OR. ( CHEAP    .EQ. 0        )  )THEN
C
C                    This list is the cheapest seen so far,
C                    possibly because it's the first one 
C                    considered.  At the moment, it's as good
C                    a candidate for removal as any.
C 
                     CHEAP  = I
                     MINEXP = BTEXP(I)

                  END IF

               END IF
 
            END DO
 

            IF ( CHEAP .EQ. 0 ) THEN
C
C              If there are no deletable segments, the Thing To
C              Do depends on the task that was suspended before
C              entering MAKE ROOM.
C
               IF ( STACK(TOP) .EQ. 'ADD TO END' ) THEN
C   
C                 The segment meta-data from the current file cannot
C                 be buffered.  
C
C                 In the corresponding SPK and CK cases, we would
C                 search the partial list of segments from this file,
C                 then proceed to search the rest of the file and any
C                 other old files. In this case, we don't support
C                 searching unbuffered segments, so this is the 
C                 end of the line. 
C
C                 We must clean up the segment list for the current
C                 body before we signal an error. All segments from the
C                 file we're currently searching must be deleted from
C                 the list. If we delete the head node of the list, the
C                 body table pointer to the list must be updated. If
C                 the segment list becomes empty, the body must be
C                 deleted from the body table.
C
                  HEAD =   BTBEG ( BINDEX )
                  TAIL = - LNKPRV( HEAD, STPOOL )

                  NODE =   TAIL

                  DO WHILE ( NODE .GT. 0 ) 
C
C                    Let PRVNOD be the predecessor of NODE. PRVNOD may 
C                    be negative (actually, a pointer to the list tail).
C
                     PRVNOD = LNKPRV( NODE, STPOOL )

                     IF ( STHAN(NODE) .EQ. FTHAN(FINDEX) ) THEN
C
C                       This segment is from the file we were 
C                       searching when we ran out of room. Free
C                       the segment list entry at index NODE.

                        CALL LNKFSL ( NODE, NODE, STPOOL )

                        IF ( NODE .EQ. HEAD ) THEN
C
C                          We just deleted the last remaining node in
C                          the list for the current body. We can delete
C                          this body from the body table. However, the
C                          body table contains no other bodies at this
C                          point, since we would have deleted them in
C                          the attempt to make room for the current
C                          body. So we don't need to compress the body
C                          table; we just indicate that it's empty.
C
                           NBT = 0

                        END IF
C
C                       This is the end of the block that handles the
C                       head node case.
C
                     END IF
C
C                    This is the end of the block that handles the 
C                    matching file case.
C
C                    Process the previous node. If the node is
C                    non-positive, the loop will terminate.
C
                     NODE = PRVNOD

                  END DO
C
C                 The segment table entries for the current body that
C                 are associated with the current file have been
C                 deleted.
C
C                 Make sure that a new search is started before this
C                 routine is called again.
C
                  STATUS = 'HOPELESS'
                  TOP    = 0
C
C                 It's finally time to signal the error.
C
                  CALL SETMSG ( 'ZZDSKSNS ran out of segment table '
     .            //            'room while trying to append to the '
     .            //            'tail of the segment list for body '
     .            //            '#. Current state is ADD TO END.'    )
                  CALL ERRINT ( '#', SVBODY                          )
                  CALL SIGERR ( 'SPICE(BUFFEROVERFLOW)'              )
                  CALL CHKOUT ( 'ZZDSKSNS'                           )
                  RETURN 

               ELSE
C
C                 STACK(TOP) is set to 'ADD TO FRONT'.  
C
C                 If there is no room left in the table in the middle
C                 of an attempt to add to the front of the list, just
C                 start from scratch by effectively initiating a 'NEW
C                 BODY' task.
C
C                 Return the current list to the segment table pool.
C                 Note this list is non-empty.
C
                  P    =   BTBEG ( BINDEX )
                  TAIL = - LNKPRV( P, STPOOL )

                  CALL LNKFSL ( P, TAIL, STPOOL )
C
C                 Re-initialize the table for this body, and
C                 initiate an 'OLD FILES' search, just as in 'NEW
C                 BODY'.
C
                  BTEXP( BINDEX ) = 0
                  BTHFS( BINDEX ) = FTNUM( NFT )
                  BTLFS( BINDEX ) = FTNUM( NFT ) + 1

                  STATUS = 'OLD FILES'

               END IF
 
C
C              Unwind the stack; we've set the target states already.
C
               TOP = 0
 
            ELSE 
C
C              Return this cheapest list to the segment pool.  This
C              list could be empty.
C
               HEAD = BTBEG( CHEAP )

               IF ( HEAD .GT. 0 ) THEN

                  TAIL = - LNKPRV ( HEAD, STPOOL )

                  CALL LNKFSL ( HEAD, TAIL, STPOOL )

               END IF

C
C              Fill the deleted body's space in the table with
C              the final entry in the table.
C
               IF ( CHEAP .NE. NBT ) THEN

                  BTBOD (CHEAP) = BTBOD (NBT)
                  BTEXP (CHEAP) = BTEXP (NBT)
                  BTHFS (CHEAP) = BTHFS (NBT)
                  BTLFS (CHEAP) = BTLFS (NBT)
                  BTBEG (CHEAP) = BTBEG (NBT)

               END IF

               IF ( BINDEX .EQ. NBT ) THEN
                  BINDEX = CHEAP
               END IF

C
C              One less body now.
C
               NBT    = NBT - 1
               STATUS = 'RESUME'
 
            END IF
C
C           At this point, we either made room by freeing a non-empty
C           segment list, or we're going to re-build the list for the
C           current body, starting with the highest-priority segments.
C           In the former case, the state is 'RESUME'; in the latter,
C           it's 'OLD FILES'.
C
 

         ELSE IF ( STATUS .EQ. 'ADD TO FRONT' ) THEN
C
C           The current segment information should be linked in at
C           the head of the segment list for the current body, 
C           and the pertinent body table entry should point 
C           to the new head of the list.
C
C           The only way to get here is from the block NEW SEGMENTS
C           after suspending that task.
 
            IF ( LNKNFN(STPOOL) .EQ. 0 ) THEN

               DOING  = 'ADD TO FRONT'
               URGENT = 'MAKE ROOM'
               STATUS = 'SUSPEND'
 
            ELSE
C
C              Allocate a node and link it to the front of the list
C              for the current body.
C
               CALL LNKAN ( STPOOL, NEW )

               STHAN( NEW ) = FTHAN( FINDEX )               
C
C              Store the DLA and DSK descriptors for this segment in
C              the segment table.
C
               CALL MOVEI ( DLALDS, DLADSZ, STDLAD(1,NEW) )
               CALL MOVED ( DSKLDS, DSKDSZ, STDSKD(1,NEW) ) 
 
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZDSKSNS' )
                  RETURN
               END IF

C
C              If the current list is empty, this append operation
C              is a no-op.
C 
               CALL LNKILB ( NEW, BTBEG(BINDEX), STPOOL )
               BTBEG( BINDEX ) = NEW
 
               STATUS  = 'RESUME'

            END IF
 

         ELSE IF ( STATUS .EQ. 'ADD TO END' ) THEN
C
C           The current segment information should be linked in at
C           the tail of the segment list for the current body.
C
C           The only way to get to this task is from the OLD SEGMENTS
C           block after suspending that task.
C 
            IF ( LNKNFN(STPOOL) .EQ. 0 ) THEN
 
               DOING  = 'ADD TO END'
               URGENT = 'MAKE ROOM'
               STATUS = 'SUSPEND'
 
            ELSE
C
C              Allocate a new node in the segment table pool.
C
               CALL LNKAN ( STPOOL, NEW )
 
               STHAN( NEW ) = FTHAN( FINDEX )
 
C
C              Store the DLA and DSK descriptors for this segment in
C              the segment table.
C
               CALL MOVEI ( DLALDS, DLADSZ, STDLAD(1,NEW) )
               CALL MOVED ( DSKLDS, DSKDSZ, STDSKD(1,NEW) ) 
 
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'ZZDSKSNS' )
                  RETURN
               END IF
  
               IF ( BTBEG(BINDEX) .LE. 0 ) THEN
C
C                 This is the first node in the list for this 
C                 body.
C
                  BTBEG( BINDEX ) = NEW
 
               ELSE
C
C                 Link the new node to the tail of the list.
C
                  TAIL = - LNKPRV ( BTBEG(BINDEX), STPOOL )
                  CALL LNKILA ( TAIL, NEW, STPOOL )

               END IF
  
               STATUS = 'RESUME'
 
            END IF
 

         ELSE IF ( STATUS .EQ. 'SUSPEND' ) THEN
C
C           When a task is suspended, the current activity is placed on
C           a stack, to be restored later. Two levels are provided, 
C           since some interrupts can be interrupted by others.
C
            TOP          = TOP + 1
            STACK( TOP ) = DOING
            STATUS       = URGENT
 
         ELSE IF ( STATUS .EQ. 'RESUME' ) THEN
 
            STATUS = STACK( TOP )
            TOP    = TOP - 1

         END IF
 
 
      END DO
 
C
C     Can only get here if status is 'HOPELESS', in which case a
C     segment was not found.
C
      FOUND = .FALSE.


      CALL CHKOUT ( 'ZZDSKSNS' )
      RETURN




C$Procedure ZZDSKCHK ( DSK, check for file updates )
 
      ENTRY ZZDSKCHK ( USRCTR, UPDATE )
 
C$ Abstract
C
C     Indicate to a calling routine whether any DSK load or unload
C     operations have occurred since the local state counter had
C     a particular, caller-supplied value. 
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
C     DAS
C
C$ Keywords
C
C     TOPOGRAPHY
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
C     USRCTR    is, on input, the value of a counter maintained by 
C               a calling routine. This counter is to be compared to
C               a local counter. 
C
C$ Detailed_Output
C
C     USRCTR    is, on output, the value of the local DSK load/unload
C               counter.
C
C     UPDATE    is a logical flag that is .TRUE. if and only if the
C               local DSK load/unload counter differs from the input
C               value of USRCTR. If the counters differ, a DSK load or
C               unload call has been made since the local counter last
C               had the value provided in USRCTR on input.
C
C$ Parameters
C
C     See zzctr.inc.
C
C$ Exceptions
C
C     1)  This routine will fail if the total count of DSK load and
C         unload operations exceeds the maximum accommodated by a
C         two-integer counter. The maximum count is about 4e19,
C         presuming a two-integer counter.
C   
C$ Files
C
C     This routine does not access DSK files, but all DSK file
C     load and unload operations affect the local counter used
C     by this routine.
C
C$ Particulars
C
C     This set of routines maintains a local counter that indicates
C     the kernel load status. When a call is made to ZZDSKLSF or
C     ZZDSKUSF, the local counter is incremented. Applications can
C     compare counter values they save to the current count by
C     calling this routine. A mismatch indicates that DSK load or
C     unload calls have been made since an application's counter
C     was set.
C
C     This type of check presumes the caller's counter is initialized
C     by a call to
C    
C        ZZCTRUIN
C
C     and is updated only by calls to this routine.
C
C     This routine is used by much of the DSK subsystem to determine
C     whether locally buffered information about DSK segments remains
C     up to date.
C
C$ Examples
C
C     See usage in ZZDSKBBL.
C
C$ Restrictions
C
C     This routine assumes the number of DSK load and unload operations
C     does not exceed the maximum value of a two-integer counter.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C     B.V. Semenov   (JPL)
C     R.E. Thurman   (JPL)
C     I.M. Underwood (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB) 
C
C        Updated version info.
C
C-&
 
C$ Index_Entries
C
C     check for changes to set of loaded dsk files
C
C-&
 
 
C$ Revisions
C
C     None.
C
C-&
 

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKCHK' )

      CALL ZZCTRCHK ( DSKCTR, USRCTR, UPDATE )

      CALL CHKOUT ( 'ZZDSKCHK' )
      RETURN
      END


