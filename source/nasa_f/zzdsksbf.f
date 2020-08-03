C$Procedure ZZDSKSBF ( DSK, manage the API segment buffer )
 
      SUBROUTINE ZZDSKSBF ( BODYID, NSURF,  SRFLST, ET,  FIXFID,  
     .                      VERTEX, RAYDIR, POINT,  XPT, HANDLE,
     .                      DLADSC, DSKDSC, DC,     IC,  FOUND, 
     .                      NORMAL                              )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Manage the DSK API segment buffer data structure.
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
      INCLUDE 'zzctr.inc'

      INTEGER               BODYID
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      DOUBLE PRECISION      ET
      INTEGER               FIXFID
      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      POINT  ( 3 )
      DOUBLE PRECISION      XPT    ( 3 )
      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      DOUBLE PRECISION      DSKDSC ( * )
      DOUBLE PRECISION      DC     ( * )
      INTEGER               IC     ( * )
      LOGICAL               FOUND
      DOUBLE PRECISION      NORMAL ( 3 )
 
C$ Brief_I/O
C
C     Variable  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     BODYID     I   ZZSBFXR, ZZSBFNRM
C     NSURF      I   ZZSBFXR, ZZSBFNRM
C     SRFLST     I   ZZSBFXR, ZZSBFNRM
C     ET         I   ZZSBFXR, ZZSBFNRM
C     FIXFID     I   ZZSBFXR, ZZSBFNRM
C     VERTEX     I   ZZSBFXR
C     RAYDIR     I   ZZSBFXR
C     POINT      I   ZZSBFNRM
C     XPT        O   ZZSBFXR
C     FOUND      O   ZZSBFXR
C     HANDLE     O   ZZSBFXRI
C     DLADSC     O   ZZSBFXRI
C     DSKDSC     O   ZZSBFXRI
C     DC         O   ZZSBFXRI
C     IC         O   ZZSRFXRI
C     NORMAL     O   ZZSBFNRM
C
C$ Detailed_Input
C
C     See the entry points.
C
C$ Detailed_Output
C
C     See the entry points. 
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If this routine is called directly, it signals the error 
C         SPICE(BOGUSENTRY).
C        
C     See the entry points for descriptions of errors specific to
C     those routines.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C
C     Data structure management
C     =========================
C
C     This routine manages data structures used by SPICELIB geometry
C     APIs to perform efficient computations using DSK data. This is
C     the umbrella routine in which these structures are declared.
C
C     This routine also contains entry points that use these data
C     structures to perform geometric computations using DSK data.
C     See the section titled "Computations" below.
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
C     tables must be re-initialized. Entry points of ZZDSKSBF 
C     perform this action by calling ZZDSKSBI.
C
C     In addition to the tables maintained by this subsystem, a 
C     counter is maintained. This counter is used to determine
C     whether the locally buffered data are in sync with the 
C     information stored in the ZZDSKBSR subsystem.
C     
C
C     Computations
C     ============
C
C     This routine contains the following entry points that use
C     buffered segment data to perform computations:
C
C        ZZSBFXR:    prepare for and compute unprioritized 
C                    ray-surface intercept using DSK data.
C
C        ZZSBFXRI:   prepare for and compute unprioritized 
C                    ray-surface intercept using DSK data; return data
C                    source information such as file handle, DLA
C                    descriptor, and plate ID as well.
C
C        ZZSBFNRM:   prepare for and compute outward normal 
C                    vector at a specified surface point, 
C                    using unprioritized DSK data.
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
C-    SPICELIB Version 1.0.0, 22-FEB-2017 (NJB) 
C
C        Added FAILED calls in each entry point.
C
C        17-MAY-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     manage DSK API segment buffer
C
C-&

C
C     SPICELIB functions
C
      INTEGER               ISRCHI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C     

C
C     Sizes of source info arrays returned by 
C     ZZDSKBUX:
C
      INTEGER               MAXDC
      PARAMETER           ( MAXDC  = 1 )

      INTEGER               MAXIC
      PARAMETER           ( MAXIC  = 1 )


      INTEGER               STSIZE
      PARAMETER           ( STSIZE = 10000 )

      INTEGER               BTSIZE
      PARAMETER           ( BTSIZE = 10 )

C
C     Body table variables
C     --------------------
C
C        BTNBOD  is the number of bodies in the body table.
C
C        BTBODY  is an array of body ID codes.
C
C        BTSEGP  is an array of pointers (start indices) to entries in
C                the segment table. The Ith pointer indicates the start
C                index for entries for the Ith body.
C
C        BTSTSZ  is an array of segment table entry counts. The Ith
C                element of BTSTSZ is the number of entries in the
C                segment table for the Ith body.
C     
C
      INTEGER               BTNBOD
      INTEGER               BTBODY ( BTSIZE )
      INTEGER               BTSEGP ( BTSIZE )
      INTEGER               BTSTSZ ( BTSIZE )
      
C
C     Segment table variables
C
      INTEGER               STHAN  ( STSIZE )
      DOUBLE PRECISION      STDSCR ( DSKDSZ, STSIZE )
      INTEGER               STDLAD ( DLADSZ, STSIZE )
      INTEGER               STFREE
      DOUBLE PRECISION      STOFF  ( 3,      STSIZE )
      DOUBLE PRECISION      STCTR  ( 3,      STSIZE )
      DOUBLE PRECISION      STRAD  ( STSIZE )


C
C     Local variables
C     
      DOUBLE PRECISION      LOCDC  ( MAXDC )

      INTEGER               BIX
      INTEGER               BSRCTR ( CTRSIZ )
      INTEGER               J
      INTEGER               LOCIC  ( MAXIC )
      INTEGER               NSEG
      INTEGER               SEGIDX

      LOGICAL               FIRST
      LOGICAL               UPDATE

C
C     Saved variables
C
      SAVE                  BSRCTR
      SAVE                  BTBODY
      SAVE                  BTNBOD
      SAVE                  BTSEGP
      SAVE                  BTSTSZ
      SAVE                  FIRST
      SAVE                  STCTR
      SAVE                  STDLAD
      SAVE                  STDSCR
      SAVE                  STFREE
      SAVE                  STHAN
      SAVE                  STOFF
      SAVE                  STRAD

C
C     Initial values
C
      DATA                  BTNBOD / 0      /
      DATA                  BTBODY / BTSIZE * 0 /
      DATA                  BSRCTR / CTRSIZ * 0 /
      DATA                  FIRST  / .TRUE. /
      DATA                  STFREE / 1      /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZDSKSBF'          )
      CALL SIGERR ( 'SPICE(BOGUSENTRY)' )
      CALL CHKOUT ( 'ZZDSKSBF'          )
      RETURN




C$Procedure ZZSBFXR ( DSK, prepare and perform unprioritized intercept )
 
      ENTRY ZZSBFXR ( BODYID, NSURF,  SRFLST, ET,    
     .                FIXFID, VERTEX, RAYDIR, XPT, FOUND )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Prepare and execute unprioritized ray-surface intercept
C     computation using DSK API segment buffers.
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
C
C     INTEGER               BODYID
C     INTEGER               NSURF
C     INTEGER               SRFLST ( * )
C     DOUBLE PRECISION      ET
C     INTEGER               FIXFID
C     DOUBLE PRECISION      VERTEX ( 3 )
C     DOUBLE PRECISION      RAYDIR ( 3 )
C     DOUBLE PRECISION      XPT    ( 3 )
C     LOGICAL               FOUND
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body. 
C     NSURF      I   Number of surface IDs in list.
C     SRFLST     I   Surface ID list.
C     ET         I   Evaluation epoch, seconds past J2000 TDB.
C     FIXFID     I   Frame ID of ray and intercept.
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     XPT        O   Surface intercept point.
C     FOUND      O   Found flag. True if intercept exists.
C
C$ Detailed_Input
C
C     BODYID     is the ID code of a target body. The ray-surface
C                intercept computation is performed using data
C                that represent the surface of this body.
C
C     NSURF,
C     SRFLST     are, respectively, a count of surface IDs and
C                a list of IDs. If the list is non-empty, only
C                the indicated surfaces will be used. If the 
C                list is empty, all surfaces associated with the
C                input body ID will be considered.
C     
C     ET         is the epoch for which the computation is to be
C                performed. This epoch is used for DSK segment
C                selection; only segments containing ET in their time
C                coverage interval will be used. ET is expressed as
C                seconds past J2000 TDB.
C
C     FIXFID     is the frame ID of a body-fixed frame associated
C                with the body designated by BODYID. This frame
C                is used to express the input ray's vertex and
C                direction vector. The output intercept will be
C                expressed in this frame as well.
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector
C                of a ray. Both vectors are expressed in the frame
C                designated by FIXFID.
C
C$ Detailed_Output
C
C     XPT        is the surface intercept on the target body 
C                nearest to the ray's vertex, if the intercept
C                exists.
C
C     FOUND      is a logical flag that is set to .TRUE. if and
C                only if an intercept exists.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If this routine is called directly, it signals the error 
C         SPICE(BOGUSENTRY).
C        
C     See the entry points for descriptions of errors specific to
C     those routines.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine prepares the local buffers for a ray-surface
C     intercept computation using unprioritized DSK data. 
C     It calls ZZDSKBUX to perform the computation.
C
C$ Examples
C
C     See usage in ZZRAYSFX.
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
C-    SPICELIB Version 1.0.0, 22-FEB-2017 (NJB) 
C
C        Added FAILED calls.
C
C        12-MAY-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     prepare unprioritized ray surface intercept 
C
C-&


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZSBFXR' )


      IF ( FIRST ) THEN
C
C        Initialize BSR counter.
C
         CALL ZZCTRUIN ( BSRCTR )
         
         FIRST = .FALSE.

      END IF

C
C     See whether the state of the loaded DSK set has changed
C     since the last call.
C
      CALL ZZDSKCHK ( BSRCTR, UPDATE )

      IF ( UPDATE ) THEN
C
C        Make sure the ZZDSKBSR subsystem has completed the segment
C        list for the input body since the last time the BSR loaded
C        kernel state changed.
C     
         CALL ZZDSKBBL ( BODYID )

C
C        Initialize the local buffers. We restart from scratch
C        each time the BSR loaded kernel state changes.
C
         CALL ZZDSKSBI ( BTSIZE, STSIZE, BTBODY, BTNBOD, BTSEGP,
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )
      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSBFXR' )
         RETURN
      END IF
      
C
C     Find the index of the input body ID in the body table. If
C     we re-initialized the tables, the index will be zero.
C
      BIX = ISRCHI( BODYID, BTNBOD, BTBODY )

       
      IF ( BIX .EQ. 0 ) THEN
C
C        We don't have buffered information for this body. Update
C        the body and segment tables to store data for it.
C         
         CALL ZZDSKSBA ( BODYID,
     .                   BTSIZE, STSIZE, BTBODY, BTNBOD, BTSEGP, 
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSBFXR' )
            RETURN
         END IF
C
C        The new body's position in the body table is at the end.
C         
         BIX = BTNBOD
         
      END IF

C
C     Find the ray-surface intercept, using the buffered segment
C     data.
C     
      J    = BTSEGP(BIX) 
      NSEG = BTSTSZ(BIX)    

      CALL ZZDSKBUX ( BODYID,     NSURF,      SRFLST, ET,  FIXFID,
     .                NSEG,       STHAN(J),   STDLAD(1,J), STDSCR(1,J),
     .                STOFF(1,J), STCTR(1,J), STRAD(J),    VERTEX,      
     .                RAYDIR,     XPT,        SEGIDX,      LOCDC,
     .                LOCIC,      FOUND                               )

      CALL CHKOUT ( 'ZZSBFXR' )
      RETURN      





C$Procedure ZZSBFXRI ( DSK, unprioritized intercept with info )
 
      ENTRY ZZSBFXRI ( BODYID, NSURF,  SRFLST, ET,    
     .                 FIXFID, VERTEX, RAYDIR, XPT, HANDLE,
     .                 DLADSC, DSKDSC, DC,     IC,  FOUND  )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Prepare and execute unprioritized ray-surface intercept
C     computation using DSK API segment buffers. Return source
C     information including handle, DLA descriptor, and
C     segment-specific information such as plate ID.
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
C
C     INTEGER               BODYID
C     INTEGER               NSURF
C     INTEGER               SRFLST ( * )
C     DOUBLE PRECISION      ET
C     INTEGER               FIXFID
C     DOUBLE PRECISION      VERTEX ( 3 )
C     DOUBLE PRECISION      RAYDIR ( 3 )
C     DOUBLE PRECISION      XPT    ( 3 )
C     INTEGER               HANDLE
C     INTEGER               DLADSC ( * )
C     DOUBLE PRECISION      DSKDSC ( * )
C     DOUBLE PRECISION      DC     ( * )
C     INTEGER               IC     ( * )
C     LOGICAL               FOUND
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body. 
C     NSURF      I   Number of surface IDs in list.
C     SRFLST     I   Surface ID list.
C     ET         I   Evaluation epoch, seconds past J2000 TDB.
C     FIXFID     I   Frame ID of ray and intercept.
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     XPT        O   Surface intercept point.
C     HANDLE     O   Handle of DSK file.
C     DLADSC     O   DLA descriptor of segment.
C     DSKDSC     O   DSK descriptor of segment.
C     DC         O   Double precision component of source info.
C     IC         O   Integer component of source info.
C     FOUND      O   Found flag. True if intercept exists.
C
C$ Detailed_Input
C
C     BODYID     is the ID code of a target body. The ray-surface
C                intercept computation is performed using data
C                that represent the surface of this body.
C
C     NSURF,
C     SRFLST     are, respectively, a count of surface IDs and
C                a list of IDs. If the list is non-empty, only
C                the indicated surfaces will be used. If the 
C                list is empty, all surfaces associated with the
C                input body ID will be considered.
C     
C     ET         is the epoch for which the computation is to be
C                performed. This epoch is used for DSK segment
C                selection; only segments containing ET in their time
C                coverage interval will be used. ET is expressed as
C                seconds past J2000 TDB.
C
C     FIXFID     is the frame ID of a body-fixed frame associated
C                with the body designated by BODYID. This frame
C                is used to express the input ray's vertex and
C                direction vector. The output intercept will be
C                expressed in this frame as well.
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector
C                of a ray. Both vectors are expressed in the frame
C                designated by FIXFID.
C
C$ Detailed_Output
C
C     XPT        is the surface intercept on the target body 
C                nearest to the ray's vertex, if the intercept
C                exists.
C
C     HANDLE     is the handle of the DSK file that contributed
C                the segment in which the intercept was found.
C
C     DLADSC     is the DLA descriptor of the segment in which
C                the intercept was found. 
C
C     DSKDSC     is the DSK descriptor of the segment in which
C                the intercept was found. 
C
C     DC         is the double precision component of the source
C                data associated with the intercept.
C
C     IC         is the integer component of the source
C                data associated with the intercept.
C
C                   For type 2 segments, this component contains
C                   a plate ID in the first element.
C
C     FOUND      is a logical flag that is set to .TRUE. if and
C                only if an intercept exists. The other outputs
C                are valid if and only if FOUND is .TRUE.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If this routine is called directly, it signals the error 
C         SPICE(BOGUSENTRY).
C        
C     See the entry points for descriptions of errors specific to
C     those routines.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine prepares the local buffers for a ray-surface
C     intercept computation using unprioritized DSK data. 
C     It calls ZZDSKBUX to perform the computation.
C
C$ Examples
C
C     See usage in ZZRAYSFX.
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
C-    SPICELIB Version 1.0.0, 22-FEB-2017 (NJB) 
C
C        Added FAILED calls.
C
C        12-MAY-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     unprioritized ray surface intercept with info
C
C-&


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN  ( 'ZZSBFXRI' )


      IF ( FIRST ) THEN
C
C        Initialize BSR counter.
C
         CALL ZZCTRUIN ( BSRCTR )
         
         FIRST = .FALSE.

      END IF

C
C     See whether the state of the loaded DSK set has changed
C     since the last call.
C
      CALL ZZDSKCHK ( BSRCTR, UPDATE )

      IF ( UPDATE ) THEN
C
C        Make sure the ZZDSKBSR subsystem has completed the segment
C        list for the input body since the last time the BSR loaded
C        kernel state changed.
C     
         CALL ZZDSKBBL ( BODYID )

C
C        Initialize the local buffers. We restart from scratch
C        each time the BSR loaded kernel state changes.
C
         CALL ZZDSKSBI ( BTSIZE, STSIZE, BTBODY, BTNBOD, BTSEGP,
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )
      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSBFXRI' )
         RETURN
      END IF
      
C
C     Find the index of the input body ID in the body table. If
C     we re-initialized the tables, the index will be zero.
C
      BIX = ISRCHI( BODYID, BTNBOD, BTBODY )

       
      IF ( BIX .EQ. 0 ) THEN
C
C        We don't have buffered information for this body. Update
C        the body and segment tables to store data for it.
C         
         CALL ZZDSKSBA ( BODYID,
     .                   BTSIZE, STSIZE, BTBODY, BTNBOD, BTSEGP, 
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSBFXRI' )
            RETURN
         END IF
C
C        The new body's position in the body table is at the end.
C         
         BIX = BTNBOD
         
      END IF

C
C     Find the ray-surface intercept, using the buffered segment
C     data.
C     
      J    = BTSEGP(BIX) 
      NSEG = BTSTSZ(BIX)    

      CALL ZZDSKBUX ( BODYID,     NSURF,      SRFLST, ET,  FIXFID,
     .                NSEG,       STHAN(J),   STDLAD(1,J), STDSCR(1,J),
     .                STOFF(1,J), STCTR(1,J), STRAD(J),    VERTEX,      
     .                RAYDIR,     XPT,        SEGIDX,      DC,
     .                IC,         FOUND                                )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSBFXRI' )
         RETURN
      END IF

      IF ( FOUND ) THEN
C
C        Adjust the segment index to make it relative to the base value
C        1, instead of the current base J.
C
         SEGIDX = SEGIDX + J - 1

         HANDLE = STHAN( SEGIDX )

         CALL MOVEI ( STDLAD(1,SEGIDX), DLADSZ, DLADSC )
         CALL MOVED ( STDSCR(1,SEGIDX), DSKDSZ, DSKDSC )

      END IF

      CALL CHKOUT ( 'ZZSBFXRI' )
      RETURN      





C$Procedure ZZSBFNRM ( DSK, prepare and compute unprioritized normal )

      ENTRY ZZSBFNRM ( BODYID, NSURF, SRFLST, ET,    
     .                 FIXFID, POINT, NORMAL     )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Prepare and execute unprioritized outward surface normal
C     computation using DSK API segment buffers.
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
C
C     INTEGER               BODYID
C     INTEGER               NSURF
C     INTEGER               SRFLST ( * )
C     DOUBLE PRECISION      ET
C     INTEGER               FIXFID
C     DOUBLE PRECISION      POINT  ( 3 )
C     DOUBLE PRECISION      NORMAL ( 3 )
C 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     BODYID     I   ID code of target body. 
C     NSURF      I   Number of surface IDs in list.
C     SRFLST     I   Surface ID list.
C     ET         I   Evaluation epoch, seconds past J2000 TDB.
C     FIXFID     I   Frame ID of surface point and normal vector.
C     POINT      I   Surface point.
C     NORMAL     O   Unit length outward normal at input point.
C
C$ Detailed_Input
C
C     BODYID     is the ID code of a target body. The outward
C                normal vector computation is performed using data
C                that represent the surface of this body.
C
C     NSURF,
C     SRFLST     are, respectively, a count of surface IDs and
C                a list of IDs. If the list is non-empty, only
C                the indicated surfaces will be used. If the 
C                list is empty, all surfaces associated with the
C                input body ID will be considered.
C     
C     ET         is the epoch for which the computation is to be
C                performed. This epoch is used for DSK segment
C                selection; only segments containing ET in their time
C                coverage interval will be used. ET is expressed as
C                seconds past J2000 TDB.
C
C     FIXFID     is the frame ID of a body-fixed frame associated with
C                the body designated by BODYID. This frame is used to
C                express the input surface point. The output normal
C                vector will be expressed in this frame as well.
C
C     POINT      is a surface point on the body designated by BODYID.
C                POINT is not arbitrary: it must be within a small 
C                tolerance of the actual surface. Usually POINT is
C                obtained as the output of a SPICE computation.
C
C                POINT is expressed in the frame designated by FIXFID.
C
C$ Detailed_Output
C
C     NORMAL     is the unit-length outward surface normal vector at
C                the input POINT.
C
C                NORMAL is expressed in the frame designated by FIXFID.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If this routine is called directly, it signals the error 
C         SPICE(BOGUSENTRY).
C        
C     See the entry points for descriptions of errors specific to
C     those routines.
C
C$ Files
C
C     This routine makes use of DSK files loaded by the ZZDSKBSR
C     subsystem. 
C
C$ Particulars
C
C     This routine is meant to be used only by the DSK subsystem.
C
C     This routine prepares the local buffers for a surface normal
C     vector computation using unprioritized DSK data. It calls
C     ZZDSKBUN to perform the computation.
C
C$ Examples
C
C     See usage in DSKNRM (provisional name).
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
C-    SPICELIB Version 1.0.0, 22-FEB-2017 (NJB) 
C
C        Added FAILED calls.
C
C        12-MAY-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     prepare unprioritized normal vector computation
C
C-&


      IF ( RETURN() ) THEN
         RETURN
      END IF
 
      CALL CHKIN ( 'ZZSBFNRM' )

C
C     The following block of code is straight cut-and-paste from
C     ZZSBFXR. (We need to consider packaging this code in a 
C     utility routine.)
C
      IF ( FIRST ) THEN
C
C        Initialize BSR counter.
C
         CALL ZZCTRUIN ( BSRCTR )
         
         FIRST = .FALSE.

      END IF

C
C     See whether the state of the loaded DSK set has changed
C     since the last call.
C
      CALL ZZDSKCHK ( BSRCTR, UPDATE )

      IF ( UPDATE ) THEN
C
C        Make sure the ZZDSKBSR subsystem has completed the segment
C        list for the input body since the last time the BSR loaded
C        kernel state changed.
C     
         CALL ZZDSKBBL ( BODYID )

C
C        Initialize the local buffers. We restart from scratch
C        each time the BSR loaded kernel state changes.
C
         CALL ZZDSKSBI ( BTSIZE, STSIZE, BTBODY, BTNBOD, BTSEGP, 
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )
      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZSBFNRM' )
         RETURN
      END IF

C
C     Find the index of the input body ID in the body table. If
C     we re-initialized the tables, the index will be zero.
C
      BIX = ISRCHI( BODYID, BTNBOD, BTBODY )

      IF ( BIX .EQ. 0 ) THEN
C
C        We don't have buffered information for this body. Update
C        the body and segment tables to store data for it.
C         
         CALL ZZDSKSBA ( BODYID,
     .                   BTSIZE, STSIZE, BTBODY, BTNBOD, BTSEGP,
     .                   BTSTSZ, STHAN,  STDSCR, STDLAD, STFREE,
     .                   STOFF,  STCTR,  STRAD                  )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZSBFNRM' )
            RETURN
         END IF
C
C        The new body's position in the body table is at the end.
C         
         BIX = BTNBOD
         
      END IF

C
C     Find the outward unit normal vector, using the buffered segment
C     data.
C     
      J    = BTSEGP(BIX) 
      NSEG = BTSTSZ(BIX)      

      CALL ZZDSKBUN ( BODYID,      NSURF,      SRFLST,     ET, 
     .                FIXFID,      NSEG,       STHAN(J),   STDLAD(1,J),
     .                STDSCR(1,J), STOFF(1,J), STCTR(1,J), STRAD(J), 
     .                POINT,       NORMAL                              )

      CALL CHKOUT ( 'ZZSBFNRM' )
      RETURN      

      END

