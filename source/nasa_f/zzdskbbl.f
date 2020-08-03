C$Procedure ZZDSKBBL ( DSK, build BSR segment list )
 
      SUBROUTINE ZZDSKBBL ( BODYID )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Force the ZZDSKBSR subsystem to build a complete, current segment
C     list for a specified body.
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
C     DSK
C     BODY
C     SEGMENT
C     TOPOGRAPHY
C
C$ Declarations

      IMPLICIT NONE
      
      INCLUDE 'zzctr.inc'
      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsk02.inc'

      INTEGER               BODYID

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     BODYID     I   Body ID code.
C
C$ Detailed_Input
C
C     BODYID     is the ID code of the body for which a complete
C                segment list is to be built by the ZZDSKBSR subsystem.
C 
C$ Detailed_Output
C
C     None. See Particulars for a description of the side effects of
C     this routine.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an error occurs while building the segment list for 
C         the specified body, the error will be diagnosed by 
C         a routine in the call tree of this routine.
C          
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     Loaded DSK files are searched by ZZDSKBSR entry points as a
C     result of a call to this routine.
C
C     Kernel data are normally loaded once per program run, NOT every
C     time this routine is called.
C     
C$ Particulars
C
C     This routine forces the ZZDSKBSR subsystem to form a complete
C     segment list for a specified body. This will allow a ZZDSKBSR
C     segment search to obtain information about all loaded segments
C     for that body.
C
C     This routine enables routines that use the DSK API segment
C     buffers to accumulate complete segment data from the ZZDSKBSR
C     subsystem.
C
C     This routine may be rather slow, since it can cause a great
C     number of physical file reads to occur.  It should be called only
C     when the set of loaded DSKs is known to have changed.
C
C$ Examples
C
C     See usage in ZZDSKSBF entry points ZZSBFXR and ZZSBFNRM.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by the DSK
C     subsystem.
C      
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-JUL-2016 (NJB) 
C
C        Updated to use ZZDSKNOT. Bug fix: change of body
C        now causes an update.
C
C        05-FEB-2016 (NJB) 
C
C           Based on Version 1.0.0 25-SEP-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     force zzdskbsr to build segment list for body
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     EXTERNAL routines
C     
      EXTERNAL ZZDSKNOT

C
C     Local parameters
C     
       
C
C     Local variables
C     
      DOUBLE PRECISION      DSKDSC ( DSKDSZ )

      INTEGER               CTR    ( CTRSIZ )
      INTEGER               DLADSC ( DLADSZ )
      INTEGER               HANDLE
      INTEGER               PRVBOD

      LOGICAL               FIRST
      LOGICAL               NEWBOD
      LOGICAL               SEGFND
      LOGICAL               UPDATE

C
C     Saved variables
C
      SAVE                  CTR
      SAVE                  FIRST
      SAVE                  PRVBOD

C
C     Initial values
C
      DATA                  CTR    / CTRSIZ * 0 /
      DATA                  FIRST  / .TRUE.     /
      DATA                  PRVBOD / 0 /

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZDSKBBL' )


      IF ( FIRST ) THEN

         CALL ZZCTRUIN ( CTR )

         FIRST  = .FALSE.
         NEWBOD = .TRUE.

      ELSE
C
C        PRVBOD is the last body for which a successful list
C        build was performed.
C
         NEWBOD = BODYID .NE. PRVBOD

      END IF

C
C     See whether the state of the loaded DSK set has changed
C     since the last call.
C
      CALL ZZDSKCHK ( CTR, UPDATE )

      IF ( UPDATE .OR. NEWBOD ) THEN
C        
C        Force the BSR subsystem to build a complete segment list
C        for the current body. We'll do this by using a matching
C        function that indicates no segment matches: this forces
C        ZZDSKSNS to search all files for segments and build a 
C        complete segment search.
C
C        Start a segment search.
C
         CALL ZZDSKBSS ( BODYID )

         CALL ZZDSKSNS ( ZZDSKNOT, HANDLE, DLADSC, DSKDSC, SEGFND )
         
      
         IF ( .NOT. FAILED() ) THEN
  
            PRVBOD = BODYID

         END IF

      END IF


      CALL CHKOUT ( 'ZZDSKBBL' )
      RETURN
      END
