C$Procedure DSKCLS ( DSK, close file )
 
      SUBROUTINE DSKCLS ( HANDLE, OPTMIZ )
 
C$ Abstract
C
C     Close a DSK file.
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
C     DAS
C     DSK
C
C$ Keywords
C
C     DAS
C     DSK
C     FILES
C
C$ Declarations

      IMPLICIT NONE
      INCLUDE 'dla.inc'     

      INTEGER               HANDLE
      LOGICAL               OPTMIZ
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle assigned to the opened DSK file.
C     OPTMIZ     I   Flag indicating whether to segregate the DSK.
C
C$ Detailed_Input
C
C     HANDLE      is the DAS file handle associated with the file.
C                 The file may be open for read or write access.
C
C     OPTMIZ      is a logical flag indicating whether the DSK 
C                 should be segregated before it is closed. This
C                 option applies only to files open for write 
C                 access. The value of OPTMIZ has no effect for
C                 files opened for read access.
C
C                 See the DAS Required Reading das.req for a 
C                 discussion of segregation of DAS files.
C
C$ Detailed_Output
C
C     None. This routine operates by side effects.
C 
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If an error occurs when the file is closed, the error will be 
C        diagnosed by routines in the call tree of this routine.
C
C$ Files
C
C     See argument HANDLE.
C
C$ Particulars
C
C     This routine provides a DSK-level interface for closing DSK files.
C
C     In cases where DSKs opened for write access are to be closed
C     without segregation, this interface is slightly simpler than that
C     available at the DAS level.
C
C$ Examples
C
C     1) Close a new DSK file using DAS segregation. HANDLE
C        is the DAS file handle of the DSK. 
C
C        This is the normal choice for DSK creation.
C 
C           CALL DSKCLS ( HANDLE, .TRUE. )
C
C     2) Close a new DSK file without using DAS segregation. The 
C        close operation will be fast, but reading the file will be 
C        less efficient than if the file had been segregated.
C
C           CALL DSKCLS ( HANDLE, .TRUE. )
C
C     3) Close an existing DSK file that had been opened
C        for read access. In this case OPTMIZ is ignored:
C 
C           CALL DSKCLS ( HANDLE, .FALSE. )
C
C        or
C
C           CALL DSKCLS ( HANDLE, .TRUE. )
C     
C$ Restrictions
C
C     1) This routine should not be called by user applications
C        that have loaded a DSK file via FURNSH. Such applications
C        should call the KEEPER entry points UNLOAD or KCLEAR instead.
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
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C
C        09-OCT-2009 (NJB)
C
C           Updated header.
C
C        20-OCT-2006 (NJB)
C
C           Original DSKLIB version.
C
C-&
 
C$ Index_Entries
C
C     close a dsk file
C
C-&
 

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               METHLN
      PARAMETER           ( METHLN  = 10 )

C
C     Local variables
C
      CHARACTER*(METHLN)    METHOD

      
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKCLS' )

      IF ( OPTMIZ ) THEN
C
C        Segregate the file to enable fast read access.  This is
C        the "normal" way to close a DSK.  Segregating a large file
C        can be slow, however.
C
         CALL DASCLS ( HANDLE )

      ELSE
C
C        Close the file without first segregating it; this allows
C        the caller to close the file quickly, but results in a
C        file that will be read more slowly.
C
C        Any buffered data to be written must be explicitly flushed
C        to the file, if the file is open for write access.
C
         CALL DASHAM ( HANDLE, METHOD )
 
         IF ( METHOD .EQ. 'WRITE ' ) THEN
C
C           Write out any buffered records belonging to the
C           indicated file.
C
            CALL DASWBR ( HANDLE )

         END IF

C
C        Close the file without segregating records.
C
         CALL DASLLC ( HANDLE )
         
      END IF

      CALL CHKOUT ( 'DSKCLS' )
      RETURN
      END
