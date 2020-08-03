C$Procedure DSKOPN ( DSK, open new file )
 
      SUBROUTINE DSKOPN ( FNAME, IFNAME, NCOMCH, HANDLE )
 
C$ Abstract
C
C     Open a new DSK file for subsequent write operations.
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

      CHARACTER*(*)         FNAME
      CHARACTER*(*)         IFNAME
      INTEGER               NCOMCH
      INTEGER               HANDLE
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     FNAME      I   Name of a DSK file to be opened.
C     IFNAME     I   Internal file name.
C     NCOMCH     I   Number of comment characters to allocate.
C     HANDLE     O   Handle assigned to the opened DSK file.
C
C$ Detailed_Input
C
C     FNAME       is the name of a new DSK file to be created.  The
C                 file will be left opened for write access.
C
C     IFNAME      is the internal file name for the new file.  The name
C                 may contain as many as 60 characters.  All characters
C                 of IFNAME should be printing characters (ASCII codes
C                 32-126 decimal). This name should uniquely identify
C                 the file.
C
C     NCOMCH      is the number of comment characters to allocate.
C                 Allocating comment characters at file creation time
C                 may reduce the likelihood of having to expand the
C                 comment area later.
C
C$ Detailed_Output
C
C     HANDLE      is the file handle associated with the file. This
C                 handle is used to identify the file in subsequent
C                 calls to other DSK routines.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the input filename is blank, the error will be diagnosed
C         by routines in the call tree of this routine.  No file will
C         be created.
C
C     2)  If the specified file cannot be opened without exceeding the
C         maximum allowed number of open DAS files, the error will be
C         diagnosed by routines in the call tree of this routine.  No
C         file will be created.
C
C     3)  If the file cannot be opened properly, the error will be
C         diagnosed by routines in the call tree of this routine.  No
C         file will be created.
C
C     4)  If the initial records in the file cannot be written, the
C         error is diagnosed by routines in the call tree of this
C         routine.  No file will be created.
C
C     5)  If no logical units are available, the error will be
C         diagnosed by routines in the call tree of this routine. No
C         file will be created.
C
C     6)  If the internal file name contains nonprinting characters
C         (ASCII codes decimal 0-31 and 127-255), the error will be
C         diagnosed by routines in the call tree of this routine.  No
C         file will be created.
C
C     7)  If the number of comment characters allocated NCOMCH is
C         negative, the error will be diagnosed by routines in the call
C         tree of this routine.  No file will be created.
C
C$ Files
C
C     See argument FNAME.
C
C$ Particulars
C
C     DSK files are built using the DLA low-level format and
C     the DAS architecture; DLA files are a specialized type of DAS
C     file in which data are organized as a doubly linked list of
C     segments.  Each segment's data belong to contiguous components of
C     character, double precision, and integer type.
C
C     This routine creates a new DSK file and sets the type of the
C     file to the mnemonic code passed to it.
C
C     DSK files created by this routine have initialized file records.
C     The ID word in a DSK file record has the form
C
C        DAS/DSK
C
C     where the characters following the slash are supplied by the
C     caller of this routine.
C
C$ Examples
C
C     1)  Create a new DSK file, using an internal file name that
C         attempts to serve as an unique identifier.  No room for
C         comments will be reserved.
C
C            FNAME  =  'TEST.DSK'
C            IFNAME =  'TEST.DSK/NAIF/NJB/20-OCT-2006/14:37:00'
C            NCOMCH =   0
C 
C            CALL DSKOPN ( FNAME, IFNAME, NCOMCH, HANDLE )
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C        Corrected a few header typos.
C
C        29-APR-2010 (NJB)
C
C           Now passes NCOMCH to DLAOPN.
C
C        08-OCT-2009 (NJB)
C
C           Updated header.
C
C        20-OCT-2006 (NJB)
C
C-&
 
C$ Index_Entries
C
C     open a new dsk file
C     open a new dsk file with write access
C
C-&
 

C
C     SPICELIB functions
C
      LOGICAL               RETURN


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKOPN' )

      CALL DLAOPN ( FNAME, 'DSK', IFNAME, NCOMCH, HANDLE )

      CALL CHKOUT ( 'DSKOPN' )
      RETURN
      END


      
