C$Procedure DLAOPN ( DLA, open new file )
 
      SUBROUTINE DLAOPN ( FNAME, FTYPE, IFNAME, NCOMCH, HANDLE )
      IMPLICIT NONE
 
C$ Abstract
C
C     Open a new DLA file and set the file type.
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
C     DLA
C
C$ Keywords
C
C     DAS
C     DLA
C     FILES
C
C$ Declarations

      INCLUDE 'dla.inc'     

      CHARACTER*(*)         FNAME
      CHARACTER*(*)         FTYPE
      CHARACTER*(*)         IFNAME
      INTEGER               NCOMCH
      INTEGER               HANDLE
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     FNAME      I   Name of a DLA file to be opened.
C     FTYPE      I   Mnemonic code for type of data in the DLA file.
C     IFNAME     I   Internal file name.
C     NCOMR      I   Number of comment records to allocate.
C     HANDLE     O   Handle assigned to the opened DLA file.
C
C$ Detailed_Input
C
C     FNAME       is the name of a new DLA file to be created.  The
C                 file will be left opened for write access.
C
C     FTYPE       is a code for type of data placed into a DLA file. The
C                 non-blank part of FTYPE is used as the "file type"
C                 portion of the ID word in the DLA file.
C
C                 The first nonblank character and the three, or fewer,
C                 characters immediately following it, giving four
C                 characters, are used to represent the type of the
C                 data placed in the DLA file. This is provided as a
C                 convenience for higher level software. It is an error
C                 if this string is blank. Also, the file type may not
C                 contain any nonprinting characters. When written to
C                 the DLA file, the value for the type IS case
C                 sensitive.
C
C                 NAIF has reserved for its own use file types
C                 consisting of the upper case letters (A-Z) and the
C                 digits 0-9. NAIF recommends lower case or mixed case
C                 file types be used by all others in order to avoid
C                 any conflicts with NAIF file types.
C
C     IFNAME      is the internal file name for the new file.  The name
C                 may contain as many as 60 characters.  This name
C                 should uniquely identify the file.
C
C     NCOMR       is the number of comment records to allocate.
C                 Allocating comment records at file creation time may
C                 reduce the likelihood of having to expand the
C                 comment area later.
C
C$ Detailed_Output
C
C     HANDLE      is the file handle associated with the file. This
C                 handle is used to identify the file in subsequent
C                 calls to other DLA routines.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the input filename is blank, the error will be diagnosed by
C        routines in the call tree of this routine.  No file will be
C        created.
C
C     2) If the specified file cannot be opened without exceeding
C        the maximum allowed number of open DAS files, the error
C        will be diagnosed by routines in the call tree of this 
C        routine.  No file will be created.
C
C     3) If the file cannot be opened properly, the error will be
C        diagnosed by routines in the call tree of this routine.  No
C        file will be created.
C
C     4) If the initial records in the file cannot be written, the
C        error is diagnosed by routines in the call tree of this
C        routine.  No file will be created.
C
C     5) If no logical units are available, the error will be diagnosed
C        by routines in the call tree of this routine. No file will be
C        created.
C
C     6) If the file type is blank, the error will be diagnosed by
C        routines in the call tree of this routine.  No file will be
C        created.
C
C     7) If the file type contains nonprinting characters, decimal
C        0-31 and 127-255, the error will be diagnosed by routines in
C        the call tree of this routine.  No file will be created.
C
C     8) If the number of comment records allocated NCOMR is negative,
C        the error SPICE(BADRECORDCOUNT) will be signaled. No file will
C        be created.
C
C$ Files
C
C     See argument FNAME.
C
C$ Particulars
C
C     DLA files are built using the DAS low-level format; DLA files are
C     a specialized type of DAS file in which data are organized as a
C     doubly linked list of segments. Each segment's data belong to
C     contiguous components of character, double precision, and integer
C     type.
C
C     This routine creates a new DLA file and sets the type of the
C     file to the mnemonic code passed to it.
C
C     DLA files created by this routine have initialized file records.
C     The ID word in a DLA file record has the form
C
C        DAS/****
C
C     where the characters following the slash are supplied by the
C     caller of this routine.
C
C$ Examples
C
C     1)  Create a new DLA file, using an internal file name that
C         attempts to serve as an unique identifier, and give the file a
C         type of 'TEST'.  No room for comments will be reserved.
C
C            FNAME  =  'TEST.DLA'
C            FTYPE  =  'TEST'
C            IFNAME =  'TEST.DLA/NAIF/NJB/07-FEB-2005-02:57:00'
C            NCOMCH =   0
C 
C            CALL DLAOPN ( FNAME, FTYPE, IFNAME, NCOMCH, HANDLE )
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
C     K.R. Gehringer  (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 08-FEB-2017 (NJB)
C
C        Updated version info.
C
C        01-APR-2016 (NJB)
C
C           Changed short error message for invalid comment
C           count. Corrected reference to "DASCLU" in comments.
C
C        08-OCT-2009 (NJB)
C
C           Updated header.
C
C        09-FEB-2005 (NJB) (KRG)
C
C-&
 
C$ Index_Entries
C
C     open a new dla file
C     open a new dla file with write access
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
      LOGICAL               RETURN

C
C     Local parameters
C

C
C     Local variables
C
      INTEGER               NCOMR


      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'DLAOPN' )

C
C     Compute the number of comment records required.
C
      IF ( NCOMCH .GT. 0 ) THEN

         NCOMR = ( (NCOMCH-1) / NCHREC )  +  1

      ELSE IF ( NCOMCH .EQ. 0 ) THEN

         NCOMR = 0

      ELSE

         CALL SETMSG ( 'Requested number of comment characters ' 
     .   //            'must be non-negative but was #.'         )
         CALL ERRINT ( '#',  NCOMCH                              )
         CALL SIGERR ( 'SPICE(BADRECORDCOUNT)'                   )
         CALL CHKOUT ( 'DLAOPN'                                  )
         RETURN

      END IF

C
C     Let the DAS "open new" routine do the work.
C
      CALL DASONW ( FNAME, FTYPE, IFNAME, NCOMR, HANDLE )

C
C     Write the format version.
C
      CALL DASADI ( HANDLE, 1, FMTVER )

C
C     Initialize the forward and backward segment list pointers.
C
      CALL DASADI ( HANDLE, 1, NULPTR )
      CALL DASADI ( HANDLE, 1, NULPTR )

C
C     We leave the file open, since further writes to the file
C     should occur next.  The file will eventually be closed
C     by a call to DASCLS or DASLLC, if all goes well.
C
      CALL CHKOUT ( 'DLAOPN' )
      RETURN
      END

