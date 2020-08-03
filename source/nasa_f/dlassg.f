C$Procedure DLASSG ( DLA, same segment? )
 
      LOGICAL FUNCTION DLASSG ( HAN1, HAN2, DSC1, DSC2 )
 
C$ Abstract
C
C     Return a logical value indicating whether a two DLA
C     segments, each identified by DAS handle and DLA descriptor,
C     are in fact the same segment.
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

      IMPLICIT NONE

      INCLUDE 'dla.inc'

      INTEGER               HAN1
      INTEGER               HAN2
      INTEGER               DSC1 ( DLADSZ )
      INTEGER               DSC2 ( DLADSZ )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HAN1       I   Handle of an open DLA file.
C     HAN2       I   Handle of a second open DLA file.
C     DSC1       I   DLA descriptor of a segment in the first file.
C     DSC2       I   DLA descriptor of a segment in the second file.
C
C     The function returns .TRUE. if and only if the DLA segments
C     match.
C     
C$ Detailed_Input
C
C     HAN1        is the integer handle associated with a DLA file.
C                 The file is open for read access.
C
C     HAN2        is the integer handle associated with a second DLA
C                 file. The file is open for read access.
C
C     DSC1        is the DLA descriptor of a segment in the file
C                 associated with HAN1.
C
C     DSC2        is the DLA descriptor of a segment in the file
C                 associated with HAN2.
C
C$ Detailed_Output
C
C     The function returns .TRUE. if and only if the DLA segments
C     match. The segments are considered to match if and only if the
C     input handles match and all elements of the DLA descriptors
C     match.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If any of the inputs are invalid, this routine will
C        fail in an unspecified manner. 
C
C$ Files
C
C     See description of input arguments HAN1 and HAN2.
C
C$ Particulars
C
C     DLA files are built using the DAS low-level format; DLA files are
C     a specialized type of DAS file in which data are organized as a
C     doubly linked list of segments.  Each segment's data belong to
C     contiguous components of character, double precision, and integer
C     type.
C
C     This routine supports DLA and DSK routines by enabling
C     them to determine whether a given DLA segment matches one
C     they've previously examined. This may allow such routines
C     to avoid buffering information redundantly.
C
C$ Examples
C
C     1)  A typical use of this routine is to enable a subroutine
C         to determine whether a DLA segment identified by a
C         handle and DLA descriptor matches one seen previously.
C         The logic of such a test can be implemented as follows:
C
C
C                   SUBROUTINE SUBA ( HANDLE, DLADSC )
C                   IMPLICIT NONE
C 
C                   INCLUDE 'dla.inc'
C
C                   INTEGER               HANDLE
C                   INTEGER               DLADSC ( * )
C
C             C
C             C     SPICELIB functions
C             C             
C                   LOGICAL               DLASSG
C                   LOGICAL               FAILED
C             C
C             C     Local variables
C             C
C                   INTEGER               PRVDSC ( DLADSZ )
C                   INTEGER               PRVHAN
C                 
C             C
C             C     Saved variables
C             C
C                   SAVE                  PRVDSC
C                   SAVE                  PRVHAN
C
C             C
C             C     Initial values
C             C
C                   DATA                  PRVHAN / 0 /
C
C                   ...
C
C                   IF ( .NOT. DLASSG( HANDLE, PRVHAN,
C                  .                   DLADSC, PRVDSC ) ) THEN
C
C                      [Examine segment]
C
C                      IF ( .NOT. FAILED() ) THEN
C             C
C             C           Save values only if no error occurred.
C             C
C                         CALL MOVEI ( DLADSC, DLADSZ, PRVDSC )
C                         PRVHAN = HANDLE
C
C                      END IF
C
C                   END IF
C
C                   [Normal case]
C
C                   ...
C
C                   END
C
C
C$ Restrictions
C
C     This routine relies on uniqueness of DAS file handles.
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
C-    SPICELIB Version 1.0.0, 19-MAY-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     test DLA segments for match
C
C-&

C
C     Local variables
C
      INTEGER               I

C
C     Give the function an initial value.
C
      DLASSG = .FALSE.

C
C     If the handles don't match, we're done.
C     
      IF ( HAN1 .NE. HAN2 ) THEN
         RETURN
      END IF
C
C     Compare the DLA descriptors. All elements, including pointers,
C     must match in order to have a matching result.
C     
      DO I = 1, DLADSZ

         IF ( DSC1(I) .NE. DSC2(I) ) THEN
            RETURN
         END IF

      END DO
C
C     At this point, everything's a match.
C
      DLASSG = .TRUE.

      RETURN
      END
