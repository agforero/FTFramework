C$Procedure ATTCMP ( Attribute comparison for DSK segments )

      SUBROUTINE ATTCMP ( HAN1,   DLADS1, HAN2,  DLADS2, 
     .                    TIMTOL, MATCH,  TMATCH         )

C$ Abstract
C
C     Given handles and DLA descriptors for a pair of DSK segments,
C     indicate whether the segments have matching attributes.
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
C     DSKBRIEF
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'

      INTEGER               HAN1
      INTEGER               DLADS1 ( * )
      INTEGER               HAN2
      INTEGER               DLADS2 ( * )
      DOUBLE PRECISION      TIMTOL
      LOGICAL               MATCH
      LOGICAL               TMATCH

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     HAN1       I   Handle of DSK containing first segment.
C     DLADS1     I   DLA descriptor of first segment.
C     HAN2       I   Handle of DSK containing second segment.
C     DLADS1     I   DLA descriptor of second segment.
C     TIMTOL     I   Tolerance to use for time comparisons.
C     MATCH      O   Flag indicating whether segments match.
C     TMATCH     O   Flag indicating whether segments match 
C                    when time bound are considered.
C
C$ Detailed_Input
C
C     HAN1,
C     DLADS1     are, respectively, the DSK handle and DLA descriptor
C                of the first segment to be compared.
C
C     HAN2,
C     DLADS2     are, respectively, the DSK handle and DLA descriptor
C                of the second segment to be compared.
C
C     TIMTOL     is a tolerance value to be used for time bound 
C                comparisons. TIMTOL is a non-negative number.
C                Units are TDB seconds.
C
C$ Detailed_Output
C
C     MATCH      is a logical flag that is set to .TRUE. if and only 
C                if the DSK segments designated by the input handles and
C                DLA descriptors have matching attributes, excluding
C                time boundaries.
C
C                In order for two DSK segments to match, the following
C                segment attributes must match:
C
C                   Body
C                   Surface
C                   Frame
C                   Coordinate system
C
C                      If the coordinate system is planetodetic, the
C                      system parameters must match exactly.
C
C                   Data type
C                   Data class
C
C     TMATCH     is a logical flag that is set to .TRUE. if and only if
C                the DSK segments designated by the input handles and
C                DLA descriptors have matching time boundaries, to
C                within the tolerance specified by TIMTOL.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If a DSK data look-up fails, an error will be signaled
C         by a routine in the call tree of this routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports grouping of segments for abbreviated
C     summaries created by DSKBRIEF.
C     
C$ Examples
C
C     See usage in DSKBRIEF.
C
C$ Restrictions
C
C     1) For use only within program DSKBRIEF.
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
C-    DSKBRIEF Version 1.0.0, 15-FEB-2017 (NJB)
C     
C        Added time bounds check as a match criterion.
C
C        30-AUG-2016 (NJB)
C
C           Original version.
C-&

 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      DSKDS1 ( DSKDSZ )
      DOUBLE PRECISION      DSKDS2 ( DSKDSZ )

      LOGICAL               FAILED
      LOGICAL               RETURN

 
      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'ATTCMP' )
C
C     No match so far.
C
      MATCH = .FALSE.

C
C     Get DSK descriptors for both segments.
C
      CALL DSKGD ( HAN1, DLADS1, DSKDS1 )
      CALL DSKGD ( HAN2, DLADS2, DSKDS2 )
         
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ATTCMP' )
         RETURN
      END IF

C
C     The following must match:
C
C        Body
C        Surface
C        Frame
C        Coordinate system
C        Data type
C        Data class
C
      MATCH =       ( DSKDS1(CTRIDX) .EQ. DSKDS2(CTRIDX) )
     .        .AND. ( DSKDS1(SRFIDX) .EQ. DSKDS2(SRFIDX) )
     .        .AND. ( DSKDS1(FRMIDX) .EQ. DSKDS2(FRMIDX) )
     .        .AND. ( DSKDS1(SYSIDX) .EQ. DSKDS2(SYSIDX) )
     .        .AND. ( DSKDS1(TYPIDX) .EQ. DSKDS2(TYPIDX) )
     .        .AND. ( DSKDS1(CLSIDX) .EQ. DSKDS2(CLSIDX) )

C
C     If the coordinate system is planetodetic, the system
C     parameters must match exactly.
C
      IF ( NINT(DSKDS1(SYSIDX)) .EQ. PDTSYS ) THEN

         MATCH =       ( DSKDS1(PARIDX  ) .EQ. DSKDS2(PARIDX  ) )
     .           .AND. ( DSKDS1(PARIDX+1) .EQ. DSKDS2(PARIDX+1) )
      END IF

C
C     Compare start and stop times.
C
      TMATCH =      ( ABS(DSKDS1(BTMIDX) - DSKDS2(BTMIDX)) .LE. TIMTOL )
     .        .AND. ( ABS(DSKDS1(ETMIDX) - DSKDS2(ETMIDX)) .LE. TIMTOL )

      CALL CHKOUT ( 'ATTCMP' )
      RETURN
      END

