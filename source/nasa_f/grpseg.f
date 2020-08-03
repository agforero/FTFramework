C$Procedure GRPSEG ( Group DSK segments having shared attributes )

      SUBROUTINE GRPSEG ( UDCOMP, USETIM, TIMTOL, N,  HANDLS, 
     .                    DLADS,  MXPOOL, MXGRP,  NG, SGPTRS,
     .                    SGPOOL, SEGLST, SEGTM              )

C$ Abstract
C
C     Given arrays of handles and DLA descriptors of DSK segments,
C     and given a callback comparison function, group the segments
C     into subsets having matching attributes.
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

      INTEGER               LBPOOL
      PARAMETER           ( LBPOOL = -5 )

      EXTERNAL              UDCOMP
      LOGICAL               USETIM
      DOUBLE PRECISION      TIMTOL
      INTEGER               N
      INTEGER               HANDLS ( * )
      INTEGER               DLADS  ( DLADSZ, * )
      INTEGER               MXPOOL
      INTEGER               MXGRP
      INTEGER               NG
      INTEGER               SGPTRS ( * )
      INTEGER               SGPOOL ( 2,  LBPOOL : * )
      INTEGER               SEGLST ( * )
      LOGICAL               SEGTM  ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     UDCOMP     I   Comparison function for DSK segments.
C     USETIM     I   Time comparison flag.
C     TIMTOL     I   Time comparison tolerance.
C     N          I   Size of handle and DLA descriptor arrays.
C     HANDLS     I   Handles of DSK files.
C     DLADS      I   DLA descriptors of DSK segments.
C     MXPOOL     I   Maximum pool size.
C     MXGRP      I   Maximum group count.
C     NG         O   Number of segment groups.
C     SGPTRS     O   Pointers from groups to segment lists.
C     SGPOOL     O   Segment pool.
C     SEGLST     O   Map from pool nodes to input segment array indices.
C     SEGTM      O   Array of time match flags for segment groups.
C     
C$ Detailed_Input
C
C     UDCOMP     is a callback function which is used by this
C                routine to compare attributes of DSK segments.
C                The calling sequence of UDCOMP is
C
C                   UDCOMP ( HAN1, DLADS1, HAN2,  DLADS2, 
C                            TOL,  MATCH,  TMATCH        ) 
C
C                Inputs:
C
C                   INTEGER           HAN1         {handle of segment 1}
C                   DOUBLE PRECISION  DLADS1(*)    {DLA descriptor of 
C                                                   segment 1}
C                   INTEGER           HAN2         {handle of segment 2}
C                   DOUBLE PRECISION  DLADS2(*)    {DLA descriptor of 
C                                                   segment 2}
C                   DOUBLE PRECISION  TOL          {Time tolerance}
C
C                Outputs:
C
C                   LOGICAL           MATCH        {flag indicating
C                                                   whether segments
C                                                   match}
C
C                   LOGICAL           TMATCH       {flag indicating
C                                                   whether segment
C                                                   time bounds match}
C
C
C     USETIM     is a logical flag that indicates whether to consider
C                segment time bounds when grouping segments.
C
C     TIMTOL     is a time tolerance to use for comparing time bounds.
C                Units are TDB seconds.
C
C     N          is the number of elements in the input array HANDLS
C                and the number of DLA descriptors in the input array
C                DLADS.
C                     
C
C     HANDLS     
C     DLADS      are, respectively, arrays of handles of DSK files and
C                DLA descriptors. There are N entries in each array.
C                The Ith element of HANDLS and the Ith descriptor in
C                DLADS correspond to the Ith segment to be grouped.
C
C     MXPOOL     is the maximum number of entries that can be 
C                accommodated in the linked list pool array SGPOOL.
C
C     MXGRP      is the maximum number of entries that can be
C                placed in the array SGPTRS.
C
C
C$ Detailed_Output
C
C     NG         is the number of groups comprised by the input
C                segments.
C
C     
C     SGPTRS     is an array that maps group numbers to head nodes of
C                segment lists in SGPOOL. SGPTRS(I) is the index of the
C                head node of the segment list of the Ith group.
C                 
C                SGPTRS must be declared with size at least MXGRP.
C
C
C     SGPOOL     is a doubly linked list pool used to store segment
C                lists. Entries for a given segment group are linked
C                together. 
C
C                SGPOOL must be declared with its second dimension
C                at least equal to MXPOOL.
C
C
C     SEGLST     is an array that maps nodes of SGPOOL to entries
C                in the input arrays HANDLS and DLADS. The Ith node
C                of SGPOOL 
C
C                   SGPOOL(*,I) 
C
C                corresponds to the segment designated by 
C
C                   HANDLS(    SEGLST(I) )
C                   DLADS ( *, SEGLST(I) )
C
C
C     SEGTM      is an array of flags associated with segment groups.
C                The Ith element of SEGTIM is associated with the 
C                group starting at node SGPTRS(I). Each element of
C                SEGTM that is associated with a group is .TRUE. if
C                and only if the segments in that group have matching
C                time bounds. The elements of SEGTM that are not 
C                assocated with groups are undefined.
C
C                SEGTM must be declared with size at least MXGRP.
C
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the number of segments is not at least 1, the error
C         SPICE(INVALIDCOUNT) is signaled.
C
C     2)  If the pool size is not at least N, the error
C         SPICE(ARRAYTOOSMALL) is signaled.
C
C     3)  If the maximum group count is not at least 1, the error
C         SPICE(INVALIDSIZE) is signaled. Normally room for multiple
C         groups should be provided by the caller.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine groups segments for abbreviated summaries
C     created by DSKBRIEF. Segments that are considered by UDCOMP
C     to have matching attributes are grouped together.
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
C     DSKBRIEF Version 1.0.0, 15-FEB-2017 (NJB)
C
C-&

C
C     SPICELIB functions
C
      INTEGER               LNKTL

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C     
      INTEGER               HEAD
      INTEGER               I
      INTEGER               J
      INTEGER               K
      INTEGER               NODE
      INTEGER               TAIL

      LOGICAL               MATCH
      LOGICAL               TMATCH


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'GRPSEG' )

      IF ( N .LT. 1 ) THEN

         CALL SETMSG ( 'Number of segments is #; must be at least 1.' )
         CALL ERRINT ( '#', N                                         )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                          )
         CALL CHKOUT ( 'GRPSEG'                                       )
         RETURN

      END IF

      IF ( MXPOOL .LT. N ) THEN

         CALL SETMSG ( 'Pool size is #; must be at least #.' )
         CALL ERRINT ( '#', MXPOOL                           )
         CALL ERRINT ( '#', N                                )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                )
         CALL CHKOUT ( 'GRPSEG'                              )
         RETURN

      END IF

      IF ( MXGRP .LT. 1 ) THEN
C
C        Normally MXGRP should allow room for multiple groups. However,
C        we can't know in advance how much room is needed, except that
C        the group count can't exceed N. We simply catch clearly 
C        invalid values here.
C
         CALL SETMSG ( 'Group array size is #; must be at least 1.' )
         CALL ERRINT ( '#', MXGRP                                   )
         CALL SIGERR ( 'SPICE(INVALIDSIZE)'                         )
         CALL CHKOUT ( 'GRPSEG'                                     )
         RETURN

      END IF

C
C     Initialize the segment pool.
C
      CALL LNKINI ( MXPOOL, SGPOOL )

C
C     Initialize the time match flags.
C
      DO I = 1, MXGRP
         SEGTM(I) = .TRUE.
      END DO

C
C     The first segment is automatically in the first group.
C
      CALL LNKAN ( SGPOOL, NODE )

      NG           = 1
      SGPTRS(NG)   = NODE
      SEGLST(NODE) = 1

C
C     Now examine the other segments. If a segment matches one
C     we've already seen, the segment is added to the latter's
C     group. Otherwise, the segment becomes the first member of
C     a new group.
C
      DO I = 2, N
C
C        Compare the Ith segment to the first segment
C        of each group, until we find a match or run out
C        of groups.
C
         MATCH = .FALSE.

         J     = 0
         
         DO WHILE (  ( .NOT. MATCH ) .AND. ( J .LT. NG )  )

            J = J + 1
C
C           Compare the Ith segment against the first segment of the
C           Jth group.
C
            HEAD = SGPTRS(J)

            K    = SEGLST( HEAD ) 

            CALL UDCOMP ( HANDLS(I), DLADS(1,I), HANDLS(K), DLADS(1,K),
     .                    TIMTOL,    MATCH,      TMATCH                )

            IF ( USETIM ) THEN
C
C              Time bounds are being considered; we have a match only
C              if the time bounds and the other attributes match.
C               
               MATCH = MATCH .AND. TMATCH

            END IF


            IF ( MATCH ) THEN
C
C              Indicate whether the group to which this segment belongs
C              has inconsistent time tags.
C
               IF ( .NOT. TMATCH ) THEN

                  SEGTM(J) = .FALSE.

               END IF

            END IF


            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'GRPSEG' )
               RETURN
            END IF

         END DO

         IF ( MATCH ) THEN
C
C           The Ith segment belongs to group J. Link a node for the
C           segment to the tail of the list for this group.
C           
            TAIL = LNKTL ( HEAD, SGPOOL ) 
 
            CALL LNKAN ( SGPOOL, NODE )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'GRPSEG' )
               RETURN
            END IF

            CALL LNKILA ( TAIL, NODE, SGPOOL )

            SEGLST(NODE) = I

         ELSE
C
C           This segment is in a category of its own. Create a new
C           segment pointer for it.
C
            IF ( NG .EQ. MXGRP ) THEN

               CALL SETMSG ( 'Size of group array is #; cannot add '
     .         //            'new element.'                         )
               CALL ERRINT ( '#', MXGRP                             )
               CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                 )
               CALL CHKOUT ( 'GRPSEG'                               )
               RETURN

            END IF

 
            CALL LNKAN ( SGPOOL, NODE )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'GRPSEG' )
               RETURN
            END IF

            NG           = NG + 1

            SGPTRS(NG)   = NODE
C
C           Associate the Ith segment with the new node.
C
            SEGLST(NODE) = I

         END IF

      END DO


      CALL CHKOUT ( 'GRPSEG' )
      RETURN
      END 
