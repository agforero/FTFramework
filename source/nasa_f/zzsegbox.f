C$Procedure ZZSEGBOX (Bounding box for DSK segment volume element)
 
      SUBROUTINE ZZSEGBOX ( DSKDSC, BOXCTR, MAXR )
 
C$ Abstract
C
C     Create a bounding box for a DSK segment volume element.
C     The outputs are the box's center and radius. The center
C     is relative to the segment's frame center.
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
C     GEOMETRY
C     MATH
C
C$ Declarations

      INCLUDE 'dskdsc.inc'

      DOUBLE PRECISION      DSKDSC ( * )
      DOUBLE PRECISION      BOXCTR ( 3 )
      DOUBLE PRECISION      MAXR

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     DSKDSC     I   DSK segment descriptor.
C     CENTER     O   Center of bounding box.
C     RADIUS     O   Radius of box.
C
C$ Detailed_Input
C
C     DSKDSC     is the DSK descriptor for the segment of interest.
C
C$ Detailed_Output
C
C     CENTER     is a double precision 3-vector representing the center
C                of a box tangent to and containing the volume
C                specified by the coordinate system, coordinate
C                parameters, if any, and bounds specified in the DSK
C                segment descriptor.
C
C                The box bounds the surface specified by the input
C                descriptor; the offset between the center of the
C                segment's reference frame and the segment's central
C                body plays no role. 
C
C     RADIUS     is the radius of the sphere that circumscribes the
C                box. RADIUS is equal to the length of a line segment
C                connecting the center of the box to any corner.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If coordinate system specified in the descriptor is not
C         recognized, the error SPICE(NOTSUPPORTED) is signaled.
C
C     2)  Any errors that occur while deriving the bounding box
C         will be signaled by routines in the call tree of this
C         routine.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine consolidates the logic branch needed to compute a
C     bounding box for a DSK segment. Using this routine can eliminate
C     the need to test higher-level routines using each type of
C     coordinate system.
C
C     The box bounds the surface specified by the input descriptor; the
C     offset between the center of the segment's reference frame and
C     the segment's central body plays no role.
C
C     To compute a body's DSK bounding sphere, centered at the body,
C     and taking into account all loaded segments, use ZZDSKSPH.
C
C$ Examples
C
C     None.
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 17-JUN-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     bounding box for dsk segment volume element
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               RETURN
      
C
C     Local variables
C
      DOUBLE PRECISION      L1
      DOUBLE PRECISION      L2
      DOUBLE PRECISION      L3

      INTEGER               CORSYS


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZSEGBOX' )


      CORSYS = NINT( DSKDSC(SYSIDX) )

      IF ( CORSYS .EQ. LATSYS ) THEN

         CALL ZZLATBOX ( DSKDSC(MN1IDX), BOXCTR, L1, 
     .                   L2,             L3,     MAXR )

      ELSE IF ( CORSYS .EQ. RECSYS ) THEN

         CALL ZZRECBOX ( DSKDSC(MN1IDX), BOXCTR, L1, 
     .                   L2,             L3,     MAXR )

      ELSE IF ( CORSYS .EQ. PDTSYS ) THEN

         CALL ZZPDTBOX ( DSKDSC(MN1IDX), DSKDSC(PARIDX), BOXCTR, 
     .                   L1,             L2,             L3,     
     .                   MAXR                                   )

      ELSE
         
         CALL SETMSG ( 'Coordinate system # is not supported.' ) 
         CALL ERRINT ( '#',  CORSYS                            )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                   )
         CALL CHKOUT ( 'ZZSEGBOX'                              )
         RETURN

      END IF

      CALL CHKOUT ( 'ZZSEGBOX' )
      RETURN
      END
