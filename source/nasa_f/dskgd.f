C$Procedure DSKGD ( DSK, return DSK segment descriptor  )
 
      SUBROUTINE DSKGD ( HANDLE, DLADSC, DSKDSC )
      IMPLICIT NONE
  
C$ Abstract
C
C     Return the DSK descriptor from a DSK segment identified
C     by a DAS handle and DLA descriptor.
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
C     NAIF_IDS
C
C$ Keywords
C
C     DAS
C     DSK
C     FILES
C     TOPOGRAPHY
C
C$ Declarations

      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'

      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      DOUBLE PRECISION      DSKDSC ( * )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of a DSK file.  
C     DLADSC     I   DLA segment descriptor.
C     DSKDSC     O   DSK segment descriptor.
C
C$ Detailed_Input
C
C     HANDLE         is the handle of a DSK file that is open for
C                    read access.
C
C     DLADSC         is the DLA segment descriptor corresponding to
C                    a DSK segment.
C
C$ Detailed_Output
C
C     DSKDSC         is the DSK segment descriptor of the segment
C                    designated by the input handle and DLA descriptor.
C
C$ Parameters
C
C     See the INCLUDE file
C
C         dskdsc.inc
C
C$ Exceptions
C
C     1) If the size of the double precision component of the 
C        segment is smaller than that of a DSK descriptor, the
C        error SPICE(INVALIDFORMAT) is signaled.
C
C     2) DAS read errors will be diagnosed by routines in the
C        call tree of this routine.
C
C$ Files
C
C     See input argument HANDLE.
C
C$ Particulars
C
C     This is a convenience routine intended for use by low-level
C     routines that read DSK segments. This routine may also be called
C     by user applications that must access DSK files at the segment
C     level.
C     
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as input,
C     the compiler and supporting libraries, and the machine specific
C     arithmetic implementation.
C
C     1) Dump the DSK descriptors of a DSK file.
C 
C        Example code begins here.
C
C
C              PROGRAM EX1
C              IMPLICIT NONE
C
C              INCLUDE 'dla.inc'
C              INCLUDE 'dskdsc.inc'
C              INCLUDE 'dsk02.inc'
C
C              INTEGER               FILSIZ
C              PARAMETER           ( FILSIZ = 255 )
C
C              CHARACTER*(FILSIZ)    DSK
C
C              DOUBLE PRECISION      DSKDSC ( DSKDSZ )
C
C              INTEGER               DLADSC ( DLADSZ )
C              INTEGER               HANDLE
C              INTEGER               I
C              INTEGER               NXTDSC ( DLADSZ )
C
C              LOGICAL               FOUND
C
C
C              CALL PROMPT ( 'Enter DSK name > ', DSK )
C        C
C        C     Open the DSK file and begin a forward search
C        C     for segments.
C        C
C              CALL DASOPR ( DSK, HANDLE )
C
C              CALL DLABFS ( HANDLE, NXTDSC, FOUND )
C
C              DO WHILE ( FOUND )
C        C
C        C        Make the DLA descriptor we just fetched
C        C        the current one.
C        C
C                 CALL MOVEI ( NXTDSC, DLADSZ, DLADSC )
C
C                 CALL DSKGD ( HANDLE, DLADSC, DSKDSC )
C
C                 WRITE (*,*) 'DSK descriptor contents: '
C
C                 DO I = 1, DSKDSZ
C                    WRITE (*,*) DSKDSC(I)
C                 END DO
C        C
C        C        Find the next segment, if it exists.
C        C
C                 CALL DLAFNS ( HANDLE, DLADSC, NXTDSC, FOUND )
C
C              END DO
C
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, and when the input DSK name was 
C
C         phobos512.bds
C
C     (this file is available on the NAIF server), the output was:
C
C
C        Enter DSK name > phobos512.bds
C         DSK descriptor contents:
C           401.00000000000000
C           401.00000000000000
C          1.00000000000000000
C           2.0000000000000000
C           10021.000000000000
C          1.00000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C           0.0000000000000000
C          -3.1415926535897931
C           3.1415926535897931
C          -1.5707963267948966
C           1.5707963267948966
C           8.0496322487215526
C           13.940939832123945
C          -1577879958.8160586
C           1577880066.1839132
C
C
C     2) Again, dump the DSK descriptors of a DSK file, this time
C        interpreting the descriptor information and displaying 
C        it in a user-friendly form. This display is a simplified
C        version of that created by the utility DSKBRIEF.
C
C        This program requests the name of an optional meta-kernel.
C        The meta-kernel can be used to define surface name-ID 
C        associations. If no meta-kernel is needed, the user can
C        enter a carriage return at the prompt for this file.
C
C
C        Example code begins here.
C
C
C               PROGRAM EX2
C               IMPLICIT NONE
C
C               INCLUDE 'dla.inc'
C               INCLUDE 'dskdsc.inc'
C               INCLUDE 'dsk02.inc'
C               INCLUDE 'srftrn.inc'
C         C
C         C     SPICELIB functions
C         C
C               DOUBLE PRECISION      DPR
C         C
C         C     Local parameters
C         C
C               CHARACTER*(*)         FMT1
C               PARAMETER           ( FMT1 = '(A,2(F21.14))' )
C
C               CHARACTER*(*)         FMT2
C               PARAMETER           ( FMT2 = '(A,I3)' )
C
C               INTEGER               BDNMLN
C               PARAMETER           ( BDNMLN = 36 )
C
C               INTEGER               FILSIZ
C               PARAMETER           ( FILSIZ = 255 )
C
C               INTEGER               FRNMLN
C               PARAMETER           ( FRNMLN = 32 )
C
C               INTEGER               NAMLEN
C               PARAMETER           ( NAMLEN = 30 )
C
C               INTEGER               TIMLEN
C               PARAMETER           ( TIMLEN = 40 )
C
C               INTEGER               NSYS
C               PARAMETER           ( NSYS   = 4 )
C
C               INTEGER               NCLASS
C               PARAMETER           ( NCLASS = 2 )
C
C               INTEGER               CLNMLN
C               PARAMETER           ( CLNMLN = 25 )
C         C
C         C     Local variables
C         C
C               CHARACTER*(BDNMLN)    BODNAM
C               CHARACTER*(TIMLEN)    BTIME
C               CHARACTER*(CLNMLN)    CLSNMS ( NCLASS )
C               CHARACTER*(FILSIZ)    DSK
C               CHARACTER*(TIMLEN)    ETIME
C               CHARACTER*(FRNMLN)    FRAME
C               CHARACTER*(FILSIZ)    META
C               CHARACTER*(SFNMLN)    SRFNAM
C               CHARACTER*(NAMLEN)    SYSNAM
C               CHARACTER*(NAMLEN)    SYSNMS ( NSYS )
C
C               DOUBLE PRECISION      DSKDSC ( DSKDSZ )
C               DOUBLE PRECISION      F
C               DOUBLE PRECISION      RE
C               DOUBLE PRECISION      RP
C
C               INTEGER               BODYID
C               INTEGER               CORSYS
C               INTEGER               DCLASS
C               INTEGER               DLADSC ( DLADSZ )
C               INTEGER               DTYPE
C               INTEGER               FRAMID
C               INTEGER               HANDLE
C               INTEGER               NXTDSC ( DLADSZ )
C               INTEGER               SEGNO
C               INTEGER               SURFID
C
C               LOGICAL               FOUND
C               LOGICAL               ISNAME
C         C
C         C     Initial values
C         C
C               DATA                  CLSNMS / 'Single-valued surface',
C              .                               'General surface'       /
C
C               DATA                  SYSNMS / 'Latitudinal',
C              .                               'Cylindrical',
C              .                               'Rectangular',
C              .                               'Planetodetic' /
C
C
C               CALL PROMPT ( 'Enter DSK name         > ', DSK  )
C               CALL PROMPT ( 'Enter meta-kernel name > ', META )
C
C               IF ( META .NE. ' ' ) THEN
C                  CALL FURNSH ( META )
C               END IF
C         C
C         C     Open the DLA file and begin a forward search
C         C     for segments.
C         C
C               CALL DASOPR ( DSK, HANDLE )
C
C               SEGNO = 0
C
C               CALL DLABFS ( HANDLE, NXTDSC, FOUND )
C
C               DO WHILE ( FOUND )
C
C                  SEGNO = SEGNO + 1
C         C
C         C        Make the DLA descriptor we just fetched
C         C        the current one.
C         C
C                  CALL MOVEI ( NXTDSC, DLADSZ, DLADSC )
C
C                  CALL DSKGD ( HANDLE, DLADSC, DSKDSC )
C
C                  BODYID = NINT( DSKDSC(CTRIDX) )
C                  SURFID = NINT( DSKDSC(SRFIDX) )
C                  FRAMID = NINT( DSKDSC(FRMIDX) )
C                  DTYPE  = NINT( DSKDSC(TYPIDX) )
C                  DCLASS = NINT( DSKDSC(CLSIDX) )
C
C                  CALL BODC2S ( BODYID, BODNAM )
C                  CALL SRFC2S ( SURFID, BODYID, SRFNAM, ISNAME )
C                  CALL FRMNAM ( FRAMID, FRAME  )
C
C                  IF ( FRAME .EQ. ' ' ) THEN
C                     CALL INTSTR ( FRAMID, FRAME )
C                  END IF
C
C                  CALL ETCAL ( DSKDSC(BTMIDX), BTIME )
C                  CALL ETCAL ( DSKDSC(ETMIDX), ETIME )
C
C                  CORSYS = NINT( DSKDSC(SYSIDX) )
C
C                  SYSNAM = SYSNMS( CORSYS )
C
C                  WRITE (*,*)    '===================================='
C                  WRITE (*,FMT2) ' DSK descriptor for segment ',
C              .                  SEGNO
C                  WRITE (*,*)    '  Body:              ', BODNAM
C                  WRITE (*,*)    '  Surface:           ', SRFNAM
C                  WRITE (*,*)    '  Frame:             ', FRAME
C                  WRITE (*,*)    '  Start time (TDB):  ', BTIME
C                  WRITE (*,*)    '  Stop time  (TDB):  ', ETIME
C                  WRITE (*,*)    '  Data type:         ', DTYPE
C                  WRITE (*,*)    '  Data class:        ', DCLASS, ' ',
C              .                                        CLSNMS(DCLASS)
C                  WRITE (*,*)    '  Coordinate system: ', SYSNAM
C
C                  IF ( CORSYS .EQ. PDTSYS ) THEN
C
C                     RE = DSKDSC(PARIDX  )
C                     F  = DSKDSC(PARIDX+1)
C                     RP = RE * ( 1.D0 - F )
C
C                     WRITE (*,*) '     Equatorial radius (km): ', RE
C                     WRITE (*,*) '     Polar radius      (km): ', RP
C
C                  END IF
C
C                  WRITE (*,*) '  Segment boundaries:'
C
C                  IF ( CORSYS .EQ. LATSYS ) THEN
C
C                     WRITE (*,FMT1) '    Longitude (deg):   ',
C              .                  DPR() * DSKDSC(MN1IDX),
C              .                  DPR() * DSKDSC(MX1IDX)
C                     WRITE (*,FMT1) '    Latitude  (deg):   ',
C              .                  DPR() * DSKDSC(MN2IDX),
C              .                  DPR() * DSKDSC(MX2IDX)
C                     WRITE (*,FMT1) '    Radius     (km):   ',
C              .                          DSKDSC(MN3IDX),
C              .                          DSKDSC(MX3IDX)
C
C                  ELSE IF ( CORSYS .EQ. CYLSYS ) THEN
C
C                     CALL SETMSG ( 'Coordinate system was '
C              .      //            'Cylindrical'           )
C                     CALL SIGERR ( 'SPICE(NOTSUPPORTED)'   )
C
C
C                  ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C                     WRITE (*,FMT1) '    X-coordinate (km): ',
C              .                          DSKDSC(MN1IDX),
C              .                          DSKDSC(MX1IDX)
C                     WRITE (*,FMT1) '    Y-coordinate (km): ',
C              .                          DSKDSC(MN2IDX),
C              .                          DSKDSC(MX2IDX)
C                     WRITE (*,FMT1) '    Z-coordinate (km): ',
C              .                          DSKDSC(MN3IDX),
C              .                          DSKDSC(MX3IDX)
C
C                  ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C                     WRITE (*,FMT1) '    Longitude (deg):   ',
C              .                  DPR() * DSKDSC(MN1IDX),
C              .                  DPR() * DSKDSC(MX1IDX)
C                     WRITE (*,FMT1) '    Latitude  (deg):   ',
C              .                  DPR() * DSKDSC(MN2IDX),
C              .                  DPR() * DSKDSC(MX2IDX)
C                     WRITE (*,FMT1) '    Altitude   (km):   ',
C              .                          DSKDSC(MN3IDX),
C              .                          DSKDSC(MX3IDX)
C                  END IF
C         C
C         C        Find the next segment, if it exists.
C         C
C                  CALL DLAFNS ( HANDLE, DLADSC, NXTDSC, FOUND )
C
C               END DO
C
C               END
C
C
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, and when the input DSK name was 
C
C         phobos512.bds
C
C     (this file is available on the NAIF server), the output was:
C 
C
C      Enter DSK name         > phobos512.bds
C      Enter meta-kernel name >
C       ====================================
C       DSK descriptor for segment   1
C         Body:              PHOBOS
C         Surface:           401
C         Frame:             IAU_PHOBOS
C         Start time (TDB):  1950 JAN 01 00:00:41.183
C         Stop time  (TDB):  2050 JAN 01 00:01:06.183
C         Data type:                    2
C         Data class:                   1  Single-valued surface
C         Coordinate system: Latitudinal
C         Segment boundaries:
C          Longitude (deg):     -180.00000000000000   180.00000000000000
C          Latitude  (deg):      -90.00000000000000    90.00000000000000
C          Radius     (km):        8.04963224872155    13.94093983212395
C
C 
C
C     3) Again, dump the DSK descriptors of a DSK file, using the 
C        program from example 2, but this time reading the DSK file
C
C           phobos_3_3_3seg.bds
C
C        which can be created by running an example program from
C        DSKW02. Use the meta-kernel below to demonstrate surface
C        name-ID mapping:
C
C
C           KPL/MK
C
C           File: dskgd_ex3.tm
C
C           This meta-kernel is intended to support operation of SPICE
C           example programs. The file contents shown here should not be
C           assumed to contain adequate or correct versions of data
C           required by SPICE-based user applications.
C
C
C           \begindata
C
C           NAIF_SURFACE_NAME += ( 'Phobos example surface 1',
C                                  'Phobos example surface 2',
C                                  'Phobos example surface 3' )
C           NAIF_SURFACE_CODE += (   1,   2,   3 )
C           NAIF_SURFACE_BODY += ( 401, 401, 401 )
C
C           \begintext
C
C
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, using the example DSK named above, the output was:
C
C
C      Enter DSK name         > phobos_3_3_3seg.bds
C      Enter meta-kernel name > dskgd_ex3.tm
C       ====================================
C       DSK descriptor for segment   1
C         Body:              PHOBOS
C         Surface:           Phobos example surface 1
C         Frame:             IAU_PHOBOS
C         Start time (TDB):  1950 JAN 01 00:00:00.000
C         Stop time  (TDB):  2050 JAN 01 00:00:00.000
C         Data type:                    2
C         Data class:                   2  General surface
C         Coordinate system: Latitudinal
C         Segment boundaries:
C          Longitude (deg):     -180.00000000000000   180.00000000000000
C          Latitude  (deg):      -90.00000000000000    90.00000000000000
C          Radius     (km):        8.22529807597397    14.01176814562576
C       ====================================
C       DSK descriptor for segment   2
C         Body:              PHOBOS
C         Surface:           Phobos example surface 2
C         Frame:             IAU_PHOBOS
C         Start time (TDB):  1950 JAN 01 00:00:00.000
C         Stop time  (TDB):  2050 JAN 01 00:00:00.000
C         Data type:                    2
C         Data class:                   2  General surface
C         Coordinate system: Rectangular
C         Segment boundaries:
C          X-coordinate (km):     -1.30000000000000     1.31000000000000
C          Y-coordinate (km):     -1.21000000000000     1.20000000000000
C          Z-coordinate (km):     -9.45293235778800     9.63817977905300
C       ====================================
C       DSK descriptor for segment   3
C         Body:              PHOBOS
C         Surface:           Phobos example surface 3
C         Frame:             IAU_PHOBOS
C         Start time (TDB):  1950 JAN 01 00:00:00.000
C         Stop time  (TDB):  2050 JAN 01 00:00:00.000
C         Data type:                    2
C         Data class:                   2  General surface
C         Coordinate system: Planetodetic
C            Equatorial radius (km):    13.000000000000000
C            Polar radius      (km):    9.0999999999999996
C         Segment boundaries:
C          Longitude (deg):     -180.00000000000000   180.00000000000000
C          Latitude  (deg):      -90.00000000000000    90.00000000000000
C          Altitude   (km):       -3.72866868360370     1.37201579108146
C
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
C        Updated version info.
C
C        22-JAN-2016 (NJB)
C
C           Added new header example programs and updated existing
C           example program. Made minor changes to code to enhance
C           readability. Corrected header typo.
C
C        09-OCT-2009 (NJB)
C
C-&
 
C$ Index_Entries
C
C     return dsk segment descriptor
C
C-&
 

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local variables
C
      INTEGER               DPBASE
      INTEGER               DPSIZE
      
      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKGD' )

C
C     Fetch the base address and size of the DP component of the
C     indicated segment.
C
      DPBASE = DLADSC ( DBSIDX )
      DPSIZE = DLADSC ( DSZIDX )

C
C     If we don't have enough d.p. elements to hold a descriptor,
C     something's wrong.
C
      IF ( DPSIZE .LT. DSKDSZ ) THEN

         CALL SETMSG ( 'Size of d.p. component of segment is #; '   
     .   //            'cannot extract descriptor.  This is a '     
     .   //            'file format error which may be indicative ' 
     .   //             'of a corrupted file.'                      )
         CALL ERRINT ( '#', DPSIZE                                  )
         CALL SIGERR ( 'SPICE(INVALIDFORMAT)'                       )
         CALL CHKOUT ( 'DSKGD'                                      )
         RETURN

      END IF

C
C     Extract the descriptor.
C     
      CALL DASRDD ( HANDLE, DPBASE+1, DPBASE+DSKDSZ, DSKDSC )

      CALL CHKOUT ( 'DSKGD' )
      RETURN
      END

