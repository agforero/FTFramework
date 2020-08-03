C$Procedure    ZZWSEG02 ( MKDSK, write type 2 DSK segment )
 
      SUBROUTINE ZZWSEG02 ( INFILE, HANDLE )
  
C$ Abstract
C
C     Write a type 2 segment to a DSK file.
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
C     FILES
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'
      INCLUDE 'mkdsk02.inc'
      INCLUDE 'mkdsk.inc'

      CHARACTER*(*)         INFILE
      INTEGER               HANDLE
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFILE     I   Name of plate model input file.
C     HANDLE     I   DAS file handle of DSK.
C
C$ Detailed_Input
C
C     INFILE         is the name of the MKDSK input plate data file
C                    from which a type 2 DSK segment is to be created.
C                    See the MKDSK User's Guide for details.
C
C     HANDLE         is the handle of a DSK that is open for writing.
C                    The file is left open by this routine and must be
C                    closed by the calling application.
C
C$ Detailed_Output
C
C     None. 
C
C     This routine operates by side effects. See Particulars
C     for details.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     This routine is meant to be operated in RETURN SPICE error
C     handling mode. The caller is expected to delete the DSK file if
C     an error occurs during file creation.
C
C
C     1) If the setup file requests creation of a segment having
C        an unrecognized data type, the error SPICE(NOTSUPPORTED)
C        is signaled.
C
C     2) If a new DSK having the specified name cannot be created,
C        the error will be diagnosed by routines in the call tree
C        of this routine.
C
C     3) If an error occurs while writing comments to the DSK file,
C        the error will be diagnosed by routines in the call tree
C        of this routine.
C
C     4) If an error occurs while writing data to the DSK file,
C        the error will be diagnosed by routines in the call tree
C        of this routine.
C
C     5) If an error is present in the setup file, the error will
C        be diagnosed, if possible, by routines in the call tree
C        of this routine.
C
C$ Files
C
C     Input
C     -----
C     
C     1) This routine expects the input plate data file designated
C        by INFILE to conform to an expected format. The supported
C        formats are described in the MKDSK user's guide.
C
C     2) This routine assumes that a MKDSK setup file has already
C        been loaded into the kernel pool.
C
C     Output
C     ------
C       
C     This routine writes a type 2 DSK segment to a DSK file that
C     has been opened for write access. 
C
C$ Particulars
C
C     This routine parses an input plate data file, creates
C     a spatial index for the plate model, and writes a 
C     type 2 DSK segment to the file designated by HANDLE.
C
C$ Examples
C
C     See usage in MKDSK.
C
C$ Restrictions
C
C     This routine should be called only from within MKDSK.
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
C-    MKDSK Version 5.0.0, 24-FEB-2017 (NJB)
C
C        Added automatic voxel scale determination, for the
C        case where scales are not specified in the setup file.
C
C        Now includes mkdsk02.inc rather than dsk02.inc.
C     
C        Last update 17-AUG-2016 (NJB)
C
C           Now supports planetodetic and rectangular coordinates.
C           Now supports unit conversion for vertices.
C
C-    MKDSK Version 4.0.0, 30-DEC-2015 (NJB)
C
C        Re-written to make use of the SPICELIB DSK type 2 
C        spatial index creation routine DSKMI2.
C
C        The minimum radius of the plate set is now computed
C        rather than estimated.
C
C        Updated to use new ZZMKSPIN interface. This routine
C        no longer fills in the coarse voxel grid pointer array,
C        since ZZMKSPIN performs that task. 
C
C        Now imports parameter declarations for pointer and cell
C        array bounds from mkdsk.inc.
C
C
C-    MKDSK Version 3.0.0, 06-AUG-2014 (NJB)
C
C        No longer passes PID and VID arrays to RDFFPL.
C        No longer uses arrays
C
C          V1, V2, V3, X, Y, Z
C
C        No longer declares workspace array PNTRS, since this
C        array is no longer needed by ZZMKSPIN.
C
C-    MKDSK Version 2.0.0, 05-MAY-2014 (NJB)
C
C        Now does not compute the vertex plate map. The call
C        to DSKW02 sets the vertex-plate map array size to 0.
C
C        Last update was Version 2.0.0, 03-MAY-2014 (NJB)
C
C        Now has improved error checking. Calls ZZ* routines rather
C        than original versions.
C
C-    SPICELIB Version 1.0.0, 29-JUN-2010 (NJB)
C
C        Removed variable VERTS, which occupied 9*MAXPLT d.p.
C        numbers; in other words, over 0.5Gb of memory. 
C        Re-wrote calls to PLCALC and MKSPIN accordingly.
C
C        Changed dimension of VTXPTR from MAXNPV to MAXVRT;
C        given the current values of these parameters, namely
C        ~24e6 and ~4e6, this saves about 80Mb of memory.
C
C        Declared the arrays CELLS and PNTRS here rather than
C        separately in MKSPIN and VRTCOM; this saves about
C        3*MAXNPV integers, or about 290Mb of memory. Re-wrote
C        calls to MKSPIN and VRTCOM accordingly.
C
C
C-&
 
C$ Index_Entries
C
C     write type 2 segment to dsk file
C
C-&
 
C
C     SPICELIB functions
C
      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )


C 
C     Maximum size of vertex-plate list:
C
      INTEGER               MAXVTL
      PARAMETER           ( MAXVTL = MAXNPV + MAXVRT )

C
C     Size of integer part of spatial index:
C      
      INTEGER               SPXISZ
      PARAMETER           ( SPXISZ =   IXIFIX + MAXVXP + MXNVLS 
     .                               + MAXVRT + MAXVTL         )
C
C     Local variables
C
      CHARACTER*(LNSIZE)    AUNITS
      CHARACTER*(LNSIZE)    DUNITS
      CHARACTER*(FRNMLN)    FRAME
      
      DOUBLE PRECISION      AVPLEX
      DOUBLE PRECISION      CORPAR ( NSYPAR )
      DOUBLE PRECISION      DVAL
      DOUBLE PRECISION      EXTENT ( 2, 3 )
      DOUBLE PRECISION      FIRST
      DOUBLE PRECISION      LAST
      DOUBLE PRECISION      MNCOR1
      DOUBLE PRECISION      MNCOR2
      DOUBLE PRECISION      MNCOR3
      DOUBLE PRECISION      MXCOR1
      DOUBLE PRECISION      MXCOR2
      DOUBLE PRECISION      MXCOR3
      DOUBLE PRECISION      SCALE
      DOUBLE PRECISION      SPAIXD ( IXDFIX )
      DOUBLE PRECISION      VOXSCL
      DOUBLE PRECISION      VRTCES ( 3, MAXVRT )
      
      INTEGER               CELLS ( 2, MAXCEL )
      INTEGER               CENTID
      INTEGER               CGRSCL
      INTEGER               CORSYS
      INTEGER               DCLASS
      INTEGER               DTYPE
      INTEGER               I
      INTEGER               NP     
      INTEGER               NV    
      INTEGER               NVXPTR
      INTEGER               NVXTOT
      INTEGER               PLATES ( 3, MAXPLT )
      INTEGER               PLTTYP
      INTEGER               SPAIXI ( SPXISZ )
      INTEGER               SURFID
      INTEGER               TRGCOR
      INTEGER               TRGFIN
      INTEGER               VGREXT ( 3 )

      LOGICAL               MAKVPM

C
C     Saved variables
C
C     Save all variables to minimize stack overflow problems.
C     
      SAVE


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ZZWSEG02' )

C
C     Extract general DSK parameters from the setup file.
C
      CALL GETGEN ( SURFID,  CENTID,  FRAME,   FIRST,
     .              LAST,    DCLASS,  DTYPE,   AUNITS,
     .              DUNITS,  CORSYS,  CORPAR,  MNCOR1,  
     .              MXCOR1,  MNCOR2,  MXCOR2,  MAKVPM )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZWSEG02' )
         RETURN
      END IF

C
C     Check the data class.
C
      IF (  ( DCLASS .LT. SVFCLS ) .OR. ( DCLASS .GT. GENCLS ) ) THEN

         CALL SETMSG ( 'Data class was #. The only supported values '
     .   //            'are 1 and 2.'                                )
         CALL ERRINT ( '#',  DCLASS                                  )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                         )
         CALL CHKOUT ( 'ZZWSEG02'                                    )
         RETURN

      END IF

C
C     Verify the segment type.
C
      IF ( DTYPE .NE. 2 ) THEN

         CALL SETMSG ( 'Segment type must be 2 but actually was #.' )
         CALL ERRINT ( '#',  DTYPE                                  )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                        )
         CALL CHKOUT ( 'ZZWSEG02'                                   )
         RETURN

      END IF
 
C
C     Get type 2 specific parameters. If the coarse voxel scale
C     is set to zero, compute the scales.
C
      CALL GETP02 ( PLTTYP, VOXSCL, CGRSCL )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'ZZWSEG02' )
         RETURN
      END IF

C
C     Read the plate model input file.
C
      CALL TOSTDO ( 'Reading plate model input file...' )


      IF ( PLTTYP .LT. GRID5 ) THEN
C
C        The input file contains plates and vertices.
C
         CALL RDFFPL ( INFILE, PLTTYP, NV, VRTCES, NP, PLATES )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'ZZWSEG02' )
            RETURN
         END IF

C
C        Convert vertices to have units of km, if necessary.
C     
         IF (       ( .NOT. EQSTR( DUNITS, 'KM'         ) )
     .        .AND. ( .NOT. EQSTR( DUNITS, 'KILOMETERS' ) )  ) THEN

C
C           Let SCALE be the number of km equivalent to one input unit.
C
            CALL CONVRT ( 1.D0, DUNITS, 'KM', SCALE )
         
            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'ZZWSEG02' )
               RETURN
            END IF

            DO I = 1, NV

               CALL VSCLIP ( SCALE, VRTCES(1,I) )

            END DO

         END IF



      ELSE
C
C        The input file contains a height grid.
C
         CALL MKGRID ( INFILE, PLTTYP, AUNITS, DUNITS, CORSYS, CORPAR, 
     .                 MAXVRT, MAXPLT, NV,     VRTCES, NP,     PLATES )
C
C        The vertices output from MKGRID always have units of km.
C
      END IF


      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZWSEG02' )
         RETURN
      END IF



      IF ( CGRSCL .EQ. 0 ) THEN
C
C        Compute vertex and plate extents.
C
         CALL ZZPSXTNT ( NV, VRTCES, NP, PLATES, EXTENT, AVPLEX )

C
C        Compute target coarse and fine voxel counts.
C
         CALL ZZTRGNVX ( NP, TRGCOR, TRGFIN )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'ZZWSEG02' )
            RETURN
         END IF

C
C        Compute coarse and fine voxel scales.
C
         CALL ZZVOXSCL ( EXTENT, AVPLEX, TRGCOR, 
     .                   TRGFIN, CGRSCL, VOXSCL )
         
         IF ( FAILED() ) THEN

            CALL TOSTDO ( 'Could not generate voxel scales. '
     .      //            'Set voxel scales in setup file.'  ) 


            CALL CHKOUT( 'ZZWSEG02' )
            RETURN
         END IF

      END IF

C
C     Generate the spatial index.
C
      CALL TOSTDO ( 'Generating Spatial Index...' ) 




      CALL DSKMI2 ( NV,     VRTCES, NP,     PLATES, VOXSCL,
     .              CGRSCL, MAXCEL, MAXVXP, MXNVLS, MAKVPM,
     .              SPXISZ, CELLS,  SPAIXD, SPAIXI          ) 

      IF ( FAILED() ) THEN
         CALL CHKOUT( 'ZZWSEG02' )
         RETURN
      END IF 

      CALL MOVEI ( SPAIXI(SIVGRX), 3, VGREXT )

      NVXTOT = VGREXT(1)*VGREXT(2)*VGREXT(3)


      IF ( NVXTOT .GT. MAXVOX ) THEN

         CALL SETMSG ( 'NVXTOT (#) should be smaller than '
     .   //            'MAXVOX (#).'                        )
         CALL ERRINT ( '#',  NVXTOT                         )
         CALL ERRINT ( '#',  MAXVOX                         )
         CALL SIGERR ( 'SPICE(VOXELGRIDTOOBIG)'             )
         CALL CHKOUT ( 'ZZWSEG02'                           )
         RETURN

      END IF


      NVXPTR = SPAIXI(SIVXNP) 

      IF ( NVXPTR .GT. MAXVOX ) THEN

         CALL SETMSG ( 'NVXPTR (#) should be smaller than '
     .   //            'MAXVXP (#).'                        )
         CALL ERRINT ( '#',  NVXPTR                         )
         CALL ERRINT ( '#',  MAXVXP                         )
         CALL SIGERR ( 'SPICE(POINTERSETTOOBIG)'            )
         CALL CHKOUT ( 'ZZWSEG02'                           )
         RETURN

      END IF


C
C     Convert units of coordinate parameters, if necessary.
C
      IF ( ( CORSYS .EQ. LATSYS ) .OR. ( CORSYS .EQ. PDTSYS ) ) THEN
C
C        Convert input coordinate bounds to radians and km.
C
         CALL CONVRT ( MNCOR1, AUNITS, 'RADIANS', DVAL )
         MNCOR1 = DVAL

         CALL CONVRT ( MXCOR1, AUNITS, 'RADIANS', DVAL )
         MXCOR1 = DVAL

         CALL CONVRT ( MNCOR2, AUNITS, 'RADIANS', DVAL )
         MNCOR2 = DVAL

         CALL CONVRT ( MXCOR2, AUNITS, 'RADIANS', DVAL )
         MXCOR2 = DVAL

         
         IF ( CORSYS .EQ. PDTSYS ) THEN
C
C           Convert equatorial radius to km.
C
            CALL CONVRT ( CORPAR(PARIDX), DUNITS, 'KM', DVAL )
            CORPAR(PARIDX) = DVAL

         END IF


      ELSE IF ( CORSYS .EQ. RECSYS ) THEN

         CALL CONVRT ( MNCOR1, DUNITS, 'KM', DVAL )
         MNCOR1 = DVAL

         CALL CONVRT ( MXCOR1, DUNITS, 'KM', DVAL )
         MXCOR1 = DVAL

         CALL CONVRT ( MNCOR2, DUNITS, 'KM', DVAL )
         MNCOR2 = DVAL

         CALL CONVRT ( MXCOR2, DUNITS, 'KM', DVAL )
         MXCOR2 = DVAL

      ELSE

         CALL SETMSG ( 'Coordinate system # is not supported.' )
         CALL ERRINT ( '#', CORSYS                             )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                   )
         CALL CHKOUT ( 'ZZWSEG02'                              )
         RETURN

      END IF

C     
C     Compute bounds on the third coordinate.
C
      CALL DSKRB2 ( NV,     VRTCES, NP,     PLATES,
     .              CORSYS, CORPAR, MNCOR3, MXCOR3 )

C
C     Write a type 2 (plate model) segment to the DSK file.
C
      CALL DSKW02 ( HANDLE, 
     .              CENTID, SURFID, DCLASS, FRAME,  CORSYS,
     .              CORPAR, MNCOR1, MXCOR1, MNCOR2, MXCOR2,
     .              MNCOR3, MXCOR3, FIRST,  LAST,   NV,      
     .              VRTCES, NP,     PLATES, SPAIXD, SPAIXI ) 


      CALL CHKOUT ( 'ZZWSEG02' )
      RETURN
      END
