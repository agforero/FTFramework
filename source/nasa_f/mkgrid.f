C$Procedure MKGRID ( MKDSK: create plate set from height grid )

      SUBROUTINE MKGRID ( INFILE, PLTTYP, AUNITS, DUNITS, 
     .                    CORSYS, CORPAR, MAXNV,  MAXNP, 
     .                    NV,     VERTS,  NP,     PLATES )

C$ Abstract
C
C     Create a DSK type 2 plate set from a height grid provided
C     in a file.
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
C     MKDSK
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'
      INCLUDE 'mkdsk.inc'


      CHARACTER*(*)         INFILE
      INTEGER               PLTTYP
      CHARACTER*(*)         AUNITS
      CHARACTER*(*)         DUNITS
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      INTEGER               MAXNV
      INTEGER               MAXNP
      INTEGER               NV
      DOUBLE PRECISION      VERTS  ( 3, * )
      INTEGER               NP
      INTEGER               PLATES ( 3, * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFILE     I   Name of input file.
C     PLTTYP     I   MKDSK input file format code.
C     AUNITS     I   Angular units.
C     DUNITS     I   Distance units.
C     CORSYS     I   Coordinate system.
C     CORPAR     I   Coordinate parameters.
C     MAXNV      I   Maximum number of vertices.
C     MAXNP      I   Maximum number of plates.
C     NV         O   Number of vertices.
C     VERTS      O   Vertex array.
C     NP         O   Number of plates.
C     PLATES     O   Plate array.
C
C$ Detailed_Input
C  
C     INFILE     is the name of an input data file containing height
C                grid data.
C
C     PLTTYP     is the MKDSK code indicating the format of the input
C                data file.
C
C     AUNITS     is the name of the angular unit associated with the
C                grid coordinates, if the grid coordinate system is
C                latitudinal or planetodetic. AUNITS must be supported
C                by the SPICELIB routine CONVRT.
C
C     DUNITS     is the name of the distance unit associated with the
C                grid coordinates. DUNITS must be supported by the
C                SPICELIB routine CONVRT.
C
C     CORSYS     is a DSK subsystem code designating the coordinate
C                system of the input coordinates.
C
C     CORPAR     is an array containing parameters associated with
C                the input coordinate system. The contents of the
C                array are as described in the DSK include file
C                dskdsc.inc.
C
C     COORDS     is a pair of domain coordinates: these may be,
C               
C                   - planetocentric longitude and latitude
C
C                   - planetodetic longitude and latitude
C
C                   - X and Y
C
C                For a given coordinate system, the order of the
C                elements of COORDS is that of the coordinate names in
C                the list above.
C
C     MAXNV      I   Maximum number of vertices to return.
C
C     MAXNP      I   Maximum number of plates to return.
C     
C
C$ Detailed_Output
C
C     NV         is the number of vertices in VERTS.
C
C     VERTS      is an array of 3-dimensional vertices corresponding
C                to the height grid. 
C
C                Units are km.
C
C     NP         is the number of plates in PLATES.
C
C     PLATES     is an array of plates representing a tessellation of
C                the height grid.
C
C$ Parameters
C
C     See the MKDSK include files 
C
C        mkdsk.inc
C        mkdsk02.inc
C
C     and the DSK include file
C
C        dskdsc.inc
C
C
C$ Exceptions
C
C     1)  If either the row or column count is insufficient to define a
C         surface, the error SPICE(INVALIDCOUNT) is signaled.
C
C     2)  If longitude wrap is specified for a rectangular coordinate
C         system, the error SPICE(SPURIOUSFLAG) is signaled.
C
C     3)  If either the north or south polar cap flag is .TRUE., and
C         the coordinate system is rectangular, the error
C         SPICE(SPURIOUSFLAG) is signaled.
C
C     4)  If either the row or column step is not strictly positive,
C         the error SPICE(INVALIDSTEP) is signaled.
C
C     5)  If the height scale is not is not strictly positive,
C         the error SPICE(INVALIDSCALE) is signaled.
C
C     6)  If the coordinate system is latitudinal and the height scale
C         is negative, the error SPICE(INVALIDREFVAL) is signaled.
C     
C     7)  If the number of vertices that must be created exceeds
C         MAXNV, the error SPICE(TOOMANYVERTICES) is signaled.
C     
C     8)  If the number of plates that must be created exceeds
C         MAXNP, the error SPICE(TOOMANYPLATES) is signaled.
C
C     9)  If the input file format code is not recognized, the
C         error SPICE(NOTSUPPORTED) is signaled.
C
C    10)  If an error occurs while reading the input file, the
C         error will be diagnosed by routines in the call tree
C         of this routine.
C
C    11)  If an error occurs while processing the setup file, the error
C         will be diagnosed by routines in the call tree of this
C         routine.
C
C    12)  If an error occurs while converting the input data to a
C         vertex array, the error will be diagnosed by routines in the
C         call tree of this routine.
C    
C$ Files
C
C     The file specified by INFILE can have any of the attributes (one
C     choice from each row below):
C
C        row-major  or column-major
C        top-down   or bottom-up
C        left-right or right-left
C
C     The number of tokens per line may vary. The number need have no
C     particular relationship to the row or column dimensions of the
C     output grid.
C
C     The file must contain only tokens that can be read as double
C     precision values. No non-printing characters can be present in
C     the file.
C
C     Tokens can be delimited by blanks or commas. Tokens must not be
C     split across lines.
C
C     Blank lines are allowed; however, their use is discouraged
C     because they'll cause line numbers in diagnostic messages to
C     be out of sync with actual line numbers in the file.
C
C     The file must end with a line terminator.
C
C$ Particulars
C
C     None.
C
C$ Examples
C
C     See usage in the MKDSK routine ZZWSEG02.
C
C$ Restrictions
C
C     1) For use only within program MKDSK.
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
C-    MKDSK Version 1.0.0, 25-FEB-2017 (NJB)
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      VNORM

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C
      DOUBLE PRECISION      COLSTP
      DOUBLE PRECISION      LFTCOR
      DOUBLE PRECISION      REFVAL
      DOUBLE PRECISION      ROWSTP
      DOUBLE PRECISION      HSCALE
      DOUBLE PRECISION      SUM
      DOUBLE PRECISION      TOPCOR

      INTEGER               B
      INTEGER               I
      INTEGER               J
      INTEGER               PLTBAS
      INTEGER               NCOLS
      INTEGER               NMID
      INTEGER               NNORTH
      INTEGER               NROWS
      INTEGER               NSOUTH
      INTEGER               POLIDX
      INTEGER               REQNP
      INTEGER               REQNV

      LOGICAL               LEFTRT
      LOGICAL               MKNCAP
      LOGICAL               MKSCAP
      LOGICAL               ROWMAJ
      LOGICAL               TOPDWN
      LOGICAL               WRAP
      


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'MKGRID' )

C
C     If we support the requested grid type, process the data.
C     The type is contained in the PLTTYP argument.
C     
      IF ( PLTTYP .EQ. GRID5 ) THEN
C
C        Fetch grid parameters from the kernel pool.
C
         CALL GETG05 ( CORSYS, WRAP,   MKNCAP, MKSCAP, ROWMAJ, TOPDWN,
     .                 LEFTRT, REFVAL, HSCALE, NCOLS,  NROWS,
     .                 LFTCOR, TOPCOR, COLSTP, ROWSTP         )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'MKGRID' )
            RETURN
         END IF

C
C        For safety, check the parameters that could get us into real
C        trouble.
C
         IF ( CORSYS .EQ. RECSYS ) THEN
          
            IF ( NROWS .LT. 2 ) THEN

               CALL SETMSG ( 'Number of rows was #; must have at least '
     .         //            'two rows to create a grid using the '
     .         //            'rectangular coordinate system.'         )
               CALL ERRINT ( '#', NROWS                               )
               CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                    )
               CALL CHKOUT ( 'MKGRID'                                 )
               RETURN

            END IF

            IF ( WRAP ) THEN

               CALL SETMSG ( 'Longitude wrap is not applicable to the '
     .         //            'rectangular coordinate system.'          )
               CALL SIGERR ( 'SPICE(SPURIOUSFLAG)'                     )
               CALL CHKOUT ( 'MKGRID'                                  )
               RETURN

            END IF

            IF ( MKNCAP .OR. MKSCAP ) THEN

               CALL SETMSG ( 'Polar cap creation is not applicable '
     .         //            'to the rectangular coordinate system.' )
               CALL SIGERR ( 'SPICE(SPURIOUSFLAG)'                   )
               CALL CHKOUT ( 'MKGRID'                                )
               RETURN

            END IF


         ELSE

            IF ( MKNCAP .OR. MKSCAP ) THEN

               IF ( NROWS .LT. 1 ) THEN

                  CALL SETMSG ( 'Number of rows was #; must have at  '
     .            //            'least one row to create a grid using '
     .            //            'the # coordinate system when at least '
     .            //            'one polar cap is created.'           )
                  CALL ERRINT ( '#', NROWS                            )
                  
                  IF ( CORSYS .EQ. LATSYS ) THEN

                     CALL ERRCH ( '#', 'latitudinal' )

                  ELSE IF ( CORSYS .EQ. PDTSYS ) THEN 

                     CALL ERRCH ( '#', 'planetodetic' )

                  ELSE
                     CALL ERRINT ( '#', CORSYS )
                  END IF

                  CALL SIGERR ( 'SPICE(INVALIDCOUNT)' )
                  CALL CHKOUT ( 'MKGRID'              )
                  RETURN

               END IF

            ELSE

               IF ( NROWS .LT. 2 ) THEN

                  CALL SETMSG ( 'Number of rows was #; must have at '
     .            //            'least two rows to create a grid using '
     .            //            'the # coordinate system when at no '
     .            //            'polar caps are created.'             )
                  CALL ERRINT ( '#', NROWS                            )
                  
                  IF ( CORSYS .EQ. LATSYS ) THEN

                     CALL ERRCH ( '#', 'latitudinal' )

                  ELSE IF ( CORSYS .EQ. PDTSYS ) THEN 

                     CALL ERRCH ( '#', 'planetodetic' )

                  ELSE
                     CALL ERRINT ( '#', CORSYS )
                  END IF

                  CALL SIGERR ( 'SPICE(INVALIDCOUNT)' )
                  CALL CHKOUT ( 'MKGRID'              )
                  RETURN

               END IF

            END IF

         END IF


         IF ( NCOLS .LT. 2 ) THEN

            CALL SETMSG ( 'Number of columns was #; must have at '
     .      //            'least two columns to create a grid.'   )
            CALL ERRINT ( '#', NCOLS                              )
            CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                   ) 
            CALL CHKOUT ( 'MKGRID'                                )
            RETURN

         END IF


         IF ( COLSTP .LE. 0.D0 ) THEN

            CALL SETMSG ( 'Column step must be strictly positive but '
     .      //            'was #.'                                    )
            CALL ERRDP  ( '#', COLSTP                                 )
            CALL SIGERR ( 'SPICE(INVALIDSTEP)'                        )
            CALL CHKOUT ( 'MKGRID'                                    )
            RETURN

         END IF

         IF ( ROWSTP .LE. 0.D0 ) THEN

            CALL SETMSG ( 'Row step must be strictly positive but '
     .      //            'was #.'                                    )
            CALL ERRDP  ( '#', ROWSTP                                 )
            CALL SIGERR ( 'SPICE(INVALIDSTEP)'                        )
            CALL CHKOUT ( 'MKGRID'                                    )
            RETURN

         END IF
 

         IF ( HSCALE .LE. 0.D0 ) THEN

            CALL SETMSG ( 'Height scale must be strictly positive but '
     .      //            'was #.'                                    )
            CALL ERRDP  ( '#', HSCALE                                 )
            CALL SIGERR ( 'SPICE(INVALIDSCALE)'                       )
            CALL CHKOUT ( 'MKGRID'                                    )
            RETURN

         END IF

         IF ( CORSYS .EQ. LATSYS ) THEN

            IF ( REFVAL .LT. 0.D0 ) THEN

               CALL SETMSG ( 'For latitudinal coordinates, the '
     .         //            'height reference value must be '
     .         //            'non-negative. It was #.'        )
               CALL ERRDP  ( '#', REFVAL                      )
               CALL SIGERR ( 'SPICE(INVALIDREFVAL)'           )
               CALL CHKOUT ( 'MKGRID'                         )
               RETURN

            END IF
   
         END IF


C
C        Let REQNV and REQNP be, respectively, the numbers of 
C        vertices and plates we need to create. Make sure we can handle 
C        these number.
C
         NV    = NROWS * NCOLS 

         REQNV = NV

         REQNP = 2 * ( NROWS - 1 ) * ( NCOLS - 1 )
         
         IF ( WRAP ) THEN
            REQNP = REQNP + ( 2 * ( NROWS - 1 ) )
         END IF


         IF ( MKNCAP ) THEN

            REQNV = REQNV + 1
            
            REQNP = REQNP + NCOLS

            IF ( WRAP ) THEN
               REQNP = REQNP + 1
            END IF

         END IF

         IF ( MKSCAP ) THEN

            REQNV = REQNV + 1

            REQNP = REQNP + NCOLS

            IF ( WRAP ) THEN
               REQNP = REQNP + 1
            END IF

         END IF
         
         IF ( REQNV .GT. MAXNV ) THEN

            CALL SETMSG ( 'The number of vertices that must be '
     .      //            'created is #. The maximum allowed '
     .      //            'number is #.'                        )
            CALL ERRINT ( '#', REQNV                            )
            CALL ERRINT ( '#', MAXNV                            )
            CALL SIGERR ( 'SPICE(TOOMANYVERTICES)'              )
            CALL CHKOUT ( 'MKGRID'                              )
            RETURN

         END IF
         
         IF ( REQNP .GT. MAXNP ) THEN

            CALL SETMSG ( 'The number of plates that must be '
     .      //            'created is #. The maximum allowed '
     .      //            'number is #.'                        )
            CALL ERRINT ( '#', REQNP                            )
            CALL ERRINT ( '#', MAXNP                            )
            CALL SIGERR ( 'SPICE(TOOMANYPLATES)'                )
            CALL CHKOUT ( 'MKGRID'                              )
            RETURN

         END IF
         

C
C        Create vertices. If we're making a north polar cap, leave
C        room for it at the start of the vertex array.
C        
         IF ( MKNCAP ) THEN
            B  = 2
            NV = NV + 1
         ELSE
            B  = 1
         END IF

         CALL MKVARR ( INFILE, AUNITS, DUNITS, ROWMAJ, TOPDWN,
     .                 LEFTRT, CORSYS, CORPAR, REFVAL, HSCALE,  
     .                 NCOLS,  NROWS,  LFTCOR, TOPCOR, 
     .                 COLSTP, ROWSTP, MAXNV,  VERTS(1,B)     )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'MKGRID' )
            RETURN
         END IF

C
C        The output vertices have units of km.
C
C
C        Make plates. Fill in the polar vertices, if they're needed.
C
C        We create the plates in top-down order, so the polar caps
C        will be adjacent to the nearby non-polar plates.
C
         PLTBAS = 1
         NNORTH = 0
         NSOUTH = 0
         NP     = 0

         IF ( MKNCAP ) THEN

            POLIDX = 1

            CALL ZZCAPPLT ( NCOLS, .TRUE.,  WRAP, 
     .                      PLTBAS, POLIDX, NNORTH, PLATES )

            NP = NNORTH       
C
C           The north vertex magnitude is the average of the
C           magnitudes of the vertices in the top row.
C
            SUM = 0.D0

            DO I = 1, NCOLS
               
               J   = B - 1 + I

               SUM = SUM + VNORM( VERTS(1,J) )

            END DO

            CALL VPACK ( 0.D0, 0.D0, SUM/NCOLS, VERTS(1,1) )

         END IF

C
C        Create the non-polar grid, if we have enough rows for it.
C         
         IF ( NROWS .GT. 1 ) THEN

             CALL ZZGRDPLT ( NROWS, NCOLS, WRAP, NMID, PLATES(1,NP+1) )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'MKGRID' )
               RETURN
            END IF
            
            NP = NP + NMID


            IF ( MKNCAP ) THEN
C
C              Adjust the vertex indices in the plate set.
C
               DO I = NNORTH+1, NP

                  DO J = 1, 3
                     
                     PLATES(J,I) = PLATES(J,I) + 1
                     
                  END DO

               END DO

            END IF


         ELSE
C
C           We need to make at least one polar cap, or we won't have
C           any output.
C
            IF ( ( .NOT. MKNCAP ) .AND. ( .NOT. MKSCAP ) ) THEN

               CALL SETMSG ( 'We have only one row of data in the '
     .         //            'input grid, and no polar caps were '
     .         //            'commanded to be constructed. This '
     .         //            'gives us an empty output plate set.' )
               CALL SIGERR ( 'SPICE(NOPLATES)'                     )
               CALL CHKOUT ( 'MKGRID'                              )
               RETURN

            END IF

         END IF


         IF ( MKSCAP ) THEN

            POLIDX = NV + 1
            PLTBAS = B  - 1 + ( (NROWS-1)*NCOLS )

            CALL ZZCAPPLT ( NCOLS, .FALSE., WRAP, 
     .                      PLTBAS, POLIDX, NSOUTH, PLATES(1,NP+1) )


            NP = NP + NSOUTH

C
C           The south vertex magnitude is the average of the
C           magnitudes of the vertices in the bottom row.
C
            SUM = 0.D0

            DO I = 1, NCOLS
               
               J   = PLTBAS + I

               SUM = SUM + VNORM( VERTS(1,J) )

            END DO

            CALL VPACK ( 0.D0, 0.D0, -SUM/NCOLS, VERTS(1,POLIDX) )

            NV = NV + 1
            
         END IF

      ELSE

         CALL SETMSG ( 'Input data format type is #; only type # '
     .   //            'is supported.'                            )
         CALL ERRINT ( '#', PLTTYP                                )
         CALL ERRINT ( '#', GRID5                                 )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                      )
         CALL CHKOUT ( 'MKGRID'                                   )
         RETURN

      END IF


      CALL CHKOUT ( 'MKGRID' )
      RETURN
      END

