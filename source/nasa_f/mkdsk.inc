C$ Abstract
C
C     Include Section:  MKDSK Global Parameters
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
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    Version 4.0.0, 28-FEB-2017 (NJB)
C
C        Added declaration of version string VER. Previously
C        this declaration was located in prcinf.for.
C
C        Declarations of parameters
C
C           MAXCEL
C           MAXVXP
C           MAXNVLS
C
C        were moved to dsk02.inc.
C
C        Declarations of parameters
C
C           ZZMAXV
C           ZZMAXP
C
C        were deleted.
C     
C
C-    Version 3.0.0, 20-OCT-2015 (NJB)
C
C        Parameter MAXQ was increased from 512 to 1024.
C
C-    Version 2.0.0, 26-MAR-2015 (NJB)
C
C        Declarations were added for the parameters
C
C           MAXCEL
C           MAXVXP
C           MXNVLS
C           
C
C-    Version 1.0.0, 04-MAY-2010 (NJB)
C
C-&

C
C     MKDSK version:
C
      CHARACTER*(*)         VER
      PARAMETER           ( VER = '2.0.0, 28-FEB-2017' )

C
C     Default time format:
C
      CHARACTER*(*)         DFTFMT
      PARAMETER           ( DFTFMT =  
     .                      'YYYY-MM-DDTHR:MN:SC.######## ::TDB' )

C
C     Command line length:
C
      INTEGER               CMDLEN
      PARAMETER           ( CMDLEN = 2000 )
      
C
C     SPICELIB cell lower bound:
C
      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )
      
C
C     Maximum file name length:
C
      INTEGER               FILSIZ
      PARAMETER           ( FILSIZ = 255 )
      
C
C     Output file line length:
C
      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 255 )

C
C     Length of string for keyword value processing:
C
      INTEGER               SHRTLN
      PARAMETER           ( SHRTLN = 80 ) 


C
C     The maximum 'Q' value for Gaskell shape models. This value
C     is always a multiple of 4.
C
      INTEGER               MAXQ
      PARAMETER           ( MAXQ = 1024 )

C
C     Input format parameters
C

C
C     The height grid "plate type" is the fifth format.
C
      INTEGER               GRID5
      PARAMETER           ( GRID5 = 5 )

C
C     End Include Section:  MKDSK Global Parameters 
C
