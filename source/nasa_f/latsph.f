C$Procedure      LATSPH ( Latitudinal to spherical coordinates )
 
      SUBROUTINE LATSPH ( RADIUS, LONG, LAT,  RHO, COLAT, LONGS )
 
C$ Abstract
C
C     Convert from latitudinal coordinates to spherical coordinates.
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
C     CONVERSION,  COORDINATES
C
C$ Declarations
 
      DOUBLE PRECISION   RADIUS
      DOUBLE PRECISION   LONG
      DOUBLE PRECISION   LAT
      DOUBLE PRECISION   RHO
      DOUBLE PRECISION   COLAT
      DOUBLE PRECISION   LONGS
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     RADIUS     I   Distance of a point from the origin.
C     LONG       I   Angle of the point from the XZ plane in radians.
C     LAT        I   Angle of the point from the XY plane in radians.
C     RHO        O   Distance of the point from the origin.
C     COLAT      O   Angle of the point from positive Z axis (radians).
C     LONGS      O   Angle of the point from the XZ plane (radians).
C
C$ Detailed_Input
C
C     RADIUS     Distance of a point from the origin.
C
C     LONG       Angle of the point from the XZ plane in radians.
C
C     LAT        Angle of the point from the XY plane in radians.
C
C$ Detailed_Output
C
C     RHO        Distance of the point from the origin.
C
C     COLAT      Angle between the vector from the origin to the point
C                and the positive Z axis in radians. COLAT is computed
C                as PI/2 - LAT.
C
C     LONGS      Angle of the point from the XZ plane (radians). LONGS
C                is set equal to LONG.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     Error free.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine returns the spherical coordinates of a point
C     whose position is input in latitudinal coordinates.
C
C     Latitudinal coordinates are defined by a distance from a central
C     reference point, an angle from a reference meridian, and an angle
C     above the equator of a sphere centered at the central reference
C     point.
C
C     Spherical coordinates are defined by a distance from a central
C     reference point, an angle from a reference meridian, and an angle
C     from the z-axis.
C
C$ Examples
C
C     Co-latitude is obtained by subtracting latitude from HALFPI()
C     Radius and longitude mean the same thing in both latitudinal
C     and spherical coordinates.  The table below lists LAT
C     corresponding COLAT in terms of degrees.
C
C             LAT            COLAT
C            ------         ------
C              0             90
C             20             70
C             45             45
C            -30            120
C             90              0
C            -45            135
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
C     W.L. Taber      (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.2, 26-JUL-2016 (BVS)
C
C        Minor headers edits.
C
C-    SPICELIB Version 1.0.1, 10-MAR-1992 (WLT)
C
C        Comment section for permuted index source lines was added
C        following the header.
C
C-    SPICELIB Version 1.0.0, 31-JAN-1990 (WLT)
C
C-&
 
C$ Index_Entries
C
C     latitudinal to spherical coordinates
C
C-&
 
 
C$ Revisions
C
C-    Beta Version 1.0.1, 1-Feb-1989 (WLT)
C
C        Example section of header upgraded.
C
C-&
 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      HALFPI
 
C
C     Local Variables
C
      DOUBLE PRECISION TH
      DOUBLE PRECISION PH
 
C
C     Convert to spherical coordinates, storing the results in
C     temporary variables
C
      TH = HALFPI() - LAT
      PH = LONG
 
C
C     Move results to output variables
C
      RHO   = RADIUS
      COLAT = TH
      LONGS = PH
 
      RETURN
      END
