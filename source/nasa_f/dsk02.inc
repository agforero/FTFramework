C
C     Include file dsk02.inc
C
C     This include file declares parameters for DSK data type 2
C     (plate model).
C
C-       SPICELIB Version 1.0.0 08-FEB-2017 (NJB)
C
C          Updated version info.
C     
C           22-JAN-2016 (NJB)
C
C              Now includes spatial index parameters.
C
C           26-MAR-2015 (NJB)
C     
C              Updated to increase MAXVRT to 16000002. MAXNPV
C              has been changed to (3/2)*MAXPLT. Set MAXVOX
C              to 100000000.
C
C           13-MAY-2010 (NJB)
C
C              Updated to reflect new no-record design.
C
C           04-MAY-2010 (NJB)
C
C              Updated for new type 2 segment design. Now uses
C              a local parameter to represent DSK descriptor
C              size (NB).
C
C           13-SEP-2008 (NJB)
C
C              Updated to remove albedo information.
C              Updated to use parameter for DSK descriptor size.
C
C           27-DEC-2006 (NJB)
C
C              Updated to remove minimum and maximum radius information
C              from segment layout.  These bounds are now included
C              in the segment descriptor.
C
C           26-OCT-2006 (NJB)
C
C              Updated to remove normal, center, longest side, albedo,
C              and area keyword parameters.
C
C           04-AUG-2006 (NJB)
C
C              Updated to support coarse voxel grid.  Area data
C              have now been removed.
C
C           10-JUL-2006 (NJB)
C
C
C     Each type 2 DSK segment has integer, d.p., and character
C     components.  The segment layout in DAS address space is as
C     follows:
C
C
C        Integer layout:
C
C           +-----------------+
C           | NV              |  (# of vertices)
C           +-----------------+
C           | NP              |  (# of plates )
C           +-----------------+
C           | NVXTOT          |  (total number of voxels)
C           +-----------------+
C           | VGREXT          |  (voxel grid extents, 3 integers)   
C           +-----------------+
C           | CGRSCL          |  (coarse voxel grid scale, 1 integer)  
C           +-----------------+
C           | VOXNPT          |  (size of voxel-plate pointer list)
C           +-----------------+
C           | VOXNPL          |  (size of voxel-plate list)
C           +-----------------+
C           | VTXNPL          |  (size of vertex-plate list)
C           +-----------------+
C           | PLATES          |  (NP 3-tuples of vertex IDs)
C           +-----------------+
C           | VOXPTR          |  (voxel-plate pointer array)
C           +-----------------+
C           | VOXPLT          |  (voxel-plate list)
C           +-----------------+
C           | VTXPTR          |  (vertex-plate pointer array)
C           +-----------------+
C           | VTXPLT          |  (vertex-plate list)
C           +-----------------+
C           | CGRPTR          |  (coarse grid occupancy pointers)
C           +-----------------+
C
C
      INTEGER               IXNV
      PARAMETER           ( IXNV   =          1 )
      
      INTEGER               IXNP
      PARAMETER           ( IXNP   = IXNV   + 1 )

      INTEGER               IXNVXT 
      PARAMETER           ( IXNVXT = IXNP   + 1 )
      
      INTEGER               IXVGRX
      PARAMETER           ( IXVGRX = IXNVXT + 1 )

      INTEGER               IXCGSC
      PARAMETER           ( IXCGSC = IXVGRX + 3 )

      INTEGER               IXVXPS
      PARAMETER           ( IXVXPS = IXCGSC + 1 )

      INTEGER               IXVXLS
      PARAMETER           ( IXVXLS = IXVXPS + 1 )

      INTEGER               IXVTLS
      PARAMETER           ( IXVTLS = IXVXLS + 1 )

      INTEGER               IXPLAT
      PARAMETER           ( IXPLAT = IXVTLS + 1 )

C
C        D.p. layout:
C            
C           +-----------------+
C           | DSK descriptor  |  DSKDSZ elements
C           +-----------------+
C           | Vertex bounds   |  6 values (min/max for each component)
C           +-----------------+
C           | Voxel origin    |  3 elements
C           +-----------------+
C           | Voxel size      |  1 element
C           +-----------------+
C           | Vertices        |  3*NV elements
C           +-----------------+
C
      INTEGER               IXDSCR
      PARAMETER           ( IXDSCR = 1 )
      
C
C     This local parameter MUST be kept consistent with
C     the parameter DSKDSZ which is declared in dskdsc.inc.
C
      INTEGER               DSCSZ2
      PARAMETER           ( DSCSZ2 = 24 )

      INTEGER               IXVTBD
      PARAMETER           ( IXVTBD = IXDSCR + DSCSZ2 )

      INTEGER               IXVXOR
      PARAMETER           ( IXVXOR = IXVTBD + 6 )

      INTEGER               IXVXSZ
      PARAMETER           ( IXVXSZ = IXVXOR + 3 )

      INTEGER               IXVERT
      PARAMETER           ( IXVERT = IXVXSZ + 1 )

 

C
C     Integer item keyword parameters used by fetch routines:
C
      INTEGER               KWNV
      PARAMETER           ( KWNV   = 1 )

      INTEGER               KWNP
      PARAMETER           ( KWNP   = KWNV   + 1 )
      
      INTEGER               KWNVXT
      PARAMETER           ( KWNVXT = KWNP   + 1 )

      INTEGER               KWVGRX
      PARAMETER           ( KWVGRX = KWNVXT + 1 )

      INTEGER               KWCGSC
      PARAMETER           ( KWCGSC = KWVGRX + 1 )

      INTEGER               KWVXPS
      PARAMETER           ( KWVXPS = KWCGSC + 1 )

      INTEGER               KWVXLS
      PARAMETER           ( KWVXLS = KWVXPS + 1 )

      INTEGER               KWVTLS
      PARAMETER           ( KWVTLS = KWVXLS + 1 )

      INTEGER               KWPLAT
      PARAMETER           ( KWPLAT = KWVTLS + 1 )

      INTEGER               KWVXPT
      PARAMETER           ( KWVXPT = KWPLAT + 1 )

      INTEGER               KWVXPL
      PARAMETER           ( KWVXPL = KWVXPT + 1 )

      INTEGER               KWVTPT
      PARAMETER           ( KWVTPT = KWVXPL + 1 )

      INTEGER               KWVTPL
      PARAMETER           ( KWVTPL = KWVTPT + 1 )

      INTEGER               KWCGPT
      PARAMETER           ( KWCGPT = KWVTPL + 1 )

C
C     Double precision item keyword parameters used by fetch routines:
C
      INTEGER               KWDSC
      PARAMETER           ( KWDSC  = KWCGPT + 1 )

      INTEGER               KWVTBD
      PARAMETER           ( KWVTBD = KWDSC  + 1 )

      INTEGER               KWVXOR
      PARAMETER           ( KWVXOR = KWVTBD + 1 )

      INTEGER               KWVXSZ
      PARAMETER           ( KWVXSZ = KWVXOR + 1 )

      INTEGER               KWVERT
      PARAMETER           ( KWVERT = KWVXSZ + 1 )

C
C     The parameters below formerly were declared in pltmax.inc.
C
C     Limits on plate model capacity:
C
C     The maximum number of bodies, vertices and
C     plates in a plate model or collective thereof are
C     provided here.
C
C     These values can be used to dimension arrays, or to
C     use as limit checks.
C
C     The value of MAXPLT is determined from MAXVRT via
C     Euler's Formula for simple polyhedra having triangular
C     faces.
C
C     MAXVRT is the maximum number of vertices the triangular
C            plate model software will support.
C
      INTEGER               MAXVRT
      PARAMETER           ( MAXVRT = 16000002  )
 
C
C     MAXPLT is the maximum number of plates that the triangular
C            plate model software will support.
C
      INTEGER               MAXPLT
      PARAMETER           ( MAXPLT = 2*(MAXVRT - 2) )

C
C     MAXNPV is the maximum allowed number of vertices, not taking into
C     account shared vertices.
C
C     Note that this value is not sufficient to create a vertex-plate
C     mapping for a model of maximum plate count.
C
      INTEGER               MAXNPV
      PARAMETER           ( MAXNPV = 3*MAXPLT/2 + 1 )

C
C     MAXVOX is the maximum number of voxels.
C
      INTEGER               MAXVOX
      PARAMETER           ( MAXVOX = 100000000 )

C
C     MAXCGR is the maximum size of the coarse voxel grid.
C     
      INTEGER               MAXCGR
      PARAMETER           ( MAXCGR = 100000 )

C
C     MAXEDG is the maximum allowed number of vertex or plate 
C     neighbors a vertex may have. 
C
      INTEGER               MAXEDG
      PARAMETER           ( MAXEDG = 120 )



C     DSK type 2 spatial index parameters
C     ===================================
C
C        DSK type 2 spatial index integer component
C        ------------------------------------------
C
C           +-----------------+
C           | VGREXT          |  (voxel grid extents, 3 integers)   
C           +-----------------+
C           | CGRSCL          |  (coarse voxel grid scale, 1 integer)  
C           +-----------------+
C           | VOXNPT          |  (size of voxel-plate pointer list)
C           +-----------------+
C           | VOXNPL          |  (size of voxel-plate list)
C           +-----------------+
C           | VTXNPL          |  (size of vertex-plate list)
C           +-----------------+
C           | CGRPTR          |  (coarse grid occupancy pointers)
C           +-----------------+
C           | VOXPTR          |  (voxel-plate pointer array)
C           +-----------------+
C           | VOXPLT          |  (voxel-plate list)
C           +-----------------+
C           | VTXPTR          |  (vertex-plate pointer array)
C           +-----------------+
C           | VTXPLT          |  (vertex-plate list)
C           +-----------------+
C
C
C        Index parameters
C
C
C     Grid extent:
C
      INTEGER               SIVGRX
      PARAMETER           ( SIVGRX = 1 )
C
C     Coarse grid scale:
C
      INTEGER               SICGSC
      PARAMETER           ( SICGSC = SIVGRX + 3 )
C
C     Voxel pointer count:
C
      INTEGER               SIVXNP
      PARAMETER           ( SIVXNP = SICGSC + 1 )
C
C     Voxel-plate list count:
C
      INTEGER               SIVXNL
      PARAMETER           ( SIVXNL = SIVXNP + 1 )
C
C     Vertex-plate list count:
C
      INTEGER               SIVTNL
      PARAMETER           ( SIVTNL = SIVXNL + 1 )
C
C     Coarse grid pointers:
C
      INTEGER               SICGRD
      PARAMETER           ( SICGRD = SIVTNL + 1 )

C
C     Size of fixed-size portion of integer component:
C     
      INTEGER               IXIFIX
      PARAMETER           ( IXIFIX = MAXCGR + 7 )


C
C        DSK type 2 spatial index double precision component
C        ---------------------------------------------------
C
C           +-----------------+
C           | Vertex bounds   |  6 values (min/max for each component)
C           +-----------------+
C           | Voxel origin    |  3 elements
C           +-----------------+
C           | Voxel size      |  1 element
C           +-----------------+
C
C
C
C        Index parameters
C
C     Vertex bounds:
C
      INTEGER               SIVTBD
      PARAMETER           ( SIVTBD = 1 )
C
C     Voxel grid origin:
C
      INTEGER               SIVXOR
      PARAMETER           ( SIVXOR = SIVTBD + 6 )
C
C     Voxel size:  
C
      INTEGER               SIVXSZ
      PARAMETER           ( SIVXSZ = SIVXOR + 3 )

C
C     Size of fixed-size portion of double precision component:
C     
      INTEGER               IXDFIX
      PARAMETER           ( IXDFIX = 10 )


C
C     The limits below are used to define a suggested maximum
C     size for the integer component of the spatial index. 
C

C
C     Maximum number of entries in voxel-plate pointer array:
C               
      INTEGER               MAXVXP
      PARAMETER           ( MAXVXP = MAXPLT / 2 )

C
C     Maximum cell size:
C
      INTEGER               MAXCEL
      PARAMETER           ( MAXCEL = 60000000 )

C
C     Maximum number of entries in voxel-plate list:
C
      INTEGER               MXNVLS
      PARAMETER           ( MXNVLS = MAXCEL + MAXVXP/2 )

C
C     Spatial index integer component size:
C
      INTEGER               SPAISZ
      PARAMETER           ( SPAISZ =   IXIFIX + MAXVXP
     .                               + MXNVLS + MAXVRT  
     .                               + MAXNPV          )


C
C     End of include file dsk02.inc
C

