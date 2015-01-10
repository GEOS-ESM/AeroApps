
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.... Picks the subregion containing (x,y,z) near an icosahedral vertex.

      integer function ipick( x, y, z, ir1, ir2, ir3, ir4, ir5,
     $                                 nv1, nv2, nv3, nv4, nv5  )
      Use m_Spherical_Partition, Only : maxreg => MAXREG_TRADITIONAL
c.......................................................................
c.... Argument declarations.

      real         x, y, z
      integer      ir1, ir2, ir3, ir4, ir5
      integer      nv1, nv2, nv3, nv4, nv5

c.......................................................................
c.... Icosahedral grid definitions.

      include     "icosdef.h"

c.......................................................................
c.... Local storage.

      real         d(5)
      integer      ir(5)

c.......................................................................
c.... Statewment function definitions.

      vx(nv) = xyzicos(1,nv)
      vy(nv) = xyzicos(2,nv)
      vz(nv) = xyzicos(3,nv)

      distsq(xx,yy,zz,nv) = (xx-vx(nv))**2+(yy-vy(nv))**2+(zz-vz(nv))**2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c.......................................................................
c.... Compute square of distances from (x,y,z) to surrounding vertices
c.... of icosahedral grid.

      dist1 = distsq(x,y,z,nv1)
      dist2 = distsq(x,y,z,nv2)
      dist3 = distsq(x,y,z,nv3)
      dist4 = distsq(x,y,z,nv4)
      dist5 = distsq(x,y,z,nv5)

c.......................................................................
c.... Compute twice average of squares of distances to adjacent vertices

      d(1)  = dist5 + dist1
      d(2)  = dist1 + dist2
      d(3)  = dist2 + dist3
      d(4)  = dist3 + dist4
      d(5)  = dist4 + dist5

c.......................................................................
c.... Set face indices corresponding to the d(k)'s.
      ir(1) = ir1
      ir(2) = ir2
      ir(3) = ir3
      ir(4) = ir4
      ir(5) = ir5

c.......................................................................
c.... Pick the face corresponding to the smallest d(k).

      ipick = ir(1)
      dmini = d(1)
      do 100 i = 2, 5
         if( d(i).lt.dmini ) then
            dmini = d(i)
            ipick = ir(i)
         endif
  100 continue

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end
