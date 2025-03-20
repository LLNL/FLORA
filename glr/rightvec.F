! $Id: rightvec.m,v 1.5 2003/06/04 17:00:30 bulmer Exp $

      subroutine rightvec
            use ArraySizes 
            use Const_1 
            use Const_1a 
            use Fstor 
            use Matrix 

      implicit none

c..... calculates right hand side vector for both equations,
c..... rhs1(k)=2*a2*xr(n)-a2*xr(n-1)+b1*(xi(l)-xi(n-1)) , and
c..... rhs2(k)=2*a3*xi(n)-a2*xi(n-1)+b1*(xr(l)-xr(n-1)) .

      integer m,mdel,koff1,koff2,koff,k,mbar,m1
      integer i, j, koff_ghost

      call myzero(kxx,0.,rhs1,1)
      call myzero(kxx,0.,rhs2,1)
cjb  what is  sscal?  some system call,  not in vdf  not in undefined symbols
cjb  just some vector 0 init?
cjb   sscal is in code.map  part of libblas.a

cjb  where kxx set?:  declrd  not intlzd in vdf
ctq  kxx is defined in setup.F, after reading ix and jx from ArraySizes

      do m=1,9
          mdel=(m-1)/3
          m1=m-1
          koff1=-mdel*5+m+(mdel-1)*ix
          koff2=(m-((mdel+1)*3-1))*jx+1-2*m1+7*mdel
          koff=.5*(1+isw)*koff1+.5*(1-isw)*koff2 
c           if isw = 1, koff = koff1
c           else isw = -1, koff = koff2
          do k=1,kxx
              ! Convert flattened index k to 2D indices i,j
              i = mod(k-1,ix-2) + 2  ! i goes from 2 to ix-1
              j = (k-1)/(ix-2) + 2   ! j goes from 2 to jx-1
              
              if (k < jx .or. (kxx-k) < jx) then
                  ! These are boundary points - need proper ghost cell handling
                  if (i == 2) then
                      ! Left boundary - use fi1 condition
                      if (fi1 > 0) then
                          ! Neumann: zero gradient
                          koff_ghost = k + 1
                      else
                          ! Dirichlet: zero value
                          koff_ghost = k - 1
                      endif
                  else if (i == ix-1) then
                      ! Right boundary - use fizx condition
                      if (fizx > 0) then
                          koff_ghost = k - 1
                      else 
                          koff_ghost = k + 1
                      endif
                  endif
                  
                  if (j == 2) then
                      ! Bottom boundary - use fj1 condition
                      if (fj1 > 0) then
                          koff_ghost = k + (ix-2)
                      else
                          koff_ghost = k - (ix-2)
                      endif
                  else if (j == jx-1) then
                      ! Top boundary - use fjrx condition
                      if (fjrx > 0) then
                          koff_ghost = k - (ix-2)
                      else
                          koff_ghost = k + (ix-2)
                      endif
                  endif
                  
                  ! Check bounds before accessing array
                  if (koff_ghost >= 1 .and. koff_ghost <= kxx) then
                      rhs1(k) = rhs1(k) + 2.*a3(k,m)*xro(koff_ghost) - a2(k,m)*xroo(koff_ghost)
                      rhs2(k) = rhs2(k) + 2.*a3(k,m)*xio(koff_ghost) - a2(k,m)*xioo(koff_ghost)
                  else
                      ! For out-of-bounds ghost cells, use the boundary value
                      rhs1(k) = rhs1(k) + 2.*a3(k,m)*xro(k) - a2(k,m)*xroo(k)
                      rhs2(k) = rhs2(k) + 2.*a3(k,m)*xio(k) - a2(k,m)*xioo(k)
                  endif
              else
                  ! Interior points - use standard stencil offset
                  rhs1(k) = rhs1(k) + 2.*a3(k,m)*xro(k+koff) - a2(k,m)*xroo(k+koff)
                  rhs2(k) = rhs2(k) + 2.*a3(k,m)*xio(k+koff) - a2(k,m)*xioo(k+koff)
              endif
          enddo
          if (m == 2 .or. m == 5 .or. m == 8) then
              do k=1,kxx
                  mbar=m-1-(m/4)*2
                  if (k < jx .or. (kxx-k) < jx) then
                  rhs1(k)=rhs1(k)+fac1*b1(k,mbar)*(xiol(k)-xioo(k))
                  rhs2(k)=rhs2(k)+fac2*b1(k,mbar)*(xrol(k)-xroo(k))
                  else
                  rhs1(k)=rhs1(k)+fac1*b1(k,mbar)*(xiol(k+koff)-xioo(k+koff))
                  rhs2(k)=rhs2(k)+fac2*b1(k,mbar)*(xrol(k+koff)-xroo(k+koff))
                  endif
              enddo
          endif
      enddo

      return
      end ! rightvec
