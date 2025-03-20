! $Id: rightvec.m,v 1.5 2003/06/04 17:00:30 bulmer Exp $

      subroutine rightvec
            use ArraySizes 
            use Const_1 
            use Const_1a 
            use Fstor 
            use Matrix 

      implicit none

c..... calculates right hand side vector for both equations in 2D form,
c..... br = 2*a3*xr(n) - a2*xr(n-1) + b1*(xi(l)-xi(n-1))
c..... bi = 2*a3*xi(n) - a2*xi(n-1) + b1*(xr(l)-xr(n-1))
c..... where br/bi are real/imaginary components of RHS vector b in Ax = b

      integer :: m, mdel, m1, mbar
      integer :: i, j, k
      integer :: i_off, j_off
      real*8, dimension(ix,jx) :: br, bi           ! RHS vector components
      real*8, dimension(ix,jx) :: xro_ghost, xio_ghost  ! Ghost cell arrays
      real*8, dimension(ix,jx) :: xroo_ghost, xioo_ghost

      ! Initialize RHS arrays to zero
      do j=1,jx
          do i=1,ix
              br(i,j) = 0.0
              bi(i,j) = 0.0
              xro_ghost(i,j) = 0.0
              xio_ghost(i,j) = 0.0
              xroo_ghost(i,j) = 0.0
              xioo_ghost(i,j) = 0.0
          enddo
      enddo

      ! Initialize ghost arrays with interior values
      do j=2,jx-1
          do i=2,ix-1
              k = (i-2) + (j-2)*(ix-2) + 1
              xro_ghost(i,j) = xro(k)
              xio_ghost(i,j) = xio(k)
              xroo_ghost(i,j) = xroo(k)
              xioo_ghost(i,j) = xioo(k)
          enddo
      enddo

      ! Handle boundaries
      do j=2,jx-1
          ! Left boundary (i=2)
          if (fi1 > 0) then
              ! Neumann: zero gradient
              xro_ghost(1,j) = xro_ghost(2,j)
              xio_ghost(1,j) = xio_ghost(2,j)
              xroo_ghost(1,j) = xroo_ghost(2,j)
              xioo_ghost(1,j) = xioo_ghost(2,j)
          else
              ! Dirichlet: zero value
              xro_ghost(1,j) = 0.0
              xio_ghost(1,j) = 0.0
              xroo_ghost(1,j) = 0.0
              xioo_ghost(1,j) = 0.0
          endif

          ! Right boundary (i=ix-1)
          if (fizx > 0) then
              xro_ghost(ix,j) = xro_ghost(ix-1,j)
              xio_ghost(ix,j) = xio_ghost(ix-1,j)
              xroo_ghost(ix,j) = xroo_ghost(ix-1,j)
              xioo_ghost(ix,j) = xioo_ghost(ix-1,j)
          else
              xro_ghost(ix,j) = 0.0
              xio_ghost(ix,j) = 0.0
              xroo_ghost(ix,j) = 0.0
              xioo_ghost(ix,j) = 0.0
          endif
      enddo

      ! Handle top/bottom boundaries
      do i=2,ix-1
          ! Bottom boundary (j=2)
          if (fj1 > 0) then
              xro_ghost(i,1) = xro_ghost(i,2)
              xio_ghost(i,1) = xio_ghost(i,2)
              xroo_ghost(i,1) = xroo_ghost(i,2)
              xioo_ghost(i,1) = xioo_ghost(i,2)
          else
              xro_ghost(i,1) = 0.0
              xio_ghost(i,1) = 0.0
              xroo_ghost(i,1) = 0.0
              xioo_ghost(i,1) = 0.0
          endif

          ! Top boundary (j=jx-1)
          if (fjrx > 0) then
              xro_ghost(i,jx) = xro_ghost(i,jx-1)
              xio_ghost(i,jx) = xio_ghost(i,jx-1)
              xroo_ghost(i,jx) = xroo_ghost(i,jx-1)
              xioo_ghost(i,jx) = xioo_ghost(i,jx-1)
          else
              xro_ghost(i,jx) = 0.0
              xio_ghost(i,jx) = 0.0
              xroo_ghost(i,jx) = 0.0
              xioo_ghost(i,jx) = 0.0
          endif
      enddo

      ! Main computation loop
      do j=2,jx-1
          do i=2,ix-1
              do m=1,9
                  mdel=(m-1)/3
                  m1=m-1
                  ! Calculate stencil offsets
                  i_off = mod(m-1,3) - 1  ! -1, 0, or 1
                  j_off = mdel - 1        ! -1, 0, or 1

                  ! Time evolution terms
                  k = (i-2) + (j-2)*(ix-2) + 1
                  br(i,j) = br(i,j) + 2.*a3(k,m)*
     &                      xro_ghost(i+i_off,j+j_off)
     &                      - a2(k,m)*
     &                      xroo_ghost(i+i_off,j+j_off)
                  bi(i,j) = bi(i,j) + 2.*a3(k,m)*
     &                      xio_ghost(i+i_off,j+j_off)
     &                      - a2(k,m)*
     &                      xioo_ghost(i+i_off,j+j_off)

                  ! Coupling terms (only for m=2,5,8)
                  if (m == 2 .or. m == 5 .or. m == 8) then
                      mbar = m-1-(m/4)*2
                      br(i,j) = br(i,j) + fac1*b1(k,mbar)*
     &                         (xiol(k)-
     &                         xioo_ghost(i+i_off,j+j_off))
                      bi(i,j) = bi(i,j) + fac2*b1(k,mbar)*
     &                         (xrol(k)-
     &                         xroo_ghost(i+i_off,j+j_off))
                  endif
              enddo
          enddo
      enddo

      ! Flatten 2D arrays to 1D for output
      do j=2,jx-1
          do i=2,ix-1
              k = (i-2) + (j-2)*(ix-2) + 1
              rhs1(k) = br(i,j)
              rhs2(k) = bi(i,j)
          enddo
      enddo

      return
      end ! rightvec
