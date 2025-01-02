      subroutine netcdfInit
             use netcdf
             use ArraySizes
             use ncInts
             use Const_8
             use Fstor

      implicit none
      integer(4) :: ierr

      ! Create the NetCDF file
      ierr = nf90_create("flora.nc", nf90_clobber, ncid)

      ! Define dimensions
      ierr = nf90_def_dim(ncid, "z", INT(ix,4), dim_z)
      ierr = nf90_def_dim(ncid, "psi", INT(jx,4), dim_psi)
      ierr = nf90_def_dim(ncid, "t", INT(nmax,4), dim_t)
      ierr = nf90_def_dim(ncid, "kxx", INT(kxx,4), dim_kxx)

      ! Define variables
      ierr = nf90_def_var(ncid, "FLUTE3", nf90_double, var_flute3)

      ierr = nf90_def_var(ncid, "BVAC", nf90_double, dim_z, var_bvac)
      ierr = nf90_def_var(ncid, "B", nf90_double, [dim_psi, dim_z], var_b)
      ierr = nf90_def_var(ncid, "R", nf90_double, [dim_psi, dim_z], var_r)
      ierr = nf90_def_var(ncid, "RHO", nf90_double, [dim_psi, dim_z], var_rho)
      ierr = nf90_def_var(ncid, "PPERP", nf90_double, [dim_psi, dim_z], var_pperp)
      ierr = nf90_def_var(ncid, "PPAR", nf90_double, [dim_psi, dim_z], var_ppar)
      ierr = nf90_def_var(ncid, "EPSI", nf90_double, [dim_psi, dim_z], var_epsi)

      ierr = nf90_def_var(ncid, "XRO", nf90_double, [dim_kxx, dim_t], var_xro)
      ierr = nf90_def_var(ncid, "XIO", nf90_double, [dim_kxx, dim_t], var_xio)

      ! End the define mode
      ierr = nf90_enddef(ncid)

      ! Write variables that are defined AFTER glrgen and BEFORE time loop
      ierr = nf90_put_var(ncid, var_flute3, flute3)

      ierr = nf90_put_var(ncid, var_bvac, bvac)
      ierr = nf90_put_var(ncid, var_b, transpose(b))
      ierr = nf90_put_var(ncid, var_r, transpose(r))
      ierr = nf90_put_var(ncid, var_rho, transpose(rho))
      ierr = nf90_put_var(ncid, var_pperp, transpose(pperp))
      ierr = nf90_put_var(ncid, var_ppar, transpose(ppar))
      ierr = nf90_put_var(ncid, var_epsi, transpose(epsi))
      end

      subroutine netcdfWrite
             use netcdf
             use ncInts
             use Fstor
             use Const_1a
             use Const_8

      implicit none
      integer(4) :: ierr
      
      xro_t(n,:) = xro
      xio_t(n,:) = xio
      end


      subroutine netcdfClose
             use netcdf
             use ncInts
             use Fstor

      implicit none
      integer(4) :: ierr

      ! Write variables
      ierr = nf90_put_var(ncid, var_xro, transpose(xro_t))
      ierr = nf90_put_var(ncid, var_xio, transpose(xio_t))

      ! Close the NetCDF file
      ierr = nf90_close(ncid)
      if (ierr /= nf90_noerr) then
        print *, "Error closing NetCDF file: ", trim(nf90_strerror(ierr))
        stop
      end if

      print *, "NetCDF file successfully written."
      end

