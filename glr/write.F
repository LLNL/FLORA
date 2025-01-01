      subroutine netcdfInit
             use netcdf
             use ArraySizes
             use ncInts

      implicit none
      integer(4) :: ierr

      ! Create the NetCDF file
      ierr = nf90_create("flora.nc", nf90_clobber, ncid)

      ! Define dimensions
      ierr = nf90_def_dim(ncid, "kxx", INT(kxx,4), dim_kxx)
      ierr = nf90_def_dim(ncid, "nmax", INT(nmax,4), dim_n)

      ! Define variables
      ierr = nf90_def_var(ncid, "XRO", nf90_double, [dim_kxx, dim_n], var_xro)
      !ierr = nf90_def_var(ncid, "XRO", nf90_double, [dim_n, dim_kxx], var_xro)

      ! End the define mode
      ierr = nf90_enddef(ncid)
      end

      subroutine netcdfWrite
             use netcdf
             use ncInts
             use Fstor
             use Const_1a

      implicit none
      integer(4) :: ierr
      
      xro_t(n,:) = xro
      end


      subroutine netcdfClose
             use netcdf
             use ncInts
             use Fstor

      implicit none
      integer(4) :: ierr

      ! Write variables
      ierr = nf90_put_var(ncid, var_xro, transpose(xro_t))
      !ierr = nf90_put_var(ncid, var_xro, xro_t)

      ! Close the NetCDF file
      ierr = nf90_close(ncid)
      if (ierr /= nf90_noerr) then
        print *, "Error closing NetCDF file: ", trim(nf90_strerror(ierr))
        stop
      end if

      print *, "NetCDF file successfully written."
      end

