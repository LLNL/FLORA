
      subroutine writeCDF
              use netcdf
              use Fstor
              use ArraySizes
              implicit none
            
              ! NetCDF file ID and other variables
              integer(4) :: ncid, dimid, varid, ierr
              integer(4) :: dim_kxx
              integer(4) :: var_xro
              integer(4), parameter :: nx = 10, ny = 10
              integer :: i, j
              real :: data_array(nx, ny)
            
              ! Create the NetCDF file
              ierr = nf90_create("flora.nc", nf90_clobber, ncid)
              if (ierr /= nf90_noerr) then
                print *, "Error creating NetCDF file: ", trim(nf90_strerror(ierr))
                stop
              end if

             ! Define dimensions
             ierr = nf90_def_dim(ncid, "x", nx, dimid)
             if (ierr /= nf90_noerr) then
               print *, "Error defining dimension: ", trim(nf90_strerror(ierr))
               stop
             end if
           
             ierr = nf90_def_dim(ncid, "y", ny, dimid)
             if (ierr /= nf90_noerr) then
               print *, "Error defining dimension: ", trim(nf90_strerror(ierr))
               stop
             end if
             ierr = nf90_def_dim(ncid, "kxx", INT(kxx,4), dim_kxx)
           
             ! Define the variable
             ierr = nf90_def_var(ncid, "temperature", nf90_real, [dimid, dimid], varid)
             if (ierr /= nf90_noerr) then
               print *, "Error defining variable: ", trim(nf90_strerror(ierr))
               stop
             end if

             ierr = nf90_def_var(ncid, "XRO", nf90_double, [dim_kxx], var_xro)
           
             ! End the define mode
             ierr = nf90_enddef(ncid)
             if (ierr /= nf90_noerr) then
               print *, "Error ending define mode: ", trim(nf90_strerror(ierr))
               stop
             end if

             ! Create a 10x10 array of data to write
             do i = 1, nx
               do j = 1, ny
                 data_array(i, j) = real(i + j)  ! Example data
               end do
             end do

             ! Write data to the variable
             ierr = nf90_put_var(ncid, varid, data_array)
             if (ierr /= nf90_noerr) then
               print *, "Error writing data: ", trim(nf90_strerror(ierr))
               stop
             end if
             print *, xro
             ierr = nf90_put_var(ncid, var_xro, xro)

             ! Close the NetCDF file
             ierr = nf90_close(ncid)
             if (ierr /= nf90_noerr) then
               print *, "Error closing NetCDF file: ", trim(nf90_strerror(ierr))
               stop
             end if

             print *, "NetCDF file successfully written."
      end
