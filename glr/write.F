      subroutine netcdfInit
             use netcdf
             use ArraySizes
             use ncInts
             use Const_1
             use Const_1a
             use Const_2
             use Const_3
             use Const_8
             use Fstor

      implicit none
      integer(4) :: ierr


      ! Create the NetCDF file
      if (input_file .ne. "in_flora.nml") then
        call create_output_filename(input_file, output_file)
      end if
      ierr = nf90_create(output_file, nf90_clobber, ncid)

      ! Define dimensions
      ierr = nf90_def_dim(ncid, "z", INT(ix,4), dim_z)
      ierr = nf90_def_dim(ncid, "psi", INT(jx,4), dim_psi)
      ierr = nf90_def_dim(ncid, "t", INT(nmax,4), dim_t)
      ierr = nf90_def_dim(ncid, "kxx", INT(kxx,4), dim_kxx)

      !!! Define variables

      ! inputs
      ierr = nf90_def_var(ncid, "dt", nf90_double, var_dt)
      ierr = nf90_def_var(ncid, "mm", nf90_double, var_mm)
      ierr = nf90_def_var(ncid, "fi1", nf90_double, var_fi1)
      ierr = nf90_def_var(ncid, "fj1", nf90_double, var_fj1)
      ierr = nf90_def_var(ncid, "fizx", nf90_double, var_fizx)
      ierr = nf90_def_var(ncid, "fjrx", nf90_double, var_fjrx)
      ierr = nf90_def_var(ncid, "sf6", nf90_double, var_sf6)
      ierr = nf90_def_var(ncid, "sf8", nf90_double, var_sf8)
      ierr = nf90_def_var(ncid, "ex0", nf90_double, var_ex0)
      ierr = nf90_def_var(ncid, "ex1", nf90_double, var_ex1)
      ierr = nf90_def_var(ncid, "lmax", nf90_double, var_lmax)
      ierr = nf90_def_var(ncid, "nmax", nf90_double, var_nmax)
      ierr = nf90_def_var(ncid, "kzs", nf90_double, var_kzs)
      ierr = nf90_def_var(ncid, "bias", nf90_double, var_bias)
      ierr = nf90_def_var(ncid, "cold", nf90_double, var_cold)
      ierr = nf90_def_var(ncid, "dip", nf90_double, var_dip)

      ierr = nf90_def_var(ncid, "ncenter", nf90_double, var_ncenter)
      ierr = nf90_def_var(ncid, "nsloshin", nf90_double, var_nsloshin)
      ierr = nf90_def_var(ncid, "var_betcent", nf90_double, var_betcent)
      ierr = nf90_def_var(ncid, "var_betcene", nf90_double, var_betcene)
      ierr = nf90_def_var(ncid, "var_betpas1", nf90_double, var_betpas1)
      ierr = nf90_def_var(ncid, "var_betslsh", nf90_double, var_betslsh)
      ierr = nf90_def_var(ncid, "var_betslse", nf90_double, var_betslse)
      ierr = nf90_def_var(ncid, "var_ppas2", nf90_double, var_ppas2)
      ierr = nf90_def_var(ncid, "var_ppas3", nf90_double, var_ppas3)
      ierr = nf90_def_var(ncid, "var_dpas1", nf90_double, var_dpas1)
      ierr = nf90_def_var(ncid, "var_d1trap", nf90_double, var_d1trap)

      ! initial conditions
      ierr = nf90_def_var(ncid, "BVAC", nf90_double, dim_z, var_bvac)
      ierr = nf90_def_var(ncid, "B", nf90_double, [dim_psi, dim_z], var_b)
      ierr = nf90_def_var(ncid, "R", nf90_double, [dim_psi, dim_z], var_r)
      ierr = nf90_def_var(ncid, "RHO", nf90_double, [dim_psi, dim_z], var_rho)
      ierr = nf90_def_var(ncid, "PPERP", nf90_double, [dim_psi, dim_z], var_pperp)
      ierr = nf90_def_var(ncid, "PPAR", nf90_double, [dim_psi, dim_z], var_ppar)
      ierr = nf90_def_var(ncid, "EPSI", nf90_double, [dim_psi, dim_z], var_epsi)

      ! results
      ierr = nf90_def_var(ncid, "XRO", nf90_double, [dim_kxx, dim_t], var_xro)
      ierr = nf90_def_var(ncid, "XIO", nf90_double, [dim_kxx, dim_t], var_xio)
      ierr = nf90_def_var(ncid, "FLUTE3", nf90_double, var_flute3)

      ! End the define mode
      ierr = nf90_enddef(ncid)

      ! Write variables that are defined AFTER glrgen and BEFORE time loop
      ierr = nf90_put_var(ncid, var_dt, dt)
      ierr = nf90_put_var(ncid, var_flute3, flute3)
      ierr = nf90_put_var(ncid, var_mm, mm)
      ierr = nf90_put_var(ncid, var_fi1, fi1)
      ierr = nf90_put_var(ncid, var_fj1, fj1)
      ierr = nf90_put_var(ncid, var_fizx, fizx)
      ierr = nf90_put_var(ncid, var_fjrx, fjrx)
      ierr = nf90_put_var(ncid, var_sf6, sf6)
      ierr = nf90_put_var(ncid, var_sf8, sf8)
      ierr = nf90_put_var(ncid, var_ex0, ex0)
      ierr = nf90_put_var(ncid, var_ex1, ex1)
      ierr = nf90_put_var(ncid, var_lmax, lmax)
      ierr = nf90_put_var(ncid, var_nmax, nmax)
      ierr = nf90_put_var(ncid, var_kzs, kzs)
      ierr = nf90_put_var(ncid, var_bias, bias)
      ierr = nf90_put_var(ncid, var_cold, cold)
      ierr = nf90_put_var(ncid, var_dip, dip)

      ierr = nf90_put_var(ncid, var_ncenter, ncenter)
      ierr = nf90_put_var(ncid, var_nsloshin, nsloshin)
      ierr = nf90_put_var(ncid, var_betcent, betcent)
      ierr = nf90_put_var(ncid, var_betcene, betcene)
      ierr = nf90_put_var(ncid, var_betpas1, betpas1)
      ierr = nf90_put_var(ncid, var_betslsh, betslsh)
      ierr = nf90_put_var(ncid, var_betslse, betslse)
      ierr = nf90_put_var(ncid, var_dpas1, dpas1)
      ierr = nf90_put_var(ncid, var_ppas2, ppas2)
      ierr = nf90_put_var(ncid, var_ppas3, ppas3)
      ierr = nf90_put_var(ncid, var_d1trap, d1trap)

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

      ! Subroutine to generate output file name based on input
      subroutine create_output_filename(input, output)

        implicit none
        character(len=256), intent(in) :: input
        character(len=256), intent(out) :: output
        integer :: len_input
    
        len_input = len_trim(input)
        
        ! Example: Append .nc to the input file name
        write(output, '(A)') input(1:len_input)

        ! Replace any file extension with '.nc' (if needed)
        if (len_input > 4 .and. input(len_input-3:len_input) == '.nml') then
           output = input(1:len_input-4) // '.nc'
        else
           output = input(1:len_input) // '.nc'
        end if

        print *, "FILE IO:"
        print *, "    input: ", trim(input)
        print *, "    output: ", trim(output)

      end subroutine create_output_filename


