      program Flora

          use ArraySizes
          use Const_1
          use Const_1a
          use Const_2
          use Const_3
          use Const_6
          use Const_8
          use Coils
          use TMInput
          implicit none
          print *, "Hi! my name is FLORA."

        
          ! Read the namelist from a file
          open(unit=10, file='in_flora.nml', status='old')
          read(10, nml=flora_ArraySizes)
          read(10, nml=flora_Coils)
          read(10, nml=flora_Const_1)
          read(10, nml=flora_Const_1a)
          read(10, nml=flora_Const_2)
          read(10, nml=flora_Const_3)
          read(10, nml=flora_Const_6)
          read(10, nml=flora_Const_8)
          read(10, nml=flora_TMInput)


          call glrgen
          do while (n < nmax)
            !print *, "n", n
            call glrexe
          end do
          call glrfin

      end program Flora
      
