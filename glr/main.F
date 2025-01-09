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
          use ncInts
          implicit none
          print *, "Hi! my name is FLORA."

          ! Get the input file name from command line argument
          call get_command_argument(1, input_file)
          print *, "Input file is: ", input_file

          ! Read the namelist from a file
          open(unit=10, file=input_file, status='old')
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
      
