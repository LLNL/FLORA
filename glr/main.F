      program Flora

          use Const_2
          use Const_3
          use Const_6
          use Coils
          use TMInput
          implicit none
          print *, "Hi! my name is FLORA."

        
          ! Read the namelist from a file
          open(unit=10, file='flora_in.nml', status='old')
          read(10, nml=flora_Coils)
          read(10, nml=flora_Const_2)
          read(10, nml=flora_Const_3)
          read(10, nml=flora_Const_6)
          read(10, nml=flora_TMInput)

          call glrgen
          print *, "Finished GLRGEN"

      end program Flora
      
