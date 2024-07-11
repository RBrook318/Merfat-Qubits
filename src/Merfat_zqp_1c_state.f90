
!_______________________________________________________________________________________________________________________________________________!
!                                                                                                                                               !
!                     2-  subroutine zqp_1_contact                                                                                             !
!                                                                                                                                               !
! 1- Variables and Arrays:                                                                                                                      !
! i_db:                    Loop variable for iterating over dynamic basis states (up to n_db_1c).                                               !
! i_df:                    Loop variable for iterating over dynamic features.                                                                   !
! state_number:            Variable to store the state number.                                                                                  !
! one_contact_states:      Array to store zero-contact states.                                                                                  !
! state:                   Array to store the current state.                                                                                    !
!                                                                                                                                               !
! 2- File I/O                                                                                                                                   !
! Opens 'one_zero.out' for reading (unit=298).                                                                                                  !
! Opens '1c_states.out' for writing (unit=299). The status='replace' ensures that if the file already exists, it will be overwritten. !
!                                                                                                                                               !
! 3- Loop Over Basis States: ! Make sure the loop runs only up to n_db_1c                                                                       !
! The loop iterates over dynamic basis states, ensuring it runs up to n_db_1c.                                                                  !
!                                                                                                                                               !
! 4- Read State Information:                                                                                                                    !
! Reads the state number and the corresponding dynamic features from 'one_zero.out'.                                                            !
!                                                                                                                                               !
! 5- Check for one-Contact State:                                                                                                               !
! Calls the has_consecutive_ones function to check if the state has consecutive ones.                                                           !
! If it doesn't have consecutive ones, writes the state to '1c_states.out'.                                                                     !
!                                                                                                                                               !
! 6- File Closure                                                                                                                               !
! Closes the input and output files.                                                                                                            !
!                                                                                                                                               !
! 7- Function has_consecutive_ones                                                                                                              !
! A logical function that checks if the array arr has consecutive ones. If yes, it returns .true.; otherwise, it returns .false.                !
! The primary aim of the function is to provide a logical indication of whether there are consecutive ones in the array.                        !
! It is used in the context of checking if a state has one contact (where 1 should be repeated consecutively in the state one time.             !
!_______________________________________________________________________________________________________________________________________________!


module Merfat_zqp_1c_state
  use Merfat_arr_prm
  use Merfat_zqp_fullbasis
 
  implicit none

contains

subroutine zqp_1_contact

    integer :: i_db, i_df, state_number
    integer :: one_contact_states(n_db_1c)
    integer, dimension(n_df) :: state

    ! File I/O
    open(file='output/one_zero.out', unit=298, status='unknown')
    open(file='output/1c_states.out', unit=299, status='unknown')

    ! Loop Over Basis States: 
    ! Make sure the loop runs for all n_db
    do i_db = 1, n_db  

      ! Read State Information:
      ! Read the state number
      read(298, *) state_number

      ! Read the state
      do i_df = 1, n_df
        read(298, *) state(i_df)
      end do

      ! Check if it is a one-contact state
      if (has_one_consecutive_ones(state)) then
      
        ! Write the state to one_contact_states.out
        write(299, '(I4)') state_number
        do i_df = 1, n_df
          write(299, '(I1)') state(i_df)
        end do
        write(299, *)
      end if
    end do

    ! File Closure
    close(unit=298)
    close(unit=299)

  contains

    ! Function has_one_consecutive_ones
    logical function has_one_consecutive_ones(arr)
      integer, dimension(:) :: arr
      integer :: i
      integer :: consecutive_ones

      consecutive_ones = 0

      do i = 1, size(arr) - 1
        if (arr(i) == 1 .and. arr(i + 1) == 1) then
          consecutive_ones = consecutive_ones + 1
          if (consecutive_ones > 1) then
            has_one_consecutive_ones = .false.
            return
          end if
        end if
      end do

      has_one_consecutive_ones = (consecutive_ones == 1)



    end function has_one_consecutive_ones


end subroutine zqp_1_contact

end module Merfat_zqp_1c_state
