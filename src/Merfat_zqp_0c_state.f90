
!_______________________________________________________________________________________________________________________________________________!
!                                                                                                                                               !
!                     2-  subroutine zqp_0_contact                                                                                              !
!                                                                                                                                               !
! 1- Variables and Arrays:                                                                                                                      !
! i_db:                    Loop variable for iterating over dynamic basis states (up to n_db_0c).                                               !
! i_df:                    Loop variable for iterating over dynamic features.                                                                   !
! state_number:            Variable to store the state number.                                                                                  !
! zero_contact_states:     Array to store zero-contact states.                                                                                  !
! state:                   Array to store the current state.                                                                                    !
!                                                                                                                                               !
! 2- File I/O                                                                                                                                   !
! Opens 'one_zero.out' for reading (unit=294).                                                                                                  !
! Opens '0c_states.out' for writing (unit=296). The status='replace' ensures that if the file already exists, it will be overwritten.           !
!                                                                                                                                               !
! 3- Loop Over Basis States: ! Make sure the loop runs only up to n_db_0c                                                                       !
! The loop iterates over dynamic basis states, ensuring it runs up to n_db_0c.                                                                  !
!                                                                                                                                               !
! 4- Read State Information:                                                                                                                    !
! Reads the state number and the corresponding dynamic features from 'one_zero.out'.                                                            !
!                                                                                                                                               !
! 5- Check for Zero-Contact State:                                                                                                              !
! Calls the has_consecutive_ones function to check if the state has consecutive ones.                                                           !
! If it doesn't have consecutive ones, writes the state to '0c_states.out'.                                                                     !
!                                                                                                                                               !
! 6- File Closure                                                                                                                               !
! Closes the input and output files.                                                                                                            !
!                                                                                                                                               !
! 7- Function has_consecutive_ones                                                                                                              !
! A logical function that checks if the array arr has consecutive ones. If yes, it returns .true.; otherwise, it returns .false.                !
! The primary aim of the function is to provide a logical indication of whether there are consecutive ones in the array.                        !
! It is used in the context of checking if a state has zero contact (where 1 should not be repeated consecutively in the state                  !
!_______________________________________________________________________________________________________________________________________________!

module Merfat_zqp_0c_state
  use Merfat_arr_prm
  use Merfat_zqp_fullbasis
  implicit none

contains

  subroutine zqp_0_contact
  
    integer :: i_db, i_df, state_number
    integer :: zero_contact_states(n_db_0c)
    integer, dimension(n_df) :: state
    integer :: count_selected_states

    ! Initialize the count of selected states
    count_selected_states = 0

    ! File I/O
    open(file='output/one_zero.out', unit=294, status='old')
    open(file='output/0c_states.out', unit=296, status='unknown')

    ! Loop Over Basis States:
    ! Make sure the loop runs for all n_db
  do i_db = 1, n_db

      ! Read State Information:
      ! Read the state number
      read(294, *) state_number

      ! Read the state
            do i_df = 1, n_df
               read(294, *) state(i_df)
            end do

      ! Check if it is a zero-contact state
      if (.not. has_consecutive_ones(state)) then

        ! Write the state to zero_contact_states.out
        write(296, '(I4)') state_number

            do i_df = 1, n_df
              write(296, '(I1)') state(i_df)
            end do

        write(296, *)

        ! Increment the count of selected states
        count_selected_states = count_selected_states + 1

        ! Check if the desired number of states is reached
        if (count_selected_states == n_db_0c) exit
      end if

 end do

    ! File Closure
    close(unit=294)
    close(unit=296)

  contains

    ! Function has_consecutive_ones
    logical function has_consecutive_ones(arr)
      integer, dimension(:) :: arr
      integer :: i

      has_consecutive_ones = .false.

      do i = 1, size(arr) - 1
        if (arr(i) == 1 .and. arr(i + 1) == 1) then
          has_consecutive_ones = .true.
          return
        end if
      end do


    end function has_consecutive_ones

  end subroutine zqp_0_contact

end module Merfat_zqp_0c_state



