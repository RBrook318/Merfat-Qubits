!___________________________________________________________________________________________________________________________________________!
!                                                                                                                                           !
!                     1-  subroutine setup_zqp_full_basis                                                                                   !
!                                                                                                                                           !
! this subroutine generates all possible combinations of 1 and 0 for a specified number of degrees of freedom and lines,                    !
! writes these combinations to an output file, and calculates and writes the corresponding zqp lines to another output file.                !
! The intermediate arrays (n_new, n_old, dummy_int) are used for temporary storage during the process.                                      !
!                                                                                                                                           !
! 1- Variables and Arrays:                                                                                                                  !
! n_new(:):                Allocatable array that holds the current set of combinations of 1 and 0.                                         !
! n_old(:):                Allocatable array that holds the previous set of combinations of 1 and 0.                                        !
! dummy_int(n_df + 1):     Array used to read values from an input file (one_zero.out).                                                     !
! i_db, i_df:              Loop indices for lines and degrees of freedom, respectively.                                                     !
! i_basis:                 Unused variable.                                                                                                 !
! i_0, i_1:                Loop indices for generating combinations.                                                                        !
! i_db_a_0, i_db_a_1:      Intermediate variables for generating combinations.                                                              !
! i_db_old, i_df_old:      Loop indices for the old set of combinations.                                                                    !
! a0_in, a1_in:            Real variables for storing coefficients in zqp lines.                                                            !
! n_df_a, n_db_a:          Variables for the current number of degrees of freedom and lines.                                                !
! n_df_old, n_db_old:      Variables for the old number of degrees of freedom and line                                                      !
!                                                                                                                                           !
! 2- (Combination Generation:)                                                                                                              !
! The subroutine starts by generating all possible combinations of 1 and 0 for the initial number of degrees of freedom (n_df_a)            !
!  and lines (n_db_a). This initial set is stored in the n_new array                                                                        !
!                                                                                                                                           !
! 3- Combination Update:                                                                                                                    !
! A loop is used to generate combinations for increasing values of n_df_a until it reaches the desired number of degrees of freedom (n_df). !
! The n_old array is used to store the previous set of combinations.                                                                        !
!                                                                                                                                           !
! 4- Combination Filling:                                                                                                                   !
! The new array (n_new) is filled up with zero and one combinations based on the previous set (n_old).                                      !
!                                                                                                                                           !
! 5- Deallocation:                                                                                                                          !
!The previous set (n_old) is deallocated after the loop.                                                                                    !
!                                                                                                                                           !
! 6- Output to File:                                                                                                                        !
! The combinations are written to an output file (one_zero.out), where each line corresponds to a different line (i_db), and columns        !
! represent the degrees of freedom (i_df).                                                                                                  !
!                                                                                                                                           !
! 7- File Operations:                                                                                                                       !
! Files are opened for reading (unit=295, input file) and writing (unit=290, output file).                                                  !
!                                                                                                                                           !
! 8- Read from Input File:                                                                                                                  !
! Reads values from the input file (one_zero.out) into the dummy_int array.                                                                 !
!                                                                                                                                           !
! 9- Full Basis Header:                                                                                                                     !
! Writes the full basis header to the output file (full_basis), indicating the line number.                                                 !
!                                                                                                                                           !
! 10- ZQP Lines:                                                                                                                            !
! Writes the zqp lines to the output file, calculated based on the coefficients read from dummy_int.                                        !
!                                                                                                                                           !
! 11- File Closure and Deallocation:                                                                                                        !
! Closes the input and output files and deallocates the n_new array.                                                                        !
!                                                                                                                                           !                                                                             
! 12- allocate(n_new) and deallocate(n_new) means:                                                                                          !
! A dynamic resizing process for the array n_new to accommodate increasing combinations of 0s and 1s until n_df_a reaches the value of n_df.!
! Let's break down the sequence of operations:                                                                                              !
!         1- Initially, n_df_a = 1 and n_db_a = 2. Memory is allocated for the array n_new with the size n_df_a * n_db_a, which is 1 * 2 =2 !
!            Then, values are assigned to the first two elements of n_new.                                                                  !
!         2- Inside the loop, the variables n_df_old and n_db_old are used to store the original dimensions of the array n_new before resizing. 
!            This is done to retain the original size for later reference.                                                                   !
!         3- n_df_a is incremented by 1 (n_df_a = n_df_a + 1). This increases the number of rows in the array n_new, allowing for additional !
!            combinations of 0s and 1s to be accommodated.                                                                                   !
!         4- n_db_a is doubled (n_db_a = 2 * n_db_a). This doubles the number of columns in the array n_new, expanding it horizontally to    ! 
!            accommodate more combinations.                                                                                                  !
! Memory for the array n_new is deallocated (deallocate(n_new)). This step releases the memory allocated for the previous size of n_new.     !
! Memory for the array n_new is allocated again with the updated size (allocate(n_new(n_df_a * n_db_a))). This reallocation step creates a   !
! new array n_new with the expanded dimensions to accommodate the increased combinations of 0s and 1s.                                       ! 
! In summary, the sequence of operations reallocates memory for the array n_new with an updated size to accommodate the dynamically changing !
! requirements of the program. The deallocation followed by reallocation ensures that memory is managed efficiently as the size of n_new     !
! changes during program execution.                                                                                                          !
!                                                                                                                                            !
! For example:                                                                                                                               !
! in the one_zero.out file, line 49 is:                                                                                                      !
!  49                                                                                                                                        !
!  0                                                                                                                                         !
!  1                                                                                                                                         !
!  1                                                                                                                                         !
!  0                                                                                                                                         !
!  0                                                                                                                                         !
!  0                                                                                                                                         ! 
!  0                                                                                                                                         !
! Full basis for line:  49                                                                                                                   !
! zqp(1,1,49)=(  1.00000000,0.0d0)                                                                                                           ! 
! zqp(2,1,49)=(  0.00000000,0.0d0)                                                                                                           ! 
! zqp(1,2,49)=(  0.00000000,0.0d0)                                                                                                           ! 
! zqp(2,2,49)=(  1.00000000,0.0d0)                                                                                                           ! 
! zqp(1,3,49)=(  0.00000000,0.0d0)                                                                                                           ! 
! zqp(2,3,49)=(  1.00000000,0.0d0)                                                                                                           ! 
! zqp(1,4,49)=(  1.00000000,0.0d0)                                                                                                           !
! zqp(2,4,49)=(  0.00000000,0.0d0)                                                                                                           ! 
! zqp(1,5,49)=(  1.00000000,0.0d0)                                                                                                           ! 
! zqp(2,5,49)=(  0.00000000,0.0d0)                                                                                                           ! 
! zqp(1,6,49)=(  1.00000000,0.0d0)                                                                                                           ! 
! zqp(2,6,49)=(  0.00000000,0.0d0)                                                                                                           ! 
! zqp(1,7,49)=(  1.00000000,0.0d0)                                                                                                           ! 
! zqp(2,7,49)=(  0.00000000,0.0d0)                                                                                                           !
!                                                                                                                                            !
! which means:                                                                                                                               !
!   49  0  1  1  0  0  0  0                                                                                                                  !
! the output will be written in the full_basis..out file as:                                                                                 !
! zqpf(1,7,49)=(1.0d0,0.0d0) (0.0d0,0.0d0) (0.0d0,0.0d0) (1.0d0,0.0d0) (1.0d0,0.0d0) (1.0d0,0.0d0) (1.0d0,0.0d0)                             !
! zqpf(2,7,49)=(0.0d0,0.0d0) (1.0d0,0.0d0) (1.0d0,0.0d0) (0.0d0,0.0d0) (0.0d0,0.0d0) (0.0d0,0.0d0) (0.0d0,0.0d0)                             !
!____________________________________________________________________________________________________________________________________________!

! The subroutine (( setup_zqp_full_basis )) performs the following tasks:
   ! 1- Generates all possible sets of 1s and 0s (combinations) and writes them to an output file named one_zero.out.
   ! 2- Reads the combinations from one_zero.out, calculates corresponding values, and writes the full basis to an output file named full_basis.out.                                                                 

module Merfat_zqp_fullbasis
  use Merfat_arr_prm
  ! use Merfat_one_zero_state
  implicit none

contains

subroutine setup_zqp_full_basis

    integer, allocatable :: n_new(:), n_old(:)
    integer :: dummy_int(n_df + 1)   ! dummy_int is an array containing values read from the file "one_zero.out", 
    integer :: i_db, i_df, i_basis
    integer :: i_0, i_1, i_db_a_0, i_db_a_1, i_db_old, i_df_old
    real(kind=8):: a0_in, a1_in
    integer :: n_df_a, n_db_a, n_df_old, n_db_old
    integer :: io   ! Declare io variable

    open(unit=11121, file='print/setup_zqp_full_basis.out', status='unknown',access= 'append')
  ! ________________________________________ 

    !  allocate(n_new()) and deallocate(n_new()): 
    ! Generatying all possible sets of 1 and 0 (Combination Generation)
    ! initializes n_df_a to 1 and n_db_a to 2,
    n_df_a = 1
    n_db_a = 2

    ! allocates memory for the array n_new with the size n_df_a * n_db_a, 
    ! and then assigns values to the first two elements of n_new.
    allocate(n_new(n_df_a * n_db_a))
    n_new(1) = 0
    n_new(2) = 1
!     
    ! Combination Update: This loop effectively updates the size of n_new to accommodate the increasing combinations of 0s and 1s until n_df_a reaches the value of n_df.
    do while (n_df_a < n_df)
      allocate(n_old(n_df_a * n_db_a))
      n_old = n_new


      ! It updates the variables
      ! These mathematical operations are part of the process of dynamically resizing the array n_new to accommodate additional elements as needed. 
      ! The values stored in n_df_old and n_db_old serve to retain the original dimensions of the array before resizing, which may be useful for subsequent calculations or operations.
      n_df_old = n_df_a
      n_db_old = n_db_a
      n_df_a = n_df_a + 1     ! the number of rows in the array n_new is being increased.
      n_db_a = 2 * n_db_a     ! doubles the number of columns in the array n_new, expanding the array horizontally.


      ! It deallocates the memory for n_new.
      ! It allocates memory for n_new with the updated size.
      deallocate(n_new)    ! الغاء التخصيص 
      allocate(n_new(n_df_a * n_db_a))   ! يقوم بتخصيص ذاكرة لـ ٠٠٠ بالحجم المحدث.

! ________________________________________  

  ! Combination Filling [ n_new() ]:
  ! Now fill up new array n_new with combinations of zero state and 1 state based on the values stored in the n_old array
      
  ! iterates over each column index in the old array n_old
      do i_db_old = 1, n_db_old
        i_db_a_0 = 2 * i_db_old - 1
        i_db_a_1 = 2 * i_db_old
  ! iterates over each row index in the old array n_old
        do i_df_old = 1, n_df_old
          i_0 = (i_db_a_0 - 1) * n_df_a + i_df_old
          i_1 = (i_db_a_1 - 1) * n_df_a + i_df_old
  ! Copying Values from n_old to n_new
          n_new(i_0) = n_old((i_db_old - 1) * n_df_old + i_df_old)
          n_new(i_1) = n_old((i_db_old - 1) * n_df_old + i_df_old)
        end do
  ! After copying the existing values, two additional values are set in the new array n_new:
        n_new(i_db_a_0 * n_df_a) = 0    ! This sets the value at index i_db_a_0 * n_df_a to 0, representing the zero state.
        n_new(i_db_a_1 * n_df_a) = 1    ! This sets the value at index i_db_a_1 * n_df_a to 1, representing the one state.
      end do
  ! Repeat for Each Column: This process is repeated for each column in the old array n_old, ensuring that the new array 
  ! n_new is filled with the appropriate combinations of zero and one states for all columns.

  ! Deallocation:
      deallocate(n_old)
    end do
! ________________________________________  

  ! Output to File:
    open(file='output/one_zero.out', unit=294)
    do i_db = 1, n_db
      write(294, '(I4)') i_db
      do i_df = 1, n_df
        write(294, '(I2)') n_new((i_db - 1) * n_df + i_df)
      end do
      write(294, *)
    end do
    close(unit=294)

! ________________________________________
  ! File Operations:

  open(unit=295, file='output/one_zero.out')
  open(unit=290, file='output/full_basis.out')
! ________________________________________

  ! Read from Input File:
    write(*,*) 'the file reads correctly from (one_zero.out) file: '
    write(*,*)
  ! Loop over each quantum state، iterates over each row in the file
  open(unit=998, file='print/zqp_fullbasis.out',status='replace')
  do i_db = 1, n_db

              ! Read the next line from the file
              read(295, *, iostat = io) dummy_int
              ! Check if the end of the file is reached
              if (io /= 0) exit ! Exit the loop if end of file is reached

    ! do i_df = 1, n_df

      ! read(295, *) dummy_int 
      ! Write the quantum state
       write(998,*) 'dummy_int =', dummy_int
       write(998,*)
! ________________________________________  

  ! Full Basis Header:
  ! Write the full basis header

      write(290,*) 
      write(290, '(A,I4)') 'Full basis for n_db:', i_db

      ! Loop over each quantum state in the current line
      do i_df = 1, n_df

  ! Calculate a1_in and a0_in
      ! This line retrieves the value stored in dummy_int at the index i_df + 1. 
      ! Since array indices typically start from 1 in Fortran, i_df + 1 allows accessing each element of dummy_int sequentially.
      ! a0_in and a1_in equals 1.0
        a1_in = dummy_int(i_df + 1) 
        if (a1_in**2 > 1.0d0) then
            print*, "Error: a1_in**2 exceeds 1.0, a1_in=", a1_in
        end if
      ! Normalisation: because we dealing with orthogonal vectors the sum of the squares of the coefficients should equal 1.0:[ a1_in**2 + a0_in**2 = 1.0 ] 
      ! Ensuring that the sum of their squares equals 1.0 ensures that the vector they represent has a magnitude of 1, making it a unit vector.
        a0_in = sqrt(1.0d0 - a1_in**2)
        if (a0_in < 0.0d0) then
            print*, "Error: sqrt argument is negative, a0_in=", a0_in
            stop
        end if  


      write(998,*) 'a1_in =', a1_in, 'a0_in =', a0_in
      write(998,*) 
! ________________________________________     
      
! Write zqp(1, i_df, i_db)
      write(290, '(A,I0,A,I0,A,F6.1,A,F11.8)') &
        'zqp(1,', i_df, ',', i_db, ')=(', a0_in, 'd0,0.0d0) '

        zqp(1,i_df, i_db)=a0_in

! Write zqp(2, i_df, i_db)
      write(290, '(A,I0,A,I0,A,F6.1,A,F11.8)') &
        'zqp(2,', i_df, ',', i_db, ')=(', a1_in, 'd0,0.0d0) '

          zqp(2,i_df, i_db)=a1_in

! Write empty line
        write(998,*)
    

! File Closure and Deallocation:
    end do

  end do
! ________________________________________   
  close(unit=998)
  close(unit=295)
  close(unit=290)


end subroutine setup_zqp_full_basis
end module Merfat_zqp_fullbasis



