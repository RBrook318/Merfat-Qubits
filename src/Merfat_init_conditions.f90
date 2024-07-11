
!======================================================================================================================!
! Subroutine setup_d_big   if they are non_orthogonal 
!======================================================================================================================!
! This subroutine sets up the amplitudes of CSs (configuration state functions) and solves the system of linear equations
! to compute the coefficients of CSs.

!======================================================================================================================!
! Input Variables:
!======================================================================================================================!
! - rhs: A complex array of dimension (n_db) representing the right-hand side of the linear equation system.
! - a_eigen_value_zinit: A complex array of dimension (n_db) representing eigenvalues of zinit.
! - zz_qp_i: A complex array of dimension (2, n_df) representing quantum parameters for configuration state i.
! - zz_qp_j: A complex array of dimension (2, n_df) representing quantum parameters for configuration state j.
! - a_config_0: A complex array of dimension (nconfig) representing configuration state coefficients.
! - alhs: A complex array of dimension (n_db, n_db) representing the left-hand side of the linear equation system.
! - lhs_old: A complex array of dimension (n_db, n_db) representing the previous left-hand side.
! - z_orb: A complex array of dimension (2, n_df, ndbmax) representing orbital parameters.
! - eigen_vectors: A complex allocatable array representing eigen vectors.
! - a_big: A complex array of dimension (n_db) representing the coefficients of CSs.
! - d_big: A complex array of dimension (n_db) representing the coefficients of CSs.
! - o_ij, o_ji: Complex variables representing overlaps between quantum parameters.
! - work: A complex array used for workspace.

!======================================================================================================================!
! Local Variables:
!======================================================================================================================!
! - s_s: A real variable representing the sum.
! - anorm_check: A real variable representing the norm check.
! - LULU: A character variable.
! - IPIV: An integer array representing pivot indices.
! - nsmall, ninp: Integer variables representing dimensions.
! - info: An integer representing the status of the linear solver.
! - i_db, j_db, i_df, i12, i: Integer loop variables.
! - io_stat: An integer representing the I/O status.

!======================================================================================================================!
! Initialization:
!======================================================================================================================!
! - The right-hand side (rhs) is computed based on overlaps between quantum parameters (zz_qp_i, zinit).
! - The left-hand side (alhs) is computed based on overlaps between quantum parameters (zz_qp_i, zz_qp_j).
! - The result is output to a file named "alhs.out".

!======================================================================================================================!
! Main Loop:
!======================================================================================================================!
! - The system of linear equations is solved using zgesv subroutine.
! - The computed coefficients of CSs (a_big) are checked for normalization.
! - Projection on eigen vectors is performed to obtain configuration state coefficients (a_config_0).

!=============================================== 1 ====================================================================!

module Merfat_init_conditions

    use Merfat_arr_prm
    use Merfat_overlap
    use Merfat_input
    ! use Merfat_overlap_zinit_zi
    implicit none

contains

subroutine setup_d_big_orthigonal

    complex(kind=16), dimension(2, n_df) :: zz_qp_i
    complex(kind=16), dimension(2, n_df) :: zinit
    complex(kind=16), dimension(2, i_df, i_db) :: z
    complex(kind=16) :: ovlp_ij  
    integer :: IOSTAT, io_stat, test
    logical :: file_status
    character(len=1024) :: line
    ! Open the file for reading zinit
    OPEN(UNIT=1292, FILE='output/input_zinit.out', STATUS='OLD', ACTION='READ', IOSTAT=io_stat)

    ! Check for errors in opening the file
    IF (io_stat /= 0) THEN
        PRINT *, 'Error: Unable to open file output/input_zinit.out'
        PRINT *, 'I/O Status:', io_stat
        ! Handle the error as needed, e.g., stop the program, attempt to open a different file, etc.
        STOP
    END IF

    ! Skip header lines
    READ(1292, *)
    READ(1292, *)

    ! Read zinit values from the file
    ! test = 1
    do i_df = 1,n_df
        READ(1292, '(A)') line! Read the odd/even qubit line (assuming it's a comment or label line)
        READ(1292, '(A)') line! Read the description line (assuming it's a comment or label line)
        ! write(*,*) 'i = ', i_df
        ! Read and parse zinit(1, i_df)
        READ(1292, '(A)') line
        ! write(*,*) line
        line = line(32:)
        ! write(*,*) 'chopped up line1 = ', line
        READ(line, *) zinit(1, i_df)
        ! write(*,*) 'value 1 saved ', zinit(1, i_df)
        ! Read and parse zinit(2, i_df)
        READ(1292, '(A)') line
        line = line(32:)
        ! write(*,*) 'chopped up line2 = ', line
        READ(line, *) zinit(2, i_df)
        ! write(*,*) 'value 2 saved = ', zinit(2, i_df)
    end do


    ! Check if the file is open before attempting to close it
    INQUIRE(UNIT=1292, OPENED=file_status)
    IF (file_status) THEN
        CLOSE(UNIT=1292)
    END IF

    !  Open the file for writing
    OPEN(UNIT=1220, FILE='output/d_big_out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=io_stat)

    ! Check for errors in opening the file
    IF (io_stat /= 0) THEN
        WRITE(*, *) 'Error opening file "d_big_out".'
        STOP
    END IF
    
!------------------------------------------------------

! Main loop to calculate d_big
do i_db = 1, n_db
      do i_df = 1 , n_df
                ! Construct z_i
                z_i(1, i_df) = zqp(1, i_df, i_db)
                z_i(2, i_df) = zqp(2, i_df, i_db)
                ! Construct z_j
                z_j(1, i_df) = zinit(1, i_df)
                z_j(2, i_df) = zinit(2, i_df) 
       
! Write data to file
    !    write(1220, '(A, I3, A, I3)' ) 'i_df = ', i_df, ' and i_db = ', i_db
       write(1220,*) 'i_df =', i_df, 'and', '    ', 'i_db =', i_db
       write(1220,*)
       write(1220,*)' z_i(1, i_df) =', z_i(1, i_df)
       write(1220,*)' z_i(2, i_df) =', z_i(2, i_df)
       write(1220,*)
       write(1220,*) 'z_j(1, i_df) =', z_j(1, i_df)
       write(1220,*) 'z_j(2, i_df) =', z_j(2, i_df)
       write(1220,*)
       write(1220,*) '----------------------------------------'
       end do 

    call ovlap_ij_z(z_i, z_j, ovlp_ij)
        d_big(i_db) = ovlp_ij

write(1220,*) 'd_big(i_db) =', d_big(i_db) 
end do


    ! Check if the file is open before attempting to close it,  is useful to avoid runtime errors that could occur if you try to close a file that is not open. 
    INQUIRE(UNIT=1220, OPENED=file_status)
    IF (file_status) THEN
        CLOSE(UNIT=1220)
    END IF 

end subroutine setup_d_big_orthigonal

!=============================================== 2 ====================================================================!
! Subroutine setup_zinit_from_file
!======================================================================================================================!
! This subroutine reads initial complex numbers from a file and sets them as the initial quantum parameters (zinit).

!======================================================================================================================!
! Input Variables:
!======================================================================================================================!
! - None.

!======================================================================================================================!
! Output Variables:
!======================================================================================================================!
! - zinit: A complex array of dimensions (2, n_df) representing the initial quantum parameters.

!======================================================================================================================!
! Local Variables:
!======================================================================================================================!
! - x_inp, y_inp: Complex variables to store the read values from the file.
! - anorm: Real variable to store the normalization factor.
! - i_df: An integer representing the index of the degree of freedom.

!======================================================================================================================!
! Main Goal:
!======================================================================================================================!
! The main goal of this subroutine is to read initial complex numbers from a file ('Zinit_out') and set them as the initial
! quantum parameters (zinit). These parameters are used as the starting point for quantum calculations.

!======================================================================================================================!
! Output:
!======================================================================================================================!
! - The array zinit containing the initial quantum parameters read from the file.
!======================================================================================================================!
subroutine setup_zinit_from_file(n_df, zinit)

    complex(kind=16), dimension(2, n_df) :: zinit 
    real(kind=8) :: x_inp, y_inp, anorm
    integer :: n_df
    ! integer :: i
    open(unit=556, file='print/setup_zinit_from_file.out', status='replace')
    ! Call input_zinit subroutine to get x_inp and y_inp values for each qubit
    write (556,*)
    write (556,*) '1- from call input_zinit(l-12, i) we obtain as before (l-12 = 1, 2): ' 
    call input_zinit(n_df, zinit)
    write (556,*)
    write (556,*)

    write (556,*) '2- Then the output after calculate anorm  zinit(1, i) = (x_inp / anorm),  zinit(2, i) = (y_inp / anorm) :'
    write (556,*)
         
    ! Loop over each qubit
do i = 1, n_df
        ! Read x_inp and y_inp values for the current qubit
        x_inp = real(zinit(1, i))
        y_inp = real(zinit(2, i))

        ! Calculate the normalization factor
        anorm = sqrt(abs(x_inp)**2 + abs(y_inp)**2)
        write (556,*) 'anorm = ', anorm, ' at i= ', i

        ! Check for division by zero
        if (anorm == 0.0d0) then
            print*, "Error: Division by zero"
            stop
        end if

        ! Normalize the complex number
        zinit(1, i) = (x_inp / anorm)
        zinit(2, i) = (y_inp / anorm)
    write (556,*) '(x_inp / anorm)=', zinit(1, i)
    write (556,*) '(y_inp / anorm)=', zinit(2, i)
    write (556,*)
    write (556,*)


  
end do
close(556)
end subroutine setup_zinit_from_file
end module Merfat_init_conditions
       


    

