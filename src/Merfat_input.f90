
 !_______________________________________________________________________________________________________________!
 !                                                                                                               !
 !                               subroutine input_lukin                                                          !
 !                                                                                                               !
 ! Units of Energy and Time:                                                                                     !
 !     1- Units of energy in MHz 10**6 Hz and the corresponding units of time in microseconds which 10**(-6) sec !               
 !                                                                                                               !
 ! Comments:                                                                                                     !
 !     1- The subroutine focuses on the nearest-neighbor interactions                                            !
 !     2- and includes a normalization step for the coupling matrix.                                             !
 !     3- Everything except the nearest neighbours can be neglected.                                             !
 !                                                                                                               !
 ! Description od Coupling Matrix (c_lh):                                                                        !
 !     1- Initialized to zeros with dimensions (n_df, n_df).                                                     !
 !     2- A loop calculates values for c_lh(i_df, j_df) using a specific formula involving (i_df-j_df)**6.       !
 !     3- The loop iterates over pairs of i_df and j_df such that j_df is greater than i_df to avoid redundancy. !
 !     4- The calculated values are symmetrically assigned to c_lh to ensure the matrix is symmetric.            !
 !     5- Normalization of c_lh: The entire c_lh matrix is divided by 2 after the loop.                          !
 !                                                                                                               !
 ! Calculating and Assigning Values to c_lh:                                                                     !
 !   The loop constructs a symmetric coupling matrix c_lh where each element (i, j) represents the coupling      !
 !   strength between modes or degrees of freedom i and j in the quantum system,                                 !
 !   ensuring symmetry and avoiding redundant calculations.                                                      !
 !     1- c_lh(i_df, j_df) = 2.0d0 * pinum * 24.0d0 / (i_df - j_df)**6: Calculates a value based on the          !
 !        formula 2.0 * π * 24 / (i_df - j_df)^6 and assigns it to the upper triangular part of the matrix.      !
 !     2- c_lh(j_df, i_df) = c_lh(i_df, j_df): Ensures symmetry in the matrix by assigning the calculated value  !
 !        to the corresponding lower triangular element.                                                         !
 !     3- nested loop iterates over the indices j_df from i_df + 1 to n_df. The use of i_df + 1 ensures that     !
 !        the elements below the main diagonal of the matrix are calculated, avoiding redundancy.                !
 !_______________________________________________________________________________________________________________!

!======================================================================================================================
! Subroutine input_lukin
!======================================================================================================================
! This subroutine initializes parameters required for Lukin's model computation.

!======================================================================================================================
! Parameters:
!======================================================================================================================
! - n_df: An integer parameter representing the number of degrees of freedom.
! - e_tol_min: A real parameter representing the minimum tolerance for energy calculations.
! - pinum: A real parameter representing the value of π.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - delta_lh: A real variable representing the delta term in Lukin's model.
! - omega_lh: A real variable representing the omega term in Lukin's model.
! - c_lh: A real array of dimensions (n_df, n_df) representing the coupling matrix in Lukin's model.

!======================================================================================================================
! Initialization:
!======================================================================================================================
! - Delta Term (delta_lh) is initialized to 0.0.
! - Omega Term (omega_lh) is computed using predefined formulae and divided by 2.0 for normalization.
! - Coupling Matrix (c_lh) is initialized to all zeros.

!======================================================================================================================
! Main Loop:
!======================================================================================================================
! - The loop iterates over indices i_df and j_df to populate the elements of the coupling matrix c_lh.
! - Elements of c_lh are computed based on the distance between indices i_df and j_df using a specific formula.
! - The entire matrix c_lh is normalized by dividing it by 2.0.

!=======================================================================================================================
! The resulting 7x7 matrix c_lh contains the calculated values of the coupling constants for all combinations of i_df and j_df. 
! Since c_lh is symmetric, the values of c_lh(i_df, j_df) and c_lh(j_df, i_df) are the same for corresponding elements across
! the main diagonal. Here's an example of what the resulting matrix might look like:
! c_lh = [[ c_lh(1, 1), c_lh(1, 2), c_lh(1, 3), c_lh(1, 4), c_lh(1, 5), c_lh(1, 6), c_lh(1, 7) ],
        ! [ c_lh(2, 1), c_lh(2, 2), c_lh(2, 3), c_lh(2, 4), c_lh(2, 5), c_lh(2, 6), c_lh(2, 7) ],
        ! [ c_lh(3, 1), c_lh(3, 2), c_lh(3, 3), c_lh(3, 4), c_lh(3, 5), c_lh(3, 6), c_lh(3, 7) ],
        ! [ c_lh(4, 1), c_lh(4, 2), c_lh(4, 3), c_lh(4, 4), c_lh(4, 5), c_lh(4, 6), c_lh(4, 7) ],
        ! [ c_lh(5, 1), c_lh(5, 2), c_lh(5, 3), c_lh(5, 4), c_lh(5, 5), c_lh(5, 6), c_lh(5, 7) ],
        ! [ c_lh(6, 1), c_lh(6, 2), c_lh(6, 3), c_lh(6, 4), c_lh(6, 5), c_lh(6, 6), c_lh(6, 7) ],
        ! [ c_lh(7, 1), c_lh(7, 2), c_lh(7, 3), c_lh(7, 4), c_lh(7, 5), c_lh(7, 6), c_lh(7, 7) ]]

!=======================================================================================================

 Module Merfat_input

  use Merfat_arr_prm

  implicit none

contains

!============================= 1 =======================================================================

 subroutine input_lukin

! Declarations
    integer :: i_df, j_df, iostat, io_stat
    ! integer, parameter :: n_df =7
    logical :: file_status
    ! real(kind=8), parameter :: e_tol_min = 1.0d-10
    ! real(kind=8), parameter :: pinum = 3.141592653589793238462643383279d0
    ! real(kind=8) :: delta_lh
    ! real(kind=8) :: omega_lh
    ! real(kind=8), dimension(n_df, n_df) :: c_lh

! Delta Term
    delta_lh = 0.0d0

! Omega Term
    omega_lh = 2.0d0 * pinum * 2.0d0
    omega_lh = omega_lh / 2.0d0

! Coupling Matrix
    c_lh = 0.0d0

    ! Open the file for writing
    OPEN(UNIT=1280, FILE='output/input_lukin.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=io_stat)

        ! Check for errors in opening the file
    IF (io_stat /= 0) THEN
        WRITE(*, *) 'Error opening file "input_lukin.out".'
        STOP
    END IF
              write(1280,*) 'delta_lh =', delta_lh
              write(1280,*)    
              write(1280,*) 'omega_lh =', omega_lh
              write(1280,*)
              write(1280,*) 'THERE ARE  ( 21 )   DIFFERENT c_lh(j_df, i_df):'
              write(1280,*)
              
! responsible for calculating and populating the elements of the coupling matrix c_lh
    do i_df = 1, n_df
        do j_df = i_df + 1, n_df                                                                           !  لانه في ورقة لوكين السيجما تبدا من i<j

           ! Each element of c_lh is calculated according to the following formula" ..  7*7=49 elements
            c_lh(i_df, j_df) = 2.0d0 * pinum * 24.0d0 / (i_df - j_df)**6                                    ! في صفحة ورقة لوكين c_lh = Vij = c/(alfa**6) * (j-i)**6    and    c/(alfa**6) = 2 *  pinum * 24.0d0 
            c_lh(j_df, i_df) = c_lh(i_df, j_df)                                                             ! resulting in a 7x7 matrix c_lh. The values of c_lh(i_df, j_df) will be the same as c_lh(j_df, i_df) due to symmetry.  
                
              write(1280,*)
              write(1280,*) 'i_df =', i_df
              write(1280,*) 'j_df =', j_df
              write(1280,*)
              write(1280,*) 'c_lh(i_df, j_df) =' , c_lh(i_df, j_df)
              write(1280,*) 'c_lh(j_df, i_df) =' , c_lh(j_df, i_df)
              write(1280,*)
              write(1280,*) ' ___________________________________________________' 
              write(1280,*)

        end do
    end do

! Normalization: divides each element of the matrix c_lh by 2.0d0. This operation is often used for normalization.
    c_lh = c_lh / 2.0d0
        write(1280,*)
        write(1280,*) ' There are 7*7 = 49 of c_lh which is:' 
        write(1280,*) 'c_lh =', c_lh
        write(1280,*)
        write(1280,*) ' ___________________________________________________' 

! Check if the file is open before attempting to close it,  is useful to avoid runtime errors that could occur if you try to close a file that is not open. 
    INQUIRE(UNIT=1280, OPENED=file_status)
    IF (file_status) THEN
        CLOSE(UNIT=1290)
    END IF 

 end subroutine input_lukin

!============================= 2 =======================================================================

subroutine input_zinit(n_df, zinit)
  
    integer :: n_df
    complex(kind = 16),intent(inout), dimension(2, n_df):: zinit 
    integer :: IOSTAT, io_stat
    logical :: file_status
    ! integer :: i

    ! Open the file for writing
    OPEN(UNIT=1290, FILE='output/input_zinit.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=io_stat)

        ! Check for errors in opening the file
    IF (io_stat /= 0) THEN
        WRITE(*, *) 'Error opening file "input_zinit.out".'
        STOP
    END IF

write (1290,*)
write (1290,*)' depends on i (qubit) there are 14 different zinit(2, n_df) which is :'
write (1290,*)

    ! Loop over each qubit
    do i = 1, n_df
        ! Assign values based on the qubit index
        if (mod(i, 2) == 0) then
            ! Even qubits
            zinit(1, i) = (1.0d0, 0.0d0)   ! x_inp = 1.0d0, image_part = 0.0d0
            zinit(2, i) = (0.0d0, 0.0d0)   ! y_inp = 0.0d0, image_part = 0.0d0
            write (1290,*) 'even qubit:', 'i=', i
            write (1290,*) 'zinit(1, i) =', zinit(1, i)
            write (1290,*) 'zinit(2, i) =', zinit(2, i)
            write (1290,*)
        else
            ! Odd qubits
            zinit(1, i) = (0.0d0, 0.0d0)   ! x_inp = 0.0d0, image_part = 0.0d0
            zinit(2, i) = (1.0d0, 0.0d0)   ! y_inp = 1.0d0, image_part= 0.0d0
            write (1290,*) 'Odd qubit:', 'i =', i
            write (1290,*) 'zinit(1, i) =', zinit(1, i)
            write (1290,*) 'zinit(2, i) =', zinit(2, i)
            write (1290,*)
            
           
        end if
    end do

   ! Check if the file is open before attempting to close it,  is useful to avoid runtime errors that could occur if you try to close a file that is not open. 
    INQUIRE(UNIT=1290, OPENED=file_status)
    IF (file_status) THEN
        CLOSE(UNIT=1290)
    END IF 

end subroutine input_zinit

end Module Merfat_input
 