!___________________________________________________________________________________________________________________________________________!
!                                                                                                                                           !
!                                      subroutine acf(acf_t): auto-correlation function                                                     !
!               The auto-correlation function is a mathematical tool used to measure how a signal correlates with a                         ! 
!                                       delayed version of itself over different time intervals                                             !
!                                                                                                                                           !  
!        ** ACF has an Applications in Quantum Mechanics ** :
!    Quantum Mechanics: In quantum systems, the ACF can be used to analyze the stability and coherence properties of wavefunctions or qubits over time.
!
!
!        ** Numerical Computation ** :
!    In the provided Fortran code snippet, the ACF is computed through the overlap of quantum states. The process involves:
!    1- Initializing the states.
!    2- Calculating overlaps between different states using a function (presumably ovlap_ij_z).
!    3- Summing these overlaps weighted by coefficients (a_big).
!
!    The ACF provides information about the correlation between the quantum states of the system at different points in time                ! 
!    based on the input quantum states and coefficients.                                                                                    ! 
!                                                                                                                                           !                               
!                                                       *                                                                                   !
!                              < Ψ(t)|Ψ(0)> =  segma (a)  < z| z >       , Ψ(0)= z                                                          !
!                                                   k   k    k  init               init                                                     !
!                                                                                                                                           !                                                                                                                                          !
!       1-  Main Calculation of "auto-correlation function,"                                                                                !
! For each basis function (i_db), calculates the overlap between the quantum states z_i and z_j using the subroutine overlap_ij_z.          !
!                                                                                                                                           !                                                  
!       2-  Overall Goal:                                                                                                                   !
! The overall goal of the acf subroutine is to quantify the temporal correlation between the quantum states of the system. By computing the !
! auto-correlation function, the module provides insights into how the quantum states evolve over time and how they are influenced by the   !
! system's dynamics and interactions.                                                                                                       !
!                                                                                                                                           ! 
!        3- There are two loops serve different purposes and operate on different arrays.                                                   !
!            The first loop initializes z_i for all degrees of freedom.                                                                     !
!            The second loop calculates the overlap and updates acf_t for all basis functions and degrees of freedom.                       !
!                                                                                                                                           !
! Here's a breakdown of the subroutine's steps:                                                                                             !
!-----------------------------------------------                                                                                            !                     
!                                                                                                                                           !
!        1- Variable Declarations:                                                                                                          !
!    - overlap, ovlp_ij:   Complex variables used to compute overlaps between quantum states.                                               !
!    - acf_t:              Complex variable representing the auto-correlation function.                                                     !
!    - z_i, z_j:           Complex arrays representing quantum states.                                                                      !
!    - a_big:              Real array storing coefficients for each basis function. amblitude                                               !         
!    - i_df, k, i_db:      Integer loop variables.                                                                                          !
!                                                                                                                                           !
!        2- Variable Initialization:                                                                                                        !
! Initialize the complex variable acf_t to zero.                                                                                            !
!                                                                                                                                           !
!        3- Loop Over Basis Functions and Spin Components (z_i):                                                                            !
! This loop initializes the z_i array, which represents the degrees of freedom for the initial state. It iterates over the degrees of       !
! freedom (n_df) and for each degree of freedom, it assigns the values of the zinit array (representing the initial state) to the z_i array.!
!  This loop sets up the initial conditions for the system.                                                                                 !                                                          
!                                                                                                                                           !
!        4- Loop Over Database Points (z_j):                                                                                                !
! For each database point, loop over basis functions and spin components and initialize the array z_j with values from the zqp array.       !
!                                                                                                                                           !
!        5- Overlap Calculation:                                                                                                            !
! Call the overlap subroutine to calculate the overlap (ovlp_ij) between z_i and z_j.                                                       !
!                                                                                                                                           !
!        5- Update Cross-Correlation Matrix (acf_t):                                                                                        !
! Update the cross-correlation matrix acf_t by adding the product of the overlap and the corresponding element from the a_big array.        !
!                                                                                                                                           !
!        7- Return Result:                                                                                                                  ! 
! The subroutine returns the calculated cross-correlation matrix acf_t.                                                                     !
!                                                                                                                                           !
!___________________________________________________________________________________________________________________________________________!                                                                                                                                           

! وظيفة الارتباط التلقائي هي أداة رياضية تستخدم لقياس كيفية ارتباط الإشارة بنسخة متأخرة من نفسها على فترات زمنية مختلفة.

module Merfat_ACF

  use Merfat_arr_prm
  use Merfat_step_t
  use Merfat_overlap

  implicit none

contains

!=============================================== 1 ====================================================================!

subroutine acf(acf_t)
  ! Cross correlation matrix

  ! Variable Declarations
  ! complex(kind=16) :: overlap
  ! complex(kind=16) :: ovlp_ij
  complex(kind=16) :: acf_t
  ! complex(kind=16), dimension(2, n_df) :: z_i, z_j
  complex(kind=16), dimension(n_db) :: a_big
  ! integer :: i_df, k, i_db
  integer :: io_status, io_stat
  logical :: file_status

    ! Open the file for writing
    OPEN(UNIT=1298, FILE='output/ACF.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=io_stat)

    ! Check for errors in opening the file
    IF (io_stat /= 0) THEN
        WRITE(*, *) 'Error opening file "ACF_.ut".'
        STOP
    END IF

  ! Initialize acf_t
  acf_t = (0.0d0, 0.0d0)

 ! initializes the array z_i with values from zinit 
 ! for each degree of freedom (n_df) and each spin component (k = 1, 2).
  do i_df = 1, n_df
    do k = 1, 2
      z_i(k, i_df) = zinit(k, i_df)
      write (1298,*) 'i_df =', i_df, ' and ', 'k = ', k 
      write (1298,*) 'z_i(k, i_df) =', z_i(k, i_df)
      write (1298,*) 
    end do
  end do

write (1298,*) ' the second loop' 

 ! initializes the values of z_j based on the zqp array, calculates the overlap, and updates acf_t.
 ! it calculates the overlap between the arrays z_i and z_j 
 ! and iterates for each degree of freedom (n_df) and each basis function (n_db).
  do i_db = 1, n_db
    do i_df = 1, n_df
      do k = 1, 2
        z_j(k, i_df) = zqp(k, i_df, i_db)
        write (1298,*)
        write (1298,*) 'i_df =', i_df, ' and ', 'k = ', k,  'and ', 'i_db =', i_db
        write (1298,*) 'z_j(k, i_df) =', z_j(k, i_df)
      end do
    end do

    ! Calculate overlap between z_i and z_j
    call ovlap_ij_z(z_i, z_j, ovlp_ij) 
    write (1298,*)
    write (1298,*) 'ovlp_ij =', ovlp_ij

    ! Update acf_t
    acf_t = acf_t + ovlp_ij * a_big(i_db)
    write (1298,*)
    write (1298,*) 'SO: acf_t = ' , acf_t
    write (1298,*)

  end do

    ! Check if the file is open before attempting to close it,  is useful to avoid runtime errors that could occur if you try to close a file that is not open. 
    INQUIRE(UNIT=1298, OPENED=file_status)
    IF (file_status) THEN
        CLOSE(UNIT=1298)
    END IF  

end subroutine acf

end module Merfat_ACF
