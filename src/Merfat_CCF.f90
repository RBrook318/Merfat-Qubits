


!======================================================================================================================
! 1- Subroutine: ccf
!======================================================================================================================
! Description:
!   This subroutine computes the cross-correlation matrix (CCF) based on the given quantum states.

!======================================================================================================================
! Parameters:
!======================================================================================================================
! - ccf_t: A complex array representing the cross-correlation matrix (CCF).

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - acf_t: A complex variable representing additional data (not utilized in this subroutine).
! - ovlp_ij: A complex variable representing the overlap between quantum states i and j.
! - z_i: A complex array of dimensions (2, nspinorb) representing the quantum states for index i.
! - z_j: A complex array of dimensions (2, nspinorb) representing the quantum states for index j.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
!   The main goal of this subroutine is to compute the cross-correlation matrix (CCF) based on the given quantum states
!   and store the result in the array ccf_t.
!     1- Understanding these correlations is essential for studying the behavior and dynamics of complex quantum systems.
!     2- undamental tool for studying and analyzing the correlations between quantum states, enabling deeper insights into the behavior and properties of quantum systems.

!======================================================================================================================
! Output:
!======================================================================================================================
!   The output of the ccf subroutine is the cross-correlation matrix (ccf_t). This matrix contains information 
!   about the correlation between pairs of quantum states based on their overlaps. Each element of the ccf_t 
!   matrix represents the correlation between two specific quantum states. The values in this matrix are 
!   complex numbers, indicating both the magnitude and phase of the correlation.

!   The ccf_t matrix is computed based on the overlaps between pairs of quantum states (z_i and z_j). 
!   The subroutine calculates the overlap between each pair of states and stores the result in the corresponding element of the ccf_t matrix.

!   After the subroutine executes, the ccf_t matrix contains the computed cross-correlation values, providing insights
!   into the relationships between different quantum states in the system.

!======================================================================================================================
! Algorithm:
!======================================================================================================================
! 1. Initialize the cross-correlation matrix ccf_t to zero.
! 2. Loop over each pair of quantum states represented by indices i_ccf and i_db.
!    a. Extract the quantum states z_i and z_j corresponding to indices i_ccf and i_db, respectively.
!    b. Compute the overlap between quantum states z_i and z_j and store it in ovlp_ij using a separate subroutine.
!    c. Update the corresponding element in the cross-correlation matrix ccf_t based on the computed overlap.
! 3. End of subroutine.

!======================================================================================================================
! حساب دالة الارتباط المتبادل (CCF).
! العلاقة بين الحالات الكمومية المختلفة CCF في سياق ميكانيكا الكم،  يقيس 
! zi and zj   فاذا كان لدينا مجموعتين من الاحالات ممثلة بالمتجهات 
! a-big يمكن حساب التداخل بين هذه الحالات باستخدام الضرب الداخلي  ثم يتم وزن التداخل ببعض المعاملات 
!CCF(i) = segma from j=1 to N-db [<zi|zj> (a-big(j))]


  
module Merfat_CCF

  use Merfat_arr_prm
  use Merfat_zqp_1c_state
  use Merfat_overlap
  use Merfat_overlap_zinit_zi
 

  implicit none

contains
    
subroutine ccf(ccf_t)

    ! Variable Declarations: 
    ! complex(kind=16) :: acf_t 
    ! complex(kind=16) :: ovlp_ij                             ! Overlap between quantum states
    complex(kind=16), dimension(n_ccf) :: ccf_t             ! Output cross-correlation matrix
    ! complex(kind=16), allocatable :: a_zqp_0_big(:)
    ! complex(kind=16), dimension(2, n_df) :: z_i
    ! complex(kind=16), dimension(2, n_df) :: z_j
    ! ! complex(kind=16), dimension(n_ccf, 2, n_df):: z_ccf
    ! integer :: i_ccf, i_df, i_db, k                         ! Loop indices
    integer :: io_stat 
    logical :: file_status    
   

    ! Open the file for writing
    OPEN(UNIT=1234, FILE='output/CCF.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=io_stat)

    ! Check for errors in opening the file
    IF (io_stat /= 0) THEN
        WRITE(*, *) 'Error opening file "CCF.out".'
        STOP
    END IF

    ! Initialize cross-correlation matrix
    ! ccf_t = cmplx(0.0d0, 0.0d0)
    ccf_t = (0.0d0, 0.0d0)

    ! Main Calculation:
    ! Loop through quantum states (! Loop over the cross-correlation states)
    do i_ccf = 1, n_ccf
! _____________________________________________________________

        ! Construct state vector z_i from z_ccf
        do i_df = 1, n_df
            do k = 1, 2
                z_i(k, i_df) = z_ccf(i_ccf, k, i_df)
            end do
        end do
! _______________________________________________________________

        ! Loop through a set of quantum states i_db(Calculate cross-correlation with database states)
        do i_db = 1, n_db
            ! Construct state vector z_j from zqp
            do i_df = 1, n_df
                do k = 1, 2
                    z_j(k, i_df) = zqp(k, i_df, i_db)
                end do
            end do

            ! Calculate overlap between z_i and z_j
            call ovlap_ij_z(z_i, z_j, ovlp_ij) 

            ! Update ccf_t
            ccf_t(i_ccf) = ccf_t(i_ccf) + ovlp_ij * a_big(i_db)
        end do
! _____________________________________________________________

            ! Write to the output file
            write(1234, *) "CCF(", i_ccf, ") = ", ccf_t(i_ccf)

    end do

    ! Check if the file is open before attempting to close it,  is useful to avoid runtime errors that could occur if you try to close a file that is not open. 
    INQUIRE(UNIT=1234, OPENED=file_status)
    IF (file_status) THEN
        CLOSE(UNIT=1234)
    END IF

end subroutine ccf

end module Merfat_CCF


!====================================================================================================================
! Subroutine: ccf_0
!====================================================================================================================
! Description:
!   This subroutine computes the cross-correlation matrix (CCF) based on a modified set of quantum states.
!   It is an alternative implementation to the ccf subroutine, specifically designed for a different set of states.

!====================================================================================================================
! Parameters:
!====================================================================================================================
! - ccf_t: A complex array representing the cross-correlation matrix (CCF).
! - a_zqp_0_big: A complex variable representing additional data associated with the modified quantum states.

!====================================================================================================================
! Local Variables:
!====================================================================================================================
! - acf_t: A complex variable representing additional data (not utilized in this subroutine).
! - ovlp_ij: A complex variable representing the overlap between quantum states i and j.
! - z_i: A complex array of dimensions (2, n_df) representing the quantum states for index i.
! - z_j: A complex array of dimensions (2, n_df) representing the modified quantum states.

!====================================================================================================================
! Main Goal:
!====================================================================================================================
!   The main goal of this subroutine is to compute the cross-correlation matrix (CCF) based on a modified set of 
!   quantum states and store the result in the array ccf_t. This alternative implementation is tailored for 
!   specific scenarios where the quantum states undergo modifications.

!   Understanding these correlations is essential for studying the behavior and dynamics of complex quantum systems.
!   The CCF serves as a fundamental tool for analyzing the correlations between quantum states, providing insights 
!   into the system's behavior and properties.

!====================================================================================================================
! Output:
!====================================================================================================================
!   The output of the ccf_0 subroutine is the cross-correlation matrix (ccf_t). This matrix contains information 
!   about the correlation between pairs of quantum states based on their overlaps. Each element of the ccf_t 
!   matrix represents the correlation between two specific quantum states. The values in this matrix are 
!   complex numbers, indicating both the magnitude and phase of the correlation.

!   The ccf_t matrix is computed based on the overlaps between pairs of quantum states (z_i and z_j). 
!   The subroutine calculates the overlap between each pair of states and stores the result in the corresponding 
!   element of the ccf_t matrix.

!   After the subroutine executes, the ccf_t matrix contains the computed cross-correlation values, providing insights
!   into the relationships between different quantum states in the system.

!====================================================================================================================
! Algorithm:
!====================================================================================================================
! 1. Initialize the cross-correlation matrix ccf_t to zero.
! 2. Loop over each pair of quantum states represented by indices i_ccf and i_db.
!    a. Extract the quantum states z_i and z_j corresponding to index i_ccf and modified quantum states, respectively.
!    b. Compute the overlap between quantum states z_i and z_j and store it in ovlp_ij using a separate subroutine.
!    c. Update the corresponding element in the cross-correlation matrix ccf_t based on the computed overlap.
! 3. End of subroutine.

!====================================================================================================================


! subroutine ccf_0(ccf_t)
    
!     ! Input and output parameters
!     complex(kind=16) :: acf_t 
!     complex(kind=16), allocatable :: a_zqp_0_big(:)
!     complex(kind=16), dimension(n_ccf) :: ccf_t             ! Output cross-correlation matrix
!     complex(kind=16), dimension(n_ccf,2,n_df):: z_ccf
!     complex(kind=16), dimension(2, n_df, n_db0) :: zqp_0
!     complex(kind=16), dimension(2, n_df) :: z_i
!     complex(kind=16), dimension(2, n_df) :: z_j             
!     complex(kind=16) :: ovlp_ij
!     integer :: i_ccf, i_df, i_db, k
    
!     ! Initialize the cross-correlation matrix
!     ccf_t = (0.0d0, 0.0d0)
    
!     ! Loop over each correlation index
!     do i_ccf = 1, n_ccf
!         ! Construct z_i from z_ccf
!         do i_df = 1, n_df
!             do k = 1, 2
!                 z_i(k, i_df) = z_ccf(i_ccf, k, i_df)
!             end do
!         end do
        
!         ! Loop over each basis state
!         do i_db = 1, n_db0
!             ! Construct z_j from zqp_0
!             do i_df = 1, n_df

!                 do k = 1, 2
!                     z_j(k, i_df) = zqp_0(k, i_df, i_db)
!                 end do
            
!             ! Compute overlap between z_i and z_j
!             call ovlap_ij_z(z_i, z_j, ovlp_ij)     
!             ! Update the cross-correlation matrix
!             ccf_t(i_ccf) = ccf_t(i_ccf) + ovlp_ij * a_big(i_db)(i_db)
          
!            end do 
!         end do
!     end do
! end subroutine ccf_0


! end module Merfat_CCF


