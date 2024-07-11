
 !_________________________________________________________________________________________________________!
 !                                                                                                         !
 !                                  subroutine norm(anorm)                                                 !
 !                              Calculate the norm of the system                                           !
 !                                                                                                         !
 !                               Ψ = segma   ( a | z > )                                                   !
 !                                       k,j    k   k                                                      !
 !                                                                                                         !
 !                                                     *   *                                               !
 !                              < Ψ |Ψ > =  segma    (a) (a)   < z | z >       , < z | z >  is overlap     !
 !                                              k,j   k    j      k   j             k   j                  !
 !                                                                                                         !
 !      1- calculating the norm ensures that the wavefunction is properly normalized, and it provides      !
 !         a measure of the probability distribution of a quantum state.                                   !
 !      2- This is essential for interpreting physical observables and understanding the behavior of       !
 !         quantum systems.t the wavefunction is properly normalized, which is a fundamental               !
 !      3- Normalization Condition: A wavefunction Ψ must satisfy the normalization condition:             !   
 !          ∫|Ψ|**2 dr =1                                                                                  !
 !      4- |Ψ|**2 : The squared magnitude of the wavefunction gives the probability density of finding a   !
 !         particle in a particular region of space.                                                       !
 !      5- ∫|Ψ|**2 : The integral of |Ψ|**2 over all space gives the total probability of finding the      !
 !         particle somewhere.                                                                             !
 !                                                                                                         !

!======================================================================================================================
! Subroutine norm
!======================================================================================================================
! This subroutine calculates the norm of a set of basis functions using overlaps and coefficients.

!======================================================================================================================
! Input:
!======================================================================================================================
! - zqp: A complex array of dimensions (2, n_df, n_db) representing quantum parameters.
! - nspinorb: An integer parameter representing the number of spin-orbitals.
! - n_db: An integer parameter representing the number of dynamic basis functions.

!======================================================================================================================
! Output:
!======================================================================================================================
! - anorm: A real parameter representing the calculated norm.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - ovlp_ij: A complex variable representing the overlap between two basis functions.
! - z_i, z_j: Complex arrays of dimensions (2, nspinorb) representing quantum parameters for overlap calculation.
! - a_big: A complex array of dimension (n_db) representing coefficients of basis functions.

!======================================================================================================================
! Initialization:
!======================================================================================================================
! - The norm (anorm) is initialized to 0.0.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to compute the norm of a set of basis functions, which is essential in quantum
! mechanics calculations. The norm represents the square root of the inner product of a function with itself and provides
! a measure of its "length" or "magnitude" in the function space.

!======================================================================================================================
! Main Loop:
!======================================================================================================================
! - The subroutine iterates over all pairs of basis functions (i, j) to compute their overlaps.
! - For each pair (i, j), it constructs matrices z_i and z_j containing quantum parameters for overlap calculation.
! - The overlap between z_i and z_j is calculated using the overlap subroutine.
! - The norm (anorm) is updated using the overlap and coefficients of basis functions.

!======================================================================================================================

module Merfat_norm

  use Merfat_arr_prm
  use Merfat_zqp_1c_state
  use Merfat_drivs
  use Merfat_input
  use Merfat_hamiltonian
  use Merfat_overlap
  use Merfat_init_conditions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

contains

  subroutine norm(anorm)

    ! Output parameter
    real(kind=8) :: anorm

    ! Local variables
    ! complex(kind=16) :: ovlp_ij
    ! complex(kind=16), dimension(2, n_df) :: z_i, z_j
    ! complex(kind=8),  dimension(n_db) :: a_big                     ! a_big(i) and a_big(j) are coefficients for the basis functions used in the norm calculation.
    ! complex(kind=16), dimension(2, n_df, n_db) :: zqp
    ! integer :: i, j, i_df, k

    ! Initialize the norm
    anorm =  0.0d0
  
 
! Loop over the basis functions
do i = 1, n_db
    do j = 1, n_db

        ! Loop over the qubits
         do i_df = 1, n_df
          ! Loop over the components (real and imaginary) of each spin-orbital
            do k = 1, 2

              ! Construct the matrices z_i and z_j for overlap calculation
              z_i(k, i_df) = zqp(k, i_df, i)
              z_j(k, i_df) = zqp(k, i_df, j)

            end do
         end do

        ! Call the overlap subroutine to calculate the overlap between z_i and z_j
        call ovlap_ij_z(z_i, z_j, ovlp_ij)
         
        ! ! Print intermediate values for debugging
        ! write (*, *) 
        ! write (*, *) 'i =',i
        ! write (*, *) 
        ! write (*, *) 'j =',j
        ! write (*, *) 
        ! write (*, *) 'a_big(',i,')', '=',a_big(i)
        ! write (*, *)
        ! write (*, *) 'a_big(',j,')', '=',a_big(j)
        ! write (*, *) 
        ! write (*, *) 'ovlp_ij =',ovlp_ij
        ! write (*, *) 
        ! write (*, *) '--------------------------------------------------------------------------'

        ! ! Before updating anorm, check if ovlp_ij is finite:
        ! if (.not. ieee_is_finite(real(ovlp_ij)) .or. .not. ieee_is_finite(aimag(ovlp_ij))) then
        !   write (*, *) 'Error: Non-finite overlap detected, ovlp_ij =', ovlp_ij
        !   stop
        ! end if

        ! Update the norm using the overlap and coefficients
        anorm = anorm + conjg(a_big(i)) * ovlp_ij * a_big(j)     !/f_gss(i)
    
        ! ! After updating anorm, Check if the updated norm is finite
        ! if (.not. ieee_is_finite(anorm)) then
        !   write (*, *) 'Error: Non-finite norm detected, anorm =', anorm
        !   stop
        ! end if

    end do
end do

  end subroutine norm

end module Merfat_norm
  