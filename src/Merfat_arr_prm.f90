

module Merfat_arr_prm
  IMPLICIT NONE

  PUBLIC :: n_df, ndbmax, n_db, n_db_0c, n_db_full, zqp, zinit, n_db_1c

!_____________________________________!
!                                     !
!            Constant                 !
!_____________________________________!
  ! 1- INTEGER Constant:
  INTEGER, PARAMETER :: n_df = 7                                   ! 1- The number of the degree of freedom (the maximum number of electrons times the number of molecular orbitals.)
  INTEGER, PARAMETER :: ndbmax = 2**n_df                           ! 2- The maximum number of determinants (dynamics basis set)
  INTEGER, PARAMETER :: n_db =  128                                ! 3- the number of dynamic basis functions used to represent the wavefunction in a quantum dynamics calculation.
  INTEGER, PARAMETER :: n_db_full =  128  
  INTEGER, PARAMETER :: n_db_0c =  34                              ! 4- The number of dynamic basis functions with zero contact .
  INTEGER, PARAMETER :: n_db_1c =  38                              ! 5- The number of dynamic basis functions with one contact .
  INTEGER, PARAMETER :: N = 10000                                  ! Number of time steps
  ! INTEGER, PARAMETER :: nconfig = 128
  INTEGER, PARAMETER :: nconfig = 6
  INTEGER, PARAMETER :: n_db0 = 34 
  INTEGER, PARAMETER :: ndfmax = n_df                              ! 6- The maximum number of spin orbitals that can be used in the calculation = number of degree of freedom (number of qubits0)
  INTEGER, PARAMETER :: n_ccf = 7                                  ! 7- The number of .
  INTEGER, PARAMETER :: n_integrals = 28                           ! 8- The number of integrals used in the calculation.
  INTEGER, PARAMETER :: n_dt = 10000
  INTEGER, PARAMETER :: n_step = n_dt                              ! 9- Total number of step= 10,000
  INTEGER, PARAMETER :: ntmax = 2                                  ! 10- The maximum number of time steps in a time-dependent calculation.
  INTEGER, PARAMETER :: nicmax = 1                                 ! 11- The maximum number of initial conditions in a time-dependent calculation.
  INTEGER, PARAMETER :: n_compress = 2
  INTEGER :: i_prop = 0
  INTEGER :: i_step = 0
  REAL(kind=8), PARAMETER :: dt = 0.01
  ! INTEGER, PARAMETER :: (e_tol_min =1.0d-10)

  ! 2- REAL Constant:
  REAL(kind=8), PARAMETER  :: e_tol_min = 1.0d-10 
  REAL(kind=8), PARAMETER  :: pinum = 3.141592653589793238462643383279d0
  REAL(kind=8), PARAMETER :: E = 0.0  
  REAL(kind=8) :: e_number
  REAL(kind=8) :: e_rpltn
  REAL(kind=8) :: ierr

!_____________________________________!
!                                     !
!            variable                 !
!_____________________________________!
  ! 1- Integer variable:
  integer :: n_flag
  integer :: k_df
  integer :: i, j, i_12, l
  integer :: i_db, i_df, j_db, j_df, mm_db, nn_db
  integer :: IPIV(n_db)
  integer :: nsmall, ninp, info 
  integer :: k, r
  integer :: k_db_out
  integer :: time_out
  integer :: i_config
  integer :: i_ccf 

  ! 2- Real variable:
  real(kind=16) :: t                          ! Input variable for time
  ! real(kind=8) :: dt                       ! Input parameter for time step
  real(kind=8) :: t_old                      ! Local variables for time
  real(kind=8) :: t_common                   ! Local variables for time
  real(kind=8) :: t0, t1, t2                 ! Time threshold values
  real(kind=8) :: delta                      ! Output variable for delta
  real(kind=8) :: omega                      ! Output variable for omega
  real(kind=8) :: s_lgrng_n_db 
  real(kind=8) :: anorm 
  real(kind=8) :: anorm_0                    ! Output parameter
  real(kind=8) :: i_db_k, i_dump, i_dump_k, n_dump

  ! 3- Complex variable:
  complex(kind=16) :: h_ii, h_1i  
  complex(kind=16):: ovlp_ij 
  complex(kind=16) :: o_ij, o_ij_1
  complex(kind=16) :: o_ji
  complex(kind=16) :: h    
  complex(kind=16) :: h_lh
  complex(kind=16) :: acf_t
  complex(kind=16) :: acf_t_1
  complex(kind=16) :: dvd_t
  complex(kind=16) :: dvd_tanan
  complex(kind=16) :: overlap
  complex(kind=16) :: o_z_ij
  complex(kind=16) :: o_db_ij
  complex(kind=16) :: o_db_ij_1
  complex(kind=16) :: ovlp_1d
  ! complex(kind=16) :: overlap_z_ij
!_____________________________________!
!                                     !
!            Arrays                   !
!_____________________________________!
  ! 1- complex Arrays:

  complex(kind=16),  dimension(n_db) :: a_big
  complex(kind=16),  dimension(n_db) :: d_big 
  complex(kind=16), dimension(n_db) :: d_old                                            ! Amplitude arrays
  complex(kind=16), dimension(n_df) :: rhs
  complex(kind=16), dimension(n_df) :: a2a
  complex(kind=16), dimension(n_df) :: an 
  complex(kind=16), dimension(n_df) :: anan 
  complex(kind=16), dimension(n_df) :: ovlp 
  complex(kind=16), dimension(n_db) :: dd_dt 
  complex(kind=16), dimension(n_ccf):: ccf_t                                           ! Output cross-correlation matrix 
  ! complex(kind=16), dimension(n_db0) :: a_zqp_0_big
  complex(kind=16), dimension(n_db) :: dc1, dc2, dc3, dc4
  complex(kind=16), dimension(n_db) :: dd1, dd2, dd3, dd4

  complex(kind=16), dimension(2, n_df) :: zinit  
  complex(kind=16), dimension(2, n_df) :: zz_qp_i
  complex(kind=16), dimension(2, n_df) :: zz_qp_j
  complex(kind=16), dimension(2, n_df) :: zi
  complex(kind=16), dimension(2, n_df) :: zic
  complex(kind=16), dimension(2, n_df) :: dh_dzi
  complex(kind=16), dimension(2, n_df) :: dh_dzic
  complex(kind=16), dimension(2, n_df) :: zj
  complex(kind=16), dimension(2, n_df) :: zjc 
  complex(kind=16), dimension(2, n_df) :: z_i, z_j                                      

  complex(kind=16), dimension(n_db, n_db) :: hmltn_ij
  complex(kind=16), dimension(n_db, n_db) :: ovlp_t_ij
  complex(kind=16), dimension(n_db, n_db) :: delta_2_h
  complex(kind=16), dimension(n_db, n_db) :: alhs
  complex(kind=16), dimension(nconfig, 2) :: a_config_0

  complex(kind=16), dimension(2,n_df,ndbmax) :: z_orb_n_db                             ! ndbmax=2**n_df
  complex(kind=16), dimension(2,n_df,ndbmax) :: z_orb
  complex(kind=16), dimension(2, n_df, n_db) :: dzqp1, dzqp2, dzqp3, dzqp4
  complex(kind=16), dimension(2, n_df, n_db) :: zqpout, dzqpdt
  complex(kind=16), dimension(2, n_df, n_db) :: zz
  complex(kind=16), dimension(2, n_df, n_db) :: dzzdt
  complex(kind=16), dimension(2, n_df, n_db) :: dzz
  complex(kind=16), dimension(2, n_df, n_db) :: dzz_1
  complex(kind=16), dimension(2, n_df, n_db0):: zqp_0
  complex(kind=16), dimension(2, n_df, n_db) :: zqpinit 
  complex(kind=16), dimension(n_ccf,2, n_df) :: z_ccf
  complex(kind=8),  dimension(2, n_df, n_db) :: zqp

  complex(kind=16), dimension(2, n_df, n_db) :: zqp_old
  complex(kind=16), allocatable :: a_zqp_0_big(:)
  
  ! 2- Real Arrays:
  REAL(kind=8), dimension(n_df):: delta_lh          
  REAL(kind=8), dimension(ndfmax) :: gamma1 
  REAL(kind=8), dimension(n_df):: omega_lh 
  REAL(kind=8), dimension(n_db) ::  s_s_old   
  REAL(kind=8), dimension(n_db) :: ss
  REAL(kind=8), dimension(n_db) :: dss
  REAL(kind=8), dimension(n_db) :: s_lgrng
  REAL(kind=8), dimension(n_db) :: ds_s1, ds_s2, ds_s3, ds_s4
  REAL(kind=8), dimension(n_db) :: s_s, s_sout, ds_sdt
  REAL(kind=8), dimension(n_df, n_df) :: c_lh

end module Merfat_arr_prm
 
!  The term "e_tol" = 1.0d-10 
!  likely stands for "error tolerance" or "tolerance level". In the context of numerical computations or algorithms, 
!  it represents the maximum allowable error or the threshold below which differences between values are considered negligible.

! For example, in numerical methods or simulations, when comparing two values, if the absolute difference between them
!  is less than the error tolerance (e_tol), they are considered equal for practical purposes. This is because numerical 
!  calculations often involve rounding errors or approximations, and it's essential to define a tolerance level to determine 
!  when two values can be treated as equal.                   

  

  



    