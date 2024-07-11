
!__________________________________________________________________________________________________________________________________________!
!                                                                                                                                          !
!                                         subroutine step_t  (makes time step using RK4 method)                                            !
!                                                                                                                                          !
!  1- The Runge-Kutta (RK4) method is a numerical technique used to solve ordinary differential equations (ODEs) for variables zqp and s_s !
!  over a time step dt..                                                                                                                   !
!  2- It's particularly popular for its simplicity and accuracy. The method involves taking four steps to estimate the next value of the   !
!  solution based on the current value and the derivative at various points within the time step.                                          !
!                                                                                                                                          !
!                                                                                                                                          !
! Key Steps in step_c Subroutine:                                                                                                          !
! --------------------------------                                                                                                         !
!  A-  Variable Initialization: Initialize time variables t represents the current time t_old stores the previous time and t_common is a   !
!      common time variable                                                                                                                !
!      t = t_common                                                                                                                        !
!      t_old = t                                                                                                                           !
!      t_common = t                                                                                                                        !
!  B- Copying State Variable: Copy the current state variables (zqp, s_s, d_big) to temporary variables (zqp_old, s_s_old, d_old).         !
!      zqp_old = zqp                                                                                                                       !
!      s_s_old = s_s                                                                                                                       !
!      d_old = d_big                                                                                                                       !
!  C- Numerical Integration (RK4 Stages): The subroutine proceeds through four stages (N1, N2, N3, N4) of the RK4 integration:             !
!      Stage N1: Compute derivatives using drv_t and drv_c. Update state variables using a weighted sum of increments.                     !
!      Stage N2: Similar to Stage N1 but using updated values from N1.                                                                     !
!      Stage N3: Similar to Stage N2 but using updated values from N2.                                                                     !
!      Stage N4:Similar to Stage N3 but using updated values from N3                                                                       !
!  D- Final Update and Time Update: Update the final state variables (s_s, zqp) based on increments computed in the four stages.           !
!      Update the time variables for the next time step                                                                                    !
!      ds_s4 = dt * ds_sdt                                                                                                                 !
!      s_s = s_s_old + ds_s1 / 6 + ds_s2 / 3 + ds_s3 / 3 + ds_s4 / 6                                                                       !
!      dzqp4 = dt * dzqpdt                                                                                                                 !
!      zqp = zqp_old + dzqp1 / 6 + dzqp2 / 3 + dzqp3 / 3 + dzqp4 / 6                                                                       !
!------------------------------------------------------------------------------------------------------------------------------------------!
!    This subroutine performs one step of  numerical integration using the Runge-Kutta 4th order (RK4) method.                             !
!             It calculates the evolution of variables `zqp` and `s_s` over a time step `dt`.                                              !
!             It involves solving a set of differential equations using a Runge-Kutta method.                                              !
!             The subroutine calls the `drv_t` subroutine with the current values of t, zqp, s_s, dzqpdt, ds_sdt, and n_flag as arguments. !
!             The purpose of `drv_t` is to compute the derivatives of `zqp` and `s_s` with respect to `t`.                                 ! 
!                                                                                                                                          !
! 1. Variable AND Arrays Declarations:                                                                                                     !
!      - Several complex arrays (dc1, dc2, dc3, dc4, dd1, dd2, dd3, dd4, dzqp1, dzqp2, dzqp3, dzqp4, zqp_old, dd_dt, zqpout, dzqpdt).      !
!      - Real arrays (ds_s1, ds_s2, ds_s3, ds_s4, s_sout, ds_sdt).                                                                         !
!      - Scalars (t, t_old, t_common, dt).                                                                                                 !
!      - The Variable represent derivatives or rates of change of zqp and s_s with respect to time (t):(dzqpdt, s_sout, ds_sdt, and dd_dt) !                                 
!      - A complex constant d_big                                                                                                          !
!      - Other variables (i_step, n_flag)                                                                                                  !
!                                                                                                                                          ! 
!  2. Initialization:                                                                                                                      !
!      - dCopy the current state (zqp, s_s, d_big) to temporary variables (zqp_old, s_s_old, d_old).                                       !                                                                      
!      - TSet the current time (t) to t_common.                                                                                            !
!      - Set t_old to t.                                                                                                                   !
!      - n_flag is a flag used in the calls to drv_t and drv_c.                                                                            !
!                                                                                                                                          !
!   3. Numerical Integration (Runge-Kutta):                                                                                                !
!      - The subroutine proceeds through four stages (N1, N2, N3, N4) of a Runge-Kutta integration.                                        !
!      - At each stage, drv_t and drv_c are called to compute derivatives.                                                                 !
!      - The derivatives are used to update the state variables (s_s, zqp, d_big) using a weighted sum of increments.                      !                                                                                                                                 
!                                                                                                                                          !
!    4. Final Update:                                                                                                                      !
!      - The final state variables (s_s, zqp, d_big) are updated based on the increments computed in the four stages.                      !
!                                                                                                                                          !
!    5. Time Update:                                                                                                                       !
!      - pdate the current time (t) to the next time step.                                                                                 !
!      _ Update t_common to match t. .                                                                                                     !
!                                                                                                                                          !
! 6. The subroutine calls the `drv_t` subroutine with the current values of `t`, `zqp`, `s_s`, `dzqpdt`, `ds_sdt`, and `n_flag` as         !
!    arguments. The purpose of `drv_t` is to compute the derivatives of `zqp` and `s_s` with respect to `t`.                               !
! 4. The subroutine then proceeds with four stages of the RK4 method (N1, N2, N3, N4). For each stage, it calculates intermediate          !
!    values of `ds_s`, `dzqp`, and `dd` using the derivatives from `drv_t` and updates the corresponding variables.                        !
! 5. Finally, it updates `s_s` and `d_big` with the final values computed from the four stages.                                            !
! 6. The subroutine returns after updating `t_common` with the final value of `t`.                                                         ! 
!__________________________________________________________________________________________________________________________________________! 

Module Merfat_step_t

  use Merfat_arr_prm
  use Merfat_drivs
  use Merfat_hamiltonian
  use Merfat_input
  
  implicit none

contains
!______________________________________!
!                                      !
!  Variable and Arrays Declarations    !    
!______________________________________!

subroutine step_c(i_step)

! 1- Several complex arrays
  complex(kind=16), dimension(2, n_df, n_db) :: zqp
  complex(kind=16), dimension(2, n_df, n_db) :: zqp_old     
  complex(kind=16), dimension(2, n_df, n_db) :: zqpout, dzqpdt
  complex(kind=16), dimension(2, n_df, n_db) :: dzqp1, dzqp2, dzqp3, dzqp4 

  complex(kind=16), dimension(n_db) :: d_old
  complex(kind=16), dimension(n_db) :: d_big 
  complex(kind=16), dimension(n_db) :: dc1, dc2, dc3, dc4
  complex(kind=16), dimension(n_db) :: dd1, dd2, dd3, dd4, dd_dt 
 
   
! 2- Real arrays
  real(kind=8), dimension(n_db) :: ds_s1, ds_s2, ds_s3, ds_s4 
  real(kind=8), dimension(n_db) :: s_s, ds_sdt   
  real(kind=8), dimension(n_db) :: s_s_old  
  real(kind=8), dimension(n_db) :: s_sout 

! 3- Scalars
  real(kind=16) :: t                          ! Input parameter for time
  real(kind=8)  :: dt                         ! Input parameter for time step
  real(kind=8)  :: t_old, t_common            ! Local variables for time

! 4- Other variables 
  integer :: i_step
  integer :: n_flag
  integer :: nspinorb, n_df, n_db

!______________________________________!
!             Key Steps                !
!       Variable  Initialization       !    
!______________________________________!

! Initialize ds_s1, dzqp1, dd1 based on initial derivatives
    !  ds_s1 = dt * ds_sdt
    !  dzqp1 = dt * dzqpdt
    !  dd1 = dt * dd_dt

     d_old = d_big
     s_s_old = s_s
     zqp_old = zqp

! Copy the current state (zqp, s_s, d_big) to temporary variables (zqp_old, s_s_old, d_old). Set the current time (t) to t_common. Set t_old to t.

     t = t_common
     t_old = t
	   t_common = t

!____________________________________________!
!   Numerical Integration (Runge-Kutta)      !
!                                            !
! To solve a system of ordinary differential !
!   equations (ODEs) over a specified time   !
!              interval (dt).                !
!                                            !
!      RK4 Stages (N1, N2, N3, N4)           !
! !Update state variables using RK4 method   !   
!____________________________________________!

! CALL N1
  call drv_t(t, zqp, s_s, dzqpdt, ds_sdt, n_flag)
  call drv_c(t, zqp, s_s, dzqpdt, dd_dt, n_flag)

     ds_s1 = dt * ds_sdt
     s_s = s_s_old + ds_s1 / 2
     
     dzqp1 = dt * dzqpdt
     zqp = zqp_old + dzqp1 / 2

     dd1 = dt * dd_dt
     d_big = d_old + dd1 / 2

     t = t_old + dt / 2
     t_common = t


! Update state variables using RK4 method
! CALL N2
  call drv_t(t, zqp, s_s, dzqpdt, ds_sdt, n_flag)
  call drv_c(t, zqp, s_s, dzqpdt, dd_dt, n_flag)
       
     ds_s2 = dt * ds_sdt
     s_s = s_s_old + ds_s2 / 2
 
     dzqp2 = dt * dzqpdt
     zqp = zqp_old + dzqp2 / 2                ! use a weighted average of derivatives over a time step

     dd2 = dt * dd_dt
     d_big = d_old + dd2 / 2

     t = t_old + dt / 2
     t_common = t


! Update state variables using RK4 method
! CALL N3
  call drv_t(t, zqp, s_s, dzqpdt, ds_sdt, n_flag)
  call drv_c(t, zqp, s_s, dzqpdt, dd_dt, n_flag)

  ds_s3 = dt * ds_sdt
  s_s = s_s_old + ds_s3
  
  dzqp3 = dt * dzqpdt
  zqp = zqp_old + dzqp3

  dd3 = dt * dd_dt
  d_big = d_old + dd3

  t = t_old + dt
  t_common = t


! Final Update and Time Update:
! CALL N4
  call drv_t(t, zqp, s_s, dzqpdt, ds_sdt, n_flag)
  call drv_c(t, zqp, s_s, dzqpdt, dd_dt, n_flag)

  ds_s4 = dt * ds_sdt
  s_s = s_s_old + (ds_s1 / 6) + (ds_s2 / 3) + (ds_s3 / 3) + (ds_s4 / 6)
  
  dzqp4 = dt * dzqpdt
  zqp = zqp_old + (dzqp1 / 6) + (dzqp2 / 3) + (dzqp3 / 3) + (dzqp4 / 6)

  dd4 = dt * dd_dt
  d_big = d_old + (dd1 / 6) + (dd2 / 3) + (dd3 / 3) + (dd4 / 6)

  t = t_old + dt
  t_common = t      

end subroutine step_c

end Module Merfat_step_t

