
module Merfat_calculation
  
  ! Include files
  use Merfat_arr_prm
  use Merfat_norm
  use Merfat_hamiltonian
  use Merfat_ACF 
  use Merfat_zqp_fullbasis
  use Merfat_zqp_0c_state
  use Merfat_zqp_1c_state
  use Merfat_step_t
  use Merfat_drivs
  use Merfat_input
  use Merfat_overlap_zinit_zi
  
  implicit none
contains

!_______________________________________________________________________________________________________________________
!                                                                                                                           
!                                subroutine output(dvd_t, anorm, acf_t, h, i_step)                                                                                                                             
!======================================================================================================================
! Subroutine output
!======================================================================================================================
! This subroutine is responsible for outputting various results and data during the computation.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - anorm: A real variable representing the computed norm.
! - dvd_t: A complex variable representing the computed time-dependent dipole moment.
! - ccf_t: A complex array of dimension (n_ccf) representing the computed time-dependent correlation functions.
! - acf_t: A complex variable representing the computed time-dependent autocorrelation function.
! - h: A complex variable representing the computed Hamiltonian.
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - k_db_out: An integer variable used for outputting data related to dynamic basis functions.
! - dvd_tanan: A complex variable representing the time-dependent dipole moment multiplied by its complex conjugate.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to facilitate the outputting of various results and data obtained during the
! computation process. This includes the norm, time-dependent dipole moment, time-dependent correlation functions,
! time-dependent autocorrelation function, Hamiltonian, and other relevant information.

!======================================================================================================================
! Subroutines and Functions:
!======================================================================================================================
! - output_norm: A subroutine to output the computed norm.
! - output_acf: A subroutine to output the computed time-dependent autocorrelation function.
! - output_h_ord: A subroutine to output the computed Hamiltonian.
! - output_zqp: A subroutine to output data related to dynamic basis functions.
! - output_dvd: A subroutine to output the computed time-dependent dipole moment and its complex conjugate.

!======================================================================================================================                                                                                                              
subroutine output(dvd_t, anorm, acf_t,ccf_t, h , i_step)

    ! Input variables
    real(kind=8) :: anorm
    complex(kind=16) :: dvd_t
    complex(kind=16), dimension(n_ccf) :: ccf_t
    complex(kind=16) :: acf_t, h
    integer :: i_step
    integer :: iostat
    integer :: k_db_out

    ! Dummy calls to assumed procedures
    call output_norm(anorm, i_step)     
    call output_acf(acf_t, i_step)      
    call output_h_ord(h, i_step)       
    k_db_out = 1                        
    call output_zqp(i_step, k_db_out)   
    call output_e_number(e_number, i_step)   
    call output_d_a(i_step)
    call output_ccf(ccf_t, i_step)     
    call output_dvd(i_step, dvd_t)

end subroutine output


!_____________________________________________________________________________________________________________________
!                                                                                                                             
!                                1- subroutine output_norm(anorm, i_step)                                                      
!                                                                                                                             
!======================================================================================================================
! Subroutine output_norm
!======================================================================================================================
! This subroutine is responsible for outputting the computed norm to a text file at each computation step.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - anorm: A real variable representing the computed norm.
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - time_out: A real variable representing the current time based on the computation step and time increment.
! - acf_t, h: Complex variables representing additional data (not utilized in this subroutine).
! - ierr: An integer representing the I/O status.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to output the computed norm to a text file named "output_norm.txt" at each computation
! step. This allows for the tracking and analysis of the norm over time during the computational process.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file named "output_norm.txt" containing the computation step, time, and computed norm at each step.

!======================================================================================================================                                                                                                            
subroutine output_norm(anorm, i_step)
  
  real(kind=8) :: anorm
  real(kind=8) :: time_out
  complex(kind=16) :: acf_t, h
  integer :: i_step
  integer :: iostat

  time_out = dt * i_step

  ! Open the file for writing
  open(unit=601, file='output/calcu_norm.out', status='unknown', action='write', form='formatted', iostat=iostat)
  if (iostat /= 0) then
    write(*, *) 'Error opening file calcu_norm.out'
  end if

  ! Write to the file
  write(601, '(1x, i10, 1x, e15.8, 1x, e15.8)') i_step, time_out, anorm

  ! Write an empty line
    write(601, *)

  ! Close the file
    close(unit=601)

end subroutine output_norm
 

!______________________________________________________________________________________________________________________
!                                                                                                                             
!                            2- subroutine output_acf(acf_t, i_step)                                             
!                                                                                                                                                                                                                                        
!======================================================================================================================
! Subroutine output_acf
!======================================================================================================================
! This subroutine calculates and outputs the autocorrelation function (ACF) at each computation step.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - acf_t: A complex variable representing the computed time-dependent autocorrelation function.
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - t: A real variable representing the current time based on the computation step and time increment.
! - eigen_values: A real array containing the eigenvalues.
! - acf_t_1: A complex variable representing the analytically calculated autocorrelation function.
! - a_config_0: A complex array representing the configuration state.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to calculate and output the autocorrelation function (ACF) at each computation step.
! This ACF represents the correlation between a system's property at a given time and its property at a later time.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file containing the computation step, time, real part, imaginary part, magnitude of the computed ACF, real part,
!   imaginary part, and magnitude of the analytically calculated ACF.

!======================================================================================================================
subroutine output_acf(acf_t, i_step)

  real(kind=8) :: anorm
  real(kind=16) :: t
  real(8), dimension(nconfig) :: eigen_values
  complex(kind=16) :: acf_t_1
  complex(kind=16) :: acf_t
  complex(kind=16) ::h
  complex(kind=16), dimension(nconfig, 2) :: a_config_0
  integer :: i_step
  integer :: i_config
  integer :: io_status
  integer :: iostat

  ! Open the file for writing
    open(unit=449, file='output/calcu_ACF.out', status='replace', action='write', form='formatted', iostat=iostat)
    if (iostat /= 0) then
        write(*,*) 'Error opening file "calcu_ACF.out".'
    end if

  ! Calculate acf analytically
  ! Initialize acf_t_1
  acf_t_1 = (0.0d0, 0.0d0)
  t = dt * i_step

do i_config = 1, nconfig
       acf_t_1 = acf_t_1 + conjg(a_config_0(i_config, 1)) &
       * a_config_0(i_config, 1) &
       * exp(-cmplx(0.0d0, 1.0d0, kind=16) * eigen_values(i_config) * t)
end do

  ! Write to the file
    write(449, 1) i_step, dt * i_step, &
    real(acf_t), aimag(acf_t), &
    sqrt(real(acf_t)**2 + aimag(acf_t)**2), &     ! computes the magnitude of the complex number acf_t
    real(acf_t_1), aimag(acf_t_1), &
    sqrt(real(acf_t_1)**2 + aimag(acf_t_1)**2)    ! computes the magnitude of the complex number 
    1 format(1x, i6, 7(2x, e22.15))  

    ! Write an empty line
    write(449, *)

    ! Close the file
    close(unit=449)

end subroutine output_acf


!______________________________________________________________________________________________________________________
!
!                                 3- Subroutine output_ccf
! 
!======================================================================================================================
! This subroutine calculates and outputs the cross-correlation function (CCF) at each computation step.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - ccf_t: A complex array representing the computed time-dependent cross-correlation function.
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - t: A real variable representing the current time based on the computation step and time increment.
! - eigen_values: A real array containing the eigenvalues.
! - ccf_t_1: A complex variable representing the analytically calculated cross-correlation function.
! - a_config_0: A complex array representing the configuration state.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to calculate and output the cross-correlation function (CCF) at each computation step.
! This CCF represents the correlation between the properties of two different systems at a given time.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file containing the computation step, time, real part, imaginary part, magnitude of the computed CCF, real part,
!   imaginary part, and magnitude of the analytically calculated CCF.
!======================================================================================================================
subroutine output_ccf(ccf_t, i_step)
    
    integer :: i_step
    complex(kind=16), dimension(n_ccf) :: ccf_t
    complex(kind=16), dimension(n_ccf) :: ccf_t_1  
    real(kind=8) :: dt = 0.0d0
    integer :: i_ccf
    integer :: iostat

    ! Open the file for writing
    open(unit=543, file='output/calcu_CCF.out', status='replace', action='write', form='formatted', iostat=iostat)
    if (iostat /= 0) then
        write(*,*) 'Error opening file "calcu_CCF.out".'
    end if

    ! Write the time step
    write(543, '(1x, i6, 2x, e22.15)') i_step, dt * i_step
    write(543, *)

    ! Write the real and imaginary parts and magnitude of each element of ccf_t
  do i_ccf = 1, n_ccf
      write(543, '(3x, 3(2x, e22.15))') real(ccf_t(i_ccf)), aimag(ccf_t(i_ccf)), &
                                        sqrt(real(ccf_t(i_ccf))**2 + aimag(ccf_t(i_ccf))**2)
      write(543, *)
  end do

    ! Write an empty line
    write(543, *)

    ! Close the file
    close(unit=543)

  end subroutine output_ccf


!______________________________________________________________________________________________________________________
!                                                                                                                             
!                   4- subroutine output_h_ord(h, i_step)                                               
!                                                                                                                                                                                                                                          !
!======================================================================================================================
! Subroutine output_h_ord
!======================================================================================================================
! This subroutine is responsible for outputting the ordered Hamiltonian (h_ord) at each computation step.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - h: A complex variable representing the computed Hamiltonian.
! - i_step: An integer representing the current computation step.
! - e_rpltn: A real variable representing energy resolution.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - time_out: A real variable representing the current time based on the computation step and time increment.
! - ierr: A real variable representing the I/O status.
! - i_db: An integer representing the dynamic basis function index.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to output the ordered Hamiltonian (h_ord) at each computation step.
! The ordered Hamiltonian represents the sorted energies of the system's dynamic basis functions.
! This information is crucial for understanding the energy structure and dynamics of the system being studied.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file containing the computation step, time, and ordered Hamiltonian (h_ord) for each dynamic basis function.

!======================================================================================================================
subroutine output_h_ord(h, i_step)

  complex(kind=16) :: h,  acf_t
  real(kind=8) :: e_rpltn
  real(kind=8) :: time_out
  integer :: i_step
  integer :: i_db
  integer :: iostat

  ! Open the file for writing
    open(unit=986, file='output/calcu_h_ord.out', status='replace', action='write', form='formatted', iostat=iostat)
    if (iostat /= 0) then
        write(*,*) 'Error opening file "calcu_h_ord.out".'
    end if

  ! Write to the file
time_out = dt * i_step

do i_db = 1, n_db
      write(986, '(1x, i10, 1x, e15.8)') i_step, time_out
      write(986, *) 

      write(986, '(1x, i10, 1x, e15.8)') i_step, time_out
      write(986, *) 

            call h_ord(i_db, i_db, h)
            if (i_db == 1) then
              write (*, *) 'h_ord=', real(h)
            end if

      write(986, '(2(1x, e15.8))') real(h), real(h) + e_rpltn
      write(986, *) 
end do

  ! Write an empty line
    write(986, *)

  ! Close the file
    close(unit=986)

end subroutine output_h_ord


!______________________________________________________________________________________________________________________
!                                                                                                                             
!                   5- subroutine output_zqp(i_step, k_db_out)                                              
!                                                                                                                                                                                                                                          !
!======================================================================================================================
! Subroutine output_zqp
!======================================================================================================================
! This subroutine is responsible for outputting the quantum parameters (zqp) at each computation step for a specific
! dynamic basis.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - i_step: An integer representing the current computation step.
! - k_db_out: An integer representing the index of the dynamic basis for which the quantum parameters are outputted.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - e_number: A real variable representing the summation of squared absolute values of the quantum parameters.
! - i_df: An integer representing the index of the degree of freedom.
! - acf_t, h: Complex variables (not utilized in this subroutine).

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to output the quantum parameters (zqp) for a specific dynamic basis at each
! computation step. Quantum parameters represent the complex coefficients associated with the dynamic basis functions
! in quantum mechanics simulations. Outputting these parameters provides insights into the state of the quantum system
! at different time steps.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file containing the computation step, time, and quantum parameters (zqp) for each degree of freedom
!   associated with the specified dynamic basis.
! - The total energy contribution (e_number) of the quantum parameters for the specified dynamic basis.

!======================================================================================================================
subroutine output_zqp(i_step, k_db_out)

  real(kind=8) :: anorm
  real(kind=8) :: e_number, dt
  complex(kind=16) :: acf_t, h
  integer :: i_df, n_df, i_step, k_db_out
  integer :: iostat

  ! Open the file for writing
    open(unit=136, file='output/calcu_zqp.out', status='replace', action='write', form='formatted', iostat=iostat)
    if (iostat /= 0) then
        write(*,*) 'Error opening file "calcu_zqp.out".'
    end if

  ! Write to the file
  write(136, '(1x, i6, 2x, e15.8)') i_step, dt * i_step

  e_number = 0.0d0
do i_df = 1, n_df
    e_number = e_number + abs(zqp(2, i_df, k_db_out))**2

    write(136, 2) i_df, &
      real(zqp(1, i_df, k_db_out)), &
      aimag(zqp(1, i_df, k_db_out)), &
      abs(zqp(1, i_df, k_db_out)), &
      real(zqp(2, i_df, k_db_out)), &
      aimag(zqp(2, i_df, k_db_out)), &
      abs(zqp(2, i_df, k_db_out)), &
      sqrt(abs(zqp(1, i_df, k_db_out))**2 &
      + abs(zqp(2, i_df, k_db_out))**2)
  2 format(2x, i3, 7(2x, e15.8), $)
 end do

write(136, *) 
write(136, '(2x, e12.5)') s_s(k_db_out)
write(136, *) 
write(136, '(2x, e12.5)') e_number
write(136, *) 

  ! Close the file
    close(unit=136)

end subroutine output_zqp


!______________________________________________________________________________________________________________________
!                                                                                                                             
!                   6- subroutine output_e_number(e_number, i_step)                                              
!                                                                                                                                                                                                                                          !
!======================================================================================================================
! Subroutine output_e_number
!======================================================================================================================
! This subroutine plays a pivotal role in the quantum simulation framework, tasked with computing and recording the total 
! energy associated with the quantum states. By meticulously processing the provided input variables, output_e_number 
! offers valuable insights into the energy dynamics and behavior of the simulated quantum system.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - e_number: A real variable representing the cumulative energy contribution of the quantum states.
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - None: output_e_number primarily operates on the provided input variables without requiring additional local variables 
!   for computation.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The primary objective of this subroutine is to quantify and document the total energy contribution (e_number) of the 
! quantum states throughout the simulation. By recording this energy metric at discrete time intervals, output_e_number 
! facilitates comprehensive analysis of the system's energy distribution and dynamics.

!======================================================================================================================
! Output:
!======================================================================================================================
! - File Content: The output file generated by output_e_number contains a detailed record of the computed energy numbers 
!   (e_number) at each computation step (i_step). These entries serve as crucial data points for monitoring and analyzing 
!   the evolution of the system's energy states.
! - File Format: Each entry in the output file corresponds to a specific computation step, providing a structured presentation 
!   of the system's energy dynamics for efficient analysis and interpretation.

!======================================================================================================================
subroutine output_e_number(e_number, i_step)
 
  integer :: i_step
  real(kind=8) :: e_number
  integer :: iostat

  ! Open the file for writing
  open(unit=1098,file='output/calcu_e_number.out',status='unknown',action='write',form='formatted',iostat=iostat,position='append')
  if (iostat /= 0) then
      print*, "Error opening file alcu_e_number.out"
  end if

  ! Write to the file in the desired format
  write(1098, '(i6, 1x, e12.5)') i_step, e_number

  ! Close the file
  close(1098)

end subroutine output_e_number




!______________________________________________________________________________________________________________________
!                                                                                                                             
!                   7- subroutine output_dvd(i_step, dvd_t, dvd_tanan)                                             
!                                                                                                                                                                                                                                        !
!======================================================================================================================
! Subroutine output_dvd
!======================================================================================================================
! This subroutine is responsible for outputting the time-dependent dipole moment (dvd_t) and its complex conjugate
! (dvd_tanan) at each computation step.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - dvd_t: A complex variable representing the computed time-dependent dipole moment.
! - dvd_tanan: A complex variable representing the computed time-dependent dipole moment multiplied by its complex conjugate.
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to output the time-dependent dipole moment (dvd_t) and its complex conjugate
! (dvd_tanan) at each computation step. The dipole moment is a crucial quantity in quantum mechanics, representing the
! spatial separation of positive and negative charges within a molecule or system. Analyzing its behavior over time
! provides insights into the dynamics and behavior of the system under study.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file containing the computation step, time, real part, and imaginary part of the time-dependent dipole moment
!   (dvd_t) and its complex conjugate (dvd_tanan) at each step.

!======================================================================================================================
subroutine output_dvd(i_step, dvd_t)

  real(kind=8) :: anorm
  complex(kind=16) :: dvd_t
  integer :: i_step
  integer :: iostat

  ! Open the file for writing
  open(unit=445,file='output/calcu_dvd_t.out',status='unknown',action='write',form='formatted',iostat=iostat)
  if (iostat /= 0) then
      print*, "Error opening file calcu_dvd_t.out"
  end if

  ! Write to the file
  write(445,'(1x, i6, 2x, e15.8)')  i_step, dt * i_step
  write(445, *) 
  write(445,'(1x, 4(2x, e15.8))') real(dvd_t), aimag(dvd_t)
  write(445, *) 

  ! Write an empty line
  write(445, *)

  ! Close the file
  close(445)

end subroutine output_dvd


!______________________________________________________________________________________________________________________
!                                                                                                                             
!                   7- subroutine output_d_a(i_step)                                             
!                                                                                                                                                                                                                                        !
!======================================================================================================================
! Subroutine output_d_a
!======================================================================================================================
! This subroutine is responsible for outputting the time-dependent dipole moment magnitude at each computation step.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - i_step: An integer representing the current computation step.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to output the time-dependent dipole moment magnitude at each computation step.
! The dipole moment magnitude is an essential quantity in quantum mechanics, representing the strength of the
! separation of positive and negative charges within a molecule or system. Analyzing its behavior over time
! provides insights into the dynamics and behavior of the system under study.

!======================================================================================================================
! Output:
!======================================================================================================================
! - A text file containing the computation step, time, and magnitude of the time-dependent dipole moment at each step.

!======================================================================================================================
subroutine output_d_a(i_step)

  integer :: i_step
  real(kind=8) :: anorm
  integer :: i_db
  integer :: iostat

    ! Open the file for writing
  open(unit=4545,file='output/calcu_d_a.out',status='unknown',action='write',form='formatted',iostat=iostat)
  if (iostat/= 0) then
      print*, "Error opening file calcu_d_a.out"
  end if

  ! Write the computation step and time
  write(4545, '(1x, i6, 2x, e15.8)') i_step, dt * i_step
  write(4545, *)

  ! Write the absolute values of d_big(i_db)
  do i_db = 1, n_db
    write(4545, '(1x, e12.5)') abs(d_big(i_db))
    write(4545, *)
  end do

  ! Write an empty line
  write(4545, *)

  ! Close the file
  close(4545)

    write (*,*) '---------------------------------------------'
    write (*,*) 'module Merfat_calculation is ---------- done '
    write (*,*) '---------------------------------------------'
  
end subroutine output_d_a

end module Merfat_calculation

