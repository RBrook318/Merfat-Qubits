! The specific details of the calculations involve quantum mechanical concepts such as 
! (Hamiltonians, overlaps, and complex exponentials) . 
! The aim is to provide a computationally efficient way to model the time evolution of a quantum system.

Module Merfat_drivs  

  
  use Merfat_arr_prm
  use Merfat_hamiltonian
  use Merfat_input
  use Merfat_overlap

  implicit none

contains

 !_______________________________________________________________________________________________________________________________________!
 !                                                                                                                                       !
 !                              subroutine drv_t(t, zz, ss, dzzdt, dd_dt, n_flag)                                                        !
 !                               calculates the derivatives with respect to time                                                         ! 
 !                                                                                                                                       !
 ! The drv_t subroutine within the module appears to be responsible for computing the derivatives of the quantum states with             !
 !  respect to time (dzz and dss).                                                                                                       !
 !                                                                                                                                       !
 !                                                                                                                                       !
 !           1- The main goal of this subroutine:                                                                                        !
 ! is to compute the time derivatives of the quantum states (dzz) and the coefficients (dss). These derivatives are essential            !
 ! for simulating the time evolution of the quantum system described by the input quantum states (zz) and coefficients (ss).             !
 !                                                                                                                                       !
 !                                                                                                                                       !
 !          2- Calculations:                                                                                                             !
 ! - Initialization:                                                                                                                     !
 !       Initializes the arrays dzz and dzz_1 with zero imaginary and real parts.                                                        !
 ! - Rearrangement:                                                                                                                      !
 !       Rearranges the input 3D array zz into a 2D array zi.                                                                            !
 ! - Calculation of dzz_1:                                                                                                               !
 !       Calculates the derivative of the Hamiltonian with respect to zi using the subroutine h_su2_lukin.                               !
 !       Assigns the calculated derivative values to dzz_1.                                                                              !
 ! - Calculation of dzz:                                                                                                                 !
 !       Computes dzz by taking the complex conjugate of dzz_1 and multiplying it by -i (imaginary unit).                                !
 ! - Calculation of dss:                                                                                                                 !
 !       Computes the derivative of the coefficients (dss) using the Hamiltonian and the input quantum states.                           !     
 !       Subtracts h_ii from dss.                                                                                                        !
 ! - Update of ss:                                                                                                                       !
 !       Updates the coefficients (ss) over the time step using the calculated derivatives dss (by adding `dt` multiplied by `dss`)      !
 !                                                                                                                                       !
 !                                                                                                                                       !
 !        3- Overall Goal:                                                                                                               !
 ! The subroutine aims to calculate the derivatives of the quantum states and coefficients, which are crucial for simulating the         !
 ! time evolution of the quantum system described by the Schrödinger equation. These derivatives are then used to update the             !
 ! quantum states and coefficients for the next time step in the simulation                                                              !
 !                                                                                                                                       !
 !                                                                                                                                       !
 !       4- output:                                                                                                                      !
 ! The drv_t subroutine doesn't directly produce any visible output that is printed to the screen or stored in files.                    !
 ! Instead, it modifies the values of the output arrays dzz and dss, and also updates the values of the input array ss.                  !
 !       4.1- Output Arrays:                                                                                                             !
 ! dzz: This array holds the derivatives of the quantum states with respect to time. It is a complex array of dimension (2, n_df, n_db). !
 ! dss: This array holds the derivatives of the coefficients with respect to time. It is a real array of dimension (n_db).               !
 !       4.2- Modified Input Array:                                                                                                      !
 ! ss: The subroutine updates the input coefficients array ss over the time step. After the subroutine execution, ss will contain        !
 ! updated coefficient values                                                                                                            !
 !_______________________________________________________________________________________________________________________________________!

  subroutine drv_t(t, zz, ss, dzz, dss, n_flag) 
  
! Input arguments 
      complex(kind=16), dimension(2, n_df, n_db) :: zz
      real(kind=8), dimension(n_db):: ss

! Output arguments
     complex(kind=16), dimension(2, n_df, n_db):: dzz
     complex(kind=16), dimension(2, n_df, n_db):: dzz_1
     real(kind=8), dimension(n_db):: dss

! Local variables
     complex(kind=16), dimension(2, n_df) :: zi, dh_dzi
     real(kind=16) :: t
     real(kind=8) :: h_ii, h_1i
     integer :: i, j, i_12
     integer :: i_db, i_df, j_db, j_df, mm_db, nn_db
     integer :: n_flag
     

! Initializing complex arrays with zero real and imaginary parts.
      ! Trajectories
       dzz = (0.0d0, 0.0d0)
       dzz_1 = (0.0d0, 0.0d0)

! Rearranging the 3D array zz into a 2D array zi 
! Calculate dzz_1
do i_db = 1, n_db
       ! Amplitudes
       do i_12 = 1, 2
         do i_df = 1, n_df
           zi(i_12, i_df) = zz(i_12, i_df, i_db) !(Loop to copy values from zz to zi) 
         end do
       end do
    
! calculate the derivative dh_dzi based on the rearranged zi.
         call h_su2_lukin(zi, zj, h_lh)
        
! Assign dh_dzi values to dzz_1
       do i_df = 1, n_df
         dzz_1(1, i_df, i_db) = dh_dzi(1, i_df)
         dzz_1(2, i_df, i_db) = dh_dzi(2, i_df)
       end do
end do

! Calculates the dzz array using dzz_1 and complex conjugation.
        !  dzz = -(0.0d0, 1.0d0) * conjg(dzz_1)
        dzz = 0

! Calculate dss    
do i_db = 1, n_db
! calculate h_ii.
        call h_ord(i, j, h)   
        
! Action, 
       dss(i_db) = 0.0d0
       dss(i_db)=-h_ii  

! Calculate dss for each i_12 and i_df    
      do i = 1, 2
         do i_df = 1, n_df
           dss(i_db) = dss(i_db) - imag(conjg(zz(i, i_df, i_db)) * dzz(i, i_df, i_db))
         end do
      end do

! Subtract h_ii from dss
          dss(i_db) = dss(i_db) - h_ii
end do

! Updates the state variable ss over the time step
    ss = ss + dt * dss

  end subroutine drv_t
  
 !_______________________________________________________________________________________________________________________________________!
 !                                                                                                                                       !
 !                              subroutine drv_c(t, zz, ss, dzzdt, dd_dt, n_flag)                                                        !
 !                               calculates the derivatives with respect to parameter C                                                  !    
 !                                                                                                                                       !
 ! The drv_c subroutine aims to calculate the derivatives of the quantum states (dzzdt) and coefficients (dd_dt) respect to parameter C  !
 !                                                                                                                                       ! 
 !           1- Input Arguments:                                                                                                         ! 
 ! - t: The current time.                                                                                                                !
 ! - zz: Quantum state array, a complex 3D array of shape (2, n_df, n_db).                                                               !
 ! - ss: Coefficients array, a real 1D array of size n_db.                                                                               !
 ! - s_lgrng: Additional coefficients array, a real 1D array of size n_db.                                                               !
 ! - n_flag: An integer flag.                                                                                                            !
 !                                                                                                                                       ! 
 !           2- Output Arguments:                                                                                                        !
 ! - dzzdt: Derivatives of the quantum states with respect to time, a complex 3D array of shape (2, n_df, n_db).                         ! 
 ! - dd_dt: Derivatives of the coefficients with respect to time, a complex 1D array of size n_db.                                       !
 !                                                                                                                                       ! 
 !          3- Local Variables:                                                                                                          !
 ! Various complex arrays used for intermediate calculations.                                                                            !
 ! - hmltn_ij: Matrix representing the Hamiltonian terms.                                                                                !
 ! - ovlp_ij: Matrix representing the overlap terms.                                                                                     !
 ! - delta_2_h: Matrix representing a combination of overlap, Hamiltonian, and time-dependent overlap terms.                             !
 ! - rhs: Right-hand side vector used in solving the system of linear equations.                                                         !
 ! - alhs: Left-hand side matrix used in solving the system of linear equations.                                                         !
 ! - IPIV: Integer array used in solving the system of linear equations.                                                                 !
 ! - Other integer variables for loop indices and other purposes.                                                                        !
 ! - info: Integer variable used to store the output status of the linear equation solver.                                               !  
 !                                                                                                                                       ! 
 !         4- Calculations:                                                                                                              !
 ! - The subroutine calculates various matrices (hmltn_ij, ovlp_ij, ovlp_t_ij, delta_2_h) representing Hamiltonian terms, overlap terms, !
 !   and their derivatives with respect to C.                                                                                            !
 ! - It constructs the left-hand side matrix (alhs) and the right-hand side vector (rhs) of a system of linear equations.                !
 ! - It solves the system of linear equations using the LAPACK subroutine zgesv.                                                         !
 ! - The derivatives of the quantum states (dzzdt) and coefficients (dd_dt) with respect to c are updated based on the solution of       !
 !   the linear equations.                                                                                                               !
 !_______________________________________________________________________________________________________________________________________!

  subroutine drv_c(t, zz, ss, dzzdt, dd_dt, n_flag)

    ! Define input parameters:
     complex(kind=16),dimension(2, n_df, n_db) :: zz
     complex(kind=16),dimension(2, n_df, n_db) :: dzzdt
     complex(kind=16),dimension(n_db) :: dd_dt
     REAL(kind=8), dimension(n_db) :: s_lgrng
     real(kind=8),dimension(n_db) :: ss
     real(kind=8) :: pinum
     integer :: n_flag

    ! Local variables:
     complex(kind=16), dimension(n_db, n_db) :: hmltn_ij
     complex(kind=16), dimension(n_db, n_db) :: ovlp_t_ij
     complex(kind=16), dimension(n_db, n_db) :: delta_2_h
     complex(kind=16), dimension(n_db, n_db) :: alhs 
     complex(kind=16), dimension(n_db) :: rhs 
     complex(kind=16), dimension(n_db, n_db):: ovlp_ij
     complex(kind=16) :: h_ij, o_db_ij, o_db_ij_1, o_ij, o_ij_1
     complex(kind=16) :: work(10000)
     real(kind=8) :: s_lgrng_n_db
     real(kind=16) :: t
     real(kind=8) :: h_ii

     integer, dimension(n_db):: IPIV
     integer :: i_db, j_db, i_df, j_df, mm_db, nn_db
     integer :: nsmall, ninp, info
     character(len=1) :: LULU
 
    ! Initialize the derivative array to zero
     dd_dt = (0.0d0, 0.0d0)

    ! Calculate the matrices hmltn_ij and ovlp_ij
    do i_db = 1, n_db
      do j_db = 1, n_db

        call h_ord(i, j, h)
        call overlap_ij(i_db, j_db, o_ij)

        hmltn_ij(i_db, j_db) = h_ij
        ovlp_ij(i_db, j_db) = o_ij
      end do
    end do

    ! Calculate the time derivative of the overlap matrix ovlp_t_ij
do mm_db = 1, n_db
    do nn_db = 1, n_db
        o_ij = (0.0d0, 0.0d0)
        do i_df = 1, n_df
          o_ij_1 = (1.0d0, 0.0d0)
          do j_df = 1, n_df
            if (j_df.ne.i_df) then
              o_ij_1 = o_ij_1 * (conjg(zz(1, j_df, mm_db)) * zz(1, j_df, nn_db) + &
                                 conjg(zz(2, j_df, mm_db)) * zz(2, j_df, nn_db))
            else
              o_ij_1 = o_ij_1 * (conjg(zz(1, j_df, mm_db)) * dzzdt(1, j_df, nn_db) + &
                                 conjg(zz(2, j_df, mm_db)) * dzzdt(2, j_df, nn_db))
            end if
          end do
          o_ij = o_ij + o_ij_1
        end do

        ovlp_t_ij(mm_db, nn_db) = (o_ij)
    end do
end do

    ! Calculate the delta_2_h matri
do i_db = 1, n_db
      do j_db = 1, n_db
        delta_2_h(i_db, j_db) = ovlp_ij(i_db, j_db) * s_lgrng(j_db) + &
                                hmltn_ij(i_db, j_db) - (0.0d0, 1.0d0) * ovlp_t_ij(i_db, j_db)
      end do
end do

    ! Calculate the right-hand side vector rhs
do i_db = 1, n_db
       rhs(i_db) = (0.0d0, 0.0d0)
       do j_db = 1, n_db
         rhs(i_db) = rhs(i_db) + delta_2_h(i_db, j_db) * d_big(j_db) * &
                     exp((0.0d0, 1.0d0) * (ss(j_db) - ss(i_db)))
      end do
end do

    ! Calculate the left-hand side matrix alhs
do i_db = 1, n_db
      do j_db = 1, n_db
        alhs(i_db, j_db) = ovlp_ij(i_db, j_db) * exp((0.0d0, 1.0d0) * (ss(j_db) - ss(i_db)))
      end do
end do

    ! Solve the system of linear equations for the derivatives
    nsmall = n_db
    ninp = n_db
    call zgesv(nsmall, 1, alhs, ninp, IPIV, rhs, ninp, info)
    
    ! Update the derivative array dd_dt
    dd_dt = -(0.0d0, 1.0d0) * rhs
    
  end subroutine drv_c

  
end Module Merfat_drivs
