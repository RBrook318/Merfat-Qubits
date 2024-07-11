

module Merfat_overlap
  
  ! Include files
  use Merfat_arr_prm
  use Merfat_zqp_fullbasis
  
  implicit none
contains

!_____________________________________ 1 _______________________________________

subroutine overlap_ij(i_db, j_db, ovlp_ij)

    ! complex(kind=8), dimension(2, n_df, n_db) :: zqp
    complex(kind=16), dimension(2, n_df) :: z_i
    complex(kind=16), dimension(2, n_df) :: z_j
    complex(kind=16) :: ovlp_ij
    integer :: i_df, iostat
    integer:: i, j, i_db, j_db


 do i_df = 1, n_df
    ! Construct z_i
        z_i(1, i_df) = zqp(1, i_df, i_db)
        z_i(2, i_df) = zqp(2, i_df, i_db)
                
    ! Construct z_j
        z_j(1, i_df) = zqp(1, i_df, j_db)
        z_j(2, i_df) = zqp(2, i_df, j_db)

        !  write(*,*) 'i_df =', i_df
        !  write(*,*) 'i_db =' ,i_db
        !  write(*,*) 'j_db =', j_db
        !  write(*,*)
        !  write(*,*) 'z_i(1, i_df) =', z_i(1, i_df) 
        !  write(*,*) 'z_i(2, i_df) =', z_i(2, i_df) 
        !  write(*,*) 'z_j(1, i_df) =', z_j(1, i_df)
        !  write(*,*) 'z_i(2, i_df) =', z_j(2, i_df)  
        !  write(*,*) 
 end do

    call ovlap_ij_z(z_i, z_j, ovlp_ij)
        !  write(*,*) 'i_db =' ,i_db
        !  write(*,*) 'j_db =', j_db
        !  write(*,*) 'ovlp_ij =' , ovlp_ij
        !  write(*,*)
        !  write(*,*)
        !  write(*,*)
end subroutine overlap_ij

!_____________________________________ 2 _______________________________________

subroutine ovlap_ij_z(z_i, z_j, o_z_ij)                      ! I need it in Module:   Merfat ACF - Merfat_norm - Merfat_CCF - Merfat_init_conditions
 ! verlap Computation Subroutine

    complex(kind=16), dimension(2, n_df) :: z_i
    complex(kind=16), dimension(2, n_df) :: z_j
    complex(kind=16) :: o_z_ij
    ! complex(kind=16) :: ovlp_1d
    real(kind=8), parameter :: e_tol = 1.0d-10              ! 1×10**−10  very small number - e_tol: error tolerance. it represents the maximum allowable erro
    ! integer :: i_df
 
    o_z_ij = (1.0d0, 0.0d0)


    do i_df = 1, n_df

        ovlp_1d = conjg(z_i(1, i_df)) * z_j(1, i_df) + &
                  conjg(z_i(2, i_df)) * z_j(2, i_df)
        o_z_ij = o_z_ij * ovlp_1d

    end do


    if (abs(o_z_ij) < e_tol) then
          o_z_ij = (0.0d0, 0.0d0)

    end if


end subroutine ovlap_ij_z

!_____________________________________ 2 _______________________________________

subroutine overlap_ij_db(i_db, j_db, o_db_ij)                   ! I need it in Module Merfat_driv.f90 
    
    complex(kind=16) :: o_db_ij
    ! complex(kind=16) :: o_db_ij_1
    ! complex(kind=16), dimension(2, n_df, n_db) :: zz, dzzdt
    ! integer :: i_df, j_df
    integer :: i_db, j_db
  
    ! ! Initialize o_ij
    ! o_db_ij = (0.0d0, 0.0d0)
    ! o_db_ij_1 = (1.0d0, 0.0d0)

    !        do j_df = 1, n_df
    !                o_db_ij_1 =o_db_ij_1 * (conjg(zz(1, j_df, i_db)) * zz(1, j_df, j_db) + &
    !                                        conjg(zz(2, j_df, i_db)) * zz(2, j_df, j_db))
    !        end do
                
    ! o_db_ij = o_db_ij_1
 

end subroutine overlap_ij_db

end module Merfat_overlap

















!======================================================================================================================
! 1- Subroutine overlap_ij_z
!======================================================================================================================
! This subroutine calculates the overlap between two sets of quantum parameters (z_i and z_j) for a specific dynamic basis.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - z_i: A complex array of dimensions (2, n_df) representing the quantum parameters for the first set.
! - z_j: A complex array of dimensions (2, n_df) representing the quantum parameters for the second set.

!======================================================================================================================
! Output Variables:
!======================================================================================================================
! - ovlp_ij: A complex variable representing the calculated overlap between the two sets of quantum parameters.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - ovlp_1d: A complex variable representing the overlap contribution from each degree of freedom.
! - e_tol: A real parameter representing the tolerance for considering overlap values as zero.
! - i_df: An integer representing the index of the degree of freedom.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to compute the overlap between two sets of quantum parameters (z_i and z_j) for a
! specific dynamic basis. The overlap represents the similarity or correlation between the quantum states described by
! the two sets of parameters.

!======================================================================================================================
! Output:
!======================================================================================================================
! - The calculated overlap (ovlp_ij) between the two sets of quantum parameters.
! - If the calculated overlap is below a specified tolerance (e_tol), it is considered zero to avoid numerical errors.

!======================================================================================================================




!======================================================================================================================
! 2- Subroutine ovlap_ij_z(z_i, z_j, o_z_ij) 
!======================================================================================================================
! This subroutine calculates the overlap between two sets of quantum parameters (z_i and z_j) for a specific dynamic basis.

!======================================================================================================================
! Input Variables:
!======================================================================================================================
! - z_i: A complex array of dimensions (2, n_df) representing the quantum parameters for the first set.
! - z_j: A complex array of dimensions (2, n_df) representing the quantum parameters for the second set.

!======================================================================================================================
! Output Variables:
!======================================================================================================================
! - o_z_ij: A complex variable representing the calculated overlap between the two sets of quantum parameters.

!======================================================================================================================
! Local Variables:
!======================================================================================================================
! - ovlp_1d: A complex variable representing the overlap contribution from each degree of freedom.
! - e_tol: A real parameter representing the tolerance for considering overlap values as zero.
! - i_df: An integer representing the index of the degree of freedom.

!======================================================================================================================
! Main Goal:
!======================================================================================================================
! The main goal of this subroutine is to compute the overlap between two sets of quantum parameters (z_i and z_j) for a
! specific dynamic basis. The overlap represents the similarity or correlation between the quantum states described by
! the two sets of parameters.

!======================================================================================================================
! Output:
!======================================================================================================================
! - The calculated overlap (o_z_ij) between the two sets of quantum parameters.
! - If the calculated overlap is below a specified tolerance (e_tol), it is considered zero to avoid numerical errors.

!======================================================================================================================

            


