  
module Merfat_overlap_zinit_zi

  use Merfat_arr_prm
  
  implicit none

contains


! _____________________ 1 ________________________________________

subroutine overlap_zz_ij(zz_qp_i, zz_qp_j, o_ij)                     ! I need it in Module Merfat_init_conditions.f90
    
    ! Declarations
    complex(kind=16), dimension(2, n_df) :: zz_qp_i
    complex(kind=16), dimension(2, n_df) :: zz_qp_j
    complex(kind=16) :: o_ij
    integer :: i_df, i_12

    ! Initialize overlap to (1, 0)
    o_ij = (1.0d0, 0.0d0)

    ! Calculate the overlap
    do i_df = 1, n_df
        do i_12 = 1, 2
            o_ij = o_ij * (conjg(zz_qp_i(i_12, i_df)) * zz_qp_j(i_12, i_df))
        end do
    end do
end subroutine overlap_zz_ij

! _____________________ 2  ________________________________________

subroutine overlap_zi_zinit(zz_qp_i, zinit, o_ij)                  ! I need it in Module Merfat_init_conditions.f90

    complex(kind=16), dimension(2, n_df) :: zz_qp_i
    complex(kind=16), dimension(2, n_df) :: zinit
    complex(kind=16) :: o_ij
    complex(kind=16) :: ovlp_1d
    real(kind=8), parameter :: e_tol = 1.0d-10
    integer :: i_df

    o_ij = (1.0d0, 0.0d0)

    do i_df = 1, n_df
        ovlp_1d = conjg(zz_qp_i(1, i_df)) * zinit(1, i_df) + &
                  conjg(zz_qp_i(2, i_df)) * zinit(2, i_df)
        o_ij = o_ij * ovlp_1d
    end do

    if (abs(o_ij) < e_tol) then
        o_ij = (0.0d0, 0.0d0)
    end if

end subroutine overlap_zi_zinit

! _____________________ 3  ________________________________________


end module Merfat_overlap_zinit_zi
