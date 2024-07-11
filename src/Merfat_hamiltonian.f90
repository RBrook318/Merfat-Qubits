


 !_____________________________________________________________________________________________________________________________!
 !                                                                                                                             !
 !                   1- subroutine h_su2_lukin(zi, zj, h_lh) , independ of time                                                !
 !                                                                                                                             !                                                                                                             !
 ! The h_su2_lukin subroutine calculates the Hamiltonian matrix elements for a quantum system described by the SU(2) model,    !
 ! particularly in the context of the Lukin model (representing the energy contributions of the given quantum states under the SU(2) model.). This subroutine is crucial for simulating the dynamics of quantum systems   !
 ! under the influence of various parameters and interactions.                                                                 !
 !                                                                                                                             !
 !  - subroutine named h_su2_lukin used to calculate a quantum Hamiltonian for a quantum system with multiple modes            !
 !      based on the Lukin method.                                                                                             ! 
 !  - The delta and omega terms involve contributions from individual modes, while the coupling terms account for              !
 !      interactions between different modes.                                                                                  !
 !  - One of the well-known achievements of Mikhail Lukin is Lukin's involvement in developing methods for manipulating        !
 !     individual quantum bits (qubits) using systems based on atoms or artificial atoms (quantum dots).                       !
 !                                                                                                                             !                                                                                       
 !                                                                                                                             !
 !                 1- Intent and Variable Declarations:                                                                                           !
 !   1- Input:                                                                                                                 !
 ! zi(2, n_df), zj(2, n_df): Arrays representing quantum states for two different states, each with two components             !
 ! (representing spin-up and spin-down) and n_df degrees of freedom.                                                           !
 !                                                                                                                             !
 !   2- Output:                                                                                                                !
 ! h_lh: Complex variable representing the Hamiltonian matrix element calculated based on the given quantum states.            !  
 !                                                                                                                             !
 !   3- Parameter Declarations:
 ! e_tol_min: Minimum tolerance for small values (set to 10*10^-19)
 ! n_df: Number of degrees of freedom = 7.(number of modes)
 !
 !   4- Variable Initialization:
 ! delta_lh, omega_lh: Parameters defining the Hamiltonian structure, initialized to appropriate values.
 ! 
 !   5- others: 
 ! zic(2, n_df), zjc(2, n_df), a2a(n_df), ovlp_ij are complex arrays or variables used in intermediate calculations. !
 ! k_df, i_df, j_df are integer variables used as loop indices.                                                      ! 
 !                                                                                                                             !
 !                2- Variable Initialization:                                                                                                    !
 !        1- delta_lh and omega_lh are real variables initialized with specific values.                                        !
 !        2- zic and zjc are initialized as the complex conjugates of zi and zj, respectively.                                 !
 !        3- h_lh is initialized to zero.                                                                                      !
 !                                                                                                                             !
 !                 3 - Calculations:
 ! Delta and Omega Term Calculation Loop:
 ! Computes contributions to the Hamiltonian matrix element from the "delta" and "omega" terms for each degree of freedom.
 ! Delta and Omega Term Calculation Loop:                                                                                      !
 !        1- A loop (do loop) iterates over n_df modes to i_df calculates the contributions to the Hamiltonian from the delta  !
 !           term and the omega term for each mode.                                                                            !
 !        2- The overlap term ovlp_ij is calculated, excluding the current mode (i_df).                                        !
 !        3- Delta term (delta_lh) and Omega term (omega_lh) contributions are added to the Hamiltonian.                       !
 !        4- Delta term: Adds a contribution proportional to delta_lh(i_df) * zic(2, i_df) * zj(2, i_df) * ovlp_ij.            !
 !        5- Omega term: Adds a contribution proportional to omega_lh(i_df) * (zic(2, i_df) * zj(1, i_df) + zic(1, i_df) * zj(2, i_df)) * ovlp_ij.
 !                                                                                                                             !                                                 
 ! Coupling Term Calculation Loop:                                                                                             !
 !        1- Two nested loops iterate over n_df modes (i_df and j_df) calculate the contributions to the Hamiltonian from      !
 !            coupling terms between different modes.                                                                          !                                             
 !        2- For each pair (i_df, j_df) where i_df is not equal to j_df, it calculates the overlap ovlp_ij of all terms except !
 !           i_df and j_df.                                                                                                    !
 !        3- The coupling term (c_lh) is added to the Hamiltonian based on the overlap and coupling coefficients (a2a).        !
 !          (Adds a contribution proportional to c_lh(i_df, j_df) * ovlp_ij * a2a(i_df) * a2a(j_df) to the Hamiltonian.)       !                                                                                                      
 ! Summary:                                                                                                                    !
 !        1- Overall, the subroutine computes the Hamiltonian matrix elements based on the provided quantum states,            !
 !           considering both individual contributions from each degree of freedom and coupling interactions between different !
 !           degrees of freedom. These Hamiltonian matrix elements are essential for understanding the dynamics and behavior   !
 !           of the quantum system under study.                                                                                !
 !        2- It includes terms related to delta, omega, and coupling between modes.                                            !
 !                                                                                                                             !                                                                                                                                                                                                       
 !_____________________________________________________________________________________________________________________________!

Module Merfat_hamiltonian

  use Merfat_arr_prm
  use Merfat_input
  use Merfat_zqp_fullbasis

  implicit none

contains
subroutine h_su2_lukin(zi, zj, h_lh)

! Array Declarations
  complex(kind=16), dimension(2, n_df) :: zi           ! Arrays representing the wavefunctions (Quantum states) of the system 
  complex(kind=16), dimension(2, n_df) :: zj
  ! complex(kind=16), dimension(2, n_df) :: zic
  ! complex(kind=16), dimension(2, n_df) :: zjc 
  complex(kind=16) :: h_lh                             !  Complex variable representing the Lukin Hamiltonian.
  complex(kind=16), dimension(n_df) :: a2a             ! Array used to store intermediate results in the calculation.
  complex(kind=16), dimension(n_df) :: a2b
  ! complex(kind=16) :: ovlp_ij                          ! < z_i | z_j > Complex variable representing the overlap between different modes.
  ! real(kind=8), dimension(n_df, n_df) :: c_lh

! Parameter Declarations:  
  real(kind=8), parameter :: e_tol_min = 1.0d-10       ! A small tolerance value.
  real(kind=8) :: delta_lh                             ! Declare delta_lh as scalar
  real(kind=8) :: omega_lh                             ! Declare omega_lh as scalar
  integer :: i_df, j_df, k_df, r
  integer, PARAMETER :: n_df = 7 
  REAL(kind=8), PARAMETER  :: pinum = 3.141592653589793238462643383279d0 

  delta_lh = 0.0d0 
  delta_lh = delta_lh * 2 * pinum
  omega_lh = 2.0d0 * pinum * 2.0d0            ! Initialize to the desired value
  omega_lh = omega_lh / 2.0d0                 ! Halve the initialized value 

! Variable Initialization
  zic = conjg(zi)
  zjc = conjg(zj)
  h_lh = (0.0d0, 0.0d0)

!_______________________________________________________________

! Delta and Omega Term Calculation Loop:
do i_df = 1, n_df

    ! 1- Calculate the overlap of all modes except i_df
    ovlp_ij = (1.0d0, 0.0d0)

     do k_df = 1, n_df

          if (k_df /= i_df) then
          ovlp_ij = ovlp_ij * (zic(1, k_df) * zj(1, k_df) + zic(2, k_df) * zj(2, k_df))
          end if 
     end do

! 2- Delta term: This operation represents a part of the computation of the Lukin Hamiltonian   هذا بيساوي الصفر لان الديلتا صفر 
    h_lh = h_lh - delta_lh * zic(2, i_df) * zj(2, i_df) * ovlp_ij


! 3- Omega term: This operation represents a part of the computation of the Lukin Hamiltonian  هذا مصفوفة
    h_lh = h_lh + omega_lh * (zic(2, i_df) * zj(1, i_df) + zic(1, i_df) * zj(2, i_df)) * ovlp_ij

end do

     
! 3- Coupling Term Calculation Loop:
do i_df = 1, n_df
    a2a(i_df) = zic(2, i_df) * zj(2, i_df)

end do

do i_df = 1, n_df
    do j_df = 1, n_df
      if (j_df /= i_df) then

        ! Calculate the overlap of all terms except i_df and j_df
         ovlp_ij = (1.0d0, 0.0d0)
         do k_df = 1, n_df
             if (k_df /= i_df .and. k_df /= j_df) then
             ovlp_ij = ovlp_ij * (zic(1, k_df) * zj(1, k_df) + zic(2, k_df) * zj(2, k_df))
             end if
        end do
 
        ! Coupling term calculation
        h_lh = h_lh + c_lh(i_df, j_df) * ovlp_ij * a2a(i_df) * a2a(j_df)   ! الهاملتونين هنا عباره عن مصفوفة ١٢٨ *١٢٨ كلها اصفار عدا القطر بعضه القيم صفر والبعض الاخر قيمة تساوي تقريبا ٦ 
      end if

    end do
end do

end subroutine h_su2_lukin

 
 !_____________________________________________________________________________________________________________________________!
 !                                                                                                                             !
 !                   2- subroutine dh_su2_lukin(zi, zj, h_lh) , independ of time                                               !
 !                                                                                                                             !     
 ! The dh_su2_lukin subroutine calculates the derivative of the Lukin Hamiltonian with respect to the wavefunctions zi         !
 !                                                                                                                             !
 ! Here's an explanation of its functionality:                                                                                 !
 !                                                                                                                             !
 !                                      1- Intent and Variable Declarations:                                                   !
 !    1- Input:                                                                                                                !
 ! zi(2, n_df): Arrays representing the wavefunctions of the system.                                                           !
 !                                                                                                                             !
 !    2- Output:                                                                                                               !
 ! dh_dzi(2, n_df): Derivative of the Lukin Hamiltonian with respect to zi.                                                    !
 !                                                                                                                             !
 !    3- Variable Initialization:                                                                                              ! 
 ! dh_dzi and dh_dzic are initialized to zero.                                                                                 !
 !                                                                                                                             !
 !                                                                                                                             !
 !                                        2- Calculations:                                                                     !
 !    1- Delta Term Calculation:                                                                                               !
 ! Computes the derivative of the Lukin Hamiltonian with respect to zi for the delta term.                                     !
 ! Uses the formula dh_dzi(2, i_df) = dh_dzi(2, i_df) - delta_lh(i_df) * zic(2, i_df) for each mode i_df.                      !
 !                                                                                                                             !
 !    2- Omega Term Calculation:                                                                                               !
 ! Computes the derivative of the Lukin Hamiltonian with respect to zi for the omega term.                                     !
 ! Uses the formula dh_dzi(1, i_df) = dh_dzi(1, i_df) + omega_lh(i_df) * zic(2, i_df) and dh_dzi(2, i_df) = dh_dzi(2, i_df) +  ! 
 ! omega_lh(i_df) * zic(1, i_df) for each mode i_df.                                                                           !
 !                                                                                                                             !
 !    3- Coupling Term Calculation:                                                                                            !
 ! Computes the derivative of the Lukin Hamiltonian with respect to zi for the coupling term.                                  !
 ! Uses the formula dh_dzi(2, i_df) = dh_dzi(2, i_df) + c_lh(i_df, j_df) * zic(2, i_df) * zic(2, j_df) * zi(2, j_df)           !
 ! for each pair of modes (i_df, j_df) where i_df != j_df.                                                                     !
 !                                                                                                                             !
 ! Summary:                                                                                                                    !  
 ! The dh_su2_lukin subroutine is essential for understanding how small changes in the wavefunctions affect the Lukin Hamiltonian. 
 ! It calculates the derivative of the Hamiltonian with respect to the wavefunctions, providing insights into the              !
 ! system's dynamics and behavior under perturbations.                                                                         !
 !                                                                                                                             !
 !_____________________________________________________________________________________________________________________________!                                                                                                                           !     

subroutine dh_su2_lukin(zi, dh_dzi)
  ! This subroutine computes the derivative of the Lukin Hamiltonian with respect to the wavefunctions zi

  ! Variable Declarations
  complex(kind=16), dimension(2, n_df) :: zi                                                               ! Input wavefunctions
  complex(kind=16), dimension(2, n_df) :: dh_dzi                                                           ! Derivative of the Lukin Hamiltonian with respect to zi
  complex(kind=16), dimension(2, n_df) :: dh_dzic                                                          ! Derivative of the Lukin Hamiltonian with respect to zic
  ! real(kind=8), parameter :: e_tol_min = 1.0d-10                                                         ! Tolerance parameter

  ! ! Local Variables
  ! complex(kind=16), dimension(2, n_df) :: zic                                                             ! Conjugate of zi
  ! integer :: i_df                                                                                         ! Loop variable
  ! complex(kind=16) :: ovlp_ij                                                                             ! Overlap between wavefunctions
  complex(kind=16), dimension(n_df) :: a2b                                                                  ! Temporary variables for coupling term calculation
  complex(kind=16), dimension(n_df) :: a2a
  
  delta_lh = 0.0d0 
  delta_lh = delta_lh * 2 * pinum
  omega_lh = 2.0d0 * pinum * 2.0d0                                                                         ! Initialize to the desired value
  omega_lh = omega_lh / 2.0d0                                                                              ! Halve the initialized value

  ! 1- Initialize output arrays
      dh_dzi = (0.0d0,0.0d0)                                                                               ! Initialize the derivative array
      dh_dzic = (0.0d0,0.0d0) 
      zic = conjg(zi)                                                                                      ! Compute the conjugate of z_i

  ! 2- Delta Term Calculation: d[ h_lh = h_lh - delta_lh * zic(2, i_df) * zj(2, i_df) * ovlp_ij ]/ dzi
do i_df = 1, n_df                                                                                         ! Loop over all modes                                   
     dh_dzi(2, i_df) =  dh_dzi(2, i_df) -  delta_lh(i_df) * zic(2, i_df)                                  ! Compute the derivative term
     dh_dzic(2, i_df) = dh_dzic(2, i_df) - delta_lh(i_df) * zi(2, i_df)                                   ! Compute the derivative term
end do


  ! 3- Omega Term Calculation: d[  h_lh = h_lh + omega_lh * (zic(2, i_df) * zj(1, i_df) + zic(1, i_df) * zj(2, i_df)) * ovlp_ij ] / dzic
do i_df = 1, n_df                                                               
     dh_dzi (1, i_df) = dh_dzi (1, i_df) + omega_lh(i_df) * zic(2, i_df)                                   ! Compute derivative for omega term (1st component)
     dh_dzi (2, i_df) = dh_dzi (2, i_df) + omega_lh(i_df) * zic(1, i_df)                                   ! Compute derivative for omega term (2nd component)
     dh_dzic(1, i_df) = dh_dzic(1, i_df) + omega_lh(i_df) * zi (2, i_df)
     dh_dzic(2, i_df) = dh_dzic(2, i_df) + omega_lh(i_df) * zi (1, i_df) 
end do

    
  ! 4- Coupling Term Calculation:
do i_df = 1, n_df
     a2a(i_df) = zic(2, i_df) * zi(2, i_df)                                                                ! Temporary variable for coupling term
     a2b(i_df) = zic(2, i_df) * zi(2, i_df)                                                                ! Temporary variable for coupling term
end do

do i_df = 1, n_df
     do j_df = 1, n_df
        if (j_df /= i_df) then   ! يحسبها فقط اذا كانوا ما يتساوون 

           dh_dzi(2, i_df) = dh_dzi(2, i_df) + c_lh(i_df, j_df) * zic(2, i_df) * zic(2, j_df) * zi (2, j_df)   ! Compute derivative for coupling term
           dh_dzi(2, j_df) = dh_dzi(2, j_df) + c_lh(i_df, j_df) * zic(2, i_df) * zi (2, i_df) * zic(2, j_df)   ! Compute derivative for coupling term
          
           dh_dzic(2, i_df) = dh_dzic(2, i_df) +  c_lh(i_df, j_df) * zi(2, i_df) * zic(2, j_df) * zi(2, j_df)
           dh_dzic(2, j_df) = dh_dzic(2, j_df) +  c_lh(i_df, j_df) * zic(2, i_df) * zi(2, i_df) * zi(2, j_df) 

        end if
     end do
end do
     
end subroutine dh_su2_lukin

!  الهاملتونين 
!______________________________________________________________________________________________________________________________________!
 !                                                                                                                                     !
 !                          3- subroutine h_ord(i, j, h)                                                                               !
 !                                                                                                                                     !
 !  "calculating an ordered Hamiltonian: h"  The h_ord subroutine appears to calculate the ordered                                     !
 !  Hamiltonian matrix element h based on the given indices i and j and the quantum states stored in the                               !
 !   arrays zi and zj.                                                                                                                 !
 !                                                                                                                                     !
 !      1- Input:                                                                                                                      !
 ! - i, j: Indices representing the quantum states to use for calculating the Hamiltonian matrix element.                              !
 ! - zqp(l, i_df, i): Quantum states stored in a 3D array zqp, where l represents the spin component (1 for spin-up, 2 for spin-down)  ! 
 !   , i_df represents the degree of freedom, and i (or j) represents the index of the quantum state.                                  !
 !                                                                                                                                     !  
 !      2- Output:                                                                                                                     ! 
 ! - h: Complex variable representing the ordered Hamiltonian matrix element based on the provided quantum states.                     !
 !                                                                                                                                     ! 
 !     3- Local Variables:                                                                                                             !
 ! - i_df, l: Loop indices for iterating over the degrees of freedom and spin components.                                              !
 !                                                                                                                                     !
 !     4- Populating zi and zj:                                                                                                        !
 ! The subroutine populates the arrays zi and zj with the quantum states corresponding to the indices i and j. It iterates over each   !
 ! degree of freedom (i_df) and spin component (l), retrieving the quantum states from the zqp array.                                  !
 !                                                                                                                                     !
 !    5- Calculation of Ordered Hamiltonian:                                                                                           !
 ! The actual calculation of the ordered Hamiltonian h seems to be missing from this subroutine. Instead, it prepares the necessary    !
 ! quantum states for the calculation, likely to be used in subsequent subroutines or functions.                                       !
 !                                                                                                                                     !
 !    6- Overall,                                                                                                                      !
 ! The h_ord subroutine sets up the quantum states required for calculating the ordered Hamiltonian matrix element. It prepares the    !
 ! input quantum states zi and zj based on the provided indices i and j.                                                               !
 !_____________________________________________________________________________________________________________________________________!

subroutine h_ord(i, j, h)

! Declarations
    integer :: i, j                                      ! Input indices i and j
    complex(kind=16) :: h                                ! Output variable representing the ordered Hamiltonian h
    complex(kind=16), dimension(2, n_df) :: zi, zj       ! 2D arrays with dimensions (2, n_df) to store quantum states corresponding to indices i and j, respectively.
    ! complex(kind=16), dimension(n_df) :: ovlp          ! 1D array for overlap calculations

! Local Variables
    ! integer :: i_df, l                                 ! Loop indices
    
! zqp بالبيانات من مصفوفة zi وzj يتم استخدام حلقتين متداخلتين لملء المصفوفتين 
! والتي من المفترض أن يتم تعريفها في مكان آخر في الكود. تعمل الحلقة الخارجية على 
! Populate zi and zj
do i_df = 1, n_df
        do l = 1, 2
            zi(l, i_df) = zqp(l, i_df, i)
            zj(l, i_df) = zqp(l, i_df, j)
        end do
 end do

!  بعد الانتهاء من ملء المصفوفتين يتم استدعاء الروتين لتخزين القيمة الهاملتونيه المحسوبة بواسطة الروتين الفرعي
! Call the subroutine to calculate the Hamiltonian
    call h_su2_lukin(zi, zj, h_lh)
    
end subroutine h_ord
end Module Merfat_hamiltonian
 
 


