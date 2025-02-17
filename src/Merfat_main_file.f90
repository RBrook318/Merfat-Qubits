



!***************************************************************************************************!
!                                                                                                   !
! Date : 014/05/23                                                                                  !
!                                                                                                   !
!  This Program performs simulations using the Couple coherent state for two_level system in |z>    !
!  by full variation method                                                                         !
!  The program has the following stucture in the Main module :                                      !
!  This code was adopted to simulate the dynamics of coupled qubits                                 !
!  This Program performs simulations the dynamics of coupled qubits using the Couple coherent state ! 
!  & with zero-contact for two_level system in |z>                                                  !
!                                                                                                   !
! The program has the following stucture in the Main module :                                       !
!      1) Read in run conditions. These are generated by the program run shell script               !
!      2) If run conditions allow, a basis set is generated and checked to ensure                   !
!            it is normalised                                                                       !
!      3) If run conditions allow, the basis set is then propagated in time                         !
!      4) Results of propagation are then output                                                    !
!                                                                                                   ! 
!  If propagation is allowed but basis set generation is disallowed, a valid basis                  !
!  set file, or set of basis set files, must be present in the run folder. This can                 !
!  be generatedby running the program with only basis set generation allowed.                       !
!                                                                                                   !
!  The program is compiled with the -openmp flag (for ifort - other compilers use.                  !
!                                                                                                   !
!***************************************************************************************************!
program main_file

! Use the module
  use Merfat_arr_prm 
  ! use Merfat_one_zero_state
  use Merfat_zqp_fullbasis
  use Merfat_zqp_0c_state
  use Merfat_zqp_1c_state
  use Merfat_overlap
  use Merfat_overlap_zinit_zi
  use Merfat_step_t
  use Merfat_drivs
  use Merfat_input
  use Merfat_hamiltonian
  use Merfat_norm
  use Merfat_ACF
  use Merfat_CCF
  use Merfat_calculation
  use Merfat_init_conditions

  implicit none

  ! Declare variables
  character(len=100) :: input_filename
  integer :: io_stat, iostat
  
  ! Setting up Constants
  real(kind=8), parameter :: sqrt2 = sqrt(2.0d0)
   i = 1
   j = 1
  !  pinum = 4.0d0*datan(1.0d0)
	!  sqrt2 = dsqrt(2.0d0)
   i_prop=0
   i_step=0


   ! Initialize random seed
  call random_seed()




  ! Call subroutines based on dependencies
  ! Reads input data or initializes variables based on some input file or parameters.
!__________________________________1- basis set : input_lukin ________________________________
  call input_lukin                 
  write(*,*) '     1- --------------- input_lukin..... done -------------------------------------'



! Sets up the zqp_full_basis for the calculations.
!__________________________________2- setup_zinit_from_fil  ____________________________ 
  call setup_zqp_full_basis()
  write(*,*) '     2- --------------- setup_zqp_full_basis..... done ----------------------------'




!__________________________________3- sqp_0_contact  ____________________________________ 
  call zqp_0_contact()
  write(*,*) '     3- --------------- zqp_0_contact..... done  ----------------------------------'
 



!__________________________________4- qp_1_contact  ____________________________________ 
  call zqp_1_contact()
  write(*,*) '     4- --------------- zqp_1_contact..... done -----------------------------------'




!__________________________________5- input_zinit  ______________________________________ 
  call input_zinit(n_df, zinit)
  write(*,*) '     5- --------------- input_zinit..... done -------------------------------------'




!__________________________________6- setup_zinit_from_file  ____________________________
  call setup_zinit_from_file(n_df, zinit)
  write(*,*) '     6- --------------- setup_zinit_from_file..... done ----------------------------'





!__________________________________7- overlap (ovlp_ij) ________________________________
open(unit=1003, file='output/ovlp_ij.out', status='unknown', action='write', iostat=iostat)
if (iostat /= 0) then
    print*, "Error opening file ovlp_ij"
end if
      
   do i_db = 1, n_db
      do j_db = 1, n_db
          call overlap_ij(i_db, j_db, ovlp_ij)
          write (1003, *) i_db, j_db, ovlp_ij
      end do
   end do

close(1003)
write(*,*) '     7- --------------- overlap (ovlp_ij)..... done ------------------------------------'




!__________________________________8- h_su2_lukin  ______________________________________  
open(unit=1006, file='output/h_su2_lukin.out', status='unknown', action='write', iostat=iostat)
if (iostat /= 0) then
    print*, "Error opening h_su2_lukin.out"  
end if
 
do i = 1, n_db 
 write (1006, *) '----------'
 write (1006, *) 'i=', i 
  
      do j = 1, n_db 
      write (1006, *) 'j =', j                    ! j= 1, 128  
      write (1006, *) '----------'

	            do i_df = 1, n_df  
              write (1006, *) 'i_df =', i_df

                  do r = 1, 2
                  write (1006, *) 'r =', r

	        zi(r,i_df) = zqp(r,i_df,i)  ! المفروض ١٤ زوج يكون بالتبادل ١ ثم صفر وهكذا بحيث ياخذ ١٢٨ احتمال 
          zic = conjg(zi)

	        zj(r,i_df) = zqp(r,i_df,j)  ! المفروض ١٤ زوج يكون بالتبادل ١ ثم صفر وهكذا بحيث ياخذ ١٢٨ احتمال 
          zjc = conjg(zj)

                  end do
              end do
      call h_su2_lukin(zi, zj, h_lh)

        ! Update delta_lh and omega_lh here
        delta_lh = delta_lh * 2 * pinum
        omega_lh = 2.0d0 * pinum * 2.0d0
        omega_lh = omega_lh / 2.0d0
              write (1006, *)
              write (1006, *) 'delta_lh =', delta_lh                     ! delta_lh = 0
              write (1006, *) 'omega_lh =', omega_lh                     ! omega_lh = 6.2831853071795862
              write (1006, *)
              write (1006, *) 'zi =', zi
              write (1006, *) 'zic =', zic
              write (1006, *) 'zj =', zj
              write (1006, *) 'zjc =', zjc
              write (1006, *)
              write (1006, *) 'a2a(i_df) =', a2a(i_df)
              write (1006, *) 'a2a(j_df) =', a2a(j_df)
              write (1006, *)
              write (1006, *) 'Hamiltonian (h_lh) =', h_lh              !  مصفوفه ١٢٨*١٢٨ 
              write (1006, *)
      end do
end do 
! close file: 
close(1006)
write(*,*) '     8- --------------- h_su2_lukin .......... done ------------------------------------'




!__________________________________9- dh_su2_lukin ______________________________________ 
write(*,*) 
write(*,*) 
write(*,*) '__________________________ dh_su2_lukin  ____________________________________'
write(*,*) 
write(*,*) 
! open(unit=1009, file='./dh_su2_lukin.out', status='unknown', action='write', iostat=iostat)
! if (iostat /= 0) then
!     print*, "Error opening h_su2_lukin.out"  
! end if



!  ! Write the output to the file
!   write(1009, *) 'Derivative of Lukin Hamiltonian with respect to zi:'
!   write(1009, *) 'dh_dzi:'
!   do i_df = 1, n_df
!     write(1009, *) 'i_df =', i_df, ': ', dh_dzi(:, i_df)
!     call dh_su2_lukin(zi, dh_dzi)
!   end do

!   ! Write the output to the file
!   write(1009, *) 'Derivative of Lukin Hamiltonian with respect to zic:'
!   write(1009, *) 'dh_dzic:'
!   do i_df = 1, n_df
!     write(1009, *) 'i_df =', i_df, ': ', dh_dzic(:, i_df)
!     call dh_su2_lukin(zi, dh_dzi)
!   end do

! close(1009)

!  write(*,*) '     9- --------------- dh_su2_lukin ..... done -----------------------------------'





!____________________________________10- h_ord(i, j, h) _________________________________ 
open(unit=1011, file='output/h_ord.out', status='unknown', action='write', iostat=iostat)
if (iostat /= 0) then
    print*, "Error opening h_ord.out"  
end if

do i = 1, n_db  
write (1011, *) 'i=', i                                  ! i= 1, 128 

      do j = 1, n_db 
      write (1011, *) 'j =', j                           ! j= 1, 128  
 
	        do i_df = 1, n_df
          write (1011, *) 'i_df =', i_df                 ! i_df = 7  

              do l = 1, 2
              write (1011, *) 'l =', l                   ! l = 2

               zi(l, i_df) = zqp(l, i_df, i)
               zj(l, i_df) = zqp(l, i_df, j)

              end do 
          end do 

      call h_ord(i, j, h)

              write (1011, *) 'zi = ' , zi
              write (1011, *) 'zj = ' , zj
              write (1011, *) 
              write (1011, *) '________________________' 
    
      end do 
end do   

close(1011)
write(*,*) '    10- --------------- h_ord(i, j, h) ....... done ------------------------------------'


!    a_big = d_big * exp((0.0d0,1.0d0) * s_s)    ! we do not really need this line. a_big is set in the smodule Merfat_init_conditions
call setup_d_big_orthigonal
! __________________________________11- step_c(i_step) ____________________________________  from here some thing wrong 
open(unit = 1020, file='output/step_c.out', status='unknown', action='write', iostat=iostat)
if (iostat /= 0) then
    print*, "Error opening drv_c.out"  
end if

        call norm(anorm)
        write(*,*) 'anorm =' , anorm

        call acf(acf_t)  
        write(*,*) 'acf_t =', acf_t

        call ccf(ccf_t)
        write(*,*) 'ccf_t =', ccf_t
       
!         ! call dvd(dvd_t, dvd_tanan)
!         call output(dvd_t, anorm, acf_t,ccf_t, h , i_step)

!  ! Loop over steps
!     do i_step = 1, ntmax 

!         ! Call the subroutines
!           call step_c(i_step)
!           call output(dvd_t, anorm, acf_t,ccf_t, h , i_step)

!         ! Write outputs to file                     ! اتاكد منه كامل 
!           write(1020,*) ' Step = ', i_step
!           write(1020,*) '___________________________________'
!           write(1020,*) ' t = ', t
!           write(1020,*) '___________________________________' 
!           write(1020,*) ' zqp = ', zqp
!           write(1020,*) '___________________________________' 
!           write(1020,*) ' s_s = ', s_s
!           write(1020,*) '___________________________________' 
!           write(1020,*) ' d_big = ', d_big
!           write(1020,*) '___________________________________'
!     end do

!           call  acf(acf_t)  

close(1020)
write(*,*) '--------------------------------------------------------------------------'
 stop
! write(*,*)
! write(*,*) '     11- --------------- step_c(i_step) ..... done --------------------------------'





! !___________________________12- drv_t(t, zz, ss, dzz, dss, n_flag) ____________________ 
! write(*,*) 
! write(*,*) 
! write(*,*) '_______________ drv_t(t, zz, ss, dzz, dss, n_flag)   ______________________'
! write(*,*) 
! write(*,*) 

! ! open(unit=1012, file='./drv_t.out', status='unknown', action='write', iostat=iostat)
! ! if (iostat /= 0) then
! !     print*, "Error opening drv_t.out"  
! ! end if

! !  do i_df = 1, n_df 
! !     do i_12 = 1, 2
! !         do i_db = 1, n_db

! !             zi(i_12, i_df) = zz(i_12, i_df, i_db) 
! !             ! write(1012,*) 
! !             ! write(1012,*) 'zi(i_12, i_df = ', zi(i_12, i_df)
! !         end do 
! !     end do
    
!     ! Call the subroutine drv_t
!     call drv_t(t, zz, ss, dzz, dss, n_flag) 
    

!             !  write(1012, *) ' drv_t calculatiom:'
!             !  write(1012, *)
!             !  write(1012, *) 't =' , t
!             !  write(1012, *)
!             !  write(1012, *) 'zz = ' , zz
!             !  write(1012, *)
!             !  write(1012, *) 'ss = ' , ss
!             !  write(1012, *) 
!             !  write(1012, *) "dzz = ", dzz
!             !  write(1012,*)
!             !  write(1012, *) "dss = ", dss
!             !  write(1012,*)
!             !  write(1012,*) '_________________________________________________________' 
!             !  write(1012,*) 
           
! ! end do

! ! close(1012)
 write(*,*) '--------------------------------------------------------------------------'
write(*,*)
write(*,*) '     12- drv_t(t, zz, ss, dzz, dss, n_flag) ..... done'
write(*,*)
write(*,*) '--------------------------------------------------------------------------'




!___________________________13- drv_c(t, zz, ss, dzzdt, dd_dt, n_flag)_________________ 
write(*,*) 
write(*,*) 
write(*,*) '_______________ drv_c(t, zz, ss, dzzdt, dd_dt, n_flag)  ____________________'
write(*,*) 
write(*,*) 
! open(unit=1013, file='./drv_c.out', status='unknown', action='write', iostat=iostat)
! if (iostat /= 0) then
!     print*, "Error opening drv_c.out"  
! end if

call drv_c(t, zz, ss, dzzdt, dd_dt, n_flag)
            !  write(1013,*) 
            !  write(1013,*) 
            !  write(1013,*) 
            !  write(1013,*) 
            !  write(1013,*) 
            !  write(1013,*) 
            !  write(1013,*) 
             
! close(1013)
write(*,*) '--------------------------------------------------------------------------'
write(*,*)
write(*,*) '     13- drv_c(t, zz, ss, dzzdt, dd_dt, n_flag) ..... done'
write(*,*)
write(*,*) '--------------------------------------------------------------------------'


end program main_file


      



















 











  