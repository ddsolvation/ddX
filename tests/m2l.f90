!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/m2l.f90
!! Performance test for M2L
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-11

program test_ddx_m2l
use ddx_core
implicit none

integer, parameter :: ntests=100000, pmax=19, pmax2=20**2
real(dp) :: vscales((2*pmax+1)**2), v4pi2lp1(2*pmax+1), &
    & vscales_rel((2*pmax+1)**2), vcnk((4*pmax+1)*(2*pmax+1)), &
    & m2l_ztranslate_coef(pmax+1, pmax+1, pmax+1), &
    & m2l_ztranslate_adj_coef(pmax+1, pmax+1, pmax+1), &
    & src_c(3, ntests), dst_c(3, ntests), src_r(ntests), dst_r(ntests), &
    & start_time, finish_time
real(dp), allocatable :: src_m(:, :), dst_l(:, :), work(:)
integer :: iseed(4), nsrc_m, nsrc_c, i, idist
real(dp), external :: dnrm2

allocate(src_m(pmax2, ntests))
allocate(dst_l(pmax2, ntests))
allocate(work(6*pmax2+19*pmax+8))
! Compute special FMM constants
!call ylmscale(2*pmax, vscales, v4pi2lp1, vscales_rel)
!call fmm_constants(2*pmax, pmax, pmax, vcnk, m2l_ztranslate_coef, &
!    & m2l_ztranslate_adj_coef)
! Init random seed
iseed = (/0, 0, 0, 1/)
idist = 3
! Generate src_m randomly
nsrc_m = ntests * ((pmax+1)**2)
call dlarnv(idist, iseed, nsrc_m, src_m)
! Set sources and destinations
nsrc_c = ntests * 3
call dlarnv(idist, iseed, nsrc_c, src_c)
call dlarnv(idist, iseed, nsrc_c, dst_c)
! Set radii
do i = 1, ntests
    src_r(i) = dnrm2(3, src_c(:, i)-dst_c(:, i), 1) / 3d0
    dst_r(i) = src_r(i)
end do
! Warmup
i = 1
call cpu_time(start_time)
call fmm_m2l_rotation_work(src_c(:, i)-dst_c(:, i), src_r(i), dst_r(i), &
    & pmax, pmax, vscales, m2l_ztranslate_coef, one, src_m(:, i), zero, &
    & dst_l(:, i), work)
! Check
call cpu_time(start_time)
do i = 1, ntests
    call fmm_m2l_rotation_work(src_c(:, i)-dst_c(:, i), src_r(i), dst_r(i), &
        & pmax, pmax, vscales, m2l_ztranslate_coef, one, src_m(:, 1), zero, &
        & dst_l(:, 1), work)
end do
call cpu_time(finish_time)
print *, "Total   time:", finish_time-start_time, "seconds"
print *, "Average time:", (finish_time-start_time) / ntests, "seconds"

deallocate(src_m, dst_l, work)

end program
