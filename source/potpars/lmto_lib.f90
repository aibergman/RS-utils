module lmto_lib
   implicit none
   
contains

subroutine defmom(konfig, lafm, lfree, lmx, nclass, nl, nsp, pnu, qnu, z)
   !- Returns default values for moments and pnu's for given nuclear charge
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   konfig:core configuration
   !i   lafm  :switch for antiferromagnetic structures
   !i   lfree :T: make spheres self-consistent for free atom
   !i   lmx   :lmx(j) =maximum l for atom j
   !i   nclass:number of classes, atoms in same class are symmetry-related
   !i   nl    :number of l's
   !i   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
   !i   z     :nuclear charge
   !o Outputs:
   !o   pnu   :principal quantum number
   !o   qnu   :moments
   !r Remarks:
   ! ----------------------------------------------------------------------
   implicit none
   ! Passed variables:
   integer konfig(0:3, *), lmx(*), nclass, nl, nsp
   real(kind=8) pnu(0:nl - 1, nsp, *), qnu(0:2, 0:nl - 1, nsp, *), z(*)
   logical lafm, lfree
   ! Local variables:
   integer ic, isp, iz, k0(0:3), konfdf(0:3), l, lmxc, q0(0:3)
   real(kind=8) fac, pnuat(0:3, 0:100), qnuat(0:3, 0:100), qnuc, &
      &                 dq1, dq0, tiny
   character*144 messg
   parameter(tiny=1.d-6)
   ! External calls:
   ! external atmoms, defcnf, fremom
   ! Intrinsic functions:
   intrinsic dble, idnint
   ! Data statements:
   data((pnuat(l, iz), l=0, 3), iz=0, 9)/                               &
      &             1.50D+00, 2.50D+00, 3.50D+00, 4.50D+00,           &
      &             1.76D+00, 2.30D+00, 3.16D+00, 4.11D+00,           &
      &             1.91D+00, 2.16D+00, 3.12D+00, 4.09D+00,           &
      &             2.63D+00, 2.46D+00, 3.18D+00, 4.11D+00,           &
      &             2.80D+00, 2.52D+00, 3.21D+00, 4.13D+00,           &
      &             2.77D+00, 2.64D+00, 3.21D+00, 4.12D+00,           &
      &             2.89D+00, 2.71D+00, 3.18D+00, 4.12D+00,           &
      &             2.92D+00, 2.79D+00, 3.16D+00, 4.11D+00,           &
      &             3.10D+00, 2.86D+00, 3.16D+00, 4.11D+00,           &
      &             3.16D+00, 2.89D+00, 3.15D+00, 4.10D+00/
   data((pnuat(l, iz), l=0, 3), iz=10, 19)/                             &
      &             3.17D+00, 2.92D+00, 3.12D+00, 4.09D+00,           &
      &             3.70D+00, 3.38D+00, 3.19D+00, 4.12D+00,           &
      &             3.76D+00, 3.50D+00, 3.22D+00, 4.12D+00,           &
      &             3.81D+00, 3.59D+00, 3.24D+00, 4.13D+00,           &
      &             3.88D+00, 3.68D+00, 3.23D+00, 4.13D+00,           &
      &             3.92D+00, 3.78D+00, 3.20D+00, 4.12D+00,           &
      &             3.94D+00, 3.84D+00, 3.18D+00, 4.11D+00,           &
      &             4.13D+00, 3.89D+00, 3.17D+00, 4.11D+00,           &
      &             4.15D+00, 3.92D+00, 3.12D+00, 4.09D+00,           &
      &             4.70D+00, 4.36D+00, 3.24D+00, 4.12D+00/
   data((pnuat(l, iz), l=0, 3), iz=20, 29)/                             &
      &             4.72D+00, 4.40D+00, 3.38D+00, 4.13D+00,           &
      &             4.69D+00, 4.40D+00, 3.53D+00, 4.13D+00,           &
      &             4.66D+00, 4.40D+00, 3.62D+00, 4.13D+00,           &
      &             4.65D+00, 4.40D+00, 3.68D+00, 4.13D+00,           &
      &             4.65D+00, 4.40D+00, 3.73D+00, 4.13D+00,           &
      &             4.65D+00, 4.41D+00, 3.77D+00, 4.13D+00,           &
      &             4.66D+00, 4.41D+00, 3.80D+00, 4.13D+00,           &
      &             4.66D+00, 4.41D+00, 3.83D+00, 4.13D+00,           &
      &             4.67D+00, 4.41D+00, 3.86D+00, 4.12D+00,           &
      &             4.69D+00, 4.42D+00, 3.88D+00, 4.12D+00/
   data((pnuat(l, iz), l=0, 3), iz=30, 39)/                             &
      &             4.77D+00, 4.50D+00, 3.92D+00, 4.11D+00,           &
      &             4.87D+00, 4.59D+00, 4.18D+00, 4.12D+00,           &
      &             4.91D+00, 4.70D+00, 4.19D+00, 4.12D+00,           &
      &             4.93D+00, 4.78D+00, 4.18D+00, 4.12D+00,           &
      &             4.95D+00, 4.84D+00, 4.17D+00, 4.11D+00,           &
      &             5.14D+00, 4.89D+00, 4.17D+00, 4.11D+00,           &
      &             5.16D+00, 4.92D+00, 4.12D+00, 4.09D+00,           &
      &             5.71D+00, 5.33D+00, 4.26D+00, 4.12D+00,           &
      &             5.73D+00, 5.36D+00, 4.40D+00, 4.13D+00,           &
      &             5.70D+00, 5.35D+00, 4.55D+00, 4.14D+00/
   data((pnuat(l, iz), l=0, 3), iz=40, 49)/                             &
      &             5.67D+00, 5.35D+00, 4.64D+00, 4.14D+00,           &
      &             5.66D+00, 5.35D+00, 4.70D+00, 4.14D+00,           &
      &             5.65D+00, 5.35D+00, 4.75D+00, 4.14D+00,           &
      &             5.65D+00, 5.35D+00, 4.79D+00, 4.14D+00,           &
      &             5.65D+00, 5.36D+00, 4.82D+00, 4.14D+00,           &
      &             5.65D+00, 5.36D+00, 4.85D+00, 4.13D+00,           &
      &             5.65D+00, 5.35D+00, 4.88D+00, 4.13D+00,           &
      &             5.69D+00, 5.37D+00, 4.90D+00, 4.12D+00,           &
      &             5.80D+00, 5.47D+00, 4.93D+00, 4.11D+00,           &
      &             5.88D+00, 5.58D+00, 5.17D+00, 4.12D+00/
   data((pnuat(l, iz), l=0, 3), iz=50, 59)/                             &
      &             5.91D+00, 5.70D+00, 5.18D+00, 4.12D+00,           &
      &             5.93D+00, 5.78D+00, 5.19D+00, 4.12D+00,           &
      &             5.94D+00, 5.84D+00, 5.18D+00, 4.12D+00,           &
      &             6.12D+00, 5.89D+00, 5.18D+00, 4.11D+00,           &
      &             6.15D+00, 5.92D+00, 5.13D+00, 4.09D+00,           &
      &             6.71D+00, 6.30D+00, 5.31D+00, 4.13D+00,           &
      &             6.69D+00, 6.27D+00, 5.52D+00, 4.15D+00,           &
      &             6.61D+00, 6.20D+00, 5.65D+00, 4.19D+00,           &
      &             6.56D+00, 6.17D+00, 5.66D+00, 4.28D+00,           &
      &             6.57D+00, 6.16D+00, 5.66D+00, 4.38D+00/
   data((pnuat(l, iz), l=0, 3), iz=60, 69)/                             &
      &             6.57D+00, 6.16D+00, 5.65D+00, 4.48D+00,           &
      &             6.57D+00, 6.16D+00, 5.64D+00, 4.55D+00,           &
      &             6.58D+00, 6.17D+00, 5.63D+00, 4.61D+00,           &
      &             6.61D+00, 6.18D+00, 5.62D+00, 4.65D+00,           &
      &             6.63D+00, 6.20D+00, 5.60D+00, 4.69D+00,           &
      &             6.65D+00, 6.22D+00, 5.59D+00, 4.73D+00,           &
      &             6.67D+00, 6.24D+00, 5.56D+00, 4.76D+00,           &
      &             6.69D+00, 6.26D+00, 5.54D+00, 4.78D+00,           &
      &             6.70D+00, 6.28D+00, 5.52D+00, 4.80D+00,           &
      &             6.73D+00, 6.32D+00, 5.47D+00, 4.82D+00/
   data((pnuat(l, iz), l=0, 3), iz=70, 79)/                             &
      &             6.74D+00, 6.34D+00, 5.44D+00, 4.85D+00,           &
      &             6.73D+00, 6.35D+00, 5.53D+00, 4.93D+00,           &
      &             6.72D+00, 6.36D+00, 5.61D+00, 5.14D+00,           &
      &             6.71D+00, 6.37D+00, 5.68D+00, 5.14D+00,           &
      &             6.71D+00, 6.36D+00, 5.73D+00, 5.14D+00,           &
      &             6.70D+00, 6.37D+00, 5.78D+00, 5.14D+00,           &
      &             6.70D+00, 6.37D+00, 5.81D+00, 5.14D+00,           &
      &             6.70D+00, 6.37D+00, 5.85D+00, 5.14D+00,           &
      &             6.71D+00, 6.36D+00, 5.87D+00, 5.13D+00,           &
      &             6.73D+00, 6.36D+00, 5.90D+00, 5.12D+00/
   data((pnuat(l, iz), l=0, 3), iz=80, 89)/                             &
      &             6.86D+00, 6.38D+00, 5.92D+00, 5.11D+00,           &
      &             6.90D+00, 6.59D+00, 6.14D+00, 5.12D+00,           &
      &             6.93D+00, 6.70D+00, 6.16D+00, 5.12D+00,           &
      &             6.94D+00, 6.79D+00, 6.17D+00, 5.12D+00,           &
      &             6.96D+00, 6.84D+00, 6.17D+00, 5.12D+00,           &
      &             7.15D+00, 6.89D+00, 6.18D+00, 5.12D+00,           &
      &             7.17D+00, 6.92D+00, 6.13D+00, 5.09D+00,           &
      &             7.81D+00, 6.98D+00, 6.19D+00, 5.10D+00,           &
      &             7.72D+00, 6.96D+00, 6.47D+00, 5.13D+00,           &
      &             7.67D+00, 6.96D+00, 6.58D+00, 5.16D+00/
   data((pnuat(l, iz), l=0, 3), iz=90, 100)/                            &
      &             7.62D+00, 6.95D+00, 6.61D+00, 5.22D+00,           &
      &             7.60D+00, 6.95D+00, 6.60D+00, 5.31D+00,           &
      &             7.56D+00, 6.94D+00, 6.58D+00, 5.41D+00,           &
      &             7.53D+00, 6.94D+00, 6.56D+00, 5.51D+00,           &
      &             7.51D+00, 6.94D+00, 6.54D+00, 5.58D+00,           &
      &             7.52D+00, 6.94D+00, 6.53D+00, 5.65D+00,           &
      &             7.54D+00, 6.94D+00, 6.52D+00, 5.70D+00,           &
      &             7.57D+00, 6.94D+00, 6.50D+00, 5.74D+00,           &
      &             7.60D+00, 6.94D+00, 6.48D+00, 5.77D+00,           &
      &             7.64D+00, 6.95D+00, 6.46D+00, 5.80D+00,           &
      &             7.68D+00, 6.95D+00, 6.44D+00, 5.82D+00/

   data((qnuat(l, iz), l=0, 3), iz=00, 09)/                             &
      &               .0D+00, .0D+00, .0D+00, .0D+00,           &
      &               .8D+00, .2D+00, .0D+00, .0D+00,           &
      &              2.0D+00, .0D+00, .0D+00, .0D+00,           &
      &               .5D+00, .5D+00, .0D+00, .0D+00,           &
      &               .8D+00, 1.0D+00, .2D+00, .0D+00,           &
      &               .9D+00, 1.8D+00, .3D+00, .0D+00,           &
      &              1.6D+00, 2.1D+00, .3D+00, .0D+00,           &
      &              1.9D+00, 2.9D+00, .2D+00, .0D+00,           &
      &               .1D+00, 3.8D+00, .1D+00, .0D+00,           &
      &               .0D+00, 4.9D+00, .1D+00, .0D+00/
   data((qnuat(l, iz), l=0, 3), iz=10, 19)/                             &
      &               .0D+00, 6.0D+00, .0D+00, .0D+00,           &
      &               .6D+00, .3D+00, .1D+00, .0D+00,           &
      &               .8D+00, 1.0D+00, .2D+00, .0D+00,           &
      &              1.1D+00, 1.5D+00, .4D+00, .0D+00,           &
      &              1.5D+00, 2.0D+00, .5D+00, .0D+00,           &
      &              1.8D+00, 2.7D+00, .5D+00, .0D+00,           &
      &              2.0D+00, 3.7D+00, .3D+00, .0D+00,           &
      &               .1D+00, 4.7D+00, .2D+00, .0D+00,           &
      &               .1D+00, 5.9D+00, .0D+00, .0D+00,           &
      &               .6D+00, .3D+00, .1D+00, .0D+00/
   data((qnuat(l, iz), l=0, 3), iz=20, 29)/                             &
      &               .7D+00, .6D+00, .7D+00, .0D+00,           &
      &               .7D+00, .6D+00, 1.7D+00, .0D+00,           &
      &               .7D+00, .7D+00, 2.6D+00, .0D+00,           &
      &               .6D+00, .8D+00, 3.6D+00, .0D+00,           &
      &               .6D+00, .8D+00, 4.6D+00, .0D+00,           &
      &               .6D+00, .8D+00, 5.6D+00, .0D+00,           &
      &               .6D+00, .8D+00, 6.6D+00, .0D+00,           &
      &               .6D+00, .8D+00, 7.6D+00, .0D+00,           &
      &               .6D+00, .8D+00, 8.6D+00, .0D+00,           &
      &               .7D+00, .8D+00, 9.5D+00, .0D+00/
   data((qnuat(l, iz), l=0, 3), iz=30, 39)/                             &
      &               .9D+00, 1.1D+00, 10.0D+00, .0D+00,           &
      &              1.4D+00, 1.3D+00, .3D+00, .0D+00,           &
      &              1.8D+00, 1.9D+00, .3D+00, .0D+00,           &
      &              1.8D+00, 2.8D+00, .4D+00, .0D+00,           &
      &              2.0D+00, 3.7D+00, .3D+00, .0D+00,           &
      &               .1D+00, 4.7D+00, .2D+00, .0D+00,           &
      &               .0D+00, 5.9D+00, .1D+00, .0D+00,           &
      &               .6D+00, .3D+00, .1D+00, .0D+00,           &
      &               .8D+00, .5D+00, .7D+00, .0D+00,           &
      &               .8D+00, .5D+00, 1.7D+00, .0D+00/
   data((qnuat(l, iz), l=0, 3), iz=40, 49)/                             &
      &               .7D+00, .6D+00, 2.6D+00, .1D+00,           &
      &               .7D+00, .7D+00, 3.5D+00, .1D+00,           &
      &               .7D+00, .7D+00, 4.5D+00, .1D+00,           &
      &               .6D+00, .8D+00, 5.4D+00, .2D+00,           &
      &               .6D+00, .8D+00, 6.4D+00, .2D+00,           &
      &               .6D+00, .7D+00, 7.5D+00, .2D+00,           &
      &               .7D+00, .6D+00, 8.6D+00, .1D+00,           &
      &               .7D+00, .7D+00, 9.5D+00, .1D+00,           &
      &               .9D+00, 1.0D+00, 10.0D+00, .1D+00,           &
      &              1.5D+00, 1.2D+00, .2D+00, .1D+00/
   data((qnuat(l, iz), l=0, 3), iz=50, 59)/                             &
      &              1.7D+00, 1.9D+00, .3D+00, .1D+00,           &
      &              1.8D+00, 2.7D+00, .4D+00, .1D+00,           &
      &              1.9D+00, 3.6D+00, .4D+00, .1D+00,           &
      &               .1D+00, 4.5D+00, .3D+00, .1D+00,           &
      &               .1D+00, 5.8D+00, .1D+00, .0D+00,           &
      &               .6D+00, .2D+00, .2D+00, .0D+00,           &
      &               .7D+00, .3D+00, 1.0D+00, .0D+00,           &
      &               .4D+00, .2D+00, 2.1D+00, .3D+00,           &
      &               .3D+00, .2D+00, 2.2D+00, 1.3D+00,           &
      &               .3D+00, .2D+00, 2.1D+00, 2.4D+00/
   data((qnuat(l, iz), l=0, 3), iz=60, 69)/                             &
      &               .3D+00, .2D+00, 2.0D+00, 3.5D+00,           &
      &               .3D+00, .2D+00, 1.9D+00, 4.6D+00,           &
      &               .3D+00, .2D+00, 1.8D+00, 5.7D+00,           &
      &               .4D+00, .2D+00, 1.7D+00, 6.7D+00,           &
      &               .4D+00, .2D+00, 1.7D+00, 7.7D+00,           &
      &               .5D+00, .2D+00, 1.6D+00, 8.7D+00,           &
      &               .5D+00, .3D+00, 1.5D+00, 9.7D+00,           &
      &               .6D+00, .3D+00, 1.4D+00, 10.7D+00,           &
      &               .7D+00, .3D+00, 1.3D+00, 11.7D+00,           &
      &               .8D+00, .4D+00, 1.0D+00, 12.8D+00/
   data((qnuat(l, iz), l=0, 3), iz=70, 79)/                             &
      &               .8D+00, .5D+00, 1.0D+00, 13.7D+00,           &
      &               .8D+00, .6D+00, 1.6D+00, 14.0D+00,           &
      &               .8D+00, .7D+00, 2.4D+00, .1D+00,           &
      &               .8D+00, .8D+00, 3.3D+00, .1D+00,           &
      &               .7D+00, .9D+00, 4.2D+00, .2D+00,           &
      &               .8D+00, .9D+00, 5.1D+00, .2D+00,           &
      &               .8D+00, .9D+00, 6.1D+00, .2D+00,           &
      &               .8D+00, .9D+00, 7.1D+00, .2D+00,           &
      &               .7D+00, .8D+00, 8.3D+00, .2D+00,           &
      &               .8D+00, .8D+00, 9.3D+00, .1D+00/
   data((qnuat(l, iz), l=0, 3), iz=80, 89)/                             &
      &              1.3D+00, .7D+00, 9.9D+00, .1D+00,           &
      &              1.6D+00, 1.1D+00, .2D+00, .1D+00,           &
      &              1.9D+00, 1.8D+00, .2D+00, .1D+00,           &
      &              1.9D+00, 2.7D+00, .3D+00, .1D+00,           &
      &              1.9D+00, 3.6D+00, .4D+00, .1D+00,           &
      &               .1D+00, 4.5D+00, .3D+00, .1D+00,           &
      &               .1D+00, 5.8D+00, .1D+00, .0D+00,           &
      &               .9D+00, 6.0D+00, .1D+00, .0D+00,           &
      &               .7D+00, 6.0D+00, 1.2D+00, .1D+00,           &
      &               .6D+00, 6.0D+00, 2.1D+00, .3D+00/
   data((qnuat(l, iz), l=0, 3), iz=90, 100)/                             &
      &               .6D+00, 5.9D+00, 2.5D+00, 1.0D+00,           &
      &               .6D+00, 5.9D+00, 2.4D+00, 2.1D+00,           &
      &               .6D+00, 5.8D+00, 2.3D+00, 3.3D+00,           &
      &               .5D+00, 5.8D+00, 2.1D+00, 4.6D+00,           &
      &               .5D+00, 5.8D+00, 2.0D+00, 5.7D+00,           &
      &               .5D+00, 5.8D+00, 1.9D+00, 6.8D+00,           &
      &               .5D+00, 5.8D+00, 1.8D+00, 7.9D+00,           &
      &               .5D+00, 5.8D+00, 1.7D+00, 9.0D+00,           &
      &               .5D+00, 5.9D+00, 1.5D+00, 10.1D+00,           &
      &               .5D+00, 5.9D+00, 1.4D+00, 11.2D+00,           &
      &               .6D+00, 5.9D+00, 1.2D+00, 12.3D+00/

   if (lfree) then
      !print *,'AB test', nclass, nl, nsp
      call fremom(lmx, nclass, nl, nsp, pnu, qnu, z)
      ! Add spin polarization to d-electrons (20% of n_d)
      do ic = 1, nclass
         l = 2
         qnuc = qnu(0, l, 1, ic)
         dq0 = 0.20d0*qnuc
         qnu(0, l, 1, ic) = qnuc - dq0
         qnu(0, l, 2, ic) = qnuc + dq0
      end do
      return
   end if
   call dinit(qnu, 3*nl*nsp*nclass)
   do ic = 1, nclass
      iz = idnint(z(ic))
      call defcnf(konfdf, lmxc, z(ic))
      call atmoms(k0, q0, z(ic))
      do isp = 1, nsp
         do l = 0, nl - 1
            ! print *,'>> zeroth',ic,l,q0(l),qnuat(l,iz), qnu(0, l, isp, ic)
            if (l .le. 3) then
               if (konfig(l, ic) .eq. konfdf(l)) then
                  ! print *,' >> check',l, lmx(ic),qnu(0, l, isp, ic) , dble(qnuat(l, iz))/nsp
                  if (l .le. lmx(ic)) qnu(0, l, isp, ic) = dble(qnuat(l, iz))/nsp
                  pnu(l, isp, ic) = pnuat(l, iz)
                  !print *,'>> ',ic,l,q0(l),qnuat(l,iz)
                  ! print *,'>> first',ic,l,q0(l),qnuat(l,iz), qnu(0, l, isp, ic)
                  ! ------------- special: partially filled core states
                  if (konfdf(l) .eq. k0(l) + 1 .and. q0(l) .ne. (4*l + 2)) then
                     ! print *,'>> special',l, q0(l), qnuat(l, iz)
                     qnu(0, l, isp, ic) = dble(qnuat(l, iz) - q0(l))/nsp
                  elseif (int(pnu(l, isp, ic)) .ne. konfig(l, ic)) then
                     ! print *,'Not here?',l, lmx(ic), konfig(l, ic), pnu(l, isp, ic)
                     write (messg, 400) iz, l, konfig(l, ic), pnu(l, isp, ic)
                     call errmsg(messg, 1)
                  end if
               elseif (konfig(l, ic) .eq. konfdf(l) + 1) then
                  pnu(l, isp, ic) = konfig(l, ic) + 0.5d0
                  ! print *,'>> second',l, lmx(ic), konfig(l, ic), konfdf(l)
                  if (l .le. lmx(ic)) then
                     if (konfdf(l) .eq. k0(l)) then
                        ! ----------------- Example N 2s into core
                        ! print *,'>> konfdf: ',ic,l,q0(l),qnuat(l,iz)
                        qnu(0, l, isp, ic) = dble(qnuat(l, iz) - q0(l))/nsp
                     elseif (konfdf(l) .eq. k0(l) - 1) then
                        ! ----------------- Example U 6p into core
                        ! print *,'>> konfdf-1: ',ic,l,q0(l),qnuat(l,iz), 2*(l + l + 1)/nsp
                        qnu(0, l, isp, ic) = dble(qnuat(l, iz) - 2*(l + l + 1))/nsp
                     else
                        write (messg, 402) iz, l, konfig(l, ic)
                        if (isp .eq. 1) call errmsg(messg, 4)
                     end if
                  end if
                  write (messg, 401) iz, l
                  if (isp .eq. 1) call errmsg(messg, 0)
               elseif (konfig(l, ic) .eq. konfdf(l) - 1) then
                  ! print *,'>> third',l, lmx(ic), konfig(l, ic), konfdf(l)
                  pnu(l, isp, ic) = konfig(l, ic) + 0.5d0
                  if (l .le. lmx(ic)) then
                     if (konfig(l, ic) .eq. k0(l) - 1) then
                        ! ----------------- Example Na 2p into valence
                        ! print *,'>> konfig-1: ',ic,l,q0(l),qnuat(l,iz), 2*(l + l + 1)/nsp
                        qnu(0, l, isp, ic) = dble(qnuat(l, iz) + 2*(l + l + 1))/nsp
                     elseif (konfig(l, ic) .eq. k0(l)) then
                        ! print *,'>> konfig: ',ic,l,q0(l),qnuat(l,iz)
                        ! ----------------- Example O 2s into valence
                        qnu(0, l, isp, ic) = dble(qnuat(l, iz) + q0(l))/nsp
                     else
                        write (messg, 402) iz, l, konfig(l, ic)
                        if (isp .eq. 1) call errmsg(messg, 4)
                     end if
                  end if
                  write (messg, 401) iz, l
                  if (isp .eq. 1) call errmsg(messg, 0)
               else
                  write (messg, 402) iz, l, konfig(l, ic)
                  if (isp .eq. 1) call errmsg(messg, 4)
               end if
            else
               pnu(l, isp, ic) = l + 1.5d0
            end if
         end do
         ! print *,'0 DEFMOM:',ic,isp,qnu(0,2,isp,ic)
         fac = dsum(nl, qnu(0, 0, isp, ic), 3)
         ! AB deactivate filling f-core states
         ! print *,' FAC:',fac, nl, dsum(4 - nl, qnuat(nl, iz), 1)/nsp
         ! if (fac .gt. tiny) then
         !    fac = (dsum(nl, qnu(0, 0, isp, ic), 3) +                              &
         !       &         dsum(4 - nl, qnuat(nl, iz), 1)/nsp)/fac
         !    call dscal(nl, fac, qnu(0, 0, isp, ic), 3)
         ! end if
         !         qnu(0,0,isp,ic)=idnint(100*dsum(nl,qnu(0,0,isp,ic),3))
         !         do l=1,nl-1
         !           qnu(0,0,isp,ic)=qnu(0,0,isp,ic)-idnint(100*qnu(0,l,isp,ic))
         !           qnu(0,l,isp,ic)=idnint(100*qnu(0,l,isp,ic))*0.01d0
         !         enddo
         !         qnu(0,0,isp,ic)=qnu(0,0,isp,ic)*0.01d0
         qnu(0, 0, isp, ic) = idnint(1000*dsum(nl, qnu(0, 0, isp, ic), 3))
         ! print *,'A DEFMOM:',ic,isp,qnu(0,2,isp,ic)
         do l = 1, nl - 1
            qnu(0, 0, isp, ic) = qnu(0, 0, isp, ic) - idnint(1000*qnu(0, l, isp, ic))
            qnu(0, l, isp, ic) = idnint(1000*qnu(0, l, isp, ic))*0.001d0
         end do
         qnu(0, 0, isp, ic) = qnu(0, 0, isp, ic)*0.001d0
         ! print *,'B DEFMOM:',ic,isp,qnu(0,2,isp,ic)
      end do
      ! ..AB. neutralize qnu if l<l_default is used
      !!! dq1 = sum(qnu(0, :, :, ic)) - sum(q0(:))
      !!! !dq1 = 0.0d0
      !!! !print *,'0-->', ic,dq1
      !!! !nl-1
      !!! do l = 0, 0
      !!!    do isp = 1, nsp
      !!!       !/nl
      !!!       qnu(0, l, isp, ic) = qnu(0, l, isp, ic) - dq1/nsp
      !!!       !print *,'qnu:',ic,isp,l,qnu(0,l,isp,ic)
      !!!    end do
      !!! end do
      !!! !dq1=sum(qnu(0,:,:,ic))-sum(q0(:))
      !!! !print *,'X-->', ic,dq1
      ! ----- for spin polarized calculation: different moments for spin up
      ! ----- and spin down
      if (nsp .eq. 2) then
         do l = 0, nl - 1
            ! New approach : spin polarize d-band by 20%
            qnuc = qnu(0, l, 1, ic)
            dq0 = 0.1d0 * l * qnuc
            qnu(0, l, 1, ic) = qnuc - dq0
            qnu(0, l, 2, ic) = qnuc + dq0
            ! qnuc = qnu(0, l, 1, ic)
            ! dq0 = 0.05d0
            ! if (l >= 2 .and. qnuc > 3.0d0) dq0 = 0.50d0
            ! !           dq0=0.05d0
            ! if (lafm .and. ic + ic .gt. nclass) dq0 = -dq0
            ! if (qnuc - dabs(dq0) .gt. 0.d0 .and.                             &
            !    &          qnuc + dabs(dq0) .lt. dble(4*l + 2)) then
            !    qnu(0, l, 1, ic) = qnuc - dq0
            !    qnu(0, l, 2, ic) = qnuc + dq0
            ! end if
         end do
      end if
      ! print *,'D DEFMOM:',ic,isp,qnu(0,2,isp,ic)
   end do

400 format(' DEFMOM: Z=', i3, '   L=', i1, ',  mismatch between',         &
      &       '|        default values for CONF=', i1, ' and P=', f4.2, '.$')
401 format(' DEFMOM: Z=', i3, '   L=', i1,                               &
      &       ',  CONF differs from default value,',                     &
      &       '|        check downfolding$')
402 format(' DEFMOM: Z=', i3, '   L=', i1, ' wrong CONF=', i3, '$')
end subroutine defmom

subroutine fremom(lmx, nclass, nl, nsp, pnu, qnu, z)
   !- Returns default values for moments for free atom
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   lmx   :lmx(j) =maximum l for atom j
   !i   nclass:number of classes, atoms in same class are symmetry-related
   !i   nl    :number of l's
   !i   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
   !i   z     :nuclear charge
   !o Outputs:
   !o   pnu   :principal quantum number
   !o   qnu   :moments
   !r Remarks:
   ! ----------------------------------------------------------------------
   implicit none
   ! Passed variables:
   integer lmx(*), nclass, nl, nsp
   real(kind=8) pnu(0:nl - 1, nsp, *), qnu(0:2, 0:nl - 1, nsp, *), z(*)
   ! Local variables:
   integer ic, isp, iqnuat(0:3, 0:100), iz, konfig(0:3), l, lmxc
   ! External calls:
   ! external defcnf, dinit
   ! Intrinsic functions:
   intrinsic dble, idnint
   ! Data statements:
   data((iqnuat(l, iz), l=0, 3), iz=0, 50)/                              &
      &           0, 0, 0, 0,                                               &
      &           1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 1, 0, 0,       &
      &           2, 2, 0, 0, 2, 3, 0, 0, 0, 4, 0, 0, 0, 5, 0, 0, 0, 6, 0, 0,       &
      &           1, 0, 0, 0, 2, 0, 0, 0, 2, 1, 0, 0, 2, 2, 0, 0, 2, 3, 0, 0,       &
      &           2, 4, 0, 0, 0, 5, 0, 0, 0, 6, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0,       &
      &           2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 3, 0, 1, 0, 5, 0, 2, 0, 5, 0,       &
      &           2, 0, 6, 0, 2, 0, 7, 0, 2, 0, 8, 0, 1, 0, 10, 0, 2, 0, 10, 0,      &
      &           2, 1, 0, 0, 2, 2, 0, 0, 2, 3, 0, 0, 2, 4, 0, 0, 0, 5, 0, 0,       &
      &           0, 6, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, 2, 0, 2, 0,       &
      &           1, 0, 4, 0, 1, 0, 5, 0, 2, 0, 5, 0, 1, 0, 7, 0, 1, 0, 8, 0,       &
      &           0, 0, 10, 0, 1, 0, 10, 0, 2, 0, 10, 0, 2, 1, 0, 0, 2, 2, 0, 0/
   data((iqnuat(l, iz), l=0, 3), iz=51, 100)/                            &
      &           2, 3, 0, 0, 2, 4, 0, 0, 0, 5, 0, 0, 0, 6, 0, 0, 1, 0, 0, 0,       &
      &           2, 0, 0, 0, 2, 0, 1, 0, 2, 0, 1, 1, 2, 0, 0, 3, 2, 0, 0, 4,       &
      &           2, 0, 0, 5, 2, 0, 0, 6, 2, 0, 0, 7, 2, 0, 1, 7, 2, 0, 0, 9,       &
      &           2, 0, 0, 10, 2, 0, 0, 11, 2, 0, 0, 12, 2, 0, 0, 13, 2, 0, 0, 14,      &
      &           2, 0, 1, 14, 2, 0, 2, 0, 2, 0, 3, 0, 2, 0, 4, 0, 2, 0, 5, 0,       &
      &           2, 0, 6, 0, 2, 0, 7, 0, 1, 0, 9, 0, 1, 0, 10, 0, 2, 0, 10, 0,      &
      &           2, 1, 0, 0, 2, 2, 0, 0, 2, 3, 0, 0, 2, 4, 0, 0, 0, 5, 0, 0,       &
      &           0, 6, 0, 0, 1, 6, 0, 0, 2, 6, 0, 0, 2, 6, 1, 0, 2, 6, 2, 0,       &
      &           2, 6, 1, 2, 2, 6, 1, 3, 2, 6, 1, 4, 2, 6, 0, 6, 2, 6, 0, 7,       &
      &           2, 6, 1, 7, 2, 6, 0, 9, 2, 6, 0, 10, 2, 6, 0, 11, 2, 6, 0, 12/

   call dinit(qnu, 3*nl*nsp*nclass)
   do ic = 1, nclass
      iz = idnint(z(ic))
      call defcnf(konfig, lmxc, z(ic))
      do isp = 1, nsp
         do l = 0, nl - 1
            if (l .le. 3) then
               if (l .le. lmx(ic)) qnu(0, l, isp, ic) = dble(iqnuat(l, iz))/dble(nsp)
               pnu(l, isp, ic) = konfig(l) + .5d0
            else
               pnu(l, isp, ic) = l + 1.5d0
            end if
         end do
      end do
   end do
end subroutine fremom

subroutine atmoms(k0, q0, z)
   !-
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   z     :nuclear charge
   !o Outputs:
   !o   k0    :modified core configuration
   !o   q0    :modified moments configuration
   !r Remarks:
   ! ----------------------------------------------------------------------
   implicit none
   ! Passed variables:
   integer k0(0:3), q0(0:3)
   real(kind=8) z
   ! Local variables:
   integer iz, iq0at(0:3, 0:100), l, li, na, k, rb, cs, fr
   parameter(li=3, na=11, k=19, rb=37, cs=55, fr=87)
   ! Intrinsic functions:
   intrinsic dble, idnint

   data((iq0at(l, iz), l=0, 3), iz=0, 49)                                &
      &  /0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, &
      &   2, 1, 0, 0, 2, 2, 0, 0, 2, 3, 0, 0, 2, 4, 0, 0, 2, 5, 0, 0, &
      &   2, 6, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 1, 0, 0, 2, 2, 0, 0, &
      &   2, 3, 0, 0, 2, 4, 0, 0, 2, 5, 0, 0, 2, 6, 0, 0, 1, 0, 0, 0, &
      &   2, 0, 0, 0, 2, 0, 1, 0, 2, 0, 2, 0, 2, 0, 3, 0, 1, 0, 5, 0, &
      &   2, 0, 5, 0, 2, 0, 6, 0, 2, 0, 7, 0, 2, 0, 8, 0, 1, 0, 10, 0, &
      &   2, 0, 10, 0, 2, 1, 10, 0, 2, 2, 10, 0, 2, 3, 10, 0, 2, 4, 10, 0, &
      &   2, 5, 10, 0, 2, 6, 10, 0, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, &
      &   2, 0, 2, 0, 1, 0, 4, 0, 1, 0, 5, 0, 2, 0, 5, 0, 1, 0, 7, 0, &
      &   1, 0, 8, 0, 0, 0, 10, 0, 1, 0, 10, 0, 2, 0, 10, 0, 2, 1, 10, 0/
   data((iq0at(l, iz), l=0, 3), iz=50, 99)                               &
      &  /2, 2, 10, 0, 2, 3, 10, 0, 2, 4, 10, 0, 2, 5, 10, 0, 2, 6, 10, 0, &
      &   1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, 2, 0, 1, 1, 2, 0, 1, 2, &
      &   2, 0, 1, 3, 2, 0, 1, 4, 2, 0, 1, 5, 2, 0, 0, 7, 2, 0, 1, 7, &
      &   2, 0, 1, 8, 2, 0, 1, 9, 2, 0, 1, 10, 2, 0, 1, 11, 2, 0, 1, 12, &
      &   2, 0, 0, 14, 2, 0, 1, 14, 2, 0, 2, 14, 2, 0, 3, 14, 2, 0, 4, 14, &
      &   2, 0, 5, 14, 2, 0, 6, 14, 2, 0, 7, 14, 1, 0, 9, 14, 1, 0, 10, 14, &
      &   2, 0, 10, 14, 2, 1, 10, 14, 2, 2, 10, 14, 2, 3, 10, 14, 2, 4, 10, 14, &
      &   2, 5, 10, 14, 2, 6, 10, 14, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 1, 0, &
      &   2, 0, 2, 0, 2, 0, 1, 2, 2, 0, 1, 3, 2, 0, 1, 4, 2, 0, 0, 6, &
      &   2, 0, 0, 7, 2, 0, 1, 7, 2, 0, 0, 9, 2, 0, 0, 10, 2, 0, 0, 11/

   iz = idnint(z)

   call icopy(4, iq0at(0, iz), 1, q0, 1)

   if (iz .ge. fr) then
      k0(0) = 7
      k0(1) = 7
      k0(2) = 6
      k0(3) = 5
   elseif (iz .ge. cs) then
      k0(0) = 6
      k0(1) = 6
      k0(2) = 5
      k0(3) = 4
   elseif (iz .ge. rb) then
      k0(0) = 5
      k0(1) = 5
      k0(2) = 4
      k0(3) = 4
   elseif (iz .ge. k) then
      k0(0) = 4
      k0(1) = 4
      k0(2) = 3
      k0(3) = 4
   elseif (iz .ge. na) then
      k0(0) = 3
      k0(1) = 3
      k0(2) = 3
      k0(3) = 4
   elseif (iz .ge. li) then
      k0(0) = 2
      k0(1) = 2
      k0(2) = 3
      k0(3) = 4
   else
      k0(0) = 1
      k0(1) = 2
      k0(2) = 3
      k0(3) = 4
   end if
end subroutine atmoms

subroutine defcnf(konfig, lmxc, z)
   !- Returns default values for core configuration
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   z     :nuclear charge
   !o Outputs:
   !o   konfig:core configuration
   !o   lmxc  :maximum l for a core state
   !r Remarks:
   ! ----------------------------------------------------------------------
   implicit none
   ! Passed variables:
   real(kind=8) z
   integer konfig(0:*), lmxc
   ! Local variables:
   integer iz
   integer li, o, f, na, cl, k, ga, br, rb, in, i, ce, cs, hf, tl, at, fr
   parameter(li=3, o=8, f=9, na=11, cl=17, k=19, ga=31, br=35, rb=37, in=49, &
      &           i=53, cs=55, ce=58, hf=72, tl=81, at=85, fr=87)
   ! Intrinsic functions:
   intrinsic dble, idnint

   iz = idnint(z)

   if (iz .ge. at) then
      konfig(0) = 7
   elseif (iz .ge. i) then
      konfig(0) = 6
   elseif (iz .ge. br) then
      konfig(0) = 5
   elseif (iz .ge. cl) then
      konfig(0) = 4
   elseif (iz .ge. o) then
      konfig(0) = 3
   elseif (iz .ge. li) then
      konfig(0) = 2
   else
      konfig(0) = 1
   end if

   !     if     (iz .ge. fr) then
   !       konfig(1) = 7
   !     elseif (iz .ge. cs) then
   if (iz .ge. cs) then
      konfig(1) = 6
   elseif (iz .ge. rb) then
      konfig(1) = 5
   elseif (iz .ge. k) then
      konfig(1) = 4
   elseif (iz .ge. na) then
      konfig(1) = 3
   else
      konfig(1) = 2
   end if

   if (iz .ge. tl) then
      konfig(2) = 6
   elseif (iz .ge. in) then
      konfig(2) = 5
   elseif (iz .ge. ga) then
      konfig(2) = 4
   else
      konfig(2) = 3
   end if

   if (iz .ge. hf) then
      !     if (iz .ge. ce) then
      konfig(3) = 5
   else
      konfig(3) = 4
   end if

   !     if     (iz .ge. hf) then
   if (iz .ge. ce) then
      lmxc = 3
   elseif (iz .ge. ga) then
      lmxc = 2
   elseif (iz .ge. na) then
      lmxc = 1
   elseif (iz .ge. li) then
      lmxc = 0
   else
      lmxc = -1
   end if

end subroutine defcnf

subroutine dinit(array, leng)
   !- Initializes real(kind=8) array to zero
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   leng  :length of array
   !o Outputs:
   !o   array :real(kind=8) array set to zero
   !r Remarks:
   ! ----------------------------------------------------------------------
   implicit none
   ! Passed variables:
   integer leng
   real(kind=8) array(leng)
   ! Local variables:
   integer i, m, mp1

   m = mod(leng, 5)
   if (m .ne. 0) then
      do i = 1, m
         array(i) = 0.d0
      end do
      if (leng .lt. 5) return
   end if
   mp1 = m + 1
   do i = mp1, leng, 5
      array(i) = 0.d0
      array(i + 1) = 0.d0
      array(i + 2) = 0.d0
      array(i + 3) = 0.d0
      array(i + 4) = 0.d0
   end do
end subroutine dinit

subroutine errmsg(messg,isev)
         implicit none
         integer isev
         character*(*) messg

         print *, isev, messg

end subroutine errmsg

subroutine dscal(n,da,dx,incx)
   !- Scales a vector by a constant  dx(i) -> a * dx(i)
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   n     :lenght of dx and dy
   !i   da    :constant
   !i   dx    :vector
   !i   incx  :incrementation for x
   !o Outputs:
   !o   dx    :vector
   !r Remarks:
   !r   Adapted from: jack dongarra, linpack, 3/11/78.
   ! ----------------------------------------------------------------------
         implicit none
   ! Passed variables:
         real(kind=8) da,dx(*)
         integer incx,n
   ! Local variables:
         integer i,m,mp1,nincx
   !
         if( n.le.0 .or. incx.le.0 )return
         if(incx.ne.1) then
   ! ----- code for increment not equal to 1
           nincx = n*incx
           do i = 1,nincx,incx
             dx(i) = da*dx(i)
           enddo
         else
   ! ----- code for increment equal to 1
           m = mod(n,5)
           if( m .ne. 0 ) then
             do i = 1,m
               dx(i) = da*dx(i)
             enddo
             if( n .lt. 5 ) return
           endif
           mp1 = m + 1
           do i = mp1,n,5
             dx(i) = da*dx(i)
             dx(i + 1) = da*dx(i + 1)
             dx(i + 2) = da*dx(i + 2)
             dx(i + 3) = da*dx(i + 3)
             dx(i + 4) = da*dx(i + 4)
           enddo
         endif
end subroutine dscal

real(kind=8) function dsum(n,dx,incx)
!- Takes the sum of the values
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of dx
!i   dx    :vector
!i   incx  :incrementation for x
!o Outputs:
!o   dsum  :sum of the values
!r Remarks:
!r  Adapted from: jack dongarra, linpack, 3/11/78.
! ----------------------------------------------------------------------
      implicit none
! Passed variables:
      real(kind=8) dx(*)
      integer incx,n
! Local variables:
      integer i,nincx
!
      dsum = 0d0
      if (n .le. 0) return

      nincx = n*incx
      do i = 1, nincx,incx
        dsum = dsum + dx(i)
      enddo
end function dsum


subroutine icopy(n,x,incx,y,incy)
   !- Copies an integer vector, x, to an integer vector, y.
   ! ----------------------------------------------------------------------
   !i Inputs:
   !i   n     :lenght of x
   !i   x     :vector to copy
   !i   incx  :incrementation for x
   !i   incy  :incrementation for y
   !o Outputs:
   !o   y     :result vector is stored in y
   !r Remarks:
   ! ----------------------------------------------------------------------
         implicit none
   ! Passed variables:
         integer x(*),y(*)
         integer incx,incy,n
   ! Local variables:
         integer i,ix,iy,m,mp1
   !
         if(n.le.0) return
         if(incx.ne.1.or.incy.ne.1) then
   ! ----  code for unequal increments or equal increments not equal to 1
           ix = 1
           iy = 1
           if(incx.lt.0)ix = (-n+1)*incx + 1
           if(incy.lt.0)iy = (-n+1)*incy + 1
           do i = 1,n
             y(iy) = x(ix)
             ix = ix + incx
             iy = iy + incy
           enddo
         else
   ! ----  code for both increments equal to 1
           m = mod(n,7)
           if( m .ne. 0 ) then
             do i = 1,m
               y(i) = x(i)
             enddo
             if( n .lt. 7 ) return
           endif
           mp1 = m + 1
           do i = mp1,n,7
             y(i) = x(i)
             y(i + 1) = x(i + 1)
             y(i + 2) = x(i + 2)
             y(i + 3) = x(i + 3)
             y(i + 4) = x(i + 4)
             y(i + 5) = x(i + 5)
             y(i + 6) = x(i + 6)
           enddo
           return
         endif
end subroutine icopy


subroutine defpp(avw, nclass, nl, nsp, pp, wsr, z)
!- Returns default values for potential parameters
! ----------------------------------------------------------------------
!i Inputs:
!i   avw   :average Wigner Seitz radius
!i   nclass:number of classes, atoms in same class are symmetry-related
!i   nl    :number of l's
!i   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
!i   wsr   :Wigner-Seitz sphere radius
!i   z     :nuclear charge
!o Outputs:
!o   pp    :potential parameters
!r Remarks:
!r  See Kampur-notes Eq 3.94
! ----------------------------------------------------------------------
         implicit none
! Passed variables:
         integer nclass, nl, nsp
         real(kind=8) avw, pp(6, 0:nl - 1, nsp, *), wsr(*), z(*)
! Local variables:
         integer ic, isp, iz, k, l, lmax
         real(kind=8) enu0(0:3, 0:99), c0(0:3, 0:99), sqdel0(0:3, 0:99),    &
        &                 p0(0:3, 0:99), gamma0(0:3, 0:99), denu0(0:3, 0:99),   &
        &                 dc0(0:3, 0:99), dsqdl0(0:3, 0:99), dp0(0:3, 0:99),    &
        &                 dgamm0(0:3, 0:99)
         real(kind=8) dlns, ror0, row, sqrow, tiny, wsr0
         parameter(tiny=1.d-4)
! Data statements:
         data((enu0(l, iz), l=0, 3), iz=0, 9)/                            &
        &              -1.3596, -.4964, 1.0177, 3.1411,           &
        &               -.5419, -.4935, -.3238, -.1540,           &
        &              -1.1716, -1.2022, -1.1638, -1.0357,           &
        &               -.3131, -.2686, -.2478, -.3479,           &
        &               -.3150, -.3717, -.1680, .0738,           &
        &               -.5658, -.3073, -.1407, -.1427,           &
        &               -.7479, -.4843, -.4941, -.1707,           &
        &              -1.1607, -.5568, -.7269, -.5361,           &
        &               -.6442, -.6120, -.6067, -.4557,           &
        &               -.8698, -.8105, -.8197, -.8305/
         data((enu0(l, iz), l=0, 3), iz=10, 19)/                            &
        &              -1.0470, -1.0224, -1.0295, -.9408,           &
        &               -.2991, -.2692, -.2394, -.1864,           &
        &               -.4155, -.3015, -.2744, -.2831,           &
        &               -.4975, -.3265, -.2592, -.1573,           &
        &               -.5819, -.3726, -.3407, -.1807,           &
        &               -.8303, -.4317, -.5021, -.3795,           &
        &              -1.1261, -.5056, -.6341, -.5806,           &
        &               -.6652, -.6193, -.6213, -.4771,           &
        &               -.8146, -.7878, -.7949, -.7371,           &
        &               -.2593, -.2409, -.2189, -.2221/
         data((enu0(l, iz), l=0, 3), iz=20, 29)/                            &
        &               -.3046, -.2612, -.2074, -.2032,           &
        &               -.2957, -.2204, -.1599, -.2084,           &
        &               -.2981, -.1727, -.1433, -.2029,           &
        &               -.2940, -.1507, -.1223, -.1853,           &
        &               -.3045, -.1478, -.1129, -.1637,           &
        &               -.3290, -.1570, -.1138, -.1402,           &
        &               -.3686, -.1907, -.1250, -.1159,           &
        &               -.4212, -.2425, -.1530, -.0942,           &
        &               -.4759, -.3047, -.1961, -.3377,           &
        &               -.5126, -.3560, -.2895, -.3188/
         data((enu0(l, iz), l=0, 3), iz=30, 39)/                            &
        &               -.5257, -.3531, -.6235, -.5933,           &
        &               -.5444, -.3786, -.3688, -.3404,           &
        &               -.6949, -.3904, -.4443, -.3727,           &
        &               -.9208, -.4302, -.5217, -.3822,           &
        &              -1.1533, -.4796, -.5983, -.5722,           &
        &               -.6163, -.5710, -.5693, -.4931,           &
        &               -.7425, -.7133, -.7205, -.6998,           &
        &               -.2491, -.2346, -.2123, -.2140,           &
        &               -.2910, -.2584, -.2150, -.2068,           &
        &               -.2929, -.2258, -.1871, -.1749/
         data((enu0(l, iz), l=0, 3), iz=40, 49)/                            &
        &               -.3006, -.1925, -.1819, -.1554,           &
        &               -.3134, -.1870, -.1917, -.1622,           &
        &               -.3286, -.1934, -.1986, -.1765,           &
        &               -.3644, -.2161, -.2254, -.2014,           &
        &               -.4097, -.2603, -.2587, -.2309,           &
        &               -.4675, -.3281, -.3060, -.2854,           &
        &               -.5263, -.4072, -.3652, -.3610,           &
        &               -.5371, -.4367, -.4914, -.4880,           &
        &               -.4993, -.3988, -.7923, -.6014,           &
        &               -.5340, -.3832, -.3876, -.3251/
         data((enu0(l, iz), l=0, 3), iz=50, 59)/                            &
        &               -.6668, -.3941, -.4442, -.3939,           &
        &               -.8297, -.4181, -.4880, -.4397,           &
        &              -1.0237, -.4624, -.5563, -.4897,           &
        &               -.5553, -.5197, -.5138, -.4812,           &
        &               -.6723, -.6408, -.6486, -.6339,           &
        &               -.2314, -.2204, -.2010, -.2015,           &
        &               -.2269, -.2048, -.1872, -.1789,           &
        &               -.1188, -.0864, -.1052, -.0814,           &
        &               -.0396, -.0168, -.0523, .0030,           &
        &               -.0029, .0008, -.0348, .0402/
         data((enu0(l, iz), l=0, 3), iz=60, 69)/                            &
        &                .0148, .0208, -.0171, .0677,           &
        &                .0242, .0348, -.0041, .0854,           &
        &                .0227, .0404, .0007, .0903,           &
        &               -.0028, .0261, -.0105, .0751,           &
        &               -.0354, .0029, -.0258, .0556,           &
        &               -.0760, -.0294, -.0474, .0289,           &
        &               -.1203, -.0674, -.0731, -.0020,           &
        &               -.1604, -.1036, -.0953, -.0270,           &
        &               -.1979, -.1386, -.1166, -.0492,           &
        &               -.2772, -.2197, -.1756, -.1112/
         data((enu0(l, iz), l=0, 3), iz=70, 79)/                            &
        &               -.3004, -.2429, -.1933, -.1355,           &
        &               -.3272, -.2406, -.1791, -.4279,           &
        &               -.3657, -.2108, -.1885, -.1578,           &
        &               -.3908, -.2074, -.2004, -.1619,           &
        &               -.4152, -.2181, -.2103, -.1808,           &
        &               -.4516, -.2350, -.2279, -.1980,           &
        &               -.4959, -.2775, -.2587, -.2209,           &
        &               -.5555, -.3442, -.3050, -.2729,           &
        &               -.6144, -.4293, -.3677, -.3557,           &
        &               -.6342, -.4892, -.4689, -.4723/
         data((enu0(l, iz), l=0, 3), iz=80, 89)/                            &
        &               -.5459, -.5194, -.7049, -.5861,           &
        &               -.6300, -.3768, -.4618, -.4202,           &
        &               -.7785, -.3779, -.4848, -.4558,           &
        &               -.9568, -.4013, -.4918, -.4716,           &
        &              -1.1579, -.4404, -.5342, -.4841,           &
        &               -.5412, -.5003, -.4892, -.4569,           &
        &               -.6402, -.6043, -.6120, -.5956,           &
        &               -.2129, -.8899, -.2951, -.3834,           &
        &               -.2356, -1.0272, -.2046, -.3577,           &
        &               -.2531, -1.1667, -.1855, -.2542/
         data((enu0(l, iz), l=0, 3), iz=90, 99)/                            &
        &               -.2310, -1.2322, -.1641, -.1400,           &
        &               -.1811, -1.2525, -.1548, -.0610,           &
        &               -.1596, -1.2406, -.1359, .0019,           &
        &               -.1526, -1.2325, -.1226, .0440,           &
        &               -.1477, -1.2323, -.1136, .0687,           &
        &               -.1527, -1.2525, -.1122, .0715,           &
        &               -.1653, -1.2817, -.1180, .0654,           &
        &               -.1762, -1.3169, -.1259, .0528,           &
        &               -.1875, -1.3588, -.1380, .0325,           &
        &               -.1982, -1.4072, -.1531, .0032/

         data((c0(l, iz), l=0, 3), iz=0, 9)/                            &
        &                .7660, 2.7078, 5.3999, 8.7865,           &
        &               -.5919, 2.5133, 6.7396, 12.2934,           &
        &              -1.2510, 1.1103, 3.0981, 5.5895,           &
        &               -.1579, .2720, 1.6066, 3.3175,           &
        &               -.4247, .2077, 2.1417, 4.4518,           &
        &               -.6455, .2851, 3.3866, 7.1058,           &
        &              -1.1194, -.1808, 2.8881, 6.1186,           &
        &              -1.4880, -.4444, 2.9151, 6.1706,           &
        &              -7.1141, -.6235, 2.9210, 6.1304,           &
        &               2.4119, -.8467, 2.8098, 5.7959/
         data((c0(l, iz), l=0, 3), iz=10, 19)/                            &
        &                .3857, -1.0424, 2.1745, 4.2112,           &
        &               -.2515, .2463, .9267, 2.0982,           &
        &               -.4247, .1948, 1.0836, 2.6817,           &
        &               -.6076, .1356, 1.2323, 3.2357,           &
        &               -.8408, -.0491, 1.1953, 3.3821,           &
        &              -1.0986, -.3015, 1.0449, 3.2131,           &
        &              -1.3364, -.4903, 1.0051, 3.2116,           &
        &               2.4791, -.6614, .9112, 2.8731,           &
        &                .2314, -.8134, .8977, 2.3158,           &
        &               -.2276, .1431, .2515, 1.2263/
         data((c0(l, iz), l=0, 3), iz=20, 29)/                            &
        &               -.2752, .2463, .1072, 1.7158,           &
        &               -.2069, .5147, .0973, 2.4393,           &
        &               -.1448, .7310, .0831, 3.0395,           &
        &               -.1015, .8922, .0709, 3.5378,           &
        &               -.0879, .9913, .0464, 3.9392,           &
        &               -.1003, 1.0388, .0121, 4.2559,           &
        &               -.1452, 1.0180, -.0342, 4.4426,           &
        &               -.2156, .9391, -.0947, 4.4966,           &
        &               -.2990, .8189, -.1647, 4.5115,           &
        &               -.3975, .6551, -.2783, 4.2917/
         data((c0(l, iz), l=0, 3), iz=30, 39)/                            &
        &               -.5696, .3416, -.6306, 3.7780,           &
        &               -.7485, .0385, 1.8972, 3.0111,           &
        &               -.9358, -.1432, 1.3110, 2.8454,           &
        &              -1.1367, -.3269, 1.0913, 2.6423,           &
        &              -1.3299, -.4686, .9986, 2.6443,           &
        &               1.2398, -.6129, .8240, 2.3490,           &
        &                .1390, -.7416, .7078, 1.9024,           &
        &               -.2275, .1412, .1842, 1.0171,           &
        &               -.2760, .2333, .1044, 1.3296,           &
        &               -.2282, .4818, .1074, 1.7907/
         data((c0(l, iz), l=0, 3), iz=40, 49)/                            &
        &               -.1806, .6846, .0915, 2.1576,           &
        &               -.1626, .8049, .0441, 2.4353,           &
        &               -.1500, .9008, -.0035, 2.6930,           &
        &               -.1678, .9294, -.0747, 2.8891,           &
        &               -.2067, .9056, -.1549, 3.0328,           &
        &               -.2707, .8153, -.2474, 3.1037,           &
        &               -.3470, .6728, -.3435, 3.0891,           &
        &               -.4410, .4794, -.4922, 2.9714,           &
        &               -.5795, .2020, -.8051, 2.6327,           &
        &               -.7121, -.0340, 1.5964, 2.1157/
         data((c0(l, iz), l=0, 3), iz=50, 59)/                            &
        &               -.8667, -.1850, 1.0741, 1.9939,           &
        &              -1.0286, -.3166, .9176, 1.9259,           &
        &              -1.1916, -.4460, .7897, 1.8437,           &
        &               1.9070, -.5618, .6252, 1.6887,           &
        &                .1332, -.6725, .4821, 1.4280,           &
        &               -.2102, .1431, .0713, .7631,           &
        &               -.1709, .4287, .0793, .8249,           &
        &                .0690, 1.1939, .2138, .3165,           &
        &                .2363, 1.8581, .3160, .1546,           &
        &                .2834, 2.1030, .3641, .1342/
         data((c0(l, iz), l=0, 3), iz=60, 69)/                            &
        &                .3136, 2.2457, .4070, .1376,           &
        &                .3250, 2.2875, .4408, .1410,           &
        &                .3101, 2.1904, .4602, .1350,           &
        &                .2490, 1.9132, .4528, .1094,           &
        &                .1791, 1.6639, .4381, .0816,           &
        &                .0985, 1.4217, .4135, .0482,           &
        &                .0153, 1.2031, .3827, .0121,           &
        &               -.0557, 1.0343, .3562, -.0167,           &
        &               -.1198, .8892, .3306, -.0417,           &
        &               -.2459, .5745, .2390, -.1070/
         data((c0(l, iz), l=0, 3), iz=70, 79)/                            &
        &               -.2868, .4867, .2182, -.1330,           &
        &               -.3028, .5518, .1839, -.4283,           &
        &               -.3140, .6550, .1356, 2.5813,           &
        &               -.3235, .7508, .0825, 2.5865,           &
        &               -.3368, .8209, .0258, 2.7328,           &
        &               -.3598, .8586, -.0399, 2.8642,           &
        &               -.4034, .8357, -.1248, 2.9337,           &
        &               -.4657, .7479, -.2263, 2.9364,           &
        &               -.5363, .6058, -.3371, 2.8695,           &
        &               -.6085, .4331, -.4687, 2.7456/
         data((c0(l, iz), l=0, 3), iz=80, 89)/                            &
        &               -.7098, .1598, -.7175, 2.3877,           &
        &               -.8034, -.0666, 5.9963, 1.9087,           &
        &               -.9552, -.1923, 1.2295, 1.8055,           &
        &              -1.1215, -.3121, .9512, 1.7128,           &
        &              -1.2935, -.4268, .8072, 1.6171,           &
        &                .5630, -.5382, .6300, 1.4621,           &
        &                .0183, -.6384, .4568, 1.2676,           &
        &               -.2439, -.8991, .0875, .7008,           &
        &               -.2045, -1.0628, .1326, .8936,           &
        &               -.1633, -1.2236, .1548, .6235/
         data((c0(l, iz), l=0, 3), iz=90, 99)/                            &
        &               -.0547, -1.3157, .2340, .3200,           &
        &                .0405, -1.3514, .3159, .2127,           &
        &                .1403, -1.3524, .4100, .2041,           &
        &                .2160, -1.3530, .4904, .2064,           &
        &                .2629, -1.3594, .5533, .2027,           &
        &                .2543, -1.3803, .5807, .1774,           &
        &                .2127, -1.4056, .5882, .1469,           &
        &                .1562, -1.4328, .5852, .1145,           &
        &                .0818, -1.4635, .5677, .0777,           &
        &               -.0053, -1.4985, .5363, .0349/

         data((sqdel0(l, iz), l=0, 3), iz=0, 9)/                            &
        &                .5647, .5694, .5749, .5802,           &
        &                .5101, .6226, .6905, .7495,           &
        &                .1709, .4404, .4851, .5183,           &
        &               -.3568, .3022, .3635, .4060,           &
        &               -.3931, .3297, .4162, .4611,           &
        &               -.4746, .3838, .5135, .5814,           &
        &               -.3847, .3057, .4829, .5406,           &
        &               -.3231, .2591, .4862, .5462,           &
        &              -2.0425, .2218, .4803, .5414,           &
        &                .9888, .1726, .4721, .5305/
         data((sqdel0(l, iz), l=0, 3), iz=10, 19)/                            &
        &                .4581, .0885, .4183, .4556,           &
        &                .2871, .2820, .2893, .3266,           &
        &                .3196, .3066, .3178, .3710,           &
        &                .3459, .3254, .3388, .4034,           &
        &                .3351, .3092, .3387, .4122,           &
        &                .2856, .2637, .3251, .4052,           &
        &                .2353, .2275, .3212, .4073,           &
        &               1.1180, .1749, .3027, .3828,           &
        &                .3837, .0951, .2951, .3483,           &
        &                .2405, .2408, .1784, .2584/
         data((sqdel0(l, iz), l=0, 3), iz=20, 29)/                            &
        &                .2966, .2902, .1522, .3020,           &
        &                .3618, .3494, .1521, .3540,           &
        &                .4064, .3899, .1536, .3914,           &
        &                .4373, .4183, .1534, .4193,           &
        &                .4579, .4372, .1508, .4401,           &
        &                .4707, .4490, .1464, .4557,           &
        &                .4733, .4518, .1391, .4643,           &
        &                .4663, .4464, .1292, .4663,           &
        &                .4506, .4337, .1171, .4709,           &
        &                .4253, .4124, .1014, .4595/
         data((sqdel0(l, iz), l=0, 3), iz=30, 39)/                            &
        &                .3723, .3637, .0708, .4384,           &
        &                .3119, .3063, .4533, .3916,           &
        &                .2846, .2786, .3549, .3829,           &
        &                .2432, .2448, .3327, .3703,           &
        &                .2066, .2198, .3237, .3735,           &
        &                .6939, .1751, .2970, .3509,           &
        &                .3512, .0994, .2740, .3204,           &
        &                .2276, .2336, .1689, .2382,           &
        &                .2763, .2777, .1648, .2700,           &
        &                .3353, .3335, .1759, .3070/
         data((sqdel0(l, iz), l=0, 3), iz=40, 49)/                            &
        &                .3759, .3721, .1827, .3330,           &
        &                .4015, .3959, .1826, .3515,           &
        &                .4215, .4146, .1810, .3677,           &
        &                .4327, .4245, .1752, .3796,           &
        &                .4359, .4274, .1660, .3883,           &
        &                .4284, .4210, .1523, .3932,           &
        &                .4092, .4047, .1342, .3936,           &
        &                .3784, .3772, .1104, .3896,           &
        &                .3290, .3299, .0754, .3723,           &
        &                .2789, .2804, .4232, .3325/
         data((sqdel0(l, iz), l=0, 3), iz=50, 59)/                            &
        &                .2554, .2571, .3320, .3263,           &
        &                .2339, .2388, .3156, .3225,           &
        &                .2013, .2156, .3028, .3169,           &
        &                .9495, .1804, .2769, .3031,           &
        &                .3384, .1055, .2466, .2836,           &
        &                .2180, .2261, .1457, .2087,           &
        &                .2893, .3016, .1676, .2016,           &
        &                .3671, .4493, .2076, .0991,           &
        &                .4005, .5854, .2283, .0640,           &
        &                .4115, .6350, .2357, .0569/
         data((sqdel0(l, iz), l=0, 3), iz=60, 69)/                            &
        &                .4197, .6598, .2412, .0541,           &
        &                .4254, .6616, .2451, .0520,           &
        &                .4277, .6330, .2466, .0496,           &
        &                .4238, .5696, .2439, .0460,           &
        &                .4178, .5178, .2401, .0425,           &
        &                .4087, .4724, .2348, .0387,           &
        &                .3973, .4358, .2285, .0351,           &
        &                .3860, .4102, .2229, .0319,           &
        &                .3742, .3898, .2173, .0290,           &
        &                .3380, .3459, .1994, .0231/
         data((sqdel0(l, iz), l=0, 3), iz=70, 79)/                            &
        &                .3264, .3341, .1949, .0209,           &
        &                .3442, .3510, .1941, .0178,           &
        &                .3688, .3736, .1956, .3861,           &
        &                .3896, .3938, .1961, .3690,           &
        &                .4050, .4091, .1946, .3759,           &
        &                .4155, .4196, .1910, .3832,           &
        &                .4174, .4226, .1831, .3870,           &
        &                .4093, .4166, .1700, .3873,           &
        &                .3898, .4009, .1516, .3839,           &
        &                .3607, .3764, .1287, .3779/
         data((sqdel0(l, iz), l=0, 3), iz=80, 89)/                            &
        &                .3048, .3254, .0889, .3568,           &
        &                .2445, .2654, 1.3924, .3194,           &
        &                .2231, .2472, .3618, .3136,           &
        &                .1991, .2302, .3185, .3072,           &
        &                .1678, .2102, .3027, .2995,           &
        &                .4621, .1790, .2772, .2848,           &
        &                .3099, .1108, .2440, .2698,           &
        &                .1887, -.0272, .1542, .2084,           &
        &                .2879, .0830, .1866, .2249,           &
        &                .3282, .1052, .2039, .1777/
         data((sqdel0(l, iz), l=0, 3), iz=90, 99)/                            &
        &                .3697, .1332, .2281, .1227,           &
        &                .3952, .1512, .2457, .1003,           &
        &                .4161, .1669, .2612, .0941,           &
        &                .4315, .1777, .2730, .0910,           &
        &                .4424, .1844, .2813, .0881,           &
        &                .4456, .1838, .2837, .0832,           &
        &                .4437, .1787, .2828, .0777,           &
        &                .4387, .1712, .2800, .0720,           &
        &                .4293, .1601, .2746, .0658,           &
        &                .4151, .1454, .2667, .0589/

         data((p0(l, iz), l=0, 3), iz=0, 9)/                            &
        &                .0124, .0059, .0034, .0022,           &
        &                .0536, .0051, .0016, .0007,           &
        &               1.9796, .0165, .0056, .0026,           &
        &                .1029, .1204, .0221, .0078,           &
        &                .1005, .0930, .0137, .0051,           &
        &                .0457, .0559, .0060, .0019,           &
        &                .1627, .1433, .0075, .0026,           &
        &                .3012, .2721, .0071, .0024,           &
        &                .4003, .4698, .0075, .0025,           &
        &                .2162, .9206, .0076, .0026/
         data((p0(l, iz), l=0, 3), iz=10, 19)/                            &
        &                .0369, 4.3807, .0110, .0045,           &
        &                .2626, .1202, .0602, .0196,           &
        &                .1880, .0949, .0444, .0120,           &
        &                .1552, .0812, .0367, .0090,           &
        &                .2258, .1095, .0379, .0083,           &
        &                .4627, .2282, .0449, .0088,           &
        &                .7788, .4209, .0472, .0084,           &
        &                .6083, 1.0315, .0611, .0107,           &
        &                .1292, 4.8618, .0583, .0140,           &
        &                .4669, .2068, .6795, .0525/
         data((p0(l, iz), l=0, 3), iz=20, 29)/                            &
        &                .1937, .1021, 1.2177, .0297,           &
        &                .0785, .0541, 1.0603, .0159,           &
        &                .0467, .0370, .9192, .0107,           &
        &                .0343, .0279, .8554, .0081,           &
        &                .0286, .0224, .8482, .0067,           &
        &                .0259, .0191, .8816, .0058,           &
        &                .0259, .0176, .9883, .0054,           &
        &                .0285, .0178, 1.1893, .0053,           &
        &                .0344, .0197, 1.5316, .0049,           &
        &                .0470, .0246, 2.1736, .0053/
         data((p0(l, iz), l=0, 3), iz=30, 39)/                            &
        &                .0983, .0448, 4.5631, .0062,           &
        &                .2787, .1002, .1734, .0101,           &
        &                .4469, .1617, .0344, .0112,           &
        &                .7603, .2930, .0344, .0130,           &
        &               1.0972, .4647, .0384, .0122,           &
        &                .5372, 1.0443, .0562, .0157,           &
        &                .1265, 4.7667, .0722, .0201,           &
        &                .5654, .2363, .7325, .0761,           &
        &                .2480, .1262, .8144, .0497,           &
        &                .1011, .0771, .6033, .0313/
         data((p0(l, iz), l=0, 3), iz=40, 49)/                            &
        &                .0601, .0555, .4980, .0233,           &
        &                .0453, .0417, .4767, .0191,           &
        &                .0369, .0327, .4729, .0162,           &
        &                .0333, .0266, .5078, .0143,           &
        &                .0327, .0230, .5849, .0130,           &
        &                .0361, .0219, .7453, .0122,           &
        &                .0454, .0238, 1.0638, .0117,           &
        &                .0692, .0317, 1.7434, .0116,           &
        &                .1570, .0603, 3.9216, .0133,           &
        &                .4051, .1321, .2112, .0222/
         data((p0(l, iz), l=0, 3), iz=50, 59)/                            &
        &                .6211, .2065, .0523, .0242,           &
        &                .8404, .3005, .0419, .0262,           &
        &               1.2275, .4754, .0470, .0287,           &
        &                .8561, .9527, .0699, .0347,           &
        &                .2219, 4.4991, .1074, .0378,           &
        &                .6129, .3046, 1.2439, .1630,           &
        &                .1684, .2242, .6990, .3463,           &
        &                .0601, .2854, .2985, 2.4115,           &
        &                .0453, .3037, .2045, 5.0603,           &
        &                .0412, .2952, .1801, 6.1280/
         data((p0(l, iz), l=0, 3), iz=60, 69)/                            &
        &                .0383, .2816, .1644, 6.5901,           &
        &                .0360, .2649, .1547, 6.9788,           &
        &                .0346, .2428, .1512, 7.5654,           &
        &                .0348, .2122, .1584, 8.8180,           &
        &                .0362, .1825, .1691, 10.4385,           &
        &                .0395, .1524, .1855, 12.7341,           &
        &                .0451, .1245, .2074, 15.9006,           &
        &                .0521, .1030, .2301, 19.6558,           &
        &                .0613, .0863, .2555, 24.3682,           &
        &                .1018, .0694, .3624, 42.6589/
         data((p0(l, iz), l=0, 3), iz=70, 79)/                            &
        &                .1231, .0668, .3994, 53.3326,           &
        &                .0971, .0552, .4014, 54.0611,           &
        &                .0705, .0434, .3769, .0783,           &
        &                .0553, .0350, .3613, .0219,           &
        &                .0469, .0290, .3611, .0151,           &
        &                .0422, .0248, .3750, .0130,           &
        &                .0420, .0223, .4238, .0123,           &
        &                .0466, .0218, .5307, .0122,           &
        &                .0593, .0241, .7484, .0125,           &
        &                .0887, .0308, 1.1774, .0131/
         data((p0(l, iz), l=0, 3), iz=80, 89)/                            &
        &                .2423, .0589, 2.8433, .0158,           &
        &                .6913, .1614, .8476, .0253,           &
        &                .9684, .2365, .1240, .0279,           &
        &               1.3039, .3401, .0505, .0314,           &
        &               1.8484, .5154, .0466, .0360,           &
        &                .4246, .9717, .0657, .0454,           &
        &                .1278, 4.1286, .1062, .0486,           &
        &               1.3674, 69.1159, .8715, .1355,           &
        &                .1734, 5.9895, .4264, .1605,           &
        &                .0935, 3.0451, .3019, .4389/
         data((p0(l, iz), l=0, 3), iz=90, 99)/                            &
        &                .0548, 1.5895, .1911, 1.1073,           &
        &                .0420, 1.1029, .1403, 1.6394,           &
        &                .0344, .8322, .1087, 1.7810,           &
        &                .0301, .6890, .0904, 1.8368,           &
        &                .0275, .6139, .0798, 1.9088,           &
        &                .0263, .6114, .0771, 2.1317,           &
        &                .0261, .6530, .0783, 2.4775,           &
        &                .0271, .7248, .0818, 2.9382,           &
        &                .0298, .8538, .0888, 3.6300,           &
        &                .0353, 1.0765, .1006, 4.7188/

         data((gamma0(l, iz), l=0, 3), iz=0, 9)/                            &
        &                .4396, .1197, .0612, .0381,           &
        &                .3761, .1125, .0627, .0428,           &
        &                .3313, .1209, .0667, .0448,           &
        &                .4242, .1012, .0608, .0426,           &
        &                .4046, .0941, .0588, .0408,           &
        &                .4040, .0861, .0584, .0414,           &
        &                .3671, .0730, .0593, .0413,           &
        &                .3512, .0591, .0605, .0421,           &
        &               -.2682, .0472, .0603, .0422,           &
        &                .6227, .0421, .0620, .0432/
         data((gamma0(l, iz), l=0, 3), iz=10, 19)/                            &
        &                .4584, .0471, .0656, .0445,           &
        &                .4191, .1118, .0589, .0415,           &
        &                .4129, .1064, .0566, .0412,           &
        &                .4049, .1014, .0541, .0402,           &
        &                .3869, .0939, .0529, .0399,           &
        &                .3678, .0815, .0531, .0405,           &
        &                .3610, .0695, .0539, .0414,           &
        &                .7474, .0579, .0544, .0414,           &
        &                .4689, .0559, .0616, .0440,           &
        &                .4246, .1167, .0452, .0409/
         data((gamma0(l, iz), l=0, 3), iz=20, 29)/                            &
        &                .4268, .1165, .0270, .0397,           &
        &                .4320, .1192, .0177, .0395,           &
        &                .4345, .1201, .0130, .0394,           &
        &                .4352, .1202, .0094, .0394,           &
        &                .4350, .1195, .0064, .0394,           &
        &                .4343, .1184, .0039, .0394,           &
        &                .4328, .1171, .0014, .0395,           &
        &                .4308, .1156, -.0008, .0396,           &
        &                .4281, .1141, -.0025, .0405,           &
        &                .4239, .1123, -.0034, .0406/
         data((gamma0(l, iz), l=0, 3), iz=30, 39)/                            &
        &                .4123, .1074, -.0001, .0416,           &
        &                .3896, .1004, .0816, .0406,           &
        &                .3756, .0925, .0605, .0404,           &
        &                .3656, .0825, .0582, .0402,           &
        &                .3637, .0732, .0577, .0410,           &
        &                .6011, .0627, .0568, .0410,           &
        &                .4654, .0585, .0616, .0437,           &
        &                .4254, .1199, .0458, .0403,           &
        &                .4280, .1206, .0349, .0387,           &
        &                .4338, .1247, .0289, .0376/
         data((gamma0(l, iz), l=0, 3), iz=40, 49)/                            &
        &                .4367, .1266, .0253, .0370,           &
        &                .4374, .1263, .0218, .0367,           &
        &                .4377, .1257, .0188, .0366,           &
        &                .4373, .1242, .0158, .0367,           &
        &                .4362, .1224, .0128, .0370,           &
        &                .4344, .1205, .0097, .0375,           &
        &                .4318, .1186, .0070, .0384,           &
        &                .4263, .1159, .0055, .0395,           &
        &                .4123, .1104, .0078, .0405,           &
        &                .3922, .1031, .0850, .0390/
         data((gamma0(l, iz), l=0, 3), iz=50, 59)/                            &
        &                .3805, .0960, .0637, .0389,           &
        &                .3733, .0884, .0602, .0387,           &
        &                .3698, .0798, .0590, .0386,           &
        &                .7261, .0696, .0571, .0388,           &
        &                .4727, .0620, .0607, .0424,           &
        &                .4292, .1241, .0419, .0380,           &
        &                .4373, .1362, .0357, .0286,           &
        &                .4471, .1770, .0353, .0046,           &
        &                .4513, .2149, .0360, -.0052,           &
        &                .4516, .2265, .0364, -.0080/
         data((gamma0(l, iz), l=0, 3), iz=60, 69)/                            &
        &                .4516, .2309, .0366, -.0093,           &
        &                .4512, .2287, .0367, -.0100,           &
        &                .4501, .2178, .0367, -.0105,           &
        &                .4480, .1971, .0365, -.0107,           &
        &                .4458, .1798, .0362, -.0108,           &
        &                .4435, .1645, .0359, -.0108,           &
        &                .4411, .1520, .0357, -.0107,           &
        &                .4389, .1434, .0355, -.0106,           &
        &                .4367, .1368, .0354, -.0104,           &
        &                .4316, .1266, .0351, -.0098/
         data((gamma0(l, iz), l=0, 3), iz=70, 79)/                            &
        &                .4292, .1239, .0351, -.0095,           &
        &                .4302, .1240, .0317, -.0046,           &
        &                .4319, .1241, .0285, .0443,           &
        &                .4327, .1243, .0257, .0390,           &
        &                .4330, .1241, .0230, .0380,           &
        &                .4329, .1234, .0204, .0376,           &
        &                .4319, .1222, .0175, .0374,           &
        &                .4301, .1206, .0144, .0375,           &
        &                .4271, .1189, .0112, .0380,           &
        &                .4217, .1168, .0088, .0389/
         data((gamma0(l, iz), l=0, 3), iz=80, 89)/                            &
        &                .4022, .1128, .0091, .0401,           &
        &                .3828, .1029, .3043, .0395,           &
        &                .3749, .0963, .0739, .0392,           &
        &                .3709, .0892, .0637, .0387,           &
        &                .3712, .0814, .0612, .0383,           &
        &                .5233, .0722, .0589, .0380,           &
        &                .4615, .0634, .0612, .0417,           &
        &                .4173, .0801, .0513, .0418,           &
        &                .4364, .0674, .0418, .0335,           &
        &                .4408, .0668, .0398, .0220/
         data((gamma0(l, iz), l=0, 3), iz=90, 99)/                            &
        &                .4453, .0663, .0403, .0087,           &
        &                .4474, .0663, .0414, .0018,           &
        &                .4498, .0664, .0424, -.0014,           &
        &                .4513, .0665, .0432, -.0033,           &
        &                .4519, .0664, .0438, -.0046,           &
        &                .4511, .0660, .0438, -.0057,           &
        &                .4496, .0655, .0437, -.0066,           &
        &                .4477, .0652, .0436, -.0073,           &
        &                .4453, .0649, .0434, -.0078,           &
        &                .4424, .0648, .0431, -.0082/

         data((denu0(l, iz), l=0, 3), iz=0, 9)/                            &
        &               2.7261, .9953, -2.0405, -6.2980,           &
        &               -.4162, -.5968, -.9522, -.5152,           &
        &               -.0762, .0790, -.1209, 2.0334,           &
        &               -.2424, -.3283, -.3047, .1140,           &
        &               -.6419, -.6134, -1.1608, -.8742,           &
        &               -.7828, -1.2444, -1.8967, -.6662,           &
        &              -1.2918, -.7721, -1.6924, -.2476,           &
        &              -1.4662, -.8436, -1.1434, .6321,           &
        &              -1.2257, -1.2141, -1.1800, .5585,           &
        &               -.5611, -.7718, -.7171, 1.4913/
         data((denu0(l, iz), l=0, 3), iz=10, 19)/                            &
        &                .0150, -.1370, -.0901, 1.9483,           &
        &               -.2011, -.2822, -.3279, -.0754,           &
        &               -.2825, -.5838, -.6799, -.0232,           &
        &               -.4329, -.8079, -1.1780, -.3426,           &
        &               -.9337, -.7104, -1.3921, -.2069,           &
        &              -1.2606, -.6951, -1.0516, .3430,           &
        &              -1.1844, -.7435, -.7565, .9068,           &
        &               -.6239, -.7198, -.7015, .8357,           &
        &               -.0038, -.1694, -.1231, 1.5797,           &
        &               -.2012, -.2421, -.2861, .1008/
         data((denu0(l, iz), l=0, 3), iz=20, 29)/                            &
        &               -.5293, -.6065, -.5677, -.1203,           &
        &              -1.0432, -1.1086, -.9815, -.2144,           &
        &              -1.5055, -1.5515, -1.4495, -.3054,           &
        &              -2.0521, -2.2718, -2.1209, .1535,           &
        &              -1.8654, -2.0075, -1.8751, -.3549,           &
        &              -1.8191, -2.0800, -1.9055, -.3480,           &
        &              -1.6457, -1.9683, -1.8949, -.3496,           &
        &              -1.3554, -1.7259, -1.7960, -.3411,           &
        &              -1.0048, -1.4023, -1.6084, .1989,           &
        &               -.6627, -1.1015, -1.3054, .2179/
         data((denu0(l, iz), l=0, 3), iz=30, 39)/                            &
        &               -.1063, -.8945, -1.0724, .7744,           &
        &               -.7211, -.5677, -1.2528, .2511,           &
        &              -1.0340, -.5270, -1.1074, .3354,           &
        &              -1.1441, -.6034, -.8042, .4281,           &
        &              -1.1137, -.6900, -.6416, .9238,           &
        &               -.6086, -.6958, -.7029, .8759,           &
        &               -.0088, -.1866, -.1419, 1.5078,           &
        &               -.2097, -.2338, -.2873, -.2670,           &
        &               -.5111, -.5667, -.4983, -.5455,           &
        &              -1.0217, -1.0699, -.8698, -.8888/
         data((denu0(l, iz), l=0, 3), iz=40, 49)/                            &
        &              -1.4999, -1.4821, -1.2930, -1.3415,           &
        &              -1.8043, -1.8132, -1.5671, -1.5948,           &
        &              -1.9608, -2.0020, -1.7674, -1.8029,           &
        &              -1.9186, -2.1283, -1.8111, -1.8819,           &
        &              -1.7249, -2.0240, -1.7578, -1.8646,           &
        &              -1.3444, -1.7451, -1.5665, -1.6738,           &
        &               -.8803, -1.2480, -1.3047, -1.3398,           &
        &               -.4309, -.8492, -.9856, -.9583,           &
        &                .0647, -.7604, -.8119, -.1641,           &
        &               -.6929, -.4537, -1.0772, -1.1981/
         data((denu0(l, iz), l=0, 3), iz=50, 59)/                            &
        &               -.8734, -.4302, -.9410, -1.2044,           &
        &              -1.0047, -.5410, -.7851, -1.1797,           &
        &               -.9772, -.5764, -.6087, -.9513,           &
        &               -.6839, -.6694, -.7026, -.8518,           &
        &               -.0067, -.1842, -.1427, -.2319,           &
        &               -.2312, -.2467, -.2757, -.2896,           &
        &               -.7362, -.7003, -.5867, -.6142,           &
        &              -1.5949, -1.3444, -1.1382, -1.2277,           &
        &              -1.8609, -1.6470, -1.4302, -1.5882,           &
        &              -1.9797, -1.6485, -1.4833, -1.7550/
         data((denu0(l, iz), l=0, 3), iz=60, 69)/                            &
        &              -1.9964, -1.6908, -1.5162, -1.8405,           &
        &              -2.0207, -1.7455, -1.5291, -1.8763,           &
        &              -2.1104, -1.8021, -1.5011, -1.8408,           &
        &              -2.0542, -1.8233, -1.4561, -1.7742,           &
        &              -1.9769, -1.7929, -1.3840, -1.6717,           &
        &              -1.8294, -1.7121, -1.2970, -1.5575,           &
        &              -1.6512, -1.5863, -1.1862, -1.4065,           &
        &              -1.4728, -1.4564, -1.0763, -1.2509,           &
        &              -1.2868, -1.3112, -.9759, -1.1060,           &
        &               -.8192, -.8961, -.7462, -.7972/
         data((denu0(l, iz), l=0, 3), iz=70, 79)/                            &
        &               -.6976, -.8008, -.6630, -.6749,           &
        &               -.9550, -1.0505, -.8592, -1.5693,           &
        &              -1.3618, -1.4027, -1.3376, -1.4275,           &
        &              -1.6435, -1.7728, -1.6046, -1.6900,           &
        &              -1.8054, -2.0186, -1.8208, -1.9006,           &
        &              -1.8060, -2.2329, -1.8874, -2.0042,           &
        &              -1.6703, -2.1825, -1.9025, -2.0714,           &
        &              -1.3356, -1.9205, -1.7982, -1.9719,           &
        &               -.8947, -1.4887, -1.5294, -1.6424,           &
        &               -.4280, -1.0509, -1.2085, -1.2097/
         data((denu0(l, iz), l=0, 3), iz=80, 89)/                            &
        &               -.1474, -.7758, -.6784, -.4362,           &
        &               -.7891, -.2517, -.8508, -1.0641,           &
        &               -.9070, -.4132, -.7350, -.9857,           &
        &               -.9963, -.5289, -.5920, -.9151,           &
        &               -.9771, -.5584, -.4680, -.7428,           &
        &               -.5996, -.6234, -.6992, -.8380,           &
        &                .0021, -.1977, -.1633, -.2616,           &
        &               -.0440, -.1858, .4968, .9386,           &
        &               -.3086, -.7203, -.2366, .1165,           &
        &               -.5126, -1.1706, -.4959, -.5288/
         data((denu0(l, iz), l=0, 3), iz=90, 99)/                            &
        &               -.8883, -1.7496, -.7903, -1.0075,           &
        &              -1.2383, -1.9246, -.9422, -1.2899,           &
        &              -1.1707, -2.0277, -1.0415, -1.5805,           &
        &              -1.0689, -2.0931, -1.0782, -1.7470,           &
        &              -1.0285, -2.1150, -1.0837, -1.8418,           &
        &               -.8829, -2.0376, -1.0698, -1.8235,           &
        &               -.7724, -1.9485, -.9783, -1.7789,           &
        &               -.7613, -1.7929, -.8647, -1.6536,           &
        &               -.6321, -1.6030, -.7579, -1.5117,           &
        &               -.6753, -1.3669, -.6022, -1.2667/

         data((dc0(l, iz), l=0, 3), iz=0, 9)/                            &
        &              -1.5360, -5.4295, -10.8273, -17.6168,           &
        &               -.4793, -6.0528, -14.3276, -25.5934,           &
        &                .2964, -1.8359, -5.6077, -11.3049,           &
        &               -.8550, -1.2723, -3.8203, -7.2670,           &
        &               -.8846, -1.6735, -5.0731, -9.7684,           &
        &              -1.4029, -2.6313, -7.8605, -15.3961,           &
        &               -.4540, -1.6432, -6.1280, -12.9782,           &
        &               -.1828, -1.2686, -5.8401, -12.8655,           &
        &            -195.3532, -1.1052, -5.6951, -12.6463,           &
        &             -30.8171, -.6084, -5.0974, -11.7573/
         data((dc0(l, iz), l=0, 3), iz=10, 19)/                            &
        &              -2.1156, -.0353, -3.5055, -8.3636,           &
        &               -.4231, -1.3481, -2.1971, -4.6460,           &
        &               -.5989, -1.7006, -2.7191, -6.0249,           &
        &               -.7442, -2.0074, -3.1442, -7.2618,           &
        &               -.5789, -1.7801, -2.9414, -7.4801,           &
        &               -.3111, -1.2703, -2.3640, -6.9395,           &
        &               -.1697, -.8977, -1.8776, -6.7060,           &
        &             -52.2208, -.5371, -1.2452, -5.8239,           &
        &              -2.0806, -.0307, -.6418, -4.4851,           &
        &               -.3511, -1.1528, -.5361, -2.7399/
         data((dc0(l, iz), l=0, 3), iz=20, 29)/                            &
        &               -.8523, -2.0155, -1.0823, -3.9248,           &
        &              -1.6243, -3.2989, -1.6586, -5.4637,           &
        &              -2.2830, -4.3327, -2.1528, -6.7450,           &
        &              -2.9198, -5.2625, -2.7710, -8.1763,           &
        &              -2.8410, -5.2862, -2.4416, -8.4704,           &
        &              -2.8218, -5.3385, -2.3739, -9.0182,           &
        &              -2.6171, -5.1265, -2.2459, -9.3201,           &
        &              -2.2562, -4.6836, -2.0293, -9.3598,           &
        &              -1.8004, -4.0816, -1.7382, -9.3341,           &
        &              -1.2845, -3.3493, -1.3604, -8.8421/
         data((dc0(l, iz), l=0, 3), iz=30, 39)/                            &
        &               -.6705, -2.3290, -1.0670, -7.8667,           &
        &               -.3709, -1.6056, -9.5278, -6.4354,           &
        &               -.2807, -1.3461, -3.7007, -6.1034,           &
        &               -.2203, -1.0469, -2.8502, -5.6197,           &
        &               -.1883, -.8141, -2.3660, -5.4776,           &
        &             -17.6083, -.5200, -1.6051, -4.7327,           &
        &              -1.6432, -.0343, -.6277, -3.6299,           &
        &               -.3194, -1.1923, -.5163, -2.1084,           &
        &               -.7741, -2.0383, -.9841, -2.8591,           &
        &              -1.5380, -3.4387, -1.5532, -3.7919/
         data((dc0(l, iz), l=0, 3), iz=40, 49)/                            &
        &              -2.2153, -4.6066, -2.0378, -4.4626,           &
        &              -2.5879, -5.2301, -2.2984, -4.8714,           &
        &              -2.8193, -5.6424, -2.4319, -5.1597,           &
        &              -2.7935, -5.6338, -2.3794, -5.2705,           &
        &              -2.5694, -5.3432, -2.1860, -5.2745,           &
        &              -2.1196, -4.7090, -1.8328, -5.1780,           &
        &              -1.5436, -3.8594, -1.4144, -5.0494,           &
        &               -.9495, -2.9143, -1.0028, -4.8394,           &
        &               -.4126, -1.9154, -.7741, -4.5785,           &
        &               -.2727, -1.3821, -8.9919, -3.6844/
         data((dc0(l, iz), l=0, 3), iz=50, 59)/                            &
        &               -.2098, -1.1572, -3.6551, -3.3795,           &
        &               -.1891, -.9999, -2.9113, -3.1081,           &
        &               -.1270, -.7521, -2.3290, -2.7587,           &
        &             -39.7635, -.5136, -1.6449, -2.3175,           &
        &              -1.9299, -.0156, -.4701, -1.7459,           &
        &               -.3488, -1.2776, -.4318, -1.3325,           &
        &              -1.1974, -3.1960, -1.0437, -1.0098,           &
        &              -2.9218, -9.9863, -1.8501, -1.9904,           &
        &              -3.9238, -19.2699, -2.2343, -2.3509,           &
        &              -4.1608, -23.3550, -2.3215, -2.2380/
         data((dc0(l, iz), l=0, 3), iz=60, 69)/                            &
        &              -4.3037, -25.2258, -2.3869, -2.2036,           &
        &              -4.3483, -24.8180, -2.4239, -2.1666,           &
        &              -4.2466, -21.5931, -2.4067, -2.0728,           &
        &              -3.9265, -16.2056, -2.3250, -1.9535,           &
        &              -3.5558, -12.4221, -2.2171, -1.8092,           &
        &              -3.1348, -9.5871, -2.0872, -1.6593,           &
        &              -2.6943, -7.5378, -1.9347, -1.4792,           &
        &              -2.3124, -6.2069, -1.7916, -1.3036,           &
        &              -1.9692, -5.2291, -1.6553, -1.1442,           &
        &              -1.2649, -3.5365, -1.3072, -.8192/
         data((dc0(l, iz), l=0, 3), iz=70, 79)/                            &
        &              -1.0493, -3.0926, -1.1937, -.6888,           &
        &              -1.3979, -3.7095, -1.6442, -1.5743,           &
        &              -1.8728, -4.5477, -2.1518, -10.0465,           &
        &              -2.1787, -5.1472, -2.4339, -6.4241,           &
        &              -2.3594, -5.5330, -2.5978, -6.1105,           &
        &              -2.3720, -5.6457, -2.5842, -6.0156,           &
        &              -2.2089, -5.4461, -2.4505, -5.8264,           &
        &              -1.8422, -4.8858, -2.1524, -5.4841,           &
        &              -1.3067, -4.0114, -1.6879, -5.0175,           &
        &               -.7690, -3.0572, -1.2306, -4.5467/
         data((dc0(l, iz), l=0, 3), iz=80, 89)/                            &
        &               -.1100, -1.7121, -.6333, -3.8218,           &
        &               -.1101, -1.1940, -347.5304, -3.1083,           &
        &               -.1407, -1.0492, -5.8833, -2.9303,           &
        &               -.1875, -.9127, -3.3708, -2.6968,           &
        &               -.1967, -.7079, -2.6085, -2.3803,           &
        &              -6.9829, -.4986, -1.8585, -1.8887,           &
        &              -1.2182, -.0233, -.6482, -1.3922,           &
        &               -.0744, -.1301, -.4500, -1.5207,           &
        &              -1.0348, -.5432, -1.0779, -1.3282,           &
        &              -1.7438, -.9133, -1.5618, -1.1238/
         data((dc0(l, iz), l=0, 3), iz=90, 99)/                            &
        &              -2.7817, -1.3703, -2.2084, -2.1993,           &
        &              -3.4439, -1.4814, -2.5417, -2.4201,           &
        &              -4.0352, -1.5469, -2.8230, -2.4950,           &
        &              -4.4434, -1.5626, -3.0117, -2.5075,           &
        &              -4.6699, -1.5485, -3.1211, -2.4777,           &
        &              -4.5607, -1.4667, -3.0757, -2.3320,           &
        &              -4.2602, -1.3722, -2.9606, -2.1747,           &
        &              -3.8568, -1.2395, -2.7869, -1.9567,           &
        &              -3.3590, -1.1081, -2.5682, -1.7368,           &
        &              -2.7641, -.9223, -2.2766, -1.4271/

         data((dsqdl0(l, iz), l=0, 3), iz=0, 9)/                            &
        &              -1.0001, -1.0001, -1.0000, -1.0000,           &
        &              -1.2189, -1.0047, -.9873, -1.0057,           &
        &              -2.7095, -.7806, -.8418, -1.0075,           &
        &              -1.2463, -1.0535, -.9942, -1.0084,           &
        &              -1.3664, -1.1961, -.9803, -1.0206,           &
        &              -1.4449, -1.3317, -.9830, -1.0201,           &
        &              -1.7298, -1.5948, -.8913, -1.0148,           &
        &              -2.2498, -1.8387, -.8613, -1.0065,           &
        &              21.1403, -2.0788, -.8491, -1.0009,           &
        &              -8.0270, -2.3658, -.8113, -.9969/
         data((dsqdl0(l, iz), l=0, 3), iz=10, 19)/                            &
        &              -1.3308, -3.0526, -.7586, -.9991,           &
        &              -1.3447, -1.1697, -.9148, -.9991,           &
        &              -1.4752, -1.2899, -.9297, -1.0160,           &
        &              -1.5670, -1.4117, -.9259, -1.0265,           &
        &              -1.7175, -1.5738, -.8732, -1.0282,           &
        &              -2.2122, -1.8065, -.8058, -1.0185,           &
        &              -2.8245, -2.0340, -.7026, -1.0041,           &
        &             -14.3604, -2.4318, -.5347, -.9896,           &
        &              -1.7867, -3.2182, -.4054, -.9881,           &
        &              -1.4105, -1.2632, -.4868, -.9705/
         data((dsqdl0(l, iz), l=0, 3), iz=20, 29)/                            &
        &              -1.4571, -1.3435, -1.3452, -.9859,           &
        &              -1.4199, -1.3769, -1.8809, -.9859,           &
        &              -1.4043, -1.3970, -2.1482, -.9882,           &
        &              -1.3961, -1.3845, -2.3654, -1.0512,           &
        &              -1.3951, -1.3684, -2.2866, -.9817,           &
        &              -1.3998, -1.3434, -2.2989, -.9772,           &
        &              -1.4116, -1.3270, -2.3184, -.9738,           &
        &              -1.4272, -1.3168, -2.3395, -.9713,           &
        &              -1.4444, -1.3103, -2.3634, -.9720,           &
        &              -1.4671, -1.3043, -2.4192, -.9715/
         data((dsqdl0(l, iz), l=0, 3), iz=30, 39)/                            &
        &              -1.5410, -1.3387, -2.8817, -.9802,           &
        &              -1.7402, -1.4832, -3.3868, -.9941,           &
        &              -2.0398, -1.6534, -1.1443, -.9985,           &
        &              -2.5732, -1.8630, -.9837, -.9979,           &
        &              -3.1213, -2.0600, -.9141, -.9909,           &
        &              -7.9980, -2.4268, -.7691, -.9786,           &
        &              -1.6144, -3.2046, -.5129, -.9785,           &
        &              -1.4474, -1.3148, -.6722, -.8350,           &
        &              -1.4899, -1.4008, -1.3005, -.8438,           &
        &              -1.4601, -1.5057, -1.7118, -.8347/
         data((dsqdl0(l, iz), l=0, 3), iz=40, 49)/                            &
        &              -1.4476, -1.5724, -1.9561, -.8006,           &
        &              -1.4493, -1.5536, -2.1399, -.7704,           &
        &              -1.4514, -1.5287, -2.2583, -.7338,           &
        &              -1.4609, -1.4697, -2.3516, -.6977,           &
        &              -1.4732, -1.4238, -2.4163, -.6664,           &
        &              -1.4896, -1.3836, -2.4685, -.6486,           &
        &              -1.5096, -1.3614, -2.5356, -.6544,           &
        &              -1.5458, -1.3558, -2.7005, -.6709,           &
        &              -1.6209, -1.3929, -3.2595, -.7548,           &
        &              -1.8587, -1.5592, -3.6682, -.6718/
         data((dsqdl0(l, iz), l=0, 3), iz=50, 59)/                            &
        &              -2.1528, -1.7236, -1.3823, -.6365,           &
        &              -2.5357, -1.8903, -1.1420, -.5934,           &
        &              -3.0625, -2.0749, -1.0575, -.5568,           &
        &             -14.0363, -2.3794, -.9492, -.5044,           &
        &              -1.9786, -3.1554, -.5922, -.5413,           &
        &              -1.4704, -1.4142, -.9145, -.5975,           &
        &              -1.4268, -1.9519, -1.5034, -.0964,           &
        &              -1.3076, -4.5750, -1.6828, -1.8220,           &
        &              -1.2670, -7.2177, -1.6661, -3.4223,           &
        &              -1.2562, -8.0774, -1.6271, -3.2640/
         data((dsqdl0(l, iz), l=0, 3), iz=60, 69)/                            &
        &              -1.2545, -8.3069, -1.5953, -3.1701,           &
        &              -1.2550, -8.0091, -1.5672, -3.1058,           &
        &              -1.2564, -7.0614, -1.5401, -3.0439,           &
        &              -1.2725, -5.5496, -1.5114, -3.0131,           &
        &              -1.2936, -4.3628, -1.4820, -2.9819,           &
        &              -1.3212, -3.4006, -1.4516, -2.9603,           &
        &              -1.3510, -2.6742, -1.4181, -2.9302,           &
        &              -1.3774, -2.2059, -1.3831, -2.8931,           &
        &              -1.4017, -1.8896, -1.3453, -2.8576,           &
        &              -1.4549, -1.5095, -1.2704, -2.8217/
         data((dsqdl0(l, iz), l=0, 3), iz=70, 79)/                            &
        &              -1.4733, -1.4340, -1.2281, -2.8044,           &
        &              -1.5093, -1.4797, -1.5885, -3.8592,           &
        &              -1.5379, -1.5170, -1.8711, -2.6171,           &
        &              -1.5504, -1.5161, -2.0680, -1.1786,           &
        &              -1.5616, -1.5052, -2.2125, -.9502,           &
        &              -1.5704, -1.4739, -2.3113, -.8658,           &
        &              -1.5859, -1.4491, -2.3997, -.8068,           &
        &              -1.6074, -1.4268, -2.4803, -.7550,           &
        &              -1.6335, -1.4061, -2.5578, -.7133,           &
        &              -1.6770, -1.3960, -2.7012, -.6861/
         data((dsqdl0(l, iz), l=0, 3), iz=80, 89)/                            &
        &              -1.7437, -1.3706, -3.1182, -.6984,           &
        &              -2.1971, -1.6212, -36.0359, -.6305,           &
        &              -2.5558, -1.7528, -2.5627, -.6155,           &
        &              -3.0254, -1.9114, -1.3789, -.5892,           &
        &              -3.6106, -2.0862, -1.1689, -.5561,           &
        &              -4.3356, -2.3709, -1.0343, -.4559,           &
        &              -1.4604, -3.0758, -.7158, -.4752,           &
        &              -1.5792, -5.4056, -1.2758, -1.0219,           &
        &              -1.4774, -4.2137, -1.5388, -.5464,           &
        &              -1.4663, -4.0994, -1.7363, -.5185/
         data((dsqdl0(l, iz), l=0, 3), iz=90, 99)/                            &
        &              -1.4297, -3.9076, -1.8364, -2.1284,           &
        &              -1.3894, -3.6924, -1.8008, -2.9009,           &
        &              -1.3726, -3.4987, -1.7551, -2.9833,           &
        &              -1.3689, -3.3722, -1.7139, -2.9691,           &
        &              -1.3679, -3.2941, -1.6786, -2.9480,           &
        &              -1.3730, -3.2786, -1.6396, -2.9257,           &
        &              -1.3809, -3.3152, -1.6095, -2.9250,           &
        &              -1.3880, -3.3641, -1.5731, -2.9097,           &
        &              -1.4046, -3.4482, -1.5313, -2.9099,           &
        &              -1.4200, -3.5654, -1.4803, -2.8850/

         data((dp0(l, iz), l=0, 3), iz=0, 9)/                            &
        &               3.9997, 3.9993, 3.9989, 3.9984,           &
        &               5.0814, 3.9771, 3.9128, 4.0307,           &
        &               6.3947, 2.5867, 3.0251, 4.0408,           &
        &               5.6297, 4.1936, 3.9471, 4.0461,           &
        &               6.5087, 4.7533, 3.8368, 4.1150,           &
        &               6.7930, 5.3427, 3.8579, 4.1156,           &
        &               7.0962, 6.0474, 3.3409, 4.0878,           &
        &               7.0976, 6.4378, 3.1881, 4.0402,           &
        &               2.1451, 6.4272, 3.1186, 4.0074,           &
        &              -3.6145, 6.5567, 2.8846, 3.9792/
         data((dp0(l, iz), l=0, 3), iz=10, 19)/                            &
        &              -2.2374, 6.9468, 2.5308, 3.9859,           &
        &               6.1665, 4.6869, 3.4665, 3.9860,           &
        &               6.9929, 5.3759, 3.5305, 4.0816,           &
        &               7.6316, 5.9317, 3.4901, 4.1416,           &
        &               7.8226, 6.5725, 3.2181, 4.1579,           &
        &               7.8031, 7.1385, 2.9570, 4.1133,           &
        &               7.8817, 7.3433, 2.5245, 4.0322,           &
        &              -2.5076, 7.3938, 1.7364, 3.9438,           &
        &              -5.1134, 7.6342, .6670, 3.9097,           &
        &               6.5656, 4.4355, 1.8147, 3.8201/
         data((dp0(l, iz), l=0, 3), iz=20, 29)/                            &
        &               6.9787, 4.3619, 5.1531, 3.9112,           &
        &               6.7017, 2.5368, 6.4832, 3.9206,           &
        &               6.4810, 1.8242, 6.9541, 3.9360,           &
        &               6.3549, 1.4604, 7.2343, 4.3775,           &
        &               6.3810, 2.3523, 7.0064, 3.8949,           &
        &               6.4344, 3.1076, 6.9264, 3.8657,           &
        &               6.5331, 3.9736, 6.8460, 3.8428,           &
        &               6.6340, 4.7304, 6.7601, 3.8244,           &
        &               6.7184, 5.2160, 6.6768, 3.8229,           &
        &               6.8247, 5.4348, 6.6343, 3.8172/
         data((dp0(l, iz), l=0, 3), iz=30, 39)/                            &
        &               7.5257, 5.6321, 6.9863, 3.8642,           &
        &               7.7844, 6.2467, -6.3854, 3.9471,           &
        &               7.8954, 6.8778, .6044, 3.9762,           &
        &               8.0394, 7.3515, 3.5168, 3.9784,           &
        &               8.1677, 7.5700, 3.4951, 3.9400,           &
        &              -3.9591, 7.6618, 2.8504, 3.8691,           &
        &              -4.1311, 7.8534, 1.4833, 3.8388,           &
        &               6.7654, 3.4085, 2.7153, 3.0363,           &
        &               7.1973, 2.7662, 5.1697, 3.1020,           &
        &               6.9621, -.0952, 6.4210, 3.0804/
         data((dp0(l, iz), l=0, 3), iz=40, 49)/                            &
        &               6.7250, -1.1385, 7.0036, 2.9100,           &
        &               6.6432, -.8331, 7.3626, 2.7546,           &
        &               6.6274, -.3491, 7.5309, 2.5578,           &
        &               6.6981, .8020, 7.6117, 2.3549,           &
        &               6.7982, 2.1241, 7.5935, 2.1628,           &
        &               6.9083, 3.5843, 7.5071, 2.0213,           &
        &               6.9907, 4.7924, 7.3995, 1.9902,           &
        &               7.2168, 5.5479, 7.4008, 2.0202,           &
        &               8.0982, 5.8557, 7.7747, 2.4671,           &
        &               8.0901, 6.6482, -6.3385, 2.0065/
         data((dp0(l, iz), l=0, 3), iz=50, 59)/                            &
        &               8.2715, 7.2513, -1.0787, 1.8016,           &
        &               8.4050, 7.6735, 3.2178, 1.5587,           &
        &               8.4890, 7.9341, 3.9952, 1.3395,           &
        &              -2.2034, 7.9963, 3.6829, 1.0025,           &
        &              -5.6710, 8.1207, 2.0509, 1.0446,           &
        &               6.9547, 1.2464, 3.7978, 1.3096,           &
        &               6.8102, -4.2752, 5.9210, .2033,           &
        &               4.3640, -5.4517, 6.5912, 6.1544,           &
        &               2.5225, -4.2786, 6.5484, 8.6347,           &
        &               2.3704, -3.9918, 6.4182, 8.1435/
         data((dp0(l, iz), l=0, 3), iz=60, 69)/                            &
        &               2.3409, -3.8041, 6.3166, 7.8953,           &
        &               2.5264, -3.7233, 6.2241, 7.7382,           &
        &               2.9914, -3.7582, 6.1307, 7.6033,           &
        &               3.8462, -4.0790, 6.0144, 7.5174,           &
        &               4.6611, -4.3836, 5.8988, 7.4457,           &
        &               5.4176, -4.6373, 5.7780, 7.3944,           &
        &               6.0143, -4.6455, 5.6495, 7.3405,           &
        &               6.3993, -4.2966, 5.5210, 7.2906,           &
        &               6.6654, -3.5779, 5.3823, 7.2510,           &
        &               7.0662, -.3453, 5.0930, 7.2233/
         data((dp0(l, iz), l=0, 3), iz=70, 79)/                            &
        &               7.1483, 1.3117, 4.9404, 7.2083,           &
        &               7.3339, 1.0325, 6.1597, 8.4030,           &
        &               7.3858, .8657, 6.8839, -10.3927,           &
        &               7.3783, .8681, 7.3490, -1.7534,           &
        &               7.3932, 1.2028, 7.6298, 1.8255,           &
        &               7.3995, 1.8557, 7.7837, 2.7391,           &
        &               7.4484, 2.8655, 7.8464, 2.8485,           &
        &               7.5214, 4.0716, 7.8224, 2.7018,           &
        &               7.5650, 5.0846, 7.7270, 2.4995,           &
        &               7.7712, 5.6768, 7.6821, 2.3196/
         data((dp0(l, iz), l=0, 3), iz=80, 89)/                            &
        &               8.3108, 5.6141, 7.8434, 2.3083,           &
        &               8.2074, 6.9828, -5.7222, 1.9124,           &
        &               8.4397, 7.3742, -6.5614, 1.8395,           &
        &               8.6258, 7.7618, .1511, 1.7172,           &
        &               8.7645, 8.0416, 3.5621, 1.5545,           &
        &              -5.5186, 8.1583, 3.9732, .9986,           &
        &              -1.6385, 8.2187, 2.6282, .7776,           &
        &               7.6610, 10.0797, 5.5571, 3.9508,           &
        &               7.5585, 9.1635, 6.3389, 2.0002,           &
        &               7.2407, 9.2681, 7.0794, 2.9620/
         data((dp0(l, iz), l=0, 3), iz=90, 99)/                            &
        &               6.1954, 9.3335, 7.5081, 7.0982,           &
        &               5.2712, 9.2429, 7.3877, 8.3960,           &
        &               4.2831, 9.1558, 7.2451, 8.3674,           &
        &               3.5587, 9.0509, 7.1096, 8.2260,           &
        &               3.2598, 8.9655, 6.9851, 8.0935,           &
        &               3.7114, 8.8996, 6.8201, 7.9617,           &
        &               4.5318, 8.8450, 6.7017, 7.8699,           &
        &               5.3865, 8.8003, 6.5565, 7.7771,           &
        &               6.2631, 8.7960, 6.3825, 7.7140,           &
        &               6.8626, 8.7733, 6.1775, 7.6340/

         data((dgamm0(l, iz), l=0, 3), iz=0, 9)/                            &
        &                .0000, .0000, .0000, .0000,           &
        &               -.0354, .0017, .0012, -.0002,           &
        &                .0611, .0196, .0090, -.0002,           &
        &               -.0363, .0000, .0010, -.0002,           &
        &               -.0724, -.0074, .0031, -.0006,           &
        &               -.0786, -.0204, .0031, -.0006,           &
        &               -.0968, -.0336, .0087, -.0004,           &
        &               -.0414, -.0376, .0096, -.0002,           &
        &             -20.0621, -.0171, .0101, .0000,           &
        &              -2.5509, .0035, .0114, .0001/
         data((dgamm0(l, iz), l=0, 3), iz=10, 19)/                            &
        &               -.0698, .0470, .0131, .0000,           &
        &               -.0464, -.0113, .0059, .0000,           &
        &               -.0738, -.0177, .0066, -.0004,           &
        &               -.1028, -.0267, .0087, -.0008,           &
        &               -.1187, -.0403, .0132, -.0008,           &
        &               -.0806, -.0535, .0165, -.0006,           &
        &               -.0170, -.0527, .0207, -.0001,           &
        &              -6.3029, -.0260, .0282, .0003,           &
        &               -.1619, .0303, .0298, .0004,           &
        &               -.0509, -.0250, .0309, .0009/
         data((dgamm0(l, iz), l=0, 3), iz=20, 29)/                            &
        &               -.0620, -.0372, -.0033, .0005,           &
        &               -.0570, -.0558, -.0218, .0004,           &
        &               -.0547, -.0646, -.0250, .0004,           &
        &               -.0531, -.0650, -.0217, -.0025,           &
        &               -.0522, -.0561, -.0206, .0006,           &
        &               -.0509, -.0467, -.0171, .0007,           &
        &               -.0503, -.0376, -.0126, .0008,           &
        &               -.0507, -.0296, -.0073, .0009,           &
        &               -.0522, -.0235, -.0016, .0009,           &
        &               -.0567, -.0190, .0048, .0009/
         data((dgamm0(l, iz), l=0, 3), iz=30, 39)/                            &
        &               -.0901, -.0183, .0244, .0006,           &
        &               -.1092, -.0295, -.2651, .0002,           &
        &               -.0978, -.0443, -.0125, .0001,           &
        &               -.0531, -.0549, .0048, .0001,           &
        &                .0014, -.0561, .0096, .0003,           &
        &              -2.1749, -.0356, .0177, .0006,           &
        &               -.1231, .0202, .0258, .0006,           &
        &               -.0540, -.0372, .0229, .0067,           &
        &               -.0656, -.0544, -.0043, .0067,           &
        &               -.0623, -.0901, -.0221, .0072/
         data((dgamm0(l, iz), l=0, 3), iz=40, 49)/                            &
        &               -.0612, -.1093, -.0290, .0089,           &
        &               -.0600, -.1041, -.0325, .0103,           &
        &               -.0588, -.0957, -.0330, .0117,           &
        &               -.0568, -.0779, -.0318, .0130,           &
        &               -.0549, -.0607, -.0282, .0141,           &
        &               -.0536, -.0437, -.0222, .0145,           &
        &               -.0534, -.0308, -.0136, .0139,           &
        &               -.0617, -.0226, -.0020, .0130,           &
        &               -.1061, -.0196, .0220, .0094,           &
        &               -.1084, -.0342, -.2963, .0139/
         data((dgamm0(l, iz), l=0, 3), iz=50, 59)/                            &
        &               -.1009, -.0483, -.0338, .0154,           &
        &               -.0733, -.0586, -.0061, .0171,           &
        &               -.0220, -.0633, .0025, .0181,           &
        &              -5.4968, -.0493, .0106, .0199,           &
        &               -.2014, .0055, .0239, .0172,           &
        &               -.0560, -.0613, .0146, .0158,           &
        &               -.0634, -.1861, -.0139, .0333,           &
        &               -.0869, -.6989, -.0247, .0005,           &
        &               -.1103, -1.3294, -.0242, -.0089,           &
        &               -.1116, -1.5729, -.0225, -.0004/
         data((dgamm0(l, iz), l=0, 3), iz=60, 69)/                            &
        &               -.1121, -1.6516, -.0213, .0035,           &
        &               -.1092, -1.5748, -.0201, .0058,           &
        &               -.1010, -1.3164, -.0190, .0073,           &
        &               -.0901, -.9348, -.0172, .0089,           &
        &               -.0800, -.6701, -.0156, .0098,           &
        &               -.0717, -.4743, -.0138, .0107,           &
        &               -.0653, -.3338, -.0119, .0110,           &
        &               -.0615, -.2432, -.0101, .0107,           &
        &               -.0597, -.1792, -.0081, .0102,           &
        &               -.0620, -.0897, -.0036, .0094/
         data((dgamm0(l, iz), l=0, 3), iz=70, 79)/                            &
        &               -.0638, -.0680, -.0015, .0089,           &
        &               -.0670, -.0765, -.0192, .0295,           &
        &               -.0662, -.0852, -.0271, -.1018,           &
        &               -.0651, -.0848, -.0328, -.0119,           &
        &               -.0647, -.0801, -.0355, .0016,           &
        &               -.0638, -.0702, -.0367, .0063,           &
        &               -.0638, -.0583, -.0350, .0091,           &
        &               -.0648, -.0451, -.0303, .0111,           &
        &               -.0661, -.0326, -.0224, .0124,           &
        &               -.0758, -.0239, -.0111, .0129/
         data((dgamm0(l, iz), l=0, 3), iz=80, 89)/                            &
        &               -.1179, -.0127, .0112, .0119,           &
        &               -.0893, -.0400, -16.4482, .0152,           &
        &               -.0710, -.0494, -.1494, .0158,           &
        &               -.0330, -.0587, -.0295, .0168,           &
        &                .0222, -.0643, -.0069, .0179,           &
        &               -.8464, -.0546, .0054, .0217,           &
        &               -.0837, -.0045, .0191, .0196,           &
        &               -.0846, .0805, -.0135, -.0019,           &
        &               -.0803, .0439, -.0225, .0158,           &
        &               -.0847, .0318, -.0336, .0187/
         data((dgamm0(l, iz), l=0, 3), iz=90, 99)/                            &
        &               -.0938, .0148, -.0397, -.0129,           &
        &               -.0982, .0018, -.0373, -.0214,           &
        &               -.1088, -.0098, -.0351, -.0181,           &
        &               -.1160, -.0152, -.0331, -.0148,           &
        &               -.1179, -.0176, -.0313, -.0116,           &
        &               -.1131, -.0161, -.0289, -.0083,           &
        &               -.1046, -.0108, -.0274, -.0051,           &
        &               -.0946, -.0050, -.0255, -.0026,           &
        &               -.0874, .0017, -.0232, -.0003,           &
        &               -.0788, .0119, -.0206, .0012/

         lmax = min0(nl - 1, 3)
         call dinit(pp, 6*nl*nsp*nclass)
         do ic = 1, nclass
            iz = idnint(z(ic))
            call defwsr(wsr0, z(ic))
            row = wsr(ic)/avw
            sqrow = dsqrt(row)
!       write(*,*)sqrow
            ror0 = wsr(ic)/wsr0
            if (dabs(ror0) .lt. tiny) ror0 = 1.d0
            dlns = log(ror0)
!       write(*,*)dlns
            do isp = 1, nsp
               do l = 0, lmax
                  pp(1, l, isp, ic) = enu0(l, iz) + denu0(l, iz)*dlns
                  pp(2, l, isp, ic) = c0(l, iz) + dc0(l, iz)*dlns
                  pp(3, l, isp, ic) = sqdel0(l, iz)*ror0**dsqdl0(l, iz)*sqrow**(l + l + 1)
                  pp(4, l, isp, ic) = p0(l, iz)*ror0**dp0(l, iz)
                  pp(5, l, isp, ic) = (gamma0(l, iz) + dgamm0(l, iz)*dlns)*row**(l + l + 1)
                  pp(6, l, isp, ic) = pp(5, l, isp, ic)
                  do k = 1, 6
                     pp(k, l, isp, ic) = int(pp(k, l, isp, ic)*1000)*0.001d0
                  end do
               end do
            end do
         end do
end subroutine defpp

subroutine defwsr(wsr, z)
!- Returns default value of radii for given nuclear charge
! ----------------------------------------------------------------------
!i Inputs:
!i   z     :nuclear charge
!o Outputs:
!o   wsr   :Wigner-Seitz sphere radius (in atomic units)
!r Remarks:
!r These are equilibrium average Wigner-Seitz radii for closed-packed
!r structure of pure element.
! ----------------------------------------------------------------------
         implicit none
! Passed variables:
         real(kind=8) wsr, z
! Local variables:
         integer iz
         real(kind=8) cpwsr(0:100)
         character*72 messg
! Data statements:
         data cpwsr/2.00,                                                  &
        &       1.39, 2.55, 3.04, 2.27, 1.96, 1.66, 1.90, 1.90, 2.17, 2.89,&
        &       3.76, 3.25, 2.95, 2.63, 2.56, 2.70, 2.85, 3.71, 4.66, 3.88,&
        &       3.31, 2.99, 2.76, 2.64, 2.57, 2.52, 2.52, 2.55, 2.62, 2.78,&
        &       2.75, 2.79, 2.83, 2.94, 3.13, 4.32, 4.95, 4.22, 3.61, 3.28,&
        &       3.03, 2.91, 2.82, 2.77, 2.78, 2.84, 2.95, 3.14, 3.30, 3.45,&
        &       3.30, 3.31, 3.50, 4.31, 5.30, 4.20, 3.91, 3.80, 3.75, 3.70,&
        &       3.65, 3.60, 3.55, 3.52, 3.61, 3.67, 3.70, 3.73, 3.75, 3.56,&
        &       3.44, 3.23, 3.04, 2.93, 2.86, 2.82, 2.83, 2.88, 2.98, 3.27,&
        &       3.57, 3.62, 3.37, 3.46, 3.63, 4.44, 5.81, 4.30, 3.84, 3.52,&
        &       3.32, 3.13, 3.02, 2.96, 2.93, 2.93, 2.95, 2.99, 3.05, 3.17/
!     data cpwsr/2d0,1.693d0, 2.561d0, 3.108d0, 2.622d0, 2.107d0,
!    .      2.253d0, 2.272d0, 2.293d0, 2.401d0, 2.885d0,
!    .      3.763d0, 3.268d0, 2.955d0, 2.886d0, 2.922d0,
!    .      2.976d0, 3.145d0, 3.709d0, 4.666d0, 3.890d0,
!    .      3.310d0, 2.995d0, 2.724d0, 2.661d0, 2.575d0,
!    .      2.533d0, 2.528d0, 2.558d0, 2.625d0, 2.806d0,
!    .      3.066d0, 3.126d0, 3.144d0, 3.214d0, 3.400d0,
!    .      3.998d0, 4.967d0, 4.212d0, 3.602d0, 3.265d0,
!    .      3.067d0, 2.923d0, 2.833d0, 2.784d0, 2.785d0,
!    .      2.843d0, 2.954d0, 3.166d0, 3.426d0, 3.476d0,
!    .      3.481d0, 3.528d0, 3.698d0, 4.329d0, 5.329d0,
!    .      4.284d0, 3.548d0, 3.294d0, 3.208d0, 3.145d0,
!    .      3.099d0, 3.073d0, 3.083d0, 3.105d0, 3.144d0,
!    .      3.197d0, 3.250d0, 3.308d0, 3.534d0, 3.600d0,
!    .      3.431d0, 3.228d0, 3.070d0, 2.953d0, 2.870d0,
!    .      2.827d0, 2.830d0, 2.886d0, 2.992d0, 3.260d0,
!    .      3.592d0, 3.623d0, 3.642d0, 3.694d0, 3.843d0,
!    .      4.441d0, 5.800d0, 4.298d0, 3.857d0, 3.522d0,
!    .      3.321d0, 3.182d0, 3.080d0, 3.013d0, 2.982d0,
!    .      2.979d0, 2.990d0, 3.025d0, 3.082d0, 3.170d0/

         !if (iprint() .ge. 120) write (lunit(1), 300)
         iz = nint(z)
         if (iz .lt. 0 .or. iz .gt. 100) then
            write (messg, 400) iz
            call errmsg(messg, 3)
            wsr = cpwsr(100)
         else
            wsr = cpwsr(iz)
         end if

300      format(/' DEFWSR: reading Wigner-Seitz radii')
400      format(' DEFWSR: bad nuclear charge, Z=', i3, '$')
end subroutine defwsr

end module lmto_lib