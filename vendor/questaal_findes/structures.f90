! Minimal stand-alone compatibility types for the in-memory FINDES entry.
module structures
  implicit none
  type str_ctrl
    integer :: nesabc(3) = [24,24,24]
    integer :: nspec = 0
    integer :: modep(3) = [2,2,2]
    integer, pointer :: ipc(:) => null()
    integer, pointer :: ips(:) => null()
    real(8) :: rmines = 1.2d0
    real(8) :: rmaxes = 2.5d0
    real(8) :: sclwsr = 1.0d0
    real(8) :: omax1(3) = 0.0d0
    real(8) :: omax2(3) = 0.0d0
    real(8) :: wsrmax = 0.0d0
    character(8), pointer :: clabl(:) => null()
    character(8), pointer :: spid(:) => null()
  end type str_ctrl

  type str_lat
    real(8) :: avw = 0.0d0
    real(8) :: plat(3,3) = 0.0d0
    real(8) :: plat0(3,3) = 0.0d0
    real(8) :: qlat(3,3) = 0.0d0
    real(8), pointer :: pos(:,:) => null()
  end type str_lat

  type str_site
  end type str_site
end module structures
