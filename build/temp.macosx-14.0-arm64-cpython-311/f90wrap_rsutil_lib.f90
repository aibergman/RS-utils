! Module rsutil_lib defined in file /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90

subroutine f90wrap_rsutil_lib__process_lattice_impl
    use rsutil_lib, only: process_lattice_impl
    implicit none
    
    call process_lattice_impl()
end subroutine f90wrap_rsutil_lib__process_lattice_impl

subroutine f90wrap_rsutil_lib__array__new_basis(dummy_this, nd, dtype, dshape, dloc)
    use system_data
    use rsutil_config
    use input_io
    use lattice_io
    use spheres
    use potpar_wrapper
    use rsutil_lib, only: rsutil_lib_new_basis => new_basis
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    if (allocated(rsutil_lib_new_basis)) then
        dshape(1:2) = shape(rsutil_lib_new_basis)
        dloc = loc(rsutil_lib_new_basis)
    else
        dloc = 0
    end if
end subroutine f90wrap_rsutil_lib__array__new_basis

subroutine f90wrap_rsutil_lib__array__new_label(dummy_this, nd, dtype, dshape, dloc)
    use system_data
    use rsutil_config
    use input_io
    use lattice_io
    use spheres
    use potpar_wrapper
    use rsutil_lib, only: rsutil_lib_new_label => new_label
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    if (allocated(rsutil_lib_new_label)) then
        dshape(1:2) = (/len(rsutil_lib_new_label(1)), shape(rsutil_lib_new_label)/)
        dloc = loc(rsutil_lib_new_label)
    else
        dloc = 0
    end if
end subroutine f90wrap_rsutil_lib__array__new_label

subroutine f90wrap_rsutil_lib__array__new_atomic_number(dummy_this, nd, dtype, dshape, dloc)
    use system_data
    use rsutil_config
    use input_io
    use lattice_io
    use spheres
    use potpar_wrapper
    use rsutil_lib, only: rsutil_lib_new_atomic_number => new_atomic_number
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(rsutil_lib_new_atomic_number)) then
        dshape(1:1) = shape(rsutil_lib_new_atomic_number)
        dloc = loc(rsutil_lib_new_atomic_number)
    else
        dloc = 0
    end if
end subroutine f90wrap_rsutil_lib__array__new_atomic_number

subroutine f90wrap_rsutil_lib__get__i(f90wrap_i)
    use rsutil_lib, only: rsutil_lib_i => i
    implicit none
    integer, intent(out) :: f90wrap_i
    
    f90wrap_i = rsutil_lib_i
end subroutine f90wrap_rsutil_lib__get__i

subroutine f90wrap_rsutil_lib__set__i(f90wrap_i)
    use rsutil_lib, only: rsutil_lib_i => i
    implicit none
    integer, intent(in) :: f90wrap_i
    
    rsutil_lib_i = f90wrap_i
end subroutine f90wrap_rsutil_lib__set__i

! End of module rsutil_lib defined in file /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90

