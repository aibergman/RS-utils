program process_lattice_main
  use rsutil_c_interface
  use iso_c_binding
  implicit none

  call process_lattice("periodic", 0_c_int, 1.2d0, 2.5d0)

end program process_lattice_main