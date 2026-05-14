program process_lattice_main
  use rsutil_c_interface
  use iso_c_binding
  implicit none

  call process_lattice("periodic", 0_c_int, 1.2d0, 2.5d0, -1_c_int, 0.0d0, 0.0d0, 1.0d0, &
       0_c_int, 0_c_int, 0_c_int, -1.0d0, 20_c_int, 0_c_int)

end program process_lattice_main
