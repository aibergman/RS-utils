# RS-LMTO format.
Lattice information is found in the lattice.nml. Lattice vectors in a(:,:) and basis vectors in crd(). Type is determined by izp. Cartesian coordinates (with alat=1).

&lattice
nbulk_bulk = 2
ntot = 2
nbas = 2
nrec = 2
a(:, 1) = 1.00000000, 0.00000000, 0.00000000
a(:, 2) = 0.00000000, 1.00000000, 0.00000000
a(:, 3) = 0.00000000, 0.00000000, 1.00000000
!
crd(:, 1) = 0.0, 0.0, 0.0
crd(:, 2) = 0.5, 0.5, 0.5
!
izp(1) = 1
izp(2) = 2
!
no(1) = 1
no(2) = 2
!
iu(1) = 1
iu(2) = 2
!
ib(1) = 1
ib(2) = 2
!
irec(1) = 1
irec(2) = 2
/

Exchange information comes in output files looking like below. Columns itype, jtype, r_ij j_ij |r_ij|
       1       2     -0.500000    0.223800   -0.276200      0.200485     0.613492
       1       3      0.223800   -0.276200   -0.500000      0.200398     0.613492
       1       4     -0.276200   -0.500000    0.223800      0.200433     0.613492
       1       1     -2.000000    0.000000    0.000000     -0.005870     2.000000
       1       1     -1.000000   -1.000000   -1.000000      0.000636     1.732051

DMI information follows a similar scheme:
Columns itype, jtype, r_ij, d_ij , |r_ij|
       1       2     -0.500000    0.223800   -0.276200      0.001839   -0.002613    0.003801     0.613492
       1       3      0.223800   -0.276200   -0.500000     -0.002425    0.002129    0.003485     0.613492
       1       4     -0.276200   -0.500000    0.223800      0.002352    0.001848   -0.003948     0.613492
       1       1     -2.000000    0.000000    0.000000      0.000054   -0.000041    0.000138     2.000000
       1       1     -1.000000   -1.000000   -1.000000     -0.000031   -0.000027    0.000114     1.732051
       1       1     -1.000000   -1.000000    0.000000      0.000190    0.000118   -0.000145     1.414214

The current `fable_v4` exchange writer emits `jij`, `dij`, and `aij` values
in mRy after multiplying the internal values by `1e3/(4*pi)`. The nine
`aij` values are a raw anisotropic `A` tensor, not the complete exchange
tensor. Because the writer outputs the whole Fortran array, the values are in
column-major order: flat values `a1..a9` reconstruct as
`A(row,column) = a[ row + 3*column ]`. The complete tensor is derived as
`J*I + [[0,Dz,-Dy],[-Dz,0,Dx],[Dy,-Dx,0]] + A`.
