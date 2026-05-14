# Example inputs for UppASD. The structure is:
- One main input file. Hardcoded to inpsd.dat name
- External files for atom positions, magnetic moments, and exchange interactions.

- inpsd.dat contains basis vectors and other simulation inputs, including names of the external files. A minimal working example is shown below:

simid bccFe100
ncell     10        10       10
BC        P         P         P
cell      1.00000   0.00000   0.00000
          0.00000   1.00000   0.00000
          0.00000   0.00000   1.00000
Sym       1

posfile    ./posfile
momfile    ./momfile
exchange   ./jfile
dm         ./dmfile

- posfile contains basis vectors and type distinctions.
First column is running index, second column is atom type (1,2,...), and the last three columns are coordinates (fractional or cartesian depending on flags in inpsd.dat)

1 1   0.000000  0.000000  0.000000
2 1   0.500000  0.500000  0.500000

- momfile contains magnetic moment info: 
First column is running index. Second column is NOT type but chemistry (for random alloys) and is always 1 in ordered structures. Then follows magnetic moment magnitude (in muB) and a unit vector of the starting magnetic direction.

1 1 2.23 1.0 0.0 0.0
2 1 2.23 1.0 0.0 0.0

- jfile contains Jij and rij for mapping the Heisenberg term of the spin Hamiltonian.
First column is type/site i, second is j, 3-5 is r_ij and 6 is J_ij (in mRy)

1 1  0.5       0.5       0.5       1.33767484769984
1 1  1.0       0.0       0.0       0.75703576545650
1 1  1.0       1.0       0.0       -0.05975437643846
1 1  1.5       0.5       0.5       -0.08819834160658

- dmfile (optional) contains Dij and rij. Follows the same conventions as jfile, but Dij is a vector from columns 6-8.
1 1  0.5       0.5       0.5  0.01 -0.01 0.00
1 1  1.0       0.0       0.0  0.02  0.02 0.02

jfile can be either symmetry reduced, or symmetry complete (controlled by `sym` flag in inpsd.dat)