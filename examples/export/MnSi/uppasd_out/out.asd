 Start UppASD from fortran main()
 --------------------------------------------------------------
             __  __          ___   _______    ____  ___        
            / / / /__  ___  / _ | / __/ _ \  / __/ / _ \       
           / /_/ / _ \/ _ \/ __ |_\ \/ // / / _ \_/ // /       
           \____/ .__/ .__/_/ |_/___/____/  \___(_)___/        
               /_/  /_/                                        
 --------------------------------------------------------------
                https://github.com/UppASD/UppASD               
 --------------------------------------------------------------
              Division of Materials Theory                     
              Department of Physics and Astronomy              
              Uppsala University                               
              Sweden                                           
 ---------------------Production-version-----------------------
 Git revision: v6.0.2-335-g73ea-dirty
 --------------------------------------------------------------
 RNG: Mersenne-Twister (mtprng)
 --------------------------------------------------------------
 Execution Date: 10/10/2025  Execution Time: 19:07:10:650
 --------------------------------------------------------------
 Using OpenMP with  8  threads out of  8 possible.
 --------------------------------------------------------------
 Read input file using new format
 Reading j-scaling data 
 Finished reading j-scaling data 
 jscaling_file:jscaling.dat
 jscaling_prefactor:   0.0000000000000000     
 jscaling_freq:   0.0000000000000000     
 jscaling_phase:   0.0000000000000000     
 done.
  Set up types of atoms done
  Set up chemical information of alloy done
  Set up global coordinates done
  Set up damping
  done
  Set up temperature
  done
 Read q-vectors
  -- Real-space lattice --
       1.00000     0.00000     0.00000
       0.00000     1.00000     0.00000
       0.00000     0.00000     1.00000
  -- Reciprocal lattice --
       1.00000     0.00000     0.00000
       0.00000     1.00000     0.00000
       0.00000     0.00000     1.00000
 Write input data done.
 Set up Hamiltonian
  Set up neighbour map for exchange nntype           4         146
 done
  Mount Heisenberg Hamiltonian done
  Calculate mean-field estimate of Tc:    92.2 K.
  Max eigenvalue:    11.923 meV
  Min eigenvalue:    -5.668 meV
  Set up neighbour map for Dzyaloshinskii-Moriya exchange done
  Mount Dzyaloshinskii-Moriya Hamiltonian done
 Initialize metatype 
 Initialize magnetic moments 
  Start random init spins  done
 Calculate the ferromagnetic stiffness  Tc-MFA from stiffness :      96.1  K
 done
 Calculate the DMI spiralization done
 Enter initial phase:  S 
  Performing initial phase:  1    Temp.   25.00       nsteps       1000       dt     0.1000E-14   lambda      1.000         0.000    
  IP  1   10% done. Mbar:  0.011593
  IP  1   20% done. Mbar:  0.011584
  IP  1   30% done. Mbar:  0.016482
  IP  1   40% done. Mbar:  0.020208
  IP  1   50% done. Mbar:  0.020201
  IP  1   60% done. Mbar:  0.020419
  IP  1   70% done. Mbar:  0.023565
  IP  1   80% done. Mbar:  0.021717
  IP  1   90% done. Mbar:  0.022106
  IP  1  100% done. Mbar:  0.023983
  Performing initial phase:  2    Temp.   10.00       nsteps       1000       dt     0.1000E-14   lambda      1.000         0.000    
  IP  2   10% done. Mbar:  0.024488
  IP  2   20% done. Mbar:  0.025535
  IP  2   30% done. Mbar:  0.024733
  IP  2   40% done. Mbar:  0.024465
  IP  2   50% done. Mbar:  0.024570
  IP  2   60% done. Mbar:  0.024942
  IP  2   70% done. Mbar:  0.023987
  IP  2   80% done. Mbar:  0.024216
  IP  2   90% done. Mbar:  0.024751
  IP  2  100% done. Mbar:  0.024227
  Performing initial phase:  3    Temp.   1.000       nsteps       1000       dt     0.1000E-14   lambda      1.000         0.000    
  IP  3   10% done. Mbar:  0.026393
  IP  3   20% done. Mbar:  0.027337
  IP  3   30% done. Mbar:  0.027766
  IP  3   40% done. Mbar:  0.027889
  IP  3   50% done. Mbar:  0.028171
  IP  3   60% done. Mbar:  0.028397
  IP  3   70% done. Mbar:  0.027983
  IP  3   80% done. Mbar:  0.027552
  IP  3   90% done. Mbar:  0.027178
  IP  3  100% done. Mbar:  0.027122
  Performing initial phase:  4    Temp.  0.1000       nsteps       1000       dt     0.1000E-14   lambda      1.000         0.000    
  IP  4   10% done. Mbar:  0.026934
  IP  4   20% done. Mbar:  0.026891
  IP  4   30% done. Mbar:  0.026847
  IP  4   40% done. Mbar:  0.026802
  IP  4   50% done. Mbar:  0.026698
  IP  4   60% done. Mbar:  0.026564
  IP  4   70% done. Mbar:  0.026148
  IP  4   80% done. Mbar:  0.025778
  IP  4   90% done. Mbar:  0.025229
  IP  4  100% done. Mbar:  0.024820
 Calculating LSWT magnon dispersions
 LSWT done.
 Enter measurement phase:
 ------------------------------------------
            Simulation data                
 ------------------------------------------
  Number of atoms   64000
  Number of ensembles   1
  Gradient: N
  Temperature:    0.00 K
  Number of simulation steps:     10000
  Time step:  1000.    as
  Simulation time:  10.00    ps
  Damping parameter: 0.5000    
  Sample averages: Y
  Sample moments: N
  Spin correlation: N
 ------------------------------------------
         Progress of simulation:           
  MP   5% done. Mbar:  0.023249. Ebar:   -0.884304. U: 0.66625. 
  MP  10% done. Mbar:  0.020496. Ebar:   -0.885031. U: 0.66263. 
  MP  15% done. Mbar:  0.017782. Ebar:   -0.885531. U: 0.65473. 
  MP  20% done. Mbar:  0.015343. Ebar:   -0.886004. U: 0.64126. 
  MP  25% done. Mbar:  0.013103. Ebar:   -0.886427. U: 0.62274. 
  MP  30% done. Mbar:  0.011603. Ebar:   -0.886870. U: 0.59878. 
  MP  35% done. Mbar:  0.010574. Ebar:   -0.887236. U: 0.57632. 
  MP  40% done. Mbar:  0.009585. Ebar:   -0.887561. U: 0.55348. 
  MP  45% done. Mbar:  0.008403. Ebar:   -0.887809. U: 0.52660. 
  MP  50% done. Mbar:  0.007087. Ebar:   -0.888073. U: 0.49034. 
  MP  55% done. Mbar:  0.006652. Ebar:   -0.888297. U: 0.45047. 
  MP  60% done. Mbar:  0.006430. Ebar:   -0.888459. U: 0.41377. 
  MP  65% done. Mbar:  0.006014. Ebar:   -0.888705. U: 0.37726. 
  MP  70% done. Mbar:  0.007209. Ebar:   -0.888958. U: 0.35538. 
  MP  75% done. Mbar:  0.008178. Ebar:   -0.889127. U: 0.36342. 
  MP  80% done. Mbar:  0.007573. Ebar:   -0.889296. U: 0.37578. 
  MP  85% done. Mbar:  0.005708. Ebar:   -0.889533. U: 0.36664. 
  MP  90% done. Mbar:  0.004197. Ebar:   -0.889671. U: 0.33237. 
  MP  95% done. Mbar:  0.003427. Ebar:   -0.889783. U: 0.28482. 
  MP 100% done. Mbar:  0.002996. Ebar:   -0.889895. U: 0.23115. 
 Simulation finished
 ----------------TIME CONSUMPTION REPORT----------------------
 CATEGORY                 TIME(sec)           PERCENT
  Startup                  1.98E-01               0.1
  Initial                  2.36E-01               0.1
  Measurement              7.51E-01               0.3
  Hamiltonian              2.49E+02              92.7
  Evolution                7.31E+00               2.7
  RNG                      5.98E+00               2.2
  Energy                   9.51E-01               0.4
  Moments                  8.71E-01               0.3
  PrintRestart             1.08E-02               0.0
  SpinCorr                 2.73E-01               0.1
 -------------------------------------------------------------
 Total wall time           2.69E+02      Sampled percent  98.9
 --------------MAIN PROGRAM PHASES TIME REPORT----------------
 CATEGORY                   EXECUTION TIME    TOTAL CPU TIME
  Time for initialization :           0.66              1.59
  Time for initial phase  :       11005.28            772.47
  Time for meas. phase    :       34882.30           1377.33
  Time for cleanup        :           0.01              0.01
 -------------------------------------------------------------
  Time total              :       45888.25           2151.39
  Time for one meas. iter :           3.49              0.14
 -------------------MEMORY CONSUMPTION REPORT-----------------------
 173 allocations and 169 deallocations, remaining memory(B):128
 Memory occupation peak: 448802 kB
 For the array: "dm_vect" in routine "allocate_dmhamiltoniandata"
 -----------------END MEMORY CONSUMPTION REPORT---------------------
