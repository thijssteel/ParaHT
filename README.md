# ParaHT
Implementation of ParaHT, a two-stage, parallel Hessenberg-triangular reduction algorithm

This implementation was used in the paper: Thijs Steel, Raf Vandebril; Parallel two-stage reduction to Hessenberg-triangular form

To run the code, make and edit a make.inc file (an example is provided). The just run `make` to build the test executable.

Example output is:

```
# Run the parallel algorithm on a 1000x1000 randomly generated matrix
./bin/profilereduction 1000 1 16 8 8 0
n:  1000 ;runtime:        1.208396 ;error:    3.970167E-15    3.780137E-15 ;phase 1 time:        0.494774 ;phase 2 time:        0.713622 ;r:  16 ;p:   8 ;nq:   8
# Use LAPACK to reduce the same matrix
./bin/profilereduction 1000 2
n:  1000 ;runtime:        3.635724 ;error:    3.663934E-15    3.169368E-15
```

This shows that the algorithm is approximately three times as fast as LAPACK for this test pencil.