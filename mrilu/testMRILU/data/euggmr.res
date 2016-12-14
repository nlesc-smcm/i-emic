
Enter name input file:
 
Construct Convec-Diff Eqn
Nr Intvals =    16
Stretching = 1.00E+03,  Coeff U_x = 0.00E+00,  Coeff U_y = 0.00E+00
 
Use Exact elimination.
Use Gustafsson in Multi Level Preconditioner.
CLSOnce =  T,  NLSFctr = 0.100,
EpsW = 0.050,  ElmFctr = 0.150,  RedFctr = 1.000,  SchTol = 0.00E+00,
DensLim = 0.900,  GlobFrac = 0.000,  LocFrac = 0.000,  SparsLim = 0.975,

ILUT factorisation of last sparse block.
CompFctr = 1.00E+00  ,  DropTol = 5.00E-03  ,  CPivTol = 1.00E+00
LUTol = 1.00E-10,  SingLU = T
 
Method: GMRES(m), with m = 30
Max Nr Iters =   40,  Residual <= 1.00E-06,  Decrease residual <= 1.00E-06

Number of equations:        289
Default block size:           1
Number of blocks:           289

Step  0:    145 of the    289 unknowns eliminated exactly.
Step  1:     54 of the    144 unknowns eliminated approximately.
Step  2:     33 of the     90 unknowns eliminated approximately.
Step  3:     16 of the     57 unknowns eliminated approximately.
Step  4:     11 of the     41 unknowns eliminated approximately.
Step  5:      8 of the     30 unknowns eliminated approximately.
Step  6:      7 of the     22 unknowns eliminated approximately.
Step  7:      5 of the     15 unknowns eliminated approximately.

Warning from incldup!
   After LDU-factorisation of last submatrix:
   Last diagonal element is (nearly) zero!
   Value of this element is changed into 1.00588E-03

Time ILDU factorization of last block:    0.00

Total time Preconditioner:    0.01


Number of non-zeros in the LDU-decomposition:      1240
Number of non-zeros in last block:                   74
Average number of non-zeros per row in L+D+U:         8.61

 Iter.            Residual  Reduction residual   Flops/unknown
     0        5.193809E-01        1.000000E+00           19.61
     1        6.662120E-02        1.282704E-01           41.75
     2        4.276840E-03        8.234495E-03           65.42
     3        7.595369E-04        1.462389E-03           91.10
     4        1.501883E-04        2.891679E-04          118.79
     5        1.467709E-05        2.825881E-05          148.49
     6        1.200901E-06        2.312178E-06          180.21
     7        1.284843E-07        2.473797E-07          221.09
Time Solution:    0.00


Absolute error:  2.8110400382419520E-08
