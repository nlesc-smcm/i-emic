
Enter name input file:
 
Construct Convec-Diff Eqn
Nr Intvals =    16
Stretching = 1.00E+03,  Coeff U_x = 0.00E+00,  Coeff U_y = 0.00E+00
 
Use Reverse Cuthill-McKee reordering.
Use Exact elimination.
GusFctr = 1.00E+00
CLSOnce =  T,  NLSFctr = 0.100,
EpsW = 0.050,  ElmFctr = 0.150,  RedFctr = 1.000,  SchTol = 0.00E+00,
DensLim = 0.900,  GlobFrac = 0.000,  LocFrac = 0.000,  SparsLim = 0.975,

ILUT factorisation of last sparse block.
CompFctr = 1.00E+00  ,  DropTol = 5.00E-03  ,  CPivTol = 1.00E+00
LUTol = 1.00E-10,  SingLU = T
 
Method: CG
Max Nr Iters =   40,  Residual <= 1.00E-06,  Decrease residual <= 1.00E-06

Number of equations:        289
Default block size:           1
Number of blocks:           289

Step  0:    145 of the    289 unknowns eliminated exactly.
Step  1:     54 of the    144 unknowns eliminated approximately.
Step  2:     32 of the     90 unknowns eliminated approximately.
Step  3:     16 of the     58 unknowns eliminated approximately.
Step  4:     10 of the     42 unknowns eliminated approximately.
Step  5:      8 of the     32 unknowns eliminated approximately.
Step  6:      8 of the     24 unknowns eliminated approximately.
Step  7:      6 of the     16 unknowns eliminated approximately.

Warning from incldup!
   After LDU-factorisation of last submatrix:
   Last diagonal element is (nearly) zero!
   Value of this element is changed into 1.12247E-03

Time ILDU factorization of last block:    0.00

Total time Preconditioner:    0.00


Number of non-zeros in the LDU-decomposition:      1260
Number of non-zeros in last block:                   72
Average number of non-zeros per row in L+D+U:         8.75

 Iter.            Residual  Reduction residual   Flops/unknown
     0        7.882791E-02        1.000000E+00           19.25
     1        1.785822E-02        2.265470E-01           42.98
     2        4.168559E-03        5.288177E-02           66.71
     3        4.575255E-04        5.804105E-03           90.44
     4        6.780956E-05        8.602228E-04          114.17
     5        5.917878E-06        7.507339E-05          137.90
     6        6.652382E-07        8.439120E-06          161.63
Time Solution:    0.00


Absolute error:  4.6735017619412011E-07
