
Enter name input file:
 
Construct Convec-Diff Eqn
Nr Intvals =    32
Stretching = 1.00E+03,  Coeff U_x = 0.00E+00,  Coeff U_y = 0.00E+00
 
GusFctr = 1.00E+00
CLSOnce =  T,  NLSFctr = 0.100,
EpsW = 0.050,  ElmFctr = 0.150,  RedFctr = 1.000,  SchTol = 0.00E+00,
DensLim = 0.900,  GlobFrac = 0.000,  LocFrac = 0.000,  SparsLim = 0.975,

ILUT factorisation of last sparse block.
CompFctr = 1.00E+00  ,  DropTol = 5.00E-03  ,  CPivTol = 1.00E+00
LUTol = 1.00E-10,  SingLU = T
 
Method: CG
Max Nr Iters =   40,  Residual <= 1.00E-06,  Decrease residual <= 1.00E-06

Number of equations:       1089
Default block size:           1
Number of blocks:          1089

Step  0:      0 of the   1089 unknowns eliminated exactly.
Visualize Schur-complement with:  vsm SchurComp01

Step  1:    541 of the   1089 unknowns eliminated approximately.
Visualize Schur-complement with:  vsm SchurComp02

Step  2:    222 of the    548 unknowns eliminated approximately.
Step  3:    117 of the    326 unknowns eliminated approximately.
Step  4:     67 of the    209 unknowns eliminated approximately.
Step  5:     40 of the    142 unknowns eliminated approximately.
Step  6:     28 of the    102 unknowns eliminated approximately.
Step  7:     19 of the     74 unknowns eliminated approximately.
Step  8:     14 of the     55 unknowns eliminated approximately.
Step  9:     12 of the     41 unknowns eliminated approximately.

Warning from incldup!
   After LDU-factorisation of last submatrix:
   Last diagonal element is (nearly) zero!
   Value of this element is changed into 1.73343E-03

Time ILDU factorization of last block:    0.00

Total time Preconditioner:    0.01


Number of non-zeros in the LDU-decomposition:      9084
Number of non-zeros in last block:                  304
Average number of non-zeros per row in L+D+U:         8.34

 Iter.            Residual  Reduction residual   Flops/unknown
     0        9.948607E-02        1.000000E+00           28.44
     1        1.073744E-02        1.079291E-01           65.88
     2        9.562582E-03        9.611981E-02          103.32
     3        2.535833E-03        2.548933E-02          140.76
     4        9.622758E-04        9.672468E-03          178.20
     5        1.072833E-04        1.078375E-03          215.64
     6        6.038161E-06        6.069353E-05          253.09
     7        1.070183E-06        1.075712E-05          290.53
     8        2.594839E-07        2.608243E-06          327.97
Time Solution:    0.00


Absolute error:  2.2572019760941364E-07
