
Enter name input file:
 
Construct Convec-Diff Eqn
Nr Intvals =    16
Stretching = 1.00E+03,  Coeff U_x = 0.00E+00,  Coeff U_y = 0.00E+00
 
Scale rowsums of Matrix |A| to 1.
Use Exact elimination.
GusFctr = 1.00E+00
CLSOnce =  T,  NLSFctr = 0.100,
EpsW = 0.050,  ElmFctr = 0.150,  RedFctr = 1.000,  SchTol = 0.00E+00,
DensLim = 0.200,  GlobFrac = 0.000,  LocFrac = 0.000,  SparsLim = 0.666,

ILUT factorisation of last sparse block.
CompFctr = 1.00E+00  ,  DropTol = 5.00E-03  ,  CPivTol = 1.00E+00
LUTol = 1.00E-10,  SingLU = T
 
Method: Bi-CGSTAB
Max Nr Iters =   40,  Residual <= 1.00E-06,  Decrease residual <= 1.00E-06

Number of equations:        289
Default block size:           1
Number of blocks:           289

Step  0:    145 of the    289 unknowns eliminated exactly.
Step  1:     53 of the    144 unknowns eliminated approximately.
Step  2:     31 of the     91 unknowns eliminated approximately.
Step  3:     18 of the     60 unknowns eliminated approximately.

Warning from incldup!
   After LDU-factorisation of last submatrix:
   Last diagonal element is (nearly) zero!
   Value of this element is changed into 3.39428E-01

Time ILDU factorization of last block:    0.00

Total time Preconditioner:    0.00


Number of non-zeros in the LDU-decomposition:      1347
Number of non-zeros in last block:                  470
Average number of non-zeros per row in L+D+U:         9.35

 Iter.            Residual  Reduction residual   Flops/unknown
     0        6.080789E+00        1.000000E+00           20.35
     1        2.156905E-01        3.547081E-02           67.02
     2        9.798967E-03        1.611463E-03          115.68
     3        2.033584E-03        3.344276E-04          164.35
     4        1.993114E-04        3.277723E-05          213.01
     5        1.016843E-05        1.672222E-06          261.67
     6        2.253781E-07        3.706395E-08          310.34
Time Solution:    0.00


Absolute error:  3.3788634230813575E-08
