
Enter name input file:
 
Construct Convec-Diff Eqn
Nr Intvals =     5
Stretching = 1.00E+03,  Coeff U_x = 0.00E+00,  Coeff U_y = 0.00E+00
 
Use Exact elimination.
GusFctr = 1.00E+00
CLSOnce =  T,  NLSFctr = 0.100,
EpsW = 0.500,  ElmFctr = 0.150,  RedFctr = 1.000,  SchTol = 0.00E+00,
DensLim = 0.200,  GlobFrac = 0.000,  LocFrac = 0.000,  SparsLim = 0.666,

ILU(2) factorisation of last sparse block.
CompFctr = 1.00E+00
LUTol = 1.00E-10,  SingLU = T
 
Method: CG
Max Nr Iters =   40,  Residual <= 1.00E-06,  Decrease residual <= 1.00E-06

Number of equations:         36
Default block size:           2
Number of blocks:            18

Step  0:     18 of the     36 unknowns eliminated exactly.

Warning from ilduk!
   After LDU-factorisation of last submatrix:
   Last diagonal element is (nearly) zero!
   Value of this element is changed into 1.02468E+00

Time ILDU factorization of last block:    0.00

Total time Preconditioner:    0.00


Number of non-zeros in the LDU-decomposition:       202
Number of non-zeros in last block:                  202
Average number of non-zeros per row in L+D+U:        11.22

 Iter.            Residual  Reduction residual   Flops/unknown
     0        1.396209E-02        1.000000E+00           20.89
     1        5.467583E-15        3.916021E-13           46.28
Time Solution:    0.00


Absolute error:  6.3796190552523058E-14
