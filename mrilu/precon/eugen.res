
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
 
Method: CG
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

Total time Preconditioner:    0.00


Number of non-zeros in the LDU-decomposition:      1240
Number of non-zeros in last block:                   74
Average number of non-zeros per row in L+D+U:         8.61

 Iter.            Residual  Reduction residual   Flops/unknown
     0        9.427141E-02        1.000000E+00           19.11
     1        1.646947E-02        1.747027E-01           42.70
     2        2.120506E-03        2.249363E-02           66.29
     3        2.693943E-04        2.857646E-03           89.88
     4        3.952377E-05        4.192550E-04          113.47
     5        1.001791E-05        1.062667E-04          137.07
     6        6.788474E-07        7.200988E-06          160.66
Time Solution:    0.00


Absolute error:  4.4954969514604548E-07

     Min. size free segment   =    95988206
     Max. size free segment   =    95988206
     Number of free segments  =           1
     Total free space         =    95988206
     Total number of segments =         110
     Actual buffer size used  =       11796
     Required buffer size    >=       14964
     Top of lower part        =       11790
     Bottom of higher part    =    96000000

