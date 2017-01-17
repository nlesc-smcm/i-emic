
Enter name input file:
 
Eqn from binary files (CSR)
 
Use Reverse Cuthill-McKee reordering.
Use Exact elimination.
Use Gustafsson in Multi Level Preconditioner.
CLSOnce =  F,  NLSFctr = 1.000,
EpsW = 0.800,  ElmFctr = 0.200,  RedFctr = 0.800,  SchTol = 0.00E+00,
DensLim = 0.100,  GlobFrac = 0.000,  LocFrac = 0.000,  SparsLim = 0.675,

ILUT factorisation of last sparse block.
CompFctr = 0.00E+00  ,  DropTol = 2.50E-03  ,  CPivTol = 0.00E+00
LUTol = 1.00E-10,  SingLU = F
 
Method: Bi-CGSTAB
Max Nr Iters =   40,  Residual <= 1.00E-06,  Decrease residual <= 1.00E-06

Enter name data file:

Number of equations:        675
Default block size:           3
Number of blocks:           225


Warning from  BlkCutMcK:  Matrix is reducible!
           3 of the      675 rows/columns in the first level set.

Step  0:    192 of the    675 unknowns eliminated exactly.
Step  1:    156 of the    483 unknowns eliminated approximately.
Step  2:     60 of the    327 unknowns eliminated approximately.
Time ILDU factorization of last block:    0.01

Total time Preconditioner:    0.03


Number of non-zeros in the LDU-decomposition:     17884
Number of non-zeros in last block:                14680
Average number of non-zeros per row in L+D+U:        37.03

 Iter.            Residual  Reduction residual   Flops/unknown
     0        6.420336E-01        1.000000E+00           71.83
     1        1.450601E-02        2.259385E-02          224.07
     2        6.955845E-04        1.083408E-03          379.17
     3        1.453348E-05        2.263664E-05          534.27
     4        1.662626E-06        2.589625E-06          689.37
     5        1.889399E-07        2.942836E-07          844.48
Time Solution:    0.01


Number of non-zeros in the LDU-decomposition:     17884
Number of non-zeros in last block:                14680
Average number of non-zeros per row in L+D+U:        37.03

 Iter.            Residual  Reduction residual   Flops/unknown
     0        1.889399E-07        1.000000E+00           71.83
Time Solution:    0.00


     Min. size free segment   =        1350
     Max. size free segment   =    47921534
     Number of free segments  =           2
     Total free space         =    47922884
     Total number of segments =          55
     Actual buffer size used  =       78468
     Required buffer size    >=      102658
     Top of lower part        =       78462
     Bottom of higher part    =    48000000

