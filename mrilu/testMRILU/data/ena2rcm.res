
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

Number of equations:       3267
Default block size:           3
Number of blocks:          1089


Warning from  BlkCutMcK:  Matrix is reducible!
           3 of the     3267 rows/columns in the first level set.

Step  0:   1635 of the   3267 unknowns eliminated exactly.
Step  1:    432 of the   1632 unknowns eliminated approximately.
Step  2:    201 of the   1200 unknowns eliminated approximately.
Step  3:    147 of the    999 unknowns eliminated approximately.
Step  4:    117 of the    852 unknowns eliminated approximately.
Time ILDU factorization of last block:    0.12

Total time Preconditioner:    0.27


Number of non-zeros in the LDU-decomposition:    116840
Number of non-zeros in last block:                87811
Average number of non-zeros per row in L+D+U:        71.59

 Iter.            Residual  Reduction residual   Flops/unknown
     0        1.815959E+00        1.000000E+00           86.62
     1        1.413188E-01        7.782047E-02          265.87
     2        2.082683E-02        1.146878E-02          447.11
     3        3.124940E-03        1.720821E-03          628.35
     4        4.771553E-03        2.627567E-03          809.59
     5        4.081478E-06        2.247561E-06          990.83
     6        5.377644E-07        2.961325E-07         1172.08
Time Solution:    0.10


Number of non-zeros in the LDU-decomposition:    116840
Number of non-zeros in last block:                87811
Average number of non-zeros per row in L+D+U:        71.59

 Iter.            Residual  Reduction residual   Flops/unknown
     0        5.377644E-07        1.000000E+00           86.62
Time Solution:    0.00


     Min. size free segment   =        6534
     Max. size free segment   =    47550524
     Number of free segments  =           2
     Total free space         =    47557058
     Total number of segments =          77
     Actual buffer size used  =      449478
     Required buffer size    >=      631632
     Top of lower part        =      449472
     Bottom of higher part    =    48000000

