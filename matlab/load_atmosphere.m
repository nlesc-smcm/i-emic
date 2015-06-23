beg  = importdata('atmos_beg.txt');
ico  = importdata('atmos_ico.txt');
jco  = importdata('atmos_jco.txt');

n   = numel(beg)-1;
nnz = beg(end)-1;

ivals = zeros(nnz,1);
jvals = jco;
vals  = ico;

row = 1;
idx = 1;
while row <= n
  for k = beg(row):beg(row+1)-1
	ivals(idx) = row;
	idx        = idx + 1;
  end
  row = row + 1;
end
A = sparse(ivals, jvals, vals, n, n);

rhs = importdata('atmos_rhs.txt');
frc = importdata('atmos_frc.txt');


xmin_ = 286 * pi / 180;
xmax_ = 350 * pi / 180;
ymin_ = 10  * pi / 180;
ymax_ = 74  * pi / 180;


% Constants
rhoa_   = 1.25       ; %//! atmospheric density 
hdima_  = 8400.      ; %//! atmospheric scale height 
cpa_    = 1000.      ; %//! heat capacity 
d0_     = 3.1e+06    ; %//! constant eddy diffusivity 
arad_   = 216.0      ; %//! radiative flux param A
brad_   = 1.5        ; %//! radiative flux param B
sun0_   = 1360.      ; %//! solar constant 
c0_     = 0.43       ; %//! atmospheric absorption coefficient
ce_     = 1.3e-03    ; %//! exchange coefficient 
ch_     = 0.94 * ce_ ; %//! exchange coefficient 
uw_     = 8.5        ; %//! mean atmospheric surface wind speed
t0_     = 15.0       ; %//! reference temperature
udim_   = 0.1e+00    ; %//! typical horizontal velocity of the ocean
r0dim_  = 6.37e+06   ; %//! radius of the earth


%// Filling the coefficients
muoa_ =  rhoa_ * ch_ * cpa_ * uw_;
amua_ = (arad_ + brad_ * t0_) / muoa_;
bmua_ =  brad_ / muoa_;
Ai_   =  rhoa_ * hdima_ * cpa_ * udim_ / (r0dim_ * muoa_);
Ad_   =  rhoa_ * hdima_ * cpa_ * d0_ / (muoa_ * r0dim_ * r0dim_);
As_   =  sun0_ * (1 - c0_) / (4 * muoa_);		

suna_ = @(y) (As_*(1 - .482 * (3 * (sin(y)).^2 - 1.) / 2.) *	(1 - .3));
