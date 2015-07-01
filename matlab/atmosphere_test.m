load_atmosphere

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

xmindeg_ = 286;
xmaxdeg_ = 350;
ymindeg_ = 10 ;
ymaxdeg_ = 74 ;

xmin_ = xmindeg_ * pi / 180;
xmax_ = xmaxdeg_ * pi / 180;
ymin_ = ymindeg_ * pi / 180;
ymax_ = ymaxdeg_ * pi / 180;

dy    =  (ymaxdeg_-ymindeg_)/(16);
yc    =  ymindeg_+dy/2:dy:ymaxdeg_;
dx    =  (xmaxdeg_-xmindeg_)/(16);
xc    =  xmindeg_+dx/2:dx:xmaxdeg_;
sun   =  suna_(yc);
ampl  =  0.0179624;
frc2  =  ampl*(sun - amua_);

To = reshape(t0_+otemp,16,16);
Ta = reshape(t0_+state,16,16);
FRC = reshape(frc,16,16);
Ta_ocean = reshape(t0_+tatm,16,16);

figure(2)
contourf((xc),(yc),To',15);
colorbar
title('T^o')
figure(3)
contourf((xc),(yc),Ta',15);
colorbar
title('T^a')
figure(4)
contourf((xc),(yc),FRC',15);
colorbar
title('frc')
figure(5)
contourf((xc),(yc),Ta_ocean',15);
colorbar
title('T^a_o')



%## figure(4)
%## imagesc(repmat(frc2,16,1));
%## colorbar
%## title('computed frc')

