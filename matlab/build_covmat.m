function [Q,M,V,D] = build_covmat(range, Q, mode)
  fprintf(1,'----------------------------------------------\n')
  %% - DEFINE CONSTANTS - ----------------------------------------------
  udim  = 0.1;                 %[m/s]   Velocity scale
  r0dim = 6.4e6;               %[m]     Radius of Earth
  T0    = 15;                  %[deg C] Reference temperature
  S0    = 0;                  %[psu]   Reference salinity
  RtD   = 180/pi;              %[-]     Radians to degrees
  %% - READ MASK - -----------------------------------------------------
  [n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
  surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
  landm_int  = landm(2:n+1,2:m+1,2:l+1);
  dx         = (xu(n+1)-xu(1))/n;
  dy         = (yv(m+1)-yv(1))/m;
  dz         = (zw(l+1)-zw(1))/l;

  [qz,dfzt,dfzw] = gridstretch(zw);
  %% - READ SOLUTION --------------
  dim = n*m*(l+la)*nun;
  %=======================================================================
  

  obs = numel(range);
  
  if nargin < 2 || isempty(Q)
	% data matrix
	M = read_many_forts(range, dim);

	% remove checkerboard modes in pressure
	p_idx   = 4:6:dim; % pressure indices
	[s1,s2] = checkerboard_modes(n,m,l);

	scale1   = s1(:)'*s1(:);
	scale2   = s2(:)'*s2(:);
	for i = 1:obs
	  M(i,p_idx) = M(i, p_idx) - (M(i, p_idx)*s1(:)/scale1)*s1(:)' ...
	                 - (M(i, p_idx)*s2(:)/scale2)*s2(:)';
	end
	
	% subtract mean to get centered data matrix
	mn = mean(M);
	M = M - repmat(mn,obs,1);
	Q = (1/(obs-1))*M'*M;
  else
	M = [];
  end
  sol = zeros(nun, n, m, l);
  
  if nargin < 3
	mode = 1;
  end
  
  fprintf(1,'------------- Calculating eigenvector--------\n')
  opts.issym = true;
  [V,D,flag] = eigs(Q,20,'lm',opts);  

  %plot_mode(V,mode)
end
