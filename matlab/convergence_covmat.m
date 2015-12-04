function [Q,M,V,D] = convergence_covmat(max, neigs)
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
  range_b = 1;
  range_e = 8;
  M = zeros(max,dim);
  obs = 0;
  fileID = fopen('eigenvalues.txt','w')
  while range_e <= max
	fprintf('adding observations %d to %d\n', range_b, range_e);
	range = 999 + (range_b:range_e);

	range_b = range_e + 1;
	range_e = range_e*2;

	obs_old  = obs;
	obs      = obs + numel(range);
	m_range  = 1:obs;
	m_append = obs_old+1:obs;
	M(m_append,:) = read_many_forts(range, dim);
	
	% remove checkerboard modes in pressure
	p_idx   = 4:6:dim; % pressure indices
	[s1,s2] = checkerboard_modes(n,m,l);

	scale1   = s1(:)'*s1(:);
	scale2   = s2(:)'*s2(:);
	for i = m_append
	  M(i,p_idx) = M(i, p_idx) - (M(i, p_idx)*s1(:)/scale1)*s1(:)' ...
	               - (M(i, p_idx)*s2(:)/scale2)*s2(:)';
	end
	
	% subtract mean to get centered data matrix
	fprintf(1,'subtract mean\n');
	mn = mean(M(m_range,:));
	Msubtr = M(m_range,:) - repmat(mn,obs,1);
	
	fprintf(1,'compute Q\n');
	Q = (1/(obs-1))*Msubtr'*Msubtr;
	total = trace(Q);
	
	fprintf(1,'compute eigenvectors.. \n');
	opts.tol = 1e-12;
	[V,D,flag] = eigs(Q,neigs,'lm',opts);
	fprintf(1,'compute eigenvectors.. done: %d\n', flag);
	if flag ~= 0
	   fprintf('NOT ALL EIGENVALUES CONVERGED.. EXITING\n');
	   break;
	end
	% scale with trace Q
	EVALS  = diag(D)/total;

	sumeigs = sum(diag(D));
	EVALS2  = diag(D) / sumeigs;

	fprintf(fileID, 'obs = %d, trace = %f, sum diag = %f\n', ...
			obs, total, sumeigs);
	for i = 1:neigs
	  fprintf(fileID, '%8.5f  %8.5f\n', EVALS(i), EVALS2(i));
	end
  end
  fclose(fileID);
  %fprintf(1,'save Q\n');
  %save(['M',num2str(max),'.mat'],'M');
end
