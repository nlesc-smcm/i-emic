% COPYRIGHT ERIK

function [svp1, svp2] = checkerboard_mode(nx, ny, nz)
  svp1 = zeros(nx, ny, nz);
  svp2 = zeros(nx, ny, nz);
  for k = 1:nz
	for j = 1:ny
	  for i = 1:nx
		if (mod(i+j-2,2) ~= 0)
		  svp1(i,j,k) = 1;
		else
		  svp2(i,j,k) = 1;
		end
	  end
	end
  end
end
