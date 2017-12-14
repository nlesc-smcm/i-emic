function [salfun] = salfun(xx, yy, ymin, ymax)

%****************************************************************************
%real FUNCTION salfun(xx,yy)
%  use m_par
%  use m_global ! we need the global value of ymin and ymax here!
%  implicit none
%  real    xx,yy
%  real    subtr    
    
    if (ymin >= 0.0) % Northern hemisphere
        salfun = cos(pi*(yy-ymin)/(ymax-ymin))
        % salfun=0.0
    else
        salfun = cos(pi*yy/ymax)
    end
    
end
