function [y_cma] = cma_q(y)
%%
% CMA_Q.M
% By Gustavo Leyva
% Contact: leyva004@umn.edu or gustavo.leyva.jimenez@gmail.com
% Department of Economics
% University of Minnesota
% Created on 03.19.13
% Last-modified on 12.15.14
%
% PURPOSE   Computes the centered-moving average of a provided raw quarterly
%           series y
% USAGE     function [y_cma] = cma_q(y)
% INPUTS    y      : raw quarterly series
% OUTPUTS	y_cma  : centered-moving average of y	
%
n = length(y);
z0 = [];

j = 1;
while j <= n-3;  
    z0 = [z0;( mean( packrv([y(j);y(j+1);y(j+2);y(j+3)]) ) )];
    j = j+1;
end;

y_cma = [];
nz = length(z0);
j = 1;
while j <= nz-1;  
    y_cma = [y_cma;( mean( packrv([z0(j);z0(j+1)]) ) )];
    j = j+1;
end;

end