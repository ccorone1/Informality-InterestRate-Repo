function [y_cma] = cma_m(y)
%%
% CMA_M.M
% Adapted for Monthly Data
%
% PURPOSE   Computes the centered-moving average of a provided raw monthly
%           series y (Using 2x12 Trend Cycle methodology)
% USAGE     function [y_cma] = cma_m(y)
% INPUTS    y      : raw monthly series
% OUTPUTS   y_cma  : centered-moving average of y
%
n = length(y);

% Step 1: 12-term Moving Average
% This results in a series centered at t + 6.5
z0 = [];
j = 1;
while j <= n-11;
    % Average of 12 months
    % Note: packrv isn't standard matlab, replaced with standard matrix construction
    % If packrv is essential for your specific environment, you can revert to:
    % packrv([y(j); ... y(j+11)])
    window_slice = y(j:j+11); 
    z0 = [z0; mean(window_slice)];
    j = j+1;
end;

% Step 2: 2-term Moving Average of the result
% This centers the series at t + 7 (Index 7 corresponds to the 7th month)
y_cma = [];
nz = length(z0);
j = 1;
while j <= nz-1;
    y_cma = [y_cma; mean([z0(j); z0(j+1)])];
    j = j+1;
end;

end