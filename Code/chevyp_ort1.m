function [P] = chevyp_ort1(n,t)

P = (1./sqrt((n^2-1)/12)).*(t-(n+1)/2);

end