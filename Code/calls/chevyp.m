function [P] = chevyp(k,n,t)

P = sqrt(2).*cos( k.*pi.*(t-0.5)./n );

end