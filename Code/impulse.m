function [impres] = impulse(B,a,s,p)
% B = VAR estimated coefficients in format 'forbet' - see main file
% a = impulse vector - Uhlig(1998)
% s = horizon of the impulse response function 
% p = VAR lag order

% References:
% Uhlig(1998): The robustness of identified VAR conclusion about money: A
% comment. Carnegie-Rochester Conference Series on Public Policy 49
% 245-263. North Holland.

k = size(B,1);
B = B(:,2:size(B,2));                                                      % reshaping B - without constants

if p == 1;
    gama = B; 
    aa = a;
else
    ident = [kron(eye(k),eye(p-1)) zeros(k*(p-1),k)];
    gama = vertcat(B,ident);
    aa = [a' zeros(1,k*(p-1))]';
end;

impres = zeros(s,k);
power = eye(length(aa));
for i = 1:s;
    aux = power*aa;
    impres(i,:) = aux(1:k)';
    power = power*gama;
end;