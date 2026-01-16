function GAM = momentxx(data,p)

k = size(data,2);
T = size(data,1);
T = T-p;
Y = data';
Z = zeros(k*p+1,T);

for i = 1:T;
    if p ==1;
        Z(:,i) = [1;Y(:,i+p-1)];
    elseif p == 2;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2)];
    elseif p == 3;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3)];
    elseif p == 4;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4)];
    elseif p == 5;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5)];
    elseif p == 6;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6)];
    elseif p == 7;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7)];
    elseif p == 8;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8)];
    elseif p == 9;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9)];
    elseif p == 10;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10)];
    elseif p == 11;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);Y(:,i+p-11)];
    elseif p == 12;
        Z(:,i) = [1;Y(:,i+p-1);Y(:,i+p-2);Y(:,i+p-3);Y(:,i+p-4);Y(:,i+p-5);Y(:,i+p-6);Y(:,i+p-7);Y(:,i+p-8);Y(:,i+p-9);Y(:,i+p-10);Y(:,i+p-11);Y(:,i+p-12)];
    end;
end;

GAM = Z*Z';