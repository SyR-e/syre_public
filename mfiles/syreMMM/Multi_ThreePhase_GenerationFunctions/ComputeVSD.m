function [T_VSD] = ComputeVSD(n_set)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% ------------------------Define parameters------------------------------%

l = n_set;      % number of winding sets
k = 3;      % number of phases per each winding set
n = k*l;    % total number of phases      

%% ---------------Define theta_n MATRIX with phases position--------------%

theta_n = zeros(1,n);  
% sequence is a1 b1 c1 a2 b2 c2 ...
for j=1:l
    for i=1:k
        theta_n(i+(j-1)*k) = (180/n)*(2*l*(i-1)+j-1);
    end
end

%% -----------------Define non-zero sequence subspace constants-----------%

% generate c vector with odd numbers <n
C0 =1; x=1;
while(C0<n)
    c(x) = C0;
    x=x+1;
    C0=C0+2;
end

% find zero sequence harmonics constants C_zs
x=1; i=1;
while(i*k<n)
    if(rem(i*k,2)~=0) 
        C_zs(x) = i*k;
        x=x+1;
    end 
    i=i+1;
end

% Find non-zero sequence harmonics constants C_nzs, removing the C_zs
% constants from the c vector

[tmp,~] = ismember(c,C_zs);
index_tmp = find(tmp==1);
C_nzs = c;   
C_nzs(index_tmp)=[];

%% Compute T VSD Matrix
 
sigma = 2/n;            % amplitude invariant transformation
VSD = zeros(n,n);     % initialize matrix 
row=1;

for x=1:length(C_nzs)
    c = C_nzs(x);
    VSD(row,:) = cosd(c*theta_n);
    VSD(row+1,:) = sind(c*theta_n);
    row=row+2;
end

% Separate neutral points zero sequence subspaces
k=1;
for x=row:n
    VSD(x,k) = 1;
    VSD(x,k+1) = 1;
    VSD(x,k+2) = 1;
    k=k+3;
end
 
T_VSD =sigma*VSD;

end