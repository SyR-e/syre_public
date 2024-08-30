function [Twind_mean, Twind_min, Twind_max] = solve_therm(geo,per,mat,filename)
%addpath("loc_gcd\")

[geo,mat,opt,T,P,M] = preprocess_therm(geo,per,mat,filename);
[Twind_mean,Twind_min,Twind_max] = funTerm_v1(geo,mat,opt,T,P,M);
end