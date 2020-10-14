function [a,b,c]=retta_mq2abc(m,q,v)
% 
% [a,b,c]=retta_mq2abc(m,q,v)
% 
% Converte retta  a forma esplicita y=m*x+q oppure x=v (retta verticale) a 
% forma implicita a*x+b*y+c=0

if isnan(v)
    [a,b,c]=retta_per_2pti(0,q,-q/m,0);
else
    [a,b,c]=retta_per_2pti(v,0,v,1);
end