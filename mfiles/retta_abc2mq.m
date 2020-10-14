function [m,q,v]=retta_abc2mq(a,b,c)
% 
% [m,q,v]=retta_abc2mq(a,b,c)
% 
% Converte retta da forma implicita a*x+b*y+c=0 a forma esplicita y=m*x+q
% oppure x=v (retta verticale)

if b~=0
    m=-a/b;
    q=-c/b;
    v=NaN;
else
    v=-c/a;
    m=NaN;
    q=NaN;
end