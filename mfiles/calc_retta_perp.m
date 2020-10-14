function [a2,b2,c2]=calc_retta_perp(a1,b1,c1,x0,y0)
% 
% [a2,b2,c2]=calc_retta_perp(a1,b1,c1,x0,y0)
% 
% Calcola retta a2*x+b2*y+c2=0 perpendicolare a retta a1*x+b1*y+c1 e
% passante per (x0,y0)

[m1,q1,v1]=retta_abc2mq(a1,b1,c1);

if isnan(v1)
    