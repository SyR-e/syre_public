function [x,y]=calc_int_retta_circ_gen(xc,yc,r,a,b,c)
% 
% [x,y]=calc_int_retta_circ_gen(xc,yc,r,a,b,c)
% 
% Calcola intersezioni tra la circonferenza di centro (xc,yc) e raggio r e
% la retta con equazione a*x+b*y+c=0
% Se non ci sono intersezioni, i risultati hanno anche parte immaginaria

A=1+a^2/b^2;
B=-2*xc+2*a*c/b^2+2*a*yc/b;
C=xc^2+c^2/b^2+yc^2+2*c*yc/b-r^2;

x=roots([A,B,C]);
y=-(a*x+c)/b;