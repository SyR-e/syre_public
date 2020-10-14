function [a,b,c]=retta_tg_circ_punto(xc,yc,r,x0,y0)
% 
% [a,b,c]=retta_tg_circ_punto(xc,yc,r,x0,y0)
% 
% Calcola la retta y=m*x+q tangente alla circonferenza di centro (xc,yc) e
% raggio r nel punto (x0,y0)

% controllo che (x0,y0) sia sulla circonferenza


if (0)
    hfig=figure();
    hold on
    grid on
    axis equal
    fi=0:0.001:2*pi;
    xplot=xc+r*cos(fi);
    yplot=yc+r*sin(fi);
    plot(xplot,yplot,'r');
    plot(x0,y0,'bo');
    keyboard
    close(hfig);
    clear hfig
end



tmp=(x0-xc)^2+(y0-yc)^2-r^2;
if tmp>1e-8
    error('(x0,y0) do not lay on the circumference')
end

% retta passante per centro circonferenza e punto di tangenza
[a0,b0,c0]=retta_per_2pti(xc,yc,x0,y0);

[m0,q0,v0]=retta_abc2mq(a0,b0,c0);

% retta voluta è perpendicolare a retta appena trovata e passa per il punto
% (x0,y0)

if isnan(v0)
    if m0~=0 % retta non orizzontale e non verticale
        m=-1/m0;
        q=y0-m*x0;
        v=NaN;
    else     % retta orizzontale
        m=NaN;
        q=NaN;
        v=x0;
    end
else        % retta verticale
    m=0;
    q=y0;
    v=NaN;
end

[a,b,c]=retta_mq2abc(m,q,v);
        
        

