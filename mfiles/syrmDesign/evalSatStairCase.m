function [level]=evalSatStairCase(x,y,alphad)
% 
% [level]=evalSatSine(alphad,b,kt)
% 
% Evaluation of the d-flux with a saturated cosine, useful for the evaluation
% of the flux carrier

debug=0;

alphad=alphad(1:end-1)*pi/180; % elimino dal vettore alphad l'ultimo valore (90°) e converto in rad

% uso solo la forma d'onda tra 0 e pi/2
% y=y(x>=0);
% x=x(x>=0);
% y=y(x<=pi/2);
% x=x(x<=pi/2);

if debug
    figure()
    hold on
    grid on
    box on
    xlim([0 90])
    ylim([0 1])
    plot(x*180/pi,y,'-r','LineWidth',1,'DisplayName','Saturated')
    plot(alphad(2:end)*180/pi,interp1(x,y,alphad(2:end)),'go','DisplayName','\alpha_d')
    legend('show','Location','SouthWest')
end

level=zeros(1,length(alphad)-1); % in alphad ho nlay+2 valori
for ii=1:length(level)
    level(ii)=trapz(x(x<=alphad(ii+1)),y(x<=alphad(ii+1)))-trapz(x(x<=alphad(ii)),y(x<=alphad(ii)));
end

level=level/(trapz(x(x<alphad(end)),y(x<alphad(end))));
%level=level/trapz(x,y);
