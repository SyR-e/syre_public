function [c1,Href] = MagHc_temp_fitting(mat)
temp=mat.LayerMag.temp.temp';
BrT=mat.LayerMag.temp.Br';
mur=mat.LayerMag.mu;
HcT=BrT/(mur*4*pi*10^-7);
k=temp==20;
Br20=BrT(k);
Href = HcT(k);
c1=(BrT(end)/Br20-1)/(temp(end)-20);


% figure()
% plot(temp,HcT,'rx',temp,Href*(1+c1*(temp-20)))
end

% ft = fittype(@(c1,temp) (Href*(1+c1*(temp-20)+c2(temp-20)^2)),'independent', {'temp'});
% fo = fitoptions( 'Method','NonlinearLeastSquares', 'StartPoint',[0 0]);
% fitresult=fit(temp,HcT,ft,fo);
% 
% coeff = coeffvalues(fitresult);
% c1=coeff(1);
% c2=coeff(2);


