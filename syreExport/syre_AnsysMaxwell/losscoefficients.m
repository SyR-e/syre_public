function [Kh,Kadd,Ke]=losscoefficients(alpha,beta,kh,ke,density)
warning off;
%[W/kg]->%[m/m^3]
kh=kh*density; 
Ke=ke*density;
B = [0:0.1:2];
f=round(logspace(log10(50),4,25));

%perdite isteresi calcolate con modello Syr-e
LossH= kh*f'.^alpha.*B.^beta;

for i=4:length(LossH(:,1))
    LossH(i,end-(i-9):end)=NaN(1,i-8);
end
    

[xData, yData, zData] = prepareSurfaceData( f,B, LossH);
ft = fittype( 'kadd*(x*y)^1.5+kh*x*y^2', 'independent', {'x', 'y'}, 'dependent', 'z' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';opts.Lower = [0 0];opts.Robust = 'LAR';opts.StartPoint = [0,kh];

fitresult = fit( [xData, yData], zData, ft, opts );

coeff = coeffvalues(fitresult);
Kadd=coeff(1);
Kh=coeff(2);

% %PLOT RESULT
% figure
% hold on
% for i=1:length(f)
%     LossHSyre= kh*f(i).^alpha*B.^beta;
%     LossHAnsys=Kh*f(i)*B.^2+Kadd*(f(i)*B).^1.5;
%     if i>3
%         LossHAnsys(end-(i-9):end)=NaN(1,i-8);
%         LossHSyre(end-(i-9):end)=NaN(1,i-8);
%     end
%     plot(B,LossHSyre,'rx',B,LossHAnsys,'bo')
% end
% xlabel('B[T]');
% ylabel('Losses[W/m^3]');
% tit=sprintf('Loss models: "x" - Syre ; "o" - Ansys\n f:[%s]',num2str(f));
% title(tit);
% legend('Syr-e','Ansys')
% hold off

warning on;
end



