function [fdfq] = interpFluxMapsTemperature(fdfq1,fdfq2,temp1,temp2,temp)

fdfq.Id = fdfq1.Id;
fdfq.Iq = fdfq1.Iq;
fdfq.Fd = nan(size(fdfq1.Id));
fdfq.Fq = nan(size(fdfq1.Id));
fdfq.T  = nan(size(fdfq1.Id));

for ii=1:numel(fdfq.Id)
    fdfq.Fd(ii) = interp1([temp1 temp2],[fdfq1.Fd(ii) fdfq2.Fd(ii)],temp);
    fdfq.Fq(ii) = interp1([temp1 temp2],[fdfq1.Fq(ii) fdfq2.Fq(ii)],temp);
    fdfq.T(ii)  = interp1([temp1 temp2],[fdfq1.T(ii) fdfq2.T(ii)],temp);
end

fdfq.dT   = nan(size(fdfq.Id));
fdfq.dTpp = nan(size(fdfq.Id));