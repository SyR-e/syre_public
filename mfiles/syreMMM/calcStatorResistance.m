function [Rs] = calcStatorResistance(d,FreqElet)

if strcmp(d.SkinEffModel.type,'0')
    kSkin=1;
else
    kSkin=interp1(d.SkinEffModel.f,d.SkinEffModel.k,FreqElet);
end

Rs=d.Rs.*kSkin.*(1+0.004*(d.temp-d.temp0));