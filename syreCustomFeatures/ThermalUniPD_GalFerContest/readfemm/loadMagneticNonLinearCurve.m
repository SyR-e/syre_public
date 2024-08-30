function mat = loadMagneticNonLinearCurve(mat,obj,H,B)

% CREATE_NLMAGCURVE creates equivalent nonlinear curve for time harmonic
% nonlinear problems
%
% USE:
% mat = create_nlmagcurve(H,B,obj,primal,mat,varargin)
%
% INPUTS:
% 'mat': struct of the material properties, as defined by function SETMATERIALPROP
% 'obj': object code (scalar)
% 'H', 'B': vectors describing H-B curve
%
% OUTPUTS:
% 'mat' this struct is added by the field 'nl' and 'inl' is switched to 1
% 'nl' has several subfields
%   * 'H', 'B': magnetic characteristic
%   * 'idx': indices of nonlinear cells
%   * 'NuDiff': differential nu for computation of Jacobian
%   * 'MuR': optimal mu_r for fixed-point technique
%
% NOTES:
%
% VERSION:
% Date: 19.12.2009
% Copyright(C) 2009-2016: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 19.01.2011: added METHOD as optional input
% 08.02.2011: added possibility of handling thermal dependance of BH
%             characteristic
% 11.06.2011: bug fix when thermal dependance is requested
% 14.06.2012: H,B forced to be column vectors
% 29.11.2012: uses function MU0
% 22.01.2016: new refurbished routine

for i = 1:length(obj)
    % load curve
    mat(obj(i)).MuR = struct('MuR',[],'B',B(:),'H',H(:),'NuDiff',zeros(size(H(:))));

    % NuDiff for Jacobian computation
    mat(obj(i)).MuR.NuDiff(2:end) = diff(mat(obj(i)).MuR.H)./diff(mat(obj(i)).MuR.B);
    mat(obj(i)).MuR.NuDiff(1) =mat(obj(i)).MuR.NuDiff(2);
    
    % optimal MuR for fixed-point technique
    mat(obj(i)).MuR.MuR = 2/(Mu0*(min(mat(obj(i)).MuR.NuDiff)+max(mat(obj(i)).MuR.NuDiff)));
end

% for i = 1:length(obj)
%     % load curve
%     mat.NonLinear(obj(i)).B = B(:);
%     mat.NonLinear(obj(i)).H = H(:);
%     
%     % set nonlinear material variable
%     mat.iNonLinear(obj(i)) = 1;
%     
%     % NuDiff for Jacobian computation
%     mat.NonLinear(obj(i)).NuDiff = ones(size(mat.NonLinear(obj(i)).H));
%     mat.NonLinear(obj(i)).NuDiff(2:end) = diff(mat.NonLinear(obj(i)).H)./diff(mat.NonLinear(obj(i)).B);
%     mat.NonLinear(obj(i)).NuDiff(1) = mat.NonLinear(obj(i)).NuDiff(2);
%     
%     % optimal MuR for fixed-point technique
%     MuR = 2/(Mu0*(min(mat.NonLinear(obj(i)).NuDiff)+max(mat.NonLinear(obj(i)).NuDiff)));
%     % update MuR parameter
%     mat.MuR(obj(i)) = MuR;
% end

end