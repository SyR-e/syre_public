function mat = setMaterialProperties(M)

% SETMATERIALPROP sets default material properties
%
% USE: mat = setMaterialProp(M)
%
% INPUT:
% 'M' vector with material codes
%
% OUTPUT:
% 'mat' is a struct with the folowing fields:
% - .Name: material name (default: '')
% - .Color: color of material (default: lines colormap array)
% electromagnetic properties
% - .EpsR      : relative permittivity (default: 1)
% - .MuR       : relative permeability (default: 1)
% - .Hc        : coercitive field of permanent magnets (default: 0)
% - .HcDir     : direction of magnetization (default: [0 0 0])
% - .Sigma     : electrical conductivity (default: 0 S/m)
% mechanical properties
% - .Young  : Young modulus (default: 0 Pa)
% - .Density: material density (default: 0 kg/m3)
% - .Poisson: Poisson coefficient (default: 0)
% thermal properties
% - .Lambda : thermal conductivity (default: 0 W/(m*K))
% - .SpecificHeat: specific heat (default: 0 J/(kg*K))
% - .Alpha  : thermal expansion (default: 0 1/K)
% 2D objects
% - .Delta  : thickness (m)
% 1D objects
% - .Section: cross section (m^2)
%
% VERSION:
% Date: 19.12.2009
% Copyright(C) 2009-2016: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 29.01.2010: change definition of mat.nMat
% 18.10.2013: added SECTION field to the structure (Luca Giaccone)
% 13.11.2013: default thickness (DELTA) is 1
% 18.01.2016: new names and updated Help

% extraction and sort of unduplicated material codes
matCodes = unique(M(:));

% check if material 1 is present
imat1 = find(matCodes == 1, 1);
if isempty(imat1)
    matCodes = [1; matCodes];
end
nMat = double(max(matCodes));

% initialization structure with default material properties

% name

mat = struct(...
    'Name',repmat({''},nMat,1),...
    'Color',loadColor(nMat,'line'),...
    ... % electromagnetic properties
    'EpsR',repmat({1},nMat,1),...
    'MuR',repmat({1},nMat,1),...
    'Hc',repmat({0},nMat,1),...
    'HcDir',repmat({[0 0 0]},nMat,1),...
    'Sigma',repmat({0},nMat,1),...
    ... % thermal properties
    'Lambda',repmat({0},nMat,1),...
    'SpecificHeat',repmat({0},nMat,1),...
    'Density',repmat({0},nMat,1),...
    'Alpha',repmat({0},nMat,1),...
    ... % mechanical properties
    'Young',repmat({0},nMat,1),...
    'Poisson',repmat({0},nMat,1),...
    ... % thickness
    'Delta',repmat({0},nMat,1),...
    ... % section
    'Section',repmat({0},nMat,1));
