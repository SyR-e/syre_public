function q = triMeshQuality(P,T,varargin)

% TRIMESHQUALITY calculates mesh quality index and plot histogram
%
% USE:
% q = triMeshQuality(P,T)
%
% INPUTS:
% 'P': node coordinates
% 'T': triangle connettivity matrix
% 'varargin': variable inputs
%   * 'plot': plot flag. Available 'y'/'n' (default: 'n')
%
% OUTPUTS:
% 'q': normalized shape ratio index
%
% NOTES:
% See D.A. Field, " Qualitative measures for initial meshes", International
% Journal for numerical methods in engineering, Vol 47, 2000, pp. 887-906
%
% VERSION:
% Date: 24.08.2012
% Copyright(C) 2012-2016: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 21.11.2016: new horizontal bar plot
% 29.05.2018: vertical bar plot
% 17.01.2019: added VARARGIN

% set defaults
plot = 'n';

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'plot')
        plot = varargin{i+1};
    end
end

% circumradius and inradius
[R,r] = triRadius2d(P,T);

% normalized shape ratio index
q = 2*r./R;

% picture
x = (0:.01:1)';
n = histc(q,x);

if strcmpi(plot,'y')
    h = figure;
    bar(x,n);
    figureSetup(h);
    ylabel('No. of elements');
    xlabel('Mesh quality index');
    figureSetup(h);
end

end