function colorscale = loadColor(n,method)

% LOADCOLOR set a colormap on the basis of the number passed as argument
%
% USE:
% colorscale = loadColor(n,method)
%
% INPUT:
% 'n': number of required colors
% 'method': 'color'/'grey' (optional, default 'color')
%
% OUTPUT:
% 'colorscale': RGB color scale
%
% VERSION:
% Date: 19.12.2009
% Copyright(C) 2009-2019: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 25.01.2019: added new default colormap based on command LINE

if nargin == 1
    method = 'line';
end

if strcmpi(method,'line')
    base = [1 1 1;
        lines(7)];
elseif strcmpi(method,'color')
    base = [1 1 1   %1 White
        1 0 0       %2 Red
        0 1 0       %3 Green
        0 0 1       %4 Blue
        1 1 0       %5 Yellow
        1 0 1       %6 Magenta
        0 1 1];     %7 Cyan
elseif strcmpi(method,'grey')
    base = [1 1 1];
end

mult = floor(double(n)/size(base,1))+1;
scale = linspace(1,1/mult,mult);
nc = size(base,1);
nb = length(scale);

scale = repmat(reshape(repmat(scale,nc,1),nb*nc,1),1,3);
colorscale = repmat(base,nb,1).*scale;