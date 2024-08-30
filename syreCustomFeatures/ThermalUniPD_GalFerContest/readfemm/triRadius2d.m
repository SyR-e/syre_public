function [R,r] = triRadius2d(P,T)

% TRIRADIUS2D calculate circumradii and inradii of triangles
% 
% USE:
% [R,r] = triRadius2d(P,T)
% 
% INPUT:
% 'P': node coordinates
% 'T': triangle connettivity matrix
% 
% OUTPUT:
% 'R': vector of circumradii
% 'r': vector of inradii
%
% NOTE:
% See http://en.wikipedia.org/wiki/Triangle
%
% VERSION:
% Date: 24.08.2012
% Copyright(C) 2012: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
%

% check input
if size(T,2) ~= 3
    error('Only triangles allowed');
end

% vectors from first node node
a = sqrt(sum((P(T(:,2),:)-P(T(:,1),:)).^2,2));
b = sqrt(sum((P(T(:,3),:)-P(T(:,1),:)).^2,2));
c = sqrt(sum((P(T(:,3),:)-P(T(:,2),:)).^2,2));

a2b2c2 = a.^2.*b.^2.*c.^2;
abc = a+b+c;
abc3 = (-a+b+c).*(a-b+c).*(a+b-c);

R = sqrt(a2b2c2./(abc.*abc3));
r = sqrt(abc3./(4*abc));