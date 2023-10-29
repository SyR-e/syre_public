function [r] = slidingBCr(location,state)

xy = location.x+j*location.y;
% theta = angle(xy);
u = state.u;

theta = reshape(angle(xy),[1 1 numel(xy)]);

r = zeros(2,length(theta));
