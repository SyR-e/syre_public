function [h] = slidingBCh(location,state)

xy = location.x+j*location.y;
% theta = angle(xy);
u = state.u;

theta = reshape(angle(xy),[1 1 numel(xy)]);

h = zeros(2,2,length(theta));
h(1,1,:) = +sin(theta);
h(1,2,:) = -cos(theta);
h(2,1,:) = +sin(theta);
h(2,2,:) = -cos(theta);

