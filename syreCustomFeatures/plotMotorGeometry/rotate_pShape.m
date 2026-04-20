function [ps] = rotate_pShape(ps,angle)

[xTmp,yTmp] = rot_point(ps.Vertices(:,1),ps.Vertices(:,2),angle);

ps.Vertices = [xTmp,yTmp];