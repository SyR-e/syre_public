function [Mat]=checkPlotMatrix(Mat,tol)
%
% [Mat]=checkPlotMatrix(Mat,tol)
% 
% Check rotor and stator geometry matrix to avoid the error that plot the
% entire circumference either a short arc.
% Two points are coincident if the distance between them is equal or
% smaller than tol

for ii=1:length(Mat(:,1))
    if Mat(ii,7)==0  %line
        d=sqrt((Mat(ii,1)-Mat(ii,3))^2+(Mat(ii,2)-Mat(ii,4))^2);
        if d<=tol
            Mat(ii,3)=Mat(ii,1);
            Mat(ii,4)=Mat(ii,2);
        end
    else
        d=sqrt((Mat(ii,3)-Mat(ii,5))^2+(Mat(ii,4)-Mat(ii,6))^2);
        if d<=tol
            Mat(ii,5)=Mat(ii,3);
            Mat(ii,6)=Mat(ii,4);
        end
    end
end
        