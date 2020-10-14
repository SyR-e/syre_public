function [MatrixWin] = MatrixWin()

%%======== MATRIX OF WINDINGS =============================================
%% This function take the vectors of koil and transform it in a matrix,
%% which is used in Syre
%%=========================================================================


thisfilepath = fileparts(which('syre.png'));
addpath (fullfile (thisfilepath,'koil'));

run koil.m

m = 3; % number of phases
[~,Q] = size(ka); % size of ka
MatrixWin = zeros(2,Q); % matrix of Syre
k = zeros(3,Q); % union of matrices of koil
k = [ka;kb;kc]; % union of matrices of koil
c = 0; % flag true or false for build MatrixWin
for i = 1 : 1 : Q
    for j = 1 : 1 : m
        if j == 1
            if abs(k(j,i)) == 1
                MatrixWin(1,i) = j*sign((k(j,i)));
                MatrixWin(2,i) = j*sign((k(j,i)));
            elseif abs(k(j,i)) == 0.5
                MatrixWin(1,i) = j*sign((k(j,i)));
            end
        elseif j == 3
            if abs(k(j,i)) == 1
                MatrixWin(1,i) = j*sign((k(j,i)));
                MatrixWin(2,i) = j*sign((k(j,i)));
            elseif abs(k(j,i)) == 0.5
                MatrixWin(2,i) = j*sign((k(j,i)));
            end
        else
            if abs(k(j,i)) == 1
                MatrixWin(1,i) = j*sign((k(j,i)));
                MatrixWin(2,i) = j*sign((k(j,i)));
            elseif abs(k(j,i)) == 0.5
                if abs(k(j+1,i)) ~= 0
                    MatrixWin(1,i) = j*sign((k(j,i)));
                else
                    MatrixWin(2,i) = j*sign((k(j,i)));
                end
            end    
        end
    end
end




