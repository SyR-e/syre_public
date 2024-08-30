function [primal,mat,X,P,T,M] = readFemm(fileName)

% READFEMM reads FEMM file with solution
%
% USE:
% [primal,mat,X] = readFemm(fileName)
%
% INPUT:
% 'fileName': name of FEMM solutionfile (with extension)
%
% OUTPUT:
% 'primal': struct with mesh information, created by CREATEPRIMAL2D
% 'mat': struct with material information
% 'X': vector with nodal solution
% 'P': node coordinates
% 'T': tetrahedra/hexahedra connectivity matrix
% 'M': vector with material codes
%
% NOTE:
%
% VERSION:
% Date: 29.04.2015
% Copyright(C) 2015-2024: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 19.06.2018: reading region IDs
% 01.08.2023: file in final version
% 13.02.2024: added P,T,M as outputs

% fprintf('reading mesh data from FEMM solution file...  ');
% tic;

% open file
fid = fopen(fileName,'r');
if fid == -1
    error('Could not open FEMM file\n');
end
frewind(fid);

% string to find
frequencyString = '[Frequency]';
solutionString = '[Solution]';
unitsString = '[LengthUnits]';
blockString = '[NumBlockLabels]';
propsString = '[BlockProps]';
BlockName = '<BlockName>';
Mur = '<Mu_x>';
Hc = '<H_c>';
HcAngle = '<H_cAngle>';
Sigma = '<Sigma>';
BHpoints = '<BHPoints>';

% read frequency
while true
    oneLine = textscan(fid,'%s %s %f\n',1);
    if strcmpi(oneLine{1},frequencyString)
        frequency = oneLine{3};
        break
    end
end

% read units
while true
    oneLine = textscan(fid,'%s %s %s\n',1);
    if strcmpi(oneLine{1},unitsString)
        if strcmpi(oneLine{3},'inches')
            scale = 2.54e-2;
        elseif strcmpi(oneLine{3},'millimeters')
            scale = 1e-3;
        elseif strcmpi(oneLine{3},'centimeters')
            scale = 1e-2;
        elseif strcmpi(oneLine{3},'mils')
            scale = 2.54e-5;
        elseif strcmpi(oneLine{3},'microns')
            scale = 1e-6;
        elseif strcmpi(oneLine{3},'meters')
            scale = 1;
        end
        break
    end
end

% material properties
while true
    oneLine = textscan(fid,'%s %s %d\n',1);
    if strcmpi (oneLine{1}, propsString)
        N = oneLine{3};
        mat = setMaterialProperties(N);
        break
    end
end

for i = 1:N
    while true
        oneLine = textscan(fid,'%s %s "%s" \n',1);
        if strcmpi (oneLine{1}, BlockName)
            mat(i).Name = oneLine{3}{1}(1:end-1);
            break
        end
    end
    while true
        oneLine = textscan(fid,'%s %s %f\n',1);
        if strcmpi (oneLine{1}, Mur)
            mat(i).MuR = oneLine{3};
            break
        end
    end
    while true
        oneLine = textscan(fid,'%s %s %f\n',1);
        if strcmpi (oneLine{1}, Hc)
            mat(i).Hc = oneLine{3};
            break
        end
    end
    while true
        oneLine = textscan(fid,'%s %s %f\n',1);
        if strcmpi (oneLine{1}, HcAngle)
            mat(i).HcAngle = oneLine{3};
            break
        end
    end
    while true
        oneLine = textscan(fid,'%s %s %f\n',1);
        if strcmpi (oneLine{1}, Sigma)
            mat(i).Sigma = 1e6*oneLine{3};
            break
        end
    end
    % BHpoints
    while true
        oneLine = textscan(fid,'%s %s %d\n',1);
        if strcmpi (oneLine{1}, BHpoints)
            if oneLine{3} > 0
                L = oneLine{3};
                for j = 1:L
                    data = textscan(fid,'%f %f\n',L);
                    B(j) = data{1};
                    H(j) = data{2};
                end
                mat = loadMagneticNonLinearCurve(mat,i,H,B);
            end
            break
        end
    end
end

% lettura NumBlockLabels
while true
    oneLine = textscan(fid,'%s %s %d\n',1);
    if strcmpi (oneLine{1}, blockString)
        Nblock = oneLine{3};
        break
    end
end
data = textscan(fid,'%f %f %f %f %f %f %f %f %f',Nblock);
K = data{3};

% go to solution part
while true
    oneLine = textscan(fid,'%s\n',1);
    if strcmpi(oneLine{1},solutionString)
        oneLine = textscan(fid,'%d\n',1);
        nPt = oneLine{1};
        break
    end
end

if frequency == 0
    data = textscan(fid,'%f %f %f %f%[^\n]',nPt);
    P = scale*[data{1:2} zeros(nPt,1)];
    X = data{3};
else
    data = textscan(fid,'%f %f %f %f%[^\n]',nPt);
    P = scale*[data{1:2} zeros(nPt,1)];
    X = data{3}+1j*data{4};
end
% number of triangles
nTri = fscanf(fid,'%i\n', 1)';

data = textscan(fid,'%d %d %d %d %d %d %d %d%[^\n]',nTri);
T = 1+[data{1:3}];
M = 1+data{4};

% remapping
M = K(M);

fclose(fid);

primal = createPrimal2d(P,T,M);

% % load triangles and object codes
% data = fscanf(fid, '%i %i %i %i', [4,nTri])';
% M = int32(data(:,4))+1; (blocco di ogni triangolo)
% T = int32(data(:,1:3))+1; (coordinate nodi triangolo)
% K = NumBlockLabels


% fprintf('done %6.2f sec\n',toc);