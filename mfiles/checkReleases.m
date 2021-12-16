% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function out = checkReleases()

out=1;

clc

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
disp('Check minimum release requirements')
disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

% Matlab
disp('Check Matlab release...')
disp('   minimum requirement      : 2016b (9.1)')
disp('   suggested release        : 2021b (9.11)')

tmp = ver('matlab');
vMatlab = eval(tmp.Version);

disp(['   installed Matlab release : ' tmp.Release(3:end-1) ' (' tmp.Version ')'])

if vMatlab>=9.1
    disp('(v)Matlab release OK')
else
    disp('(x)Matlab release too old. Some compatibility issues may arise')
    out = 0;
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

% FEMM

disp('Check FEMM version...')
disp('   minimum requirement      : 25Feb2018')

if exist('mo_getgapb')
    disp('(v)FEMM version OK')
else
    disp('(x)FEMM version too old. Please update FEMM')
    out = 0;
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

% Parallel computing toolbox
disp('Parallel Computing Toolbox check...')
ppState = parallelComputingCheck(0);

if ppState<0
    disp('(x)Parallel Computing Toolbox not available')
else
    disp('(v)Parallel Computing Toolbox available')
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

% CF Toolbox
disp('Curve Fitting Toolbox check...')
flag=0;

test = ver();
for ii=1:length(test)
    if strcmp(test(ii).Name,'Curve Fitting Toolbox')
        flag=1;
    end
end

if ~flag
    disp('(x)Curve Fitting Toolbox not available')
else
    disp('(v)Curve Fitting Toolbox available')
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

% PDE Toolbox
disp('PDE Toolbox check...')
flag=0;

test = ver();
for ii=1:length(test)
    if strcmp(test(ii).Name,'Partial Differential Equation Toolbox')
        flag=1;
    end
end

if ~flag
    disp('(x)PDE Toolbox not available')
else
    disp('(v)PDE Toolbox available')
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

if out
    disp('Minimum requirements fulfilled')
else
    disp('Minimum requirements NOT fulfilled!')
end

disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')

if nargout==0
    clear out
end

