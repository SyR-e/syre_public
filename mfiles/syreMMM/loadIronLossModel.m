% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IronLossModel] = loadIronLossModel(filename)

if strcmp(filename,'0')
    IronLossModel.type='0';
else
    load(filename)
    if exist('Pfes_c','var') % model from MagNet maps in idiq
        IronLossModel.type='map';
%         IronLossModel.Pfe_h=Pfes_h+Pfer_h;
%         IronLossModel.Pfe_c=Pfes_c+Pfer_c;
        IronLossModel.Pfes_h = Pfes_h;
        IronLossModel.Pfes_c = Pfes_c;
        IronLossModel.Pfer_h = Pfer_h;
        IronLossModel.Pfer_c = Pfer_c;
        IronLossModel.Ppm    = Ppm;
        IronLossModel.n0     = velDim;
        IronLossModel.Id     = Id;
        IronLossModel.Iq     = Iq;
        
        if exist('geo','var') % from FEA simulations
            p = geo.p;
            IronLossModel.f0    = IronLossModel.n0/60*p;
            IronLossModel.expC  = 2;
            IronLossModel.expH  = mat.Stator.alpha;
            IronLossModel.expPM = 2;
            IronLossModel.segPM = 1;
        elseif exist('factors','var') % from MMM GUI
            IronLossModel.f0    = factors.f0;
            IronLossModel.expC  = factors.expC;
            IronLossModel.expH  = factors.expH;
            IronLossModel.expPM = factors.expPM;
            IronLossModel.segPM = 1;
        else
            prompt={'Pole pair number',...
                'Exponent for eddy current loss', ...
                'Exponent for hysteresis loss', ...
                'Exponent for PM loss', ...
                'Segmentation factor for PM loss'};
            name='Input';
            numlines=1;
            defaultanswer={'2','2','1.2','2','1'};
            setup=inputdlg(prompt,name,numlines,defaultanswer);
            
            p=eval(setup{1});
            IronLossModel.f0=IronLossModel.n0/60*p;
            IronLossModel.expC=eval(setup{2});
            IronLossModel.expH=eval(setup{3});
            IronLossModel.expPM=eval(setup{4});
            IronLossModel.segPM=eval(setup{5});
            
%             % check axis style
%             if exist('dataSet','var')
%                 if ~isfield(dataSet,'axisType')
%                     if strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Vtype')
%                         dataSet.axisType = 'PM';
%                     else
%                         dataSet.axisType = 'SR';
%                     end
%                 end
%                 if ~strcmp(dataSet.axisType,motorModel.data.axisType)
%                     tmp.FluxMap_dq.Id    = Id;
%                     tmp.FluxMap_dq.Iq    = Iq;
%                     tmp.FluxMap_dq.Fd    = Fd;
%                     tmp.FluxMap_dq.Fq    = Fq;
%                     tmp.IronPMLossMap_dq = IronLossModel;
%                     if strcmp(dataSet.axisType,'SR')
%                         tmp = MMM_sr2pm(tmp);
%                     else
%                         tmp = MMM_pm2sr(tmp);
%                     end
%                     IronLossModel = tmp.IronPMLossMap_dq;
%                 end
%             end
        end
        
    elseif exist('IronLossModel','var') % model from fit
        if ~isfield(IronLossModel,'type')
            if isfield(IronLossModel,'param')
                IronLossModel.type='fitIM1';
            else
                error('Model not valid')
            end
        end
    elseif exist('out','var') %single point, to verify
        IronLossModel.type='point';
        IronLossModel.Pfe0=out.Pfes_h+out.Pfes_c+out.Pfer_h+out.Pfer_c;
        IronLossModel.F0=abs(out.Fd+j*out.Fq);
        IronLossModel.n0=out.velDim;
        
        if exist('geo','var')
            p = geo.p;
            IronLossModel.f0      = IronLossModel.n0/60*p;
            IronLossModel.expFlux = mat.Stator.beta;
            IronLossModel.expFreq = mat.Stator.alpha;
        else
            prompt={'Pole pair number',...
                'Exponent for flux linkage', ...
                'Exponent for frequency'};
            name='Input';
            numlines=1;
            defaultanswer={'2','2','1.8'};
            setup=inputdlg(prompt,name,numlines,defaultanswer);
            
            p=eval(setup{1});
            IronLossModel.f0=IronLossModel.n0*p/60;
            IronLossModel.expFlux=eval(setup{2});
            IronLossModel.expFreq=eval(setup{3});
        end
    else
        IronLossModel = [];
    end
end