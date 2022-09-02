% Copyright 2022
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

function [f] = defineCmatrix(region,y,dataForCmaxtrix)

E_Fe     = dataForCmaxtrix.E_Fe;
E_PM     = dataForCmaxtrix.E_PM;
nu       = dataForCmaxtrix.nu;
psMagnet = dataForCmaxtrix.psMagnet;
E_sleeve = dataForCmaxtrix.E_sleeve;
psSleeve = dataForCmaxtrix.psSleeve;


xy = (region.x)+j*(region.y);

f = zeros(10,length(xy));
if ~isempty(psMagnet)
    for ii=1:length(xy)
        psTest = nsidedpoly(20,'Center',[real(xy(ii)) imag(xy(ii))],'Radius',eps*1e10);
        [psInt] = intersect(psMagnet,psTest);
        if psInt.NumRegions==0
            [psInt] = intersect(psSleeve,psTest);
            if psInt.NumRegions==0
                % iron section
                f(:,ii) = [
                    E_Fe/(1-nu^2)
                    0
                    E_Fe/(2*(1+nu))
                    0
                    E_Fe/(2*(1+nu))
                    E_Fe*nu/(1-nu^2)
                    0
                    E_Fe/(2*(1+nu))
                    0
                    E_Fe/(1-nu^2)];
            else
                % sleeve section
                f(:,ii) = [
                    E_sleeve/(1-nu^2)
                    0
                    E_sleeve/(2*(1+nu))
                    0
                    E_sleeve/(2*(1+nu))
                    E_sleeve*nu/(1-nu^2)
                    0
                    E_sleeve/(2*(1+nu))
                    0
                    E_sleeve/(1-nu^2)];
            end
        else
            % PM section --> no structural and set as gum
            f(:,ii) = [
                E_PM/(1-nu^2)
                0
                E_PM/(2*(1+nu))
                0
                E_PM/(2*(1+nu))
                E_PM*nu/(1-nu^2)
                0
                E_PM/(2*(1+nu))
                0
                E_PM/(1-nu^2)];
        end
    end
end


