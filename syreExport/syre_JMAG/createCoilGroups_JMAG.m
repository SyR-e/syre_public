    %  Copyright 2024
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

function createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, coil_groupname, coilType, Ph, color)
        [row, col] = find(geo.win.avv == Ph);
        for slot_ID = 1:length(row)
            if ~isempty(row)
                P = JDesigner.CreatePoint(CoilCGX_reshape(row(slot_ID), col(slot_ID)), ...
                                          CoilCGY_reshape(row(slot_ID), col(slot_ID)), 0);
                part = model.GetPartByPosition(P);
                part.SetName(strcat(coilType, num2str(slot_ID)));
                part.SetColor(color);
                model.GetGroupList().AddPartToGroup(strcat(coil_groupname, '_', coilType(1)), ...
                                                    strcat(coilType, num2str(slot_ID)));
            end
        end
 end



