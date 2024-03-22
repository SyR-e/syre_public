    function createCoilGroups(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, coil_groupname, coilType, Ph, color)
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



