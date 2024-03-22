function createSet (model, groupname,partname,MagnetCG)

        set = model.GetSetList().CreatePartSet(groupname);
        set.SetName(groupname);
        set.SetMatcherType('selection');
        set.ClearParts();
        sel = set.GetSelection();
        if strcmp (partname, 'Magnet')
for ID=1:1:length(MagnetCG)
    sel.SelectPart(strcat('Magnet',num2str(ID)));
end
        else
        sel.SelectPart(partname);
        end
        set.AddSelected(sel);
end

