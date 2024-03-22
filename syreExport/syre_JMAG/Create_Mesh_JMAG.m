function Create_Mesh_JMAG(Mesh,partname,groupname,meshsize)
        Part_mesh = Mesh.CreateCondition('Part', partname);
        Part_mesh.SetName(partname);
        Part_mesh.SetValue('Size', meshsize);
        Part_mesh.ClearParts();
        sel = Part_mesh.GetSelection();
        sel.SelectPart(groupname);
        Part_mesh.AddSelected(sel);
end