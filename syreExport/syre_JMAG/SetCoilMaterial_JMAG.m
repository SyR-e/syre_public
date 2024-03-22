function SetCoilMaterial(study,coil_mat,eddycurrent_c,Ph)
        study.SetMaterialByName(strcat('Coil',Ph), coil_mat);
        study.GetMaterial(strcat('Coil',Ph)).SetValue('Laminated', 0);
        study.GetMaterial(strcat('Coil',Ph)).SetValue('EddyCurrentCalculation', eddycurrent_c);
        study.GetMaterial(strcat('Coil',Ph)).SetValue('UserConductivityType', 0);%'0=use material resistivity 1=Electric conductivity 2=Electric resistivity 3=Temperature dependant conductivity 4=Temperature dependant resistivity
end
% study.GetMaterial(strcat(coil_groupname,'_U')).SetValue('UserResistivityValue', coil_conductivity);
