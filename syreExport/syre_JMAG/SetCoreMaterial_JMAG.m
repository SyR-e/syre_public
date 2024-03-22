
function SetCoreMaterial_JMAG(study,core_groupname,core_mat,laminated_type,laminated_factor,hysteresisloop_index,eddycurrent_index)

    study.SetMaterialByName(core_groupname, core_mat);
    study.GetMaterial(core_groupname).SetValue('Laminated', laminated_type);
    study.GetMaterial(core_groupname).SetValue('LaminationFactor', laminated_factor);
    study.GetMaterial(core_groupname).SetValue('MagnetizationCorrection', 100);
if hysteresisloop_index==1
    study.GetMaterial(core_groupname).SetValue('UseMaterialHysteresisLoop', hysteresisloop_index);
end
    study.GetMaterial(core_groupname).SetValue('EddyCurrentCalculation', eddycurrent_index);
if eddycurrent_index==0
    study.GetMaterial(core_groupname).SetValue('ConductivityType', 0);
    study.GetMaterial(core_groupname).SetValue('UserConductivityType', 0);
end
end

