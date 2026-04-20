function FemmProblem = XFEMM_initialize(geo,mat)

FemmProblem = newproblem_mfemm('planar',...
    'LengthUnits','Millimiters',...
    'Precision',1e-8,...
    'Depth',geo.l,...
    'MinAngle',15,...
    'SmartMesh',false);


% iron (stator)
indMat = 1;
FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.Stator.MatName);
FemmProblem.Materials(indMat).BHPoints = mat.Stator.BH;
FemmProblem.Materials(indMat).Mu_x     = 1;
FemmProblem.Materials(indMat).Mu_y     = 1;
FemmProblem.Materials(indMat).LamFill  = geo.stackingFactor;

indMat = indMat+1;

% iron (rotor)
if ~strcmp(mat.Stator.MatName,mat.Rotor.MatName)
    FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.Rotor.MatName);
    FemmProblem.Materials(indMat).BHPoints = mat.Rotor.BH;
    FemmProblem.Materials(indMat).Mu_x     = 1;
    FemmProblem.Materials(indMat).Mu_y     = 1;
    FemmProblem.Materials(indMat).LamFill  = geo.stackingFactor;

    indMat = indMat+1;
end

% iron (shaft)
if ~strcmp(mat.Stator.MatName,mat.Shaft.MatName)
    if ~strcmp(mat.Shaft.MatName,'Air')
        FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.Shaft.MatName);
        FemmProblem.Materials(indMat).BHPoints = mat.Shaft.BH;
        FemmProblem.Materials(indMat).Mu_x     = 1;
        FemmProblem.Materials(indMat).Mu_y     = 1;
        FemmProblem.Materials(indMat).LamFill  = geo.stackingFactor;
    
        indMat = indMat+1;
    end
end

% conductor (stator)
FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.SlotCond.MatName);
FemmProblem.Materials(indMat).Sigma = mat.SlotCond.sigma/1e6;

indMat = indMat+1;

if ~strcmp(mat.SlotCond.MatName,mat.BarCond.MatName)
    FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.BarCond.MatName);
    FemmProblem.Materials(indMat).Sigma = mat.BarCond.sigma/1e6;

    indMat = indMat+1;
end

% permanent magnet
FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.LayerMag.MatName);
FemmProblem.Materials(indMat).Mu_x  = mat.LayerMag.mu;
FemmProblem.Materials(indMat).Mu_y  = mat.LayerMag.mu;
FemmProblem.Materials(indMat).H_c   = mat.LayerMag.Hc;
FemmProblem.Materials(indMat).Sigma = mat.LayerMag.sigmaPM/1e6;

indMat = indMat+1;

% sleeve
FemmProblem.Materials(indMat) = newmaterial_mfemm(mat.Sleeve.MatName);

% air
FemmProblem.Materials(indMat) = newmaterial_mfemm('Air');
FemmProblem.Materials(indMat).Mu_x    = 1;
FemmProblem.Materials(indMat).Mu_y    = 1;
FemmProblem.Materials(indMat).LamFill = 0;







