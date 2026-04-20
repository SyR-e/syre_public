function FemmProblem = modifyboundprop_mfemm(FemmProblem,BCname,BCprop,value)

foundBC = false;

for ii=1:numel(FemmProblem.BoundaryProps)
    if strcmp(FemmProblem.BoundaryProps(ii).Name, BCname)
        FemmProblem.BoundaryProps(ii).(BCprop) = value;
        foundBC = true;
    end
end

if ~foundBC
    warning('Boundary condition "%s" not found.', BCname);
end