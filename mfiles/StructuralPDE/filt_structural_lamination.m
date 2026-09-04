function [sVonMises_filt,R_filt,psIron] = filt_structural_lamination(structModel,R,sVonMises,data4GeoMat)
% Set sVonMises to NaN for nodes outside the iron region.
% Builds psIron = fullModel \ union(psMagnet,psSleeve,psShaft,psResin).
% Does NOT modify the mesh. Returns R_filt identical to R.
%
% Outputs:
%   sVonMises_filt - sVonMises with NaN outside psIron
%   R_filt         - identical copy of R
%   psIron         - resulting polyshape (iron region)

sVonMises_filt = sVonMises;
R_filt = R;
psIron = polyshape();  % default empty

% --- Extract node coordinates ---
nodes = [];
% try
%     mesh = R.Mesh;
%     if isfield(mesh,'Nodes')
%         nodes = mesh.Nodes;
%     end
% catch
%     % try R.Nodes
%     try
%         if isfield(R,'Nodes')
%             nodes = R.Nodes;
%         end
%     catch
%     end
% end
nodes = R.Mesh.Nodes;

if isempty(nodes)
    error('Could not locate node coordinates in R. Expected R.Mesh.Nodes or R.Nodes.');
end

% Normalize to N x 2
if size(nodes,1) == 2
    pts = nodes';
elseif size(nodes,2) == 2
    pts = nodes;
else
    error('Unexpected Nodes shape. Expected 2xN or Nx2.');
end
xq = pts(:,1);
yq = pts(:,2);

% --- Collect exclusion polyshapes (to subtract from full model) ---
excludePolys = polyshape.empty;
if ~isempty(data4GeoMat) && isstruct(data4GeoMat)
    fieldsExclude = {'psMagnet','psSleeve','psShaft','psResin'};
    for k = 1:numel(fieldsExclude)
        f = fieldsExclude{k};
        if isfield(data4GeoMat,f)
            p = data4GeoMat.(f);
            if isempty(p)
                continue;
            end
            if isa(p,'polyshape')
                excludePolys(end+1) = p; %#ok<AGROW>
            elseif isnumeric(p)
                if size(p,2) == 2
                    excludePolys(end+1) = polyshape(p(:,1), p(:,2)); %#ok<AGROW>
                elseif size(p,1) == 2
                    excludePolys(end+1) = polyshape(p(1,:)', p(2,:)'); %#ok<AGROW>
                end
            elseif iscell(p)
                if numel(p) >= 2 && isnumeric(p{1}) && isnumeric(p{2})
                    excludePolys(end+1) = polyshape(p{1}(:), p{2}(:)); %#ok<AGROW>
                elseif numel(p) == 1 && isnumeric(p{1})
                    arr = p{1};
                    if size(arr,2) == 2
                        excludePolys(end+1) = polyshape(arr(:,1), arr(:,2)); %#ok<AGROW>
                    end
                end
            end
        end
    end
end

% --- Build fullModel polyshape ---
fullPoly = polyshape.empty;
% Prefer explicit full model if provided
if ~isempty(data4GeoMat) && isstruct(data4GeoMat)
    if isfield(data4GeoMat,'psModel') && isa(data4GeoMat.psModel,'polyshape')
        fullPoly = data4GeoMat.psModel;
    elseif isfield(data4GeoMat,'psFull') && isa(data4GeoMat.psFull,'polyshape')
        fullPoly = data4GeoMat.psFull;
    end
end
% If not available, use convex hull of mesh nodes as conservative full region
if isempty(fullPoly)
    try
        k = convhull(xq,yq);
        fullPoly = polyshape(xq(k), yq(k));
    catch
        % fallback bounding box polygon
        xmin = min(xq); xmax = max(xq);
        ymin = min(yq); ymax = max(yq);
        fullPoly = polyshape([xmin xmax xmax xmin],[ymin ymin ymax ymax]);
    end
end

% --- Compute union of exclusions ---
if isempty(excludePolys)
    unionExcl = polyshape(); % empty
else
    unionExcl = excludePolys(1);
    for k = 2:numel(excludePolys)
        unionExcl = union(unionExcl, excludePolys(k));
    end
end

% --- psIron is full minus exclusions (use subtract) ---
if isempty(unionExcl) || isempty(unionExcl.Vertices)
    psIron = fullPoly;
else
    % prefer subtract; polyshape.subtract performs fullPoly - unionExcl
    try
        psIron = subtract(fullPoly, unionExcl);
    catch
        % fallback: compute difference by splitting into pieces (approx)
        psIron = fullPoly;
        try
            % if unionExcl is array, subtract each
            for k = 1:numel(excludePolys)
                psIron = subtract(psIron, excludePolys(k));
            end
        catch
            % leave psIron as fullPoly if subtract fails
            psIron = fullPoly;
        end
    end
end

% If psIron ended empty, nothing to keep
if isempty(psIron) || isempty(psIron.Vertices)
    sVonMises_filt(:) = NaN;
    return;
end

% --- Mark nodes inside or on border of psIron as kept ---
% inpolygon returns in (in or on) and on. We want inside OR on border => use in flag
[xv, yv] = boundary(psIron); % returns boundary vertices (works R2020a+)
% If boundary fails, fallback to Vertices
if isempty(xv)
    V = psIron.Vertices;
    xv = V(:,1); yv = V(:,2);
end

[in, on] = inpolygon(xq, yq, xv, yv);
keep = in; % in is true for points inside or on boundary

% Set NaN for nodes outside psIron
sVonMises_filt(~keep) = NaN;

end
