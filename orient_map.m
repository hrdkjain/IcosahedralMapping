function [map] = orient_map(sph_verts, verts, faces, filePath)
% c_faces: face connectivity 
% sph_verts: spherically parameterized mesh vertices
% verts: input mesh vertices

% Compute BB of the verts
 [xmin,ymin,zmin,xmax,ymax,zmax] = bbox(verts);
 
% Select points closest to BB's North and East
% North and West are defined seperately for each dataset
if contains(filePath, "airplane") && contains(filePath, "ShapeNetCore.v2")
    BBCorner_north= [(xmin+xmax)/2 (ymax+ymin)/2 zmin];
    BBCorner_east = [(xmin+xmax)/2 ymax zmin];
elseif contains(filePath, "chair") && contains(filePath, "ShapeNetCore.v2")
    BBCorner_north= [xmax ymin zmax];
    BBCorner_east = [xmax ymax zmax];
elseif contains(filePath, "table") && contains(filePath, "ShapeNetCore.v2")
    BBCorner_north = [xmax ymin zmax];
    BBCorner_east = [xmax ymax zmax];
elseif contains(filePath, "chair") && contains(filePath, "trail")
    BBCorner_north = [xmax ymin zmax];
    BBCorner_east = [xmax ymax zmax];
elseif contains(filePath, "ModelNet10")
    BBCorner_north = [xmax ymax zmax];
    BBCorner_east = [xmax ymin zmax];
end

% draw a line from north to east
BBCorner_northEastLine = computeLinePoints(BBCorner_north, BBCorner_east, 1);

DT =  delaunayTriangulation(verts);
id_northEastLine = [];
for i = 1:size(BBCorner_northEastLine,1)
    id_northEastLine = [id_northEastLine; nearestNeighbor(DT,BBCorner_northEastLine(i,:))];  
end

% Make north as positive Z-axis
R = vrrotvec2mat(vrrotvec([0,0,1], sph_verts(id_northEastLine(1),:)));
map = sph_verts * R;

% Make east on XZ plane
% compute the Z rotation of the east point
angleZ = atan2d(map(id_northEastLine(end),2), map(id_northEastLine(end),1));
% Should only be done by rotation around Z-axis
R = rotz(angleZ);
map = map * R;

% Alternatively, more points on the BBCorner_northEastLine can be used for 
% alignment to the XZ Plane. More robust ??

% For visualization
% C = 10*ones(size(sph_verts));
% C(id_northEastLine,:) = 255;
% write_off("tmp.off", map, faces, [], C);
end

function out = computeLinePoints(A, B, n_divisions)
out = [];
for f = 0:n_divisions
    out = [out; (B*f + A*(n_divisions-f))/n_divisions];
end
end

