function [ area ] = polygon_area( v , normal)
% assumes normal has norm == 1

normal_abs = abs(normal);
b = 3; indices = [1 2];
if normal_abs(1) > normal_abs(2)
    if normal_abs(1) > normal_abs(3)
        b = 1; indices = [2 3];
    end
elseif normal_abs(2) > normal_abs(3)
    b = 2; indices = [3 1];                
end           

% all intermediate vertices
area = v(2:end-1,indices(1))' * (v(3:end,indices(2))-v(1:end-2,indices(2)));

% first vertex
area = area + v(1,indices(1))' * (v(2,indices(2))-v(end,indices(2)));

%last vertex
area = area + v(end,indices(1))' * (v(1,indices(2))-v(end-1,indices(2)));

%area normalization
area = abs(area / (2*normal(b)));
end

