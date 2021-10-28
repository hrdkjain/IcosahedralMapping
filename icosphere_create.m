function [ v_out, f_out, planar_mapping_out] = icosphere_create( subdivisions )

% create vertices
t = (1.0 + sqrt(5.0)) / 2.0;
v = [...
    -1,  t,  0; 1,  t,  0; -1, -t,  0;  1, -t,  0; ...
     0, -1,  t; 0,  1,  t;  0, -1, -t;  0,  1, -t; ...
     t,  0, -1; t,  0,  1; -t,  0, -1; -t,  0,  1;];
v = v ./ sqrt(sum(v.^2,2));

% rotate v5 around x axis so it is at [0,0,1]
angle = atan2(v(5,2),v(5,3));
rx = [1 0 0;0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
v = (rx*v')';

% zero out values very close to zero
tol = 1.e-10;
v(abs(v)<tol) = 0;

% reorder vertices for easier map and face creation
map = [5 10 6 12 3 4 2 1 11 7 9 8];
v= v(map,:);


% figure;
% hold on;
% scatter3(v(:,1),v(:,2),v(:,3));
% text(v(:,1),v(:,2),v(:,3), cellstr('v'+string(1:size(v,1))'))
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');

% create base faces
f = [...
    1 2 3; 2 7  3; 3 7   8; 7  12  8;   % map 1
    1 3 4; 3 8  4; 4 8   9; 8  12  9;   % map 2
    1 4 5; 4 9  5; 5 9  10; 9  12 10;   % map 3
    1 5 6; 5 10 6; 6 10 11; 10 12 11;   % map 4
    1 6 2; 6 11 2; 2 11  7; 11 12  7;   % map 5
    ];


% create subdivisons and planar mapping
num_vertical_vertices = 2^subdivisions+1;
num_horizontal_vertices = 2^(subdivisions+1)+1;
vert_mat = zeros(num_vertical_vertices,num_horizontal_vertices,5,3);
planar_mapping_out = zeros(num_vertical_vertices,num_horizontal_vertices,5);


% fill maps
for i = 1:5
    f_idx = (i-1)*4+1;  
    
    % left quad
    % fill corners of quad
    vert_mat(1,1,i,:) = v(f(f_idx,1),:);
    vert_mat(num_vertical_vertices,1,i,:) = v(f(f_idx,2),:);
    vert_mat(1,num_vertical_vertices,i,:) = v(f(f_idx,3),:);
    
    f_idx = f_idx+1;
    vert_mat(num_vertical_vertices,num_vertical_vertices,i,:) = v(f(f_idx,2),:);
    
    % interpolate other values of quad
    vert_mat(:,:,i,:) = fill_quad_recursively(squeeze(vert_mat(:,:,i,:)), 1, num_vertical_vertices, 1, num_vertical_vertices);
    
    % right quad
    f_idx = f_idx+1;
    
    % fill corners of quad
    vert_mat(1,num_vertical_vertices,i,:) = v(f(f_idx,1),:);
    vert_mat(num_vertical_vertices,num_vertical_vertices,i,:) = v(f(f_idx,2),:);
    vert_mat(1,end,i,:) = v(f(f_idx,3),:);
    
    f_idx = f_idx+1;
    vert_mat(end,end,i,:) = v(f(f_idx,2),:);
    
    % interpolate other values of quad
    vert_mat(:,:,i,:) = fill_quad_recursively(squeeze(vert_mat(:,:,i,:)), 1, num_vertical_vertices, num_vertical_vertices, num_horizontal_vertices);

end

% flip map dimension, so map 5 is top map when ordering on plane
vert_mat = flip(vert_mat,3);

% extract vertices and planar mapping
v_out = reshape(permute(vert_mat(2:end,1:end-1,:,:),[2,1,3,4]), [], 3);
ind = 1:(num_vertical_vertices-1)*5*(num_horizontal_vertices-1);
planar_mapping_out(2:end,1:end-1, :) = permute(reshape(ind, num_horizontal_vertices-1, num_vertical_vertices-1,5),[2,1,3]);

% add vertices at north/south pole 
v_out(end+1,:) = vert_mat(1,1,1,:);
planar_mapping_out(1,1,:) = size(v_out,1);
v_out(end+1,:) = vert_mat(end,end,1,:);
planar_mapping_out(end,end,:) = size(v_out,1);

% fill in boundaries of the planar mapping
% next = [2:5 1];
prev = [5 1:4];
% top left
planar_mapping_out(1, 2:num_vertical_vertices,:) = planar_mapping_out(2:num_vertical_vertices,1,prev);
% top right
planar_mapping_out(1, num_vertical_vertices+1:end,:) = planar_mapping_out(end,2:num_vertical_vertices,prev);
% right 
planar_mapping_out(2:num_vertical_vertices,end,:) = planar_mapping_out(end,num_vertical_vertices+1:end,prev);  

% create faces
v_idx_1 = reshape(planar_mapping_out(1:end-1,1:end-1,:),[],1);
v_idx_2 = reshape(planar_mapping_out(2:end,1:end-1,:),[],1);
v_idx_3 = reshape(planar_mapping_out(1:end-1,2:end,:),[],1);
v_idx_4 = reshape(planar_mapping_out(2:end,2:end,:),[],1);

f_out = [v_idx_1 v_idx_2 v_idx_3; v_idx_2 v_idx_4 v_idx_3];

% delete boundary in planar mapping and reshape as required
planar_mapping_out = reshape(permute(planar_mapping_out(2:end,1:end-1,:),[2,1,3]),num_horizontal_vertices-1, []).';

end

function [vert_mat] = fill_quad_recursively(vert_mat, r_min, r_max, c_min, c_max)
    if r_min >= r_max-1 || c_min >= c_max-1
        return;
    end
    
    r_center = ceil((r_min+r_max)/2); 
    c_center = ceil((c_min+c_max)/2);   
    
    % diagonal edge
    vert_mat(r_center,c_center,:) = (vert_mat(r_min,c_max,:) + vert_mat(r_max,c_min,:)) / 2;
    vert_mat(r_center,c_center,:) = vert_mat(r_center,c_center,:) / sqrt(sum(vert_mat(r_center,c_center,:).^2));

    % left edge
    vert_mat(r_center,c_min,:) = (vert_mat(r_min,c_min,:) + vert_mat(r_max,c_min,:)) / 2; 
    vert_mat(r_center,c_min,:) = vert_mat(r_center,c_min,:) / sqrt(sum(vert_mat(r_center,c_min,:).^2));
    
    % right edge
    vert_mat(r_center,c_max,:) = (vert_mat(r_min,c_max,:) + vert_mat(r_max,c_max,:)) / 2;
    vert_mat(r_center,c_max,:) = vert_mat(r_center,c_max,:) / sqrt(sum(vert_mat(r_center,c_max,:).^2));
    
    % top edge
    vert_mat(r_min,c_center,:) = (vert_mat(r_min,c_min,:) + vert_mat(r_min,c_max,:)) / 2;
    vert_mat(r_min,c_center,:) = vert_mat(r_min,c_center,:) / sqrt(sum(vert_mat(r_min,c_center,:).^2));
    
    % bottom edge
    vert_mat(r_max,c_center,:) = (vert_mat(r_max,c_min,:) + vert_mat(r_max,c_max,:)) / 2;  
    vert_mat(r_max,c_center,:) = vert_mat(r_max,c_center,:) / sqrt(sum(vert_mat(r_max,c_center,:).^2));
    
    % top left quad
    vert_mat = fill_quad_recursively(vert_mat, r_min, r_center, c_min, c_center);
    
    % top right quad
    vert_mat = fill_quad_recursively(vert_mat, r_min, r_center, c_center, c_max);
    
    % bottom left quad
    vert_mat = fill_quad_recursively(vert_mat, r_center, r_max, c_min, c_center);
    
    % bottom right quad
    vert_mat = fill_quad_recursively(vert_mat, r_center, r_max, c_center, c_max);       
end

