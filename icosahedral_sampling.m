function [ v_ico, f_ico, planar_mapping, sampled_pos, sparse_indices, sparse_weights] = icosahedral_sampling( subdivisions, c_faces, c_sph_verts, c_verts, c_face_normals)



%% create icosahedron with wanted subdivisions
[v_ico, f_ico, planar_mapping] = icosphere_create(subdivisions);
f_ico_cell = num2cell(f_ico,2);
[v_dual, f_dual] = unitsphere_createDualMesh(v_ico, f_ico_cell);


% ignore faces that have no or nearly no area.
% this shouldn't happen in a spherical parametrization, but it does currently
v12 = c_sph_verts(c_faces(:, 2),:) - c_sph_verts(c_faces(:, 1),:);
v13 = c_sph_verts(c_faces(:, 3),:) - c_sph_verts(c_faces(:, 1),:);
tmp_area = 0.5 * sqrt(sum(cross_fast(v12, v13).^2, 2));
c_faces_filtered = c_faces(tmp_area > 1e-7, :);


%% get latitude & longitude bounds of all meshes involved

bounds_latlong  = get_latlong_bounds( c_faces_filtered, c_sph_verts );
bounds_latlong_ico = get_latlong_bounds( f_dual, v_dual );

% for each face of the ico sphere get candidate faces of target meshes
candidate_faces = cell(1, length(f_dual));

% normal icosphere faces
for f_idx = 1:length(bounds_latlong_ico.I_normal)
    
    % normal faces
    if pi-bounds_latlong_ico.bounds_theta_normal(f_idx,1) < bounds_latlong_ico.bounds_theta_normal(f_idx,2)
        indices = binarySearchGreater(bounds_latlong.bounds_theta_normal_max_sorted, bounds_latlong_ico.bounds_theta_normal(f_idx,1)) + 1;
        II = bounds_latlong.I_theta_normal_max_sorted(indices:end);
        Icand_normal = II(bounds_latlong.bounds_theta_normal(II,1) < bounds_latlong_ico.bounds_theta_normal(f_idx,2));
    else
        indices = binarySearchLess(bounds_latlong.bounds_theta_normal_min_sorted, bounds_latlong_ico.bounds_theta_normal(f_idx,2)) - 1;
        II = bounds_latlong.I_theta_normal_min_sorted(1:indices);
        Icand_normal = II(bounds_latlong.bounds_theta_normal(II,2) > bounds_latlong_ico.bounds_theta_normal(f_idx,1));
    end
    Icand_normal2 = bounds_latlong.bounds_phi_normal(Icand_normal,1) < bounds_latlong_ico.bounds_phi_normal(f_idx,2) &  bounds_latlong.bounds_phi_normal(Icand_normal,2) > bounds_latlong_ico.bounds_phi_normal(f_idx,1);
    
    % jump faces
    Icand_jump = find(bounds_latlong.bounds_theta_jump(:,1) < bounds_latlong_ico.bounds_theta_normal(f_idx,2) & bounds_latlong.bounds_theta_jump(:,2) > bounds_latlong_ico.bounds_theta_normal(f_idx,1));
    Icand_jump2 = bounds_latlong.bounds_phi_jump(Icand_jump,2) > bounds_latlong_ico.bounds_phi_normal(f_idx,1) | bounds_latlong.bounds_phi_jump(Icand_jump,1) < bounds_latlong_ico.bounds_phi_normal(f_idx,2);
    
    % pole faces
    Ipole = [];
    if bounds_latlong.bounds_theta_northpole(2) >= bounds_latlong_ico.bounds_theta_normal(f_idx,1)
        Ipole = bounds_latlong.I_northpole;
    elseif bounds_latlong.bounds_theta_southpole(1) <= bounds_latlong_ico.bounds_theta_normal(f_idx,2)
        Ipole = bounds_latlong.I_southpole;
    end
    
    Icand = [bounds_latlong.I_normal(Icand_normal(Icand_normal2)); bounds_latlong.I_jump(Icand_jump(Icand_jump2)); Ipole];
    candidate_faces{bounds_latlong_ico.I_normal(f_idx)} = Icand;
end

% jump icosphere faces
for f_idx = 1:length(bounds_latlong_ico.I_jump)
    
    % normal faces
    Icand_normal = find(bounds_latlong.bounds_theta_normal(:,1) < bounds_latlong_ico.bounds_theta_jump(f_idx,2) & bounds_latlong.bounds_theta_normal(:,2) > bounds_latlong_ico.bounds_theta_jump(f_idx,1));
    Icand_normal2 = bounds_latlong.bounds_phi_normal(Icand_normal,2) > bounds_latlong_ico.bounds_phi_jump(f_idx,1) | bounds_latlong.bounds_phi_normal(Icand_normal,1) < bounds_latlong_ico.bounds_phi_jump(f_idx,2);
    
    % jump faces
    Icand_jump = find(bounds_latlong.bounds_theta_jump(:,1) < bounds_latlong_ico.bounds_theta_jump(f_idx,2) & bounds_latlong.bounds_theta_jump(:,2) > bounds_latlong_ico.bounds_theta_jump(f_idx,1));
    
    % pole faces
    Ipole = [];
    if bounds_latlong.bounds_theta_northpole(2) >= bounds_latlong_ico.bounds_theta_jump(f_idx,1)
        Ipole = bounds_latlong.I_northpole;
    elseif bounds_latlong.bounds_theta_southpole(1) <= bounds_latlong_ico.bounds_theta_jump(f_idx,2)
        Ipole = bounds_latlong.I_southpole;
    end
    
    Icand = [bounds_latlong.I_normal(Icand_normal(Icand_normal2)); bounds_latlong.I_jump(Icand_jump); Ipole];
    candidate_faces{bounds_latlong_ico.I_jump(f_idx)} = Icand;
end

% pole faces
Icand = [bounds_latlong.I_normal(bounds_latlong.bounds_theta_normal(:,1) <= bounds_latlong_ico.bounds_theta_northpole(2)); bounds_latlong.I_jump(bounds_latlong.bounds_theta_jump(:,1) <= bounds_latlong_ico.bounds_theta_northpole(2)); bounds_latlong.I_northpole];
candidate_faces{bounds_latlong_ico.I_northpole} = Icand;

Icand = [bounds_latlong.I_normal(bounds_latlong.bounds_theta_normal(:,2) >= bounds_latlong_ico.bounds_theta_southpole(1)); bounds_latlong.I_jump(bounds_latlong.bounds_theta_jump(:,2) >= bounds_latlong_ico.bounds_theta_southpole(1)); bounds_latlong.I_southpole];
candidate_faces{bounds_latlong_ico.I_southpole} = Icand;


%% initialize outputs and helper arrays

% sampled_pos = zeros(length(f_dual), 3);
% sampled_nor = zeros(length(f_dual), 3);
sampled_area = zeros(length(f_dual), 1);

f_normals = cross_fast(c_verts(c_faces_filtered(:,2),:)-c_verts(c_faces_filtered(:,1),:),c_verts(c_faces_filtered(:,3),:)-c_verts(c_faces_filtered(:,1),:));
f_normals = f_normals ./ sqrt(sum(f_normals.^2,2));

% inv_weights = spalloc(length(c_verts), length(f_dual),12*length(c_verts));
weights = spalloc(length(f_dual), length(c_verts), 12*length(f_dual));

%% do icosahedral sampling

for f_idx = 1:length(f_dual)
    
    v_ico1 = v_dual(f_dual{f_idx},:);
    normal_ico1 = cross_fast(v_ico1, v_ico1([2:end,1],:));
    
    Icand = candidate_faces{f_idx};
    for f_cand_idx = 1:size(Icand,1)
        %f_cand_idx = 12; %7 %12
        vert_indices = c_faces_filtered(Icand(f_cand_idx),:);
        v_edges1 = c_sph_verts(vert_indices,:);
        
        % get clipped polygon
        cl_points = clip_convex_polygon(v_edges1, normal_ico1);
        
        if ~isempty(cl_points)
            % barycentric coordinates
            %[a1, a2, a3] = compute_spherical_barycentric(v_edges1(1,:), v_edges1(2,:), v_edges1(3,:), clipped_points(1:points_cnt,:));
            [a1, a2, a3] = compute_spherical_barycentric_geom(v_edges1(1,:), v_edges1(2,:), v_edges1(3,:), cl_points);
            
            % flat 3d polygon within the mesh triangle
            PV = c_verts(vert_indices(1),:) .* a1 + c_verts(vert_indices(2),:) .* a2 + c_verts(vert_indices(3),:) .* a3;
            
            %polygon centroid and area
            %centroid = sum(PV,1) / size(PV,1);
            area = polygon_area( PV , f_normals(Icand(f_cand_idx),:));
            
            % weights for inverse mapping
            % inv_weights(vert_indices, f_idx) = inv_weights(vert_indices, f_idx) + (area * sum([a1, a2, a3],1).')/ size(PV,1);
            weights(f_idx, vert_indices) = weights(f_idx, vert_indices) + (area * sum([a1, a2, a3],1))/ size(PV,1);
           
            % add to geometry image pixel
            %sampled_pos(f_idx, :) =  sampled_pos(f_idx, :) + centroid*area;
            %sampled_nor(f_idx, :) =  sampled_nor(f_idx, :) + c_face_normals(Icand(f_cand_idx),:)*area;
            sampled_area(f_idx) = sampled_area(f_idx) + area;
        end
        
    end
    
end

%% normalize sampled signals

%sampled_nor = sampled_nor ./ sampled_area;
%sampled_pos = sampled_pos ./ sampled_area;

% compose sampled signals
weights = weights ./ sampled_area;
[idx_sampling, idx_orig, sparse_weights] = find(weights);
sparse_indices = [idx_sampling, idx_orig];

sampled_pos = weights * c_verts;

%sum(abs(sampled_pos - sp2))
%max(abs(sampled_pos - sp2))

% %% compute inverse mapping offsets
% 
% % find vertices that are not yet assigned a weight
% idx = find(full(sum(inv_weights ~= 0, 2)) == 0);
% v_tmp = c_sph_verts(idx, :);
% v_tmp = reshape(v_tmp, [size(v_tmp, 1) 1 3]);
% 
% %find closest vertex on icosahedron and assign
% v_ico_tmp = reshape(v_ico, [1 size(v_ico, 1) 3]);
% dist = sum((v_tmp - v_ico_tmp).^2,3);
% [~, tmp_inv_idx] = min(dist, [], 2);
% full_idx = (tmp_inv_idx-1) * size(inv_weights,1) + idx;
% inv_weights(full_idx) = 1.0;
% 
% cnt_ico_vertices = full(sum(inv_weights>0,2));
% uniform_inv_weights = inv_weights > 0;
% uniform_inv_weights = uniform_inv_weights .* (1 ./cnt_ico_vertices);
% [inv_idx_orig, inv_idx_sampling, inv_sparse_weights] = find(uniform_inv_weights);
% inv_offset_exp = sampled_pos(inv_idx_sampling, :) - c_verts(inv_idx_orig, :);
% ox = sparse(inv_idx_orig, inv_idx_sampling, inv_offset_exp(:,1), size(c_sph_verts, 1), size(v_ico, 1));
% oy = sparse(inv_idx_orig, inv_idx_sampling, inv_offset_exp(:,2), size(c_sph_verts, 1), size(v_ico, 1));
% oz = sparse(inv_idx_orig, inv_idx_sampling, inv_offset_exp(:,3), size(c_sph_verts, 1), size(v_ico, 1));
% x = full(sum(uniform_inv_weights .* ox, 2));
% y = full(sum(uniform_inv_weights .* oy, 2));
% z = full(sum(uniform_inv_weights .* oz, 2));
% 
% inv_offset = cat(2, x, y, z);
% inv_sparse_indices = [inv_idx_orig, inv_idx_sampling];
% 
% %inv_weights_normalized = inv_weights ./ sum(inv_weights,1);
% %[inv_idx_orig, inv_idx_sampling] = find(inv_weights_normalized);
% %inv_offset = sampled_pos(inv_idx_sampling, :) - c_verts(inv_idx_orig, :);
% % [inv_idx_orig, inv_idx_sampling, inv_factor] = find(inv_weights_normalized);
% % inv_offset = sampled_pos(inv_idx_sampling, :) - c_verts(inv_idx_orig, :) .* inv_factor;


end


function [a1, a2, a3] = compute_spherical_barycentric_geom (p1, p2, p3, cand)
% triangle normal
%normal = cross_fast(p2-p1,p3-p1);

%project cands to triangle spanned by p1,p2,p3
%t = normal * (p1 - cand)' ./ (normal * cand');
%cand_proj = cand + t' .* cand;

%comp areas
p1c = cand - p1;
p2c = cand - p2;
p3c = cand - p3;

a1 = 0.5 * sqrt(sum(cross_fast(p2c, p3c).^2, 2));
a2 = 0.5 * sqrt(sum(cross_fast(p1c, p3c).^2, 2));
a3 = 0.5 * sqrt(sum(cross_fast(p1c, p2c).^2, 2));

a = a1+a2+a3;
a1 = a1 ./a;
a2 = a2 ./a;
a3 = a3 ./a;

%     close all;
%     figure;
%     hold on;
%     grid on;
%
%     % draw axis
%     plot3([0,0],[0,0],[-2,2],'-b')
%
%     v = [p1;p2;p3];
%     f = [1,2,3];
%
%     options.face_vertex_color = [0.5 0.5 0.5];
%     plot_mesh(v,f,options);
%     shading faceted;
%
%     plot3( ...
%         [zeros(size(cand,1),1) cand(:,1)]', ...
%         [zeros(size(cand,1),1) cand(:,2)]', ...
%         [zeros(size(cand,1),1) cand(:,3)]', '-m');
%
%     scatter3(cand_proj(:,1),cand_proj(:,2),cand_proj(:,3),'r.');
%
%     axis equal;
%
%     ax = gca;               % get the current axis
%     ax.Clipping = 'off';
end

