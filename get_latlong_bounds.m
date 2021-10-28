function [ s ] = get_latlong_bounds( f, v )

%% spherical coordinates
[theta, phi] = s2.car2sph(v(:,1),v(:,2),v(:,3));

%% latitude and longitude bounds
if ~iscell(f)
    bounds_phi = [min(phi(f),[],2) max(phi(f),[],2)];
    [theta_min, theta_max] = get_bounding_latitude(theta(f), ...
        reshape(v(f,1),[],3), reshape(v(f,2),[],3), reshape(v(f,3),[],3));
    bounds_theta = [theta_min theta_max];
else
    bounds_theta = zeros(length(f),2);
    bounds_phi = zeros(length(f),2);
    for i=1:length(f)
        bounds_phi(i,:) = [min(phi(f{i})) max(phi(f{i}))];
        [theta_min, theta_max] = get_bounding_latitude(theta(f{i})', ...
            v(f{i}, 1)', v(f{i}, 2)', v(f{i}, 3)' );
        bounds_theta(i,:) = [theta_min, theta_max];
    end
end


%% classify faces
% find faces on longitude discontinuity
is_on_jump = (bounds_phi(:,2) - bounds_phi(:,1)) > pi;
I_normal = find(~is_on_jump);
I_jump = find(is_on_jump);


% find faces covering the poles
poles = [0 0 1; 0 0 -1];
if ~iscell(f)
    for i=1:length(I_jump)
        inside = points_in_spherical_polygon(v(f(I_jump(i), :),:), poles);
        if inside(1) == true
            I_northpole = I_jump(i);
        elseif inside(2) == true
            I_southpole = I_jump(i);
        end
    end
else
    for i=1:length(I_jump)
        inside = points_in_spherical_polygon(v(f{I_jump(i)},:), poles);
        if inside(1) == true
            I_northpole = I_jump(i);
        elseif inside(2) == true
            I_southpole = I_jump(i);
        end
    end
end

%delete poles from jump
I_jump(I_jump==I_northpole | I_jump == I_southpole) = [];


%% adjust bounds according to face 'type' (normal, jump, pole)
% normal faces
bounds_phi_normal = bounds_phi(I_normal,:);
bounds_theta_normal = bounds_theta(I_normal,:);

% sort bounds for speedup
[bounds_theta_normal_min_sorted,I_theta_normal_min_sorted] = sort(bounds_theta_normal(:,1),1);
[bounds_theta_normal_max_sorted,I_theta_normal_max_sorted] = sort(bounds_theta_normal(:,2),1);

% faces on jump
bounds_theta_jump = bounds_theta(I_jump,:);
if ~iscell(f)
    jump_phi_min = phi(f(I_jump, :));
    jump_phi_max = jump_phi_min;
    jump_phi_min(jump_phi_min < pi) = inf;
    jump_phi_max(jump_phi_max >= pi) = -inf;
    bounds_phi_jump = [min(jump_phi_min,[],2) max(jump_phi_max,[],2)];
else
    jump_phi = cellfun(@(x) phi(x),f(I_jump)', 'UniformOutput', false);
    bounds_phi_jump = cell2mat(cellfun(@(x) [min(x(x>=pi)), max(x(x<pi))], jump_phi, 'UniformOutput', false));
end

% faces on poles
bounds_phi_northpole = [0 2*pi];
bounds_theta_northpole = [0 bounds_theta(I_northpole, 2)];
bounds_phi_southpole = [0 2*pi];
bounds_theta_southpole = [bounds_theta(I_southpole, 1) pi];

%% create output

% normal faces
s.I_normal = I_normal;
s.bounds_theta_normal = bounds_theta_normal;
s.bounds_phi_normal = bounds_phi_normal;

% sorted bounds for speedup
s.I_theta_normal_min_sorted = I_theta_normal_min_sorted;
s.bounds_theta_normal_min_sorted = bounds_theta_normal_min_sorted;
s.I_theta_normal_max_sorted = I_theta_normal_max_sorted;
s.bounds_theta_normal_max_sorted = bounds_theta_normal_max_sorted;

%faces on jump
s.I_jump = I_jump;
s.bounds_theta_jump = bounds_theta_jump;
s.bounds_phi_jump = bounds_phi_jump;

% faces on poles
s.I_northpole = I_northpole;
s.bounds_theta_northpole = bounds_theta_northpole;
s.bounds_phi_northpole = bounds_phi_northpole;
s.I_southpole = I_southpole;
s.bounds_theta_southpole = bounds_theta_southpole;
s.bounds_phi_southpole = bounds_phi_southpole;




function [theta_min, theta_max] = get_bounding_latitude(theta, x, y, z)
    
    % prepare great circle arc points
    v1 = cat(3,x,y,z);
    v2 = cat(3,x(:,[2:end,1]), y(:,[2:end,1]), z(:,[2:end,1]));
    
    % angle of great circle arc
    angle = dot(v1,v2,3);
    angle(angle>1) = 1;
    angle(angle<-1) = -1;
    angle = acos(angle);
    
    % point at 90 degree from v1 in direction of v2
    v3 = (v2-v1.*cos(angle)) ./ sin(angle);
    v3 =real(v3);
    % angle of northmost/southmost point of great circle
    angle_north = atan2(v3(:,:,3), v1(:,:,3));
    angle_south = angle_north + pi;
    angle_south(angle_south>pi) = angle_south(angle_south>pi) - 2*pi;
    
    % is extreme point within arc?
    In = find(sign(angle) == sign(angle_north) & abs(angle_north) < abs(angle));
    Is = find(sign(angle) == sign(angle_south) & abs(angle_south) < abs(angle));
    
    % theta min (north)
    v1z = v1(:,:,3);
    v3z = v3(:,:,3);
    z_north = v1z(In) .* cos(angle_north(In)) + v3z(In) .* sin(angle_north(In));  
    theta_arc_north = pi * ones([size(v1,1) size(v1,2)]);
    theta_arc_north(In) = acos(z_north);
    theta_min = min([theta theta_arc_north],[],2);
    
    % theta max (south)
    z_south = v1z(Is) .* cos(angle_south(Is)) + v3z(Is) .* sin(angle_south(Is));  
    theta_arc_south = zeros([size(v1,1) size(v1,2)]);
    theta_arc_south(Is) = acos(z_south);
    theta_max = max([theta theta_arc_south],[],2);
    
%     % northmost/southmost point of great circle
%     v_north = v1 .* cos(angle_north) + v3 .* sin(angle_north);    
%     v_south = v1 .* cos(angle_south) + v3 .* sin(angle_south);    
%     
%     f_idx = 18965;
%     v_idx = 2;
%     %figure;
%     hold on;
%     grid on;
%     
%     % draw axis
%     plot3([0,0],[0,0],[-2,2],'-b')
%     
%      %plot equator
%     patch('XData', [2 -2 -2 2],'YData', [2 2 -2 -2],'ZData', [0 0 0 0], 'FaceVertexCData', [0,1,0], 'FaceAlpha', 0.2);
%     
%     %plot face
%     patch('XData', v1(f_idx,:,1),'YData', v1(f_idx,:,2),'ZData', v1(f_idx,:,3), 'FaceVertexCData', [0,0,1], 'FaceAlpha', 1.0);
%     
%     %plot points
%     scatter3(v1(f_idx,v_idx,1), v1(f_idx,v_idx,2), v1(f_idx,v_idx,3))
%     text(v1(f_idx,v_idx,1), v1(f_idx,v_idx,2), v1(f_idx,v_idx,3),'v1')
%     scatter3(v2(f_idx,v_idx,1), v2(f_idx,v_idx,2), v2(f_idx,v_idx,3))
%     text(v2(f_idx,v_idx,1), v2(f_idx,v_idx,2), v2(f_idx,v_idx,3),'v2')
%     scatter3(v3(f_idx,v_idx,1), v3(f_idx,v_idx,2), v3(f_idx,v_idx,3))
%     text(v3(f_idx,v_idx,1), v3(f_idx,v_idx,2), v3(f_idx,v_idx,3),'v3')
%     
%     scatter3(v_north(f_idx,v_idx,1), v_north(f_idx,v_idx,2), v_north(f_idx,v_idx,3))
%     text(v_north(f_idx,v_idx,1), v_north(f_idx,v_idx,2), v_north(f_idx,v_idx,3),'v_north')
%     
%     scatter3(v_south(f_idx,v_idx,1), v_south(f_idx,v_idx,2), v_south(f_idx,v_idx,3))
%     text(v_south(f_idx,v_idx,1), v_south(f_idx,v_idx,2), v_south(f_idx,v_idx,3),'v_south')
%     
%     %scatter3(reshape(v_north(:,:,1),[],1), reshape(v_north(:,:,2),[],1), reshape(v_north(:,:,3),[],1),'k.')
%     %scatter3(reshape(v_south(:,:,1),[],1), reshape(v_south(:,:,2),[],1), reshape(v_south(:,:,3),[],1),'y.')
%    
%     
%     %plot great circle
%     i=0:0.001:2*pi;
%     i = [i 2*pi];
%     fi = reshape(v1(f_idx,v_idx,:),1,[]) .* cos(i)' + reshape(v3(f_idx,v_idx,:),1,[]) .* sin(i)';
%     plot3(fi(:,1),fi(:,2),fi(:,3),'-r')
%     
%     % plot rue north and south as reference
%     [~,max_idx] = max(fi(:,3));
%     scatter3(fi(max_idx,1), fi(max_idx,2), fi(max_idx,3));
%     text(fi(max_idx,1), fi(max_idx,2), fi(max_idx,3),'v_north_true');
%     
%     [~,min_idx] = min(fi(:,3))
%     scatter3(fi(min_idx,1), fi(min_idx,2), fi(min_idx,3));
%     text(fi(min_idx,1), fi(min_idx,2), fi(min_idx,3),'v_south_true');
%     
%     shading faceted;    
%     ax = gca;               % get the current axis
%     ax.Clipping = 'off';    % turn clipping off   
    
end

    function inside = points_in_spherical_polygon(polygon_points, target_points)
        normals = cross_fast(polygon_points, polygon_points([2:end,1],:));
        dist = target_points * normals';
        inside = sum(dist<0,2) == 0;
    end




end

