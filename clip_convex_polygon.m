function [ clipped_points ] = clip_convex_polygon(polygon_points, clip_normals)
        
    points_cnt = size(polygon_points,1);  
    clip_planes_cnt = size(clip_normals,1);
    total_points_cnt = points_cnt;

    %allocate points
    clipped_indices = 1:points_cnt;    
    tmp_points = polygon_points;     
    
    for i=1:clip_planes_cnt
        point_indices = clipped_indices(1:points_cnt);            
        last_points_cnt = points_cnt;
        points_cnt = 0;        

        % points within clipping plane
        dist = clip_normals(i,:) * tmp_points(point_indices,:)';
        
        prev = last_points_cnt;
        for j=1:last_points_cnt
            if dist(j) >= 0
                if dist(prev) < 0
                    % determine and store intersection point
                    d = dist(j) / (dist(j) - dist(prev));                    
                    total_points_cnt = total_points_cnt + 1;
                    tmp_points(total_points_cnt,:) = tmp_points(point_indices(j),:) + d * (tmp_points(point_indices(prev),:)- tmp_points(point_indices(j),:));
                    
                    points_cnt = points_cnt + 1;
                    clipped_indices(points_cnt) = total_points_cnt;                    
                end
                points_cnt = points_cnt + 1;
                clipped_indices(points_cnt) = point_indices(j);                    
            elseif dist(prev) >= 0
                % determine and store intersection point
                d = dist(j) / (dist(j) - dist(prev));                
                total_points_cnt = total_points_cnt + 1;
                tmp_points(total_points_cnt,:) = tmp_points(point_indices(j),:) + d * (tmp_points(point_indices(prev),:)- tmp_points(point_indices(j),:));
                    
                points_cnt = points_cnt + 1;
                clipped_indices(points_cnt) = total_points_cnt;                
            end 
            prev = j;
        end
        if points_cnt == 0
            break;
        end
    end
    
    if points_cnt > 0
        clipped_points = tmp_points(clipped_indices(1:points_cnt),:);
    else
        clipped_points = [];
    end 
end

