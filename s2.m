classdef s2

    methods(Static)
        
        
        function sample_grid = avg_sampling(grid, v_planar, f_planar, v_mesh) 
           
            
           selx = min(floor(P(1,:))):max(ceil(P(1,:)));
           sely = min(floor(P(2,:))):max(ceil(P(2,:))); 
        end
        
        function [v_out, f_out, dup_map, f_map] = splitTrianglesOnCut(v_in, f_in, I)
            % Splits triangles that are intersected by the unwrapping cut.
            % The triangles are dupicated in the plane above and below the
            % cut. Only works with triangles that cover less than half a
            % sphere.
            %
            % Input:
            % v: nv x 2 vertex coordinates on the projected plane
            % f: nf x 3 triangulations 
            
            
            %get vertices that will be duplicated           
            idx_dupl = unique(f_in(I,:));
            dupl_len =  max(size(idx_dupl));
            
            %create new vertices and append to list
            v_out = [v_in; zeros(dupl_len,2)];
            v_out_len = max(size(v_out));
            
            %create vertex mapping            
            v_in_len = max(size(v_in));            
            vertex_map = containers.Map(idx_dupl, (v_in_len+1):v_out_len);
            dup_map = [1:length(v_in) idx_dupl'];
            
            %get upper and lower vertices 
            idx_lowerHalf = idx_dupl(v_in(idx_dupl,2) < pi);
            idx_upperHalf = setdiff(idx_dupl, idx_lowerHalf);
            lower_map = containers.Map([idx_lowerHalf; idx_upperHalf], [cell2mat(vertex_map.values(num2cell(idx_lowerHalf))); idx_upperHalf]);
            upper_map = containers.Map([idx_upperHalf; idx_lowerHalf], [cell2mat(vertex_map.values(num2cell(idx_upperHalf))); idx_lowerHalf]);
            
            %create class mapping
            %class_map = [containers.Map(idx_lowerHalf, zeros(size(idx_lowerHalf))); containers.Map(idx_upperHalf, ones(size(idx_upperHalf)))]
            
            %create duplicated vertices
            v_out(cell2mat(vertex_map.values(num2cell(idx_lowerHalf))),:) = v_in(idx_lowerHalf,:) + [0,2*pi];
            v_out(cell2mat(vertex_map.values(num2cell(idx_upperHalf))),:) = v_in(idx_upperHalf,:) - [0,2*pi];
            
            %create new faces for with lower vhalf vertices
            f_out = [cell2mat(lower_map.values(num2cell(f_in(I,:))));cell2mat(upper_map.values(num2cell(f_in(I,:))))];             
            f_out = [f_in(setdiff(1:end,I),:); f_out];           
            
            %create face mapping 
            f_map = 1:size(f_in,1);
            f_map(I) = [];
            f_map(end+1:end+2*length(I)) = [I' I'];  
        end
        
        function [theta_s, phi_s] = scale2int(theta, phi)
            x = sin(theta) .* cos(phi);
            y = sin(theta) .* sin(phi);
            z = cos(theta);
        end
            
        function [x, y, z, theta, phi] = get_projection_grid(b, grid_type)
            %''' returns the spherical grid in euclidean
            %coordinates, where the sphere's center is moved
            %to (0, 0, 1)'''
            [theta, phi] = s2.meshgrid(b, grid_type);
            [x,y,z] = s2.sph2car(theta, phi);
        end
        
        function [v_out, f_out, theta, phi] = get_uv_mesh(b, grid_type)
            [theta, phi] = s2.meshgrid(b, grid_type);
            [x,y,z] = s2.sph2car(theta, phi);
            
            if strcmp('Driscoll-Healy', grid_type)
                % pole vertices are the same for all phi - yields triangles             
                % create vertices
                v_northpole = [x(1,1), y(1,1), z(1,1)];
                %v_southpole = [x(end,end), y(end,end), z(end,end)];
                
                % delete northpole
                x = x(:,2:end);
                y = y(:,2:end);
                z = z(:,2:end);
                
                % create all vertices                
                v_out = [ x(:), y(:), z(:); v_northpole]; %; v_southpole];
                lat_vertices = size(x,1);
                long_vertices = size(x,2);
                idx_northpole = size(v_out,1);%-1;
                %idx_southpole = size(v_out,1);
                
                % create faces
                f_out = {};
                
                % create north pole faces
                I = sub2ind([lat_vertices, long_vertices],1:lat_vertices,ones(1,lat_vertices));
                f_out = [f_out; num2cell([I' circshift(I',-1) ones(lat_vertices,1)*idx_northpole],2)];
                
                % create quads
                [A1,B1] = meshgrid(1:lat_vertices,1:long_vertices-1);
                [A2,B2] = meshgrid(1:lat_vertices,2:long_vertices);
                [A3,B3] = meshgrid([2:lat_vertices, 1],2:long_vertices);
                [A4,B4] = meshgrid([2:lat_vertices, 1],1:long_vertices-1);                
                f_out = [f_out; num2cell(sub2ind([lat_vertices, long_vertices],[A1(:) A2(:) A3(:) A4(:)],[B1(:) B2(:) B3(:) B4(:)]),2)];
                
                
                % create south pole face
                I = sub2ind([lat_vertices, long_vertices],1:lat_vertices, ones(1,lat_vertices)*long_vertices);
                f_out = [f_out; num2cell(reverse(I),2)];
                %I = sub2ind([lat_vertices, long_vertices],1:lat_vertices,ones(1,lat_vertices)*long_vertices);
                %f_out = [f_out; num2cell([I' ones(lat_vertices,1)*idx_southpole circshift(I',-1)], 2)];
            
            elseif strcmp('SOFT', grid_type)
                %no vertex at poles - polygon at poles
                
                % create all vertices                
                v_out = [ x(:), y(:), z(:)];
                lat_vertices = size(x,1);
                long_vertices = size(x,2);
                
                % create faces
                f_out = {};
                
                % north pole face
                I = sub2ind([lat_vertices, long_vertices],1:lat_vertices, ones(1,lat_vertices));
                f_out = [f_out; num2cell(I,2)];
                                
                % create quads
                [A1,B1] = meshgrid(1:lat_vertices,1:long_vertices-1);
                [A2,B2] = meshgrid(1:lat_vertices,2:long_vertices);
                [A3,B3] = meshgrid([2:lat_vertices, 1],2:long_vertices);
                [A4,B4] = meshgrid([2:lat_vertices, 1],1:long_vertices-1);                
                f_out = [f_out; num2cell(sub2ind([lat_vertices, long_vertices],[A1(:) A2(:) A3(:) A4(:)],[B1(:) B2(:) B3(:) B4(:)]),2)];
          
                                
                % create south pole face
                I = sub2ind([lat_vertices, long_vertices],1:lat_vertices, ones(1,lat_vertices)*long_vertices);
                f_out = [f_out; num2cell(reverse(I),2)];
            elseif strcmp('Clenshaw-Curtis', grid_type)
                % pole vertices are the same for all phi - yields triangles             
                % create vertices
                v_northpole = [x(1,1), y(1,1), z(1,1)];
                v_southpole = [x(end,end), y(end,end), z(end,end)];
                
                 % delete northpole and southpole
                x = x(:,2:end-1);
                y = y(:,2:end-1);
                z = z(:,2:end-1);
                
                % create all vertices                
                v_out = [ x(:), y(:), z(:); v_northpole; v_southpole];
                lat_vertices = size(x,1);
                long_vertices = size(x,2);
                idx_northpole = size(v_out,1)-1;
                idx_southpole = size(v_out,1);
                
                % create faces
                f_out = {};
                
                % create north pole faces
                I = sub2ind([lat_vertices, long_vertices],1:lat_vertices,ones(1,lat_vertices));
                f_out = [f_out; num2cell([I' circshift(I',-1) ones(lat_vertices,1)*idx_northpole],2)];
                
                % create quads
                [A1,B1] = meshgrid(1:lat_vertices,1:long_vertices-1);
                [A2,B2] = meshgrid(1:lat_vertices,2:long_vertices);
                [A3,B3] = meshgrid([2:lat_vertices, 1],2:long_vertices);
                [A4,B4] = meshgrid([2:lat_vertices, 1],1:long_vertices-1);                
                f_out = [f_out; num2cell(sub2ind([lat_vertices, long_vertices],[A1(:) A2(:) A3(:) A4(:)],[B1(:) B2(:) B3(:) B4(:)]),2)];
                                
                % create south pole face
                I = sub2ind([lat_vertices, long_vertices],1:lat_vertices,ones(1,lat_vertices)*long_vertices);
                f_out = [f_out; num2cell([I' ones(lat_vertices,1)*idx_southpole circshift(I',-1)], 2)];          
            end      
            
        end

        function [x,y,z] = sph2car(theta, phi)
            x = sin(theta) .* cos(phi);
            y = sin(theta) .* sin(phi);
            z = cos(theta);
            x(abs(x)<1e-10) = 0;
            y(abs(y)<1e-10) = 0;
            z(abs(z)<1e-10) = 0;
        end

        function [theta, phi] = car2sph(x,y,z)
            theta = acos(z);
            phi = atan2(y,x);
            I = find (phi<0);
            phi(I) = phi(I) + 2*pi;
        end

        function [theta, phi]  = meshgrid(b, grid_type)
            %"""
            %Create a coordinate grid for the 2-sphere.
            %There are various ways to setup a grid on the sphere.
            %if grid_type == 'Driscoll-Healy', we follow the grid_type from [4], which is also used in [5]:
            %beta_j = pi j / (2 b)     for j = 0, ..., 2b - 1
            %alpha_k = pi k / b           for k = 0, ..., 2b - 1
            %if grid_type == 'SOFT', we follow the grid_type from [1] and [6]
            %beta_j = pi (2 j + 1) / (4 b)   for j = 0, ..., 2b - 1
            %alpha_k = pi k / b                for k = 0, ..., 2b - 1
            %if grid_type == 'Clenshaw-Curtis', we use the Clenshaw-Curtis grid, as defined in [2] (section 6):
            %beta_j = j pi / (2b)     for j = 0, ..., 2b
            %alpha_k = k pi / (b + 1)    for k = 0, ..., 2b + 1
            %if grid_type == 'Gauss-Legendre', we use the Gauss-Legendre grid, as defined in [2] (section 6) and [7] (eq. 2):
            %beta_j = the Gauss-Legendre nodes    for j = 0, ..., b
            %alpha_k = k pi / (b + 1),               for k = 0, ..., 2 b + 1
            %if grid_type == 'HEALPix', we use the HEALPix grid, see [2] (section 6):
            %TODO
            %if grid_type == 'equidistribution', we use the equidistribution grid, as defined in [2] (section 6):
            %TODO
            %[1] SOFT: SO(3) Fourier Transforms
            %Kostelec, Peter J & Rockmore, Daniel N.
            %[2] Fast evaluation of quadrature formulae on the sphere
            %Jens Keiner, Daniel Potts
            %[3] A Fast Algorithm for Spherical Grid Rotations and its Application to Singular Quadrature
            %Zydrunas Gimbutas Shravan Veerapaneni
            %[4] Computing Fourier transforms and convolutions on the 2-sphere
            %Driscoll, JR & Healy, DM
            %[5] Engineering Applications of Noncommutative Harmonic Analysis
            %Chrikjian, G.S. & Kyatkin, A.B.
            %[6] FFTs for the 2-Sphere ? Improvements and Variations
            %Healy, D., Rockmore, D., Kostelec, P., Moore, S
            %[7] A Fast Algorithm for Spherical Grid Rotations and its Application to Singular Quadrature
            %Zydrunas Gimbutas, Shravan Veerapaneni
            %:param b: the bandwidth / resolution
            %:return: a meshgrid on S^2
            %"""
            [l_beta, l_alpha] = s2.linspace(b, grid_type);
            [theta, phi] = meshgrid(l_beta, l_alpha);
        end


        function [l_beta, l_alpha] = linspace(b, grid_type)
            if strcmp('Driscoll-Healy', grid_type)
                l_beta = (0:(2*b-1)) * pi / (2. * b);
                l_alpha = (0:(2*b-1))  * pi / b;
            elseif strcmp('SOFT', grid_type)
                l_beta = pi * (2 * (0:(2*b-1)) + 1) / (4. * b);
                l_alpha = (0:(2*b-1)) * pi / b;
            elseif strcmp('Clenshaw-Curtis', grid_type)
                % l_beta = np.arange(2 * b + 1) * np.pi / (2 * b);
                % l_alpha = np.arange(2 * b + 2) * np.pi / (b + 1);
                % Must use np.linspace to prevent numerical errors that cause beta > pi
                l_beta = linspace(0, pi, 2 * b + 1);
                l_alpha = linspace(0, 2 * pi, 2 * b + 3);
                l_alpha = l_alpha(1:end-1);
            elseif strcmp('Gauss-Legendre', grid_type)
                %x, _ = leggauss(b + 1)  % TODO: leggauss docs state that this may not be only stable for orders > 100
                %l_beta = np.arccos(x);
                %l_alpha = np.arange(2 * b + 2) * np.pi / (b + 1);
            elseif strcmp('HEALPix', grid_type)
                %TODO: implement this here so that we don't need the dependency on healpy / healpix_compat
                %from healpix_compat import healpy_sphere_meshgrid
                %return healpy_sphere_meshgrid(b)
            elseif strcmp('equidistribution', grid_type)
                %raise NotImplementedError('Not implemented yet; see Fast evaluation of quadrature formulae on the sphere.')
            else
                %raise ValueError('Unknown grid_type:' + grid_type)
            end
        end
    end

end