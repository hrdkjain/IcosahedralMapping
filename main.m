%% Title: IcosahedralMapping
% Date: 2020
% Link: https://github.com/hrdkjain/IcosahedralMapping.git

%% Main script to create planar icosahedron image from the mesh and its spherical parameterization

if exist('../Matlab-Functions', 'dir') ~= 7
    fprintf("Error: Matlab-Functions not found, make sure Matlab-Functions folder is placed in same directory as IcosahedralMapping");
    return
end

addpath('../Matlab-Functions', '../Matlab-Functions/MeshUtils', '../Matlab-Functions/natsort', '../Matlab-Functions/ProgressBar')
warning('off')
include  % header file include.m

C = textscan(fopen(listFile,'r'),"%s %s");
inFileList1 = C{1};
inFileList2 = C{2};

inFileListLen = length(inFileList1);
hbar = ParforProgressbar(inFileListLen, 'title', listFile, 'parpool', {'local', numcores});
parfor i=1:inFileListLen
    hbar.increment();
    inFilePath1 = fullfile(inFileList1{i});
    inFilePath2 = fullfile(inFileList2{i});
    [inFileFld2,inFile2,~] = fileparts(inFilePath2);
    outIcoMesh = fullfile(strrep(inFileFld2, inFld2, outIcoFld),strcat(inFile2,'_I',num2str(subdiv),'.off'));
    outIcoMat = strrep(outIcoMesh, '.off', '.mat');
    if (saveIcoMesh && exist(outIcoMesh,'file')==2)
        continue
    end
    if (saveIcoMat && exist(outIcoMat,'file')==2)
        continue
    end
    
    [outFld,~,~] = fileparts(outIcoMesh);
    if ~exist(outFld,'dir')
        mkdir(outFld)
    end
    
    try
        [v,f] = read_mesh(inFilePath1);
        assert(sum(sum(isnan(v))) == 0,'Input is NaN');
        [map,f2] = read_mesh(inFilePath2);
        f = f(1:length(f2),:);
        assert(sum(sum(isnan(map))) == 0,'Input is NaN');
        assert(length(v)==length(map),'Size of vertices mismatch');
        assert(length(f)==length(f2),'Size of faces mismatch');
        fn = face_normals(v,f);
        if orientMap
            map = orient_map(map, v, f, inFilePath1);
        end
        [ v_ico, f_ico, planar_mapping, sampled_pos, sparse_indices, sparse_weights] = icosahedral_sampling( subdiv, f, map, v, fn);
        if saveIcoMesh
            write_off(outIcoMesh, sampled_pos, f_ico);
        end
        if saveIcoMesh
            weights_save(outIcoMat, sparse_indices, sparse_weights);
        end
        if viewMesh
            %plotting Ico Mesh
            figure;
            patch('Faces',f_ico,'Vertices',sampled_pos,'FaceColor','flat','FaceVertexCData',[0.5,0.5,0.5], 'LineWidth',0.5);
            shading faceted;
            lighting none
            ax = gca;
            ax.Clipping = 'off';
            axis equal
        end
    catch err
        printError(err);
        continue
    end
end
close(hbar);

function weights_save(file, sparse_indices, sparse_weights)
    save(file,'sparse_indices', 'sparse_weights','-v6');
end
