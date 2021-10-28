% Main script to list pair of mesh and spherically parameterized files

if exist('../Matlab-Functions', 'dir') ~= 7
    fprintf("Error: Matlab-Functions not found, make sure Matlab-Functions folder is placed in same directory as IcosahedralMapping");
    return
end

addpath('../Matlab-Functions', '../Matlab-Functions/natsort', 'Skeleton3D', 'Voxelization')
warning('off')
include  % header file include.m

list_fld1 = natsortfiles({subdir(fullfile(inFld1, ext1),flStr).name}');
if dual
    list_fld2 = natsortfiles({subdir(fullfile(inFld2, ext2),flStr).name}');
end

fid = fopen(listFile,'w+');
listed = 0;
for i=1:length(list_fld1)
    if dual
        idxAinB = AinListB(list_fld1{i},list_fld2,ext2);
        if idxAinB
            fprintf(fid,'%s %s\n',fullfile(list_fld1{i}), fullfile(list_fld2{idxAinB}));
            listed = listed+1;
        end
    else
        fprintf(fid,'%s\n',fullfile(list_fld1{i}));
        listed = listed+1;
    end
end
fprintf('Listed %d files and saved to %s\n', listed, listFile);
fclose(fid);

function idx = AinListB(A,listB,extB)
% Count the number of underscores in the two lists
nUnderscoreA = count(A,'_');
nUnderscoreB = count(listB{1},'_');

% Based on the number of underscores the type of comparison can be estimated
[~,A,~] = fileparts(A);
if nUnderscoreA == nUnderscoreB
    A = strcat(A,extB(2:end));
    idx = find(contains(listB,A));
else
    if nUnderscoreA > nUnderscoreB
        underscoreIdx = strfind(A,'_');
        A = A(1:underscoreIdx(nUnderscoreB+1)-1);
        idx = find(contains(listB,strcat(filesep,A,'.')));
    elseif  nUnderscoreA < nUnderscoreB
        idx = find(contains(listB,strcat(filesep,A,'_')));
    end
end

if size(idx,1)> 1
    error(strcat('more than one occurences of: ',A));
end

end
