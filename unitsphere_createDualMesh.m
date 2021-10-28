function [ v_out, f_out ] = unitsphere_createDualMesh( v, f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% First pass: create vertices
v_out = zeros(length(f), 3);
for i=1:length(f)
    v_out(i,:) = sum(v(f{i},:),1)/size(v(f{i},:),1);
end

% enlarge so vertices are on unit sphere
v_out = v_out ./ sqrt(sum(v_out.^2,2));

%Second pass: create faces
f_out = cell(1, size(v,1));
for i=1:size(f,1)
    tmp = f{i};
    for j=1:length(tmp)
       f_out{tmp(j)} = [f_out{tmp(j)} i]; 
    end
end

wrapN = @(x, N) (1 + mod(x-1, N));

% Third pass: order face indices counter clockwise
for i=1:length(f_out)
    ind = f_out{i};

    vals = zeros(length(ind),2);
    for j=1:length(ind)
        x = f{ind(j)};
        index = find(x==i);
        vals(j,:) = [x(wrapN(index-1,length(x))) x(wrapN(index+1,length(x)))];
    end

    ind_ordered = ind(1);
    target = vals(1,1);    
    while length(ind_ordered) < length(ind)
        tmp = find(vals(:,2)==target);
        ind_ordered = [ind_ordered ind(tmp)];
        target = vals(tmp,1);
    end
    
    f_out(i) = {ind_ordered};   
end
 
end


function [vals] = get_pred_and_suc(x, target)
    wrapN = @(x, N) (1 + mod(x-1, N));
    i = find(x==target);
    vals = [x(wrapN(i-1,length(x))) x(wrapN(i+1,length(x)))];
end
