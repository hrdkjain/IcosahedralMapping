function [index] = binarySearchLess(A, num)
left = 1;
right = length(A);
index = length(A)+1;
while left <= right
    mid = ceil((left + right) / 2);
    
    if A(mid) == num
        index = mid;
        while index > 1 && A(index-1) == num
           index = index -1; 
        end
        break;
    else if A(mid) > num
        right = mid - 1;
        index = mid;
        else
            left = mid + 1;
        end
    end
end
end


