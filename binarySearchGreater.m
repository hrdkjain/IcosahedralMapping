function [index] = binarySearchGreater(A, num)
left = 1;
right = length(A);
index = 0;
while left <= right
    mid = ceil((left + right) / 2);
    
    if A(mid) == num
        index = mid;
        while index < length(A) && A(index+1) == num
           index = index + 1; 
        end
        break;
    else if A(mid) > num
        right = mid - 1;       
        else
            left = mid + 1;
            index = mid;            
        end
    end
end
end