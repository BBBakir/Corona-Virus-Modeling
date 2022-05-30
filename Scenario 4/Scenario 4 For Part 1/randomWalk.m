function [b] = randomWalk(a)
T =20;
p = rand;
n= randi([0,3]);
b = [0 0];



%%%% right top corner
if a(1) == T && a(2) == T
    if p < 1/3
        a(1) = a(1) - n;
    elseif p > 1/3 && p < 2/3
        a(2) = a(2) - n;
    elseif p > 2/3 
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    end
%%%left bottom corner
elseif a(1) == 1 && a(2) == 1
    if p < 1/3
        a(1) = a(1) + n;
    elseif p > 1/3 && p < 2/3
        a(2) = a(2) + n;
    elseif p > 2/3 
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    end
%%right bottom corner
elseif a(1) == T && a(2) == 1
    if p < 1/3
        a(1) = a(1) - n;
    elseif p > 1/3 && p < 2/3
        a(2) = a(2) + n;
    elseif p > 2/3 
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    end
%%%left top corner
elseif a(1) == 1 && a(2) == T
    if p < 1/3
        a(1) = a(1) + n;
    elseif p > 1/3 && p < 2/3
        a(2) = a(2) - n;
    elseif p > 2/3 
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    end
%%%%left edge
elseif a(1) == 1 && (a(2) < T && a(2) > 1)
    if p < 1/5
        a(1) = a(1) + n;
    elseif p > 1/5 && p < 2/5
        a(2) = a(2) + n;
    elseif p > 2/5  && p < 3/5
        a(2) = a(2) - n;
    elseif p > 3/5  && p < 4/5
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    elseif p > 4/5 
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    end
%% right edge
elseif a(1) == T && (a(2) < T && a(2) > 1)
    if p < 1/5
        a(1) = a(1) - n;
    elseif p > 1/5 && p < 2/5
        a(2) = a(2) + n;
    elseif p > 2/5  && p < 3/5
        a(2) = a(2) - n;
    elseif p > 3/5  && p < 4/5
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    elseif p > 4/5 
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    end

    %% UP edge
elseif a(2) == T && (a(1) < T && a(1) > 1)
    if p < 1/5
        a(1) = a(1) - n;
    elseif p > 1/5 && p < 2/5
        a(1) = a(1) + n;
    elseif p > 2/5  && p < 3/5
        a(2) = a(2) - n;
    elseif p > 3/5  && p < 4/5
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    elseif p > 4/5 
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    end

      %% BOTTOM edge
elseif a(2) == 1 && (a(1) < T && a(1) > 1)
    if p < 1/5
        a(1) = a(1) - n;
    elseif p > 1/5 && p < 2/5
        a(1) = a(1) + n;
    elseif p > 2/5  && p < 3/5
        a(2) = a(2) + n;
    elseif p > 3/5  && p < 4/5
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    elseif p > 4/5 
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    end

elseif a(1)>1 && a(1)< T && a(2)> 1 && a(2) <T
    if p<1/8 %left
        a(1) = a(1) - n;
    elseif p>=1/8 && p <2/8  % left & up
        if a(1)-1 < n
            n = a(1)-1;
        elseif T-a(2) < n
            n = T-a(2);
        end
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    elseif p >= 2/8 && p< 3/8 % up
        a(2) = a(2) +n;
    elseif  p >= 3/8 && p< 4/8 % up & right
        if T-a(1) < n
            n = T-a(1);
        elseif T-a(2) < n
            n = T-a(2);
        end
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    elseif p >= 4/8 && p< 5/8 % right
        a(1) = a(1) + n;
    elseif  p >= 5/8 && p< 6/8 % down & right
        if a(2)-1 < n
            n = a(2)-1;
        elseif T-a(1) < n
            n = T-a(1);
        end
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    elseif p >= 6/8 && p< 7/8 % down
        a(1) = a(1) + n;
    else % down & left
        if a(1)-1 < n
            n = a(1)-1;
        elseif a(2)-1 < n
            n = a(2)-1;
        end
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    end
end

if a(1) > T
    a(1) = T ;
elseif a(1)< 1 
    a(1) = 1;
end
if a(2) > T
    a(2) = T;
elseif a(2)< 1 
    a(2) = 1;
end
b = a;
end
