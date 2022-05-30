function [c] = isoWalk(a)
T =20;
p = rand;

n = 1;
c = [0 0];

%%%% right top corner
if a(1) == T && a(2) == T
    if p < 1/4
        a(1) = a(1) - n;
    elseif p > 1/4 && p < 2/4
        a(2) = a(2) - n;
    elseif p < 3/4 && p >= 2/4 
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end
%%%left bottom corner
elseif a(1) == 1 && a(2) == 1
    if p < 1/4
        a(1) = a(1) + n;
    elseif p > 1/4 && p < 2/4
        a(2) = a(2) + n;
    elseif p < 3/4 && p >= 2/4 
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end
%%right bottom corner
elseif a(1) == T && a(2) == 1
    if p < 1/4
        a(1) = a(1) - n;
    elseif p > 1/4 && p < 2/4
        a(2) = a(2) + n;
    elseif p < 3/4 && p >= 2/4  
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end
%%%left top corner
elseif a(1) == 1 && a(2) == T
    if p < 1/4
        a(1) = a(1) + n;
    elseif p > 1/4 && p < 2/4
        a(2) = a(2) - n;
    elseif p < 3/4 && p >= 2/4 
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end
%%%%left edge
elseif a(1) == 1 && (a(2) < T && a(2) > 1)
    if p < 1/6
        a(1) = a(1) + n;
    elseif p > 1/6 && p < 2/6
        a(2) = a(2) + n;
    elseif p > 2/6  && p < 3/6
        a(2) = a(2) - n;
    elseif p > 3/6  && p < 4/6
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    elseif p > 4/6 && p < 5/6
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end
%% right edge
elseif a(1) == T && (a(2) < T && a(2) > 1)
    if p < 1/6
        a(1) = a(1) - n;
    elseif p > 1/6 && p < 2/6
        a(2) = a(2) + n;
    elseif p > 2/6  && p < 3/6
        a(2) = a(2) - n;
    elseif p > 3/6  && p < 4/6
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    elseif p > 4/6 && p < 5/6 
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end

    %% UP edge
elseif a(2) == T && (a(1) < T && a(1) > 1)
    if p < 1/6
        a(1) = a(1) - n;
    elseif p > 1/6 && p < 2/6
        a(1) = a(1) + n;
    elseif p > 2/6  && p < 3/6
        a(2) = a(2) - n;
    elseif p > 3/6  && p < 4/6
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    elseif p > 4/6 && p < 5/6
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end

      %% BOTTOM edge
elseif a(2) == 1 && (a(1) < T && a(1) > 1)
    if p < 1/6
        a(1) = a(1) - n;
    elseif p > 1/6 && p < 2/6
        a(1) = a(1) + n;
    elseif p > 2/6  && p < 3/6
        a(2) = a(2) + n;
    elseif p > 3/6  && p < 4/6
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    elseif p > 4/6 && p < 5/6 
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    else
        a(1) = a(1);
        a(2) = a(2);
    end

elseif a(1)>1 && a(1)< T && a(2)> 1 && a(2) <T
    if p< 1/9 %left
        a(1) = a(1) - n;
    elseif p>=1/9 && p <2/9  % left & up
        a(1) = a(1) - n;
        a(2) = a(2) + n;
    elseif p >= 2/9 && p< 3/9 % up
        a(2) = a(2) +n;
    elseif  p >= 3/9 && p< 4/9 % up & right
        a(1) = a(1) + n;
        a(2) = a(2) + n;
    elseif p >= 4/9 && p< 5/9 % right
        a(1) = a(1) + n;
    elseif  p >= 5/9 && p< 6/9 % down & right
        a(1) = a(1) + n;
        a(2) = a(2) - n;
    elseif p >= 6/9 && p< 7/9 % down
        a(1) = a(1) + n;
    elseif p >= 7/9 && p< 8/9 % same
        a(1) = a(1);
        a(2) = a(2);
    else % down & left
        a(1) = a(1) - n;
        a(2) = a(2) - n;
    end
end
c =a;
end