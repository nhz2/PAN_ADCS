function [ vx ] = crs( v )
%crs cross product matrix
%   Gives the skew symmetric matrix for vector vused for calculating cross 
%   product in matrix form

if (length(v) ~= 3)
    disp('vector not length 3')
end

vx = [ 0  -v(3) v(2);
      v(3)  0   -v(1);
      -v(2) v(1)  0];

end

