function [orb_matched,dr_maxs] = orb_matched_dr(x_target,x_current,muE,nsteps,rtols,rntol)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[deltars,deltarns] = delta_r_data_hEe(x_target,x_current,muE,nsteps);

deltarx_max = max(deltars(1,:));
deltary_max = max(deltars(2,:));
deltarz_max = max(deltars(3,:));
deltarn_max = max(deltarns);
dr_maxs = [deltarx_max;deltary_max;deltarz_max;deltarn_max];

if ((deltarx_max > rtols(1)) || (deltary_max > rtols(2)) || (deltarz_max > rtols(3)))
    rmatched = 0;
else
    rmatched = 1;
end

if (deltarn_max < rntol)
    rnmatched = 1;
else
    rnmatched = 0;
end

if (rmatched && rnmatched)
    orb_matched=1;
else
    orb_matched=0;
end

end
