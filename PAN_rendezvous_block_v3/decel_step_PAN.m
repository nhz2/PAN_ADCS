function [da_step,d_step] = decel_step_PAN(along_track_sep,follower_ahead)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (along_track_sep > 250000)
    da_step = 50;
    d_step = 100000;
elseif (along_track_sep > 100000)
    da_step = 25;
    d_step = 50000;
elseif (along_track_sep > 50000)
    da_step = 10;
    d_step = 20000;
elseif (along_track_sep > 20000)
    da_step = 10;
    d_step = 20000;
elseif (along_track_sep > 2000)
    da_step = 5;
    d_step = 2000;
else
    da_step = 1;
    d_step = 0;
end

if (~follower_ahead)
    da_step = -da_step;
end

end

