function xround = round_HAH( x, d )
%% ROUND_HAH is ruounding decimals
% Example:  round_HAH(0.8236,1) --> 1
%           round_HAH(0.8236,0.1) --> 0.8
%           round_HAH(0.8236,0.01) --> 0.82
    if d <= 0
        error('ERROR: Rounding precision must be > 0')
    end
    xround = round(x/d)*d;
end

