function [ data_out ] = gray2mat_ccd( data_in )
%GRAY2MAT_CCD Summary of this function goes here
%   Detailed explanation goes here
     minn_data = min(data_in(:));
     maxn_data = max(data_in(:));
     data_out = data_in*(maxn_data-minn_data)+minn_data;
end

