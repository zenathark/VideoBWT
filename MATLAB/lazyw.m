function [ O E ] = lazyw( signal, type )
%LAZYW Summary of this function goes here
%   Detailed explanation goes here
O = dyaddown(signal,1,type);
E = dyaddown(signal,2,type);
end

