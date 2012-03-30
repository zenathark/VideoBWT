function [ sch ] = splitsubb( mtx )
%SPLITSUBB Summary of this function goes here
%   Detailed explanation goes here
%sch = (zeros([8 size(mtx)])); 
for i=7:-1:0
    sch(i+1,:,:) = bitand(mtx,2^i) > 0;
end
end

