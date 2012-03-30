function [ img ] = fusesubch( sch )
%FUSESUBCH Summary of this function goes here
%   Detailed explanation goes here
img = zeros(size(sch,2), size(sch,3));
for l=8:-1:1
    plane(:,:) = sch(l,:,:) * 2^(l-1);
    img = bitor(img,plane);
end

end

