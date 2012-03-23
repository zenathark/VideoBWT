function [ wave ] = bwtg1( signal )
%BWTG1 Summary of this function goes here
%   Detailed explanation goes here
m = (size(signal,2) - 4*ceil(size(signal,2),4));
if m>0
    signal = [signal signal(:,1:m)];
end
[sO sE] = lazyw(signal,'c');
E = bitxor(sO,sE);
O = sO;
wave = [O,E];
end

