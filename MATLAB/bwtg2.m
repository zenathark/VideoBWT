function [ wave ] = bwtg2( signal )
%BWTG2 Summary of this function goes here
%   Detailed explanation goes here
m = (size(signal,2) - 4*ceil(size(signal,2)/4));
if m>0
    signal = [signal signal(:,1:m)];
end
%Lazy Wavelet Transform
[sO sE] = lazyw(signal,'c');
%Update
E = sE;
%Predict
O = bitxor(sO,E);
%Update
E(:,1:size(sE,2)-1) = E(:,2:size(sE,2));
E(:,size(sE,2)) = sE(:,1);
%Update
E = bitxor(O,E);
%return wavelet

wave = [E,O];
end

