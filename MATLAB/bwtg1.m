function [ wave ] = bwtg1( signal )
%BWTG1 Summary of this function goes here
%   Detailed explanation goes here
%check if signal can be split it two parts of even length
m = (size(signal,2) - 4*ceil(size(signal,2)/4));
if m>0
    signal = [signal signal(:,1:m)];
end
%Lazy Wavelet Transform
[sO sE] = lazyw(signal,'c');
%Predict
E = bitxor(sO,sE);
%update*
O = sO;
%return wavelet
wave = [O,E];
end

