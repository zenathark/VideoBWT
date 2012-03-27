function [ signal ] = ibwtg1( wave )
%IBWTG1 Summary of this function goes here
%   Detailed explanation goes here
%Lazy Wavelet Transform
s = size(wave,2);
sO = wave(:,1:s/2);
sE = wave(:,s/2+1:s);
%Predict
E = sO;
%update*
O = bitxor(sO,sE);
%return wavelet
O = uint8(dyadup(O,2));
O = [zeros(size(wave,1),1) O];
E = uint8(dyadup(E,2));
E = [E zeros(size(wave,1),1)];
signal = O + E;
end

