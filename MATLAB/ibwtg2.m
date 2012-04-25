function [ signal ] = ibwtg2( wave )
%IBWTG2 Summary of this function goes here
%   Detailed explanation goes here
s = size(wave,2);
sE = wave(:,1:s/2);
sO = wave(:,s/2+1:s);
%Predict
E = bitxor(sO,sE);
%Update
E1(:,2:size(E,2)) = E(:,1:size(E,2)-1);
E1(:,1) = E(:,size(E,2));
E = E1;
%Predict
O = bitxor(sO,E);
%Update
%E = E;
%return wavelet
O = uint8(dyadup(O,2));
O = [O zeros(size(wave,1),1)];
E = uint8(dyadup(E,2));
E = [zeros(size(wave,1),1) E];
signal = O + E;
end

