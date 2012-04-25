function [ dec ] = bwt( signal, level, name )
%BWT Summary of this function goes here
%   Detailed explanation goes here
s = size(signal,1);
dec = signal;
for l=1:level
    switch name
        case 1
            dec(1:s/2^(l-1),1:s/2^(l-1)) = bwtg12d(dec(1:s/2^(l-1),1:s/2^(l-1)));
        case 2
            dec(1:s/2^(l-1),1:s/2^(l-1)) = bwtg22d(dec(1:s/2^(l-1),1:s/2^(l-1)));
    end
end
end