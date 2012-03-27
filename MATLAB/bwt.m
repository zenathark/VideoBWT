function [ dec ] = bwt( signal, level, name )
%BWT Summary of this function goes here
%   Detailed explanation goes here
s = size(signal,1);
dec = signal;
for l=1:level
    switch name
        case 1
            dec(1:s/l,1:s/l) = bwtg12d(dec(1:s/l,1:s/l));
    end
end

end