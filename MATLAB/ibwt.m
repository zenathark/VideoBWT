function [ im ] = ibwt( wave, level, name )
%IBWT Summary of this function goes here
%   Detailed explanation goes here
s = size(wave,1);
im = wave;
for l=level:-1:1
    switch name
        case 1
            im(1:s/l,1:s/l) = ibwtg12d(im(1:s/l,1:s/l));
        case 2
            im(1:s/l,1:s/l) = ibwtg22d(im(1:s/l,1:s/l));
    end
end

end

