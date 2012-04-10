function [ xwave ] = predictsingle( wave, level )
%PREDICTSINGLE Summary of this function goes here
%   Detailed explanation goes here
if level > 1
    xwave = wave;
    height = size(wave,1);
    width = size(wave,2);
%HORIZONTAL------------------------
    for l=2:level
        for y=1:(height/2^l)
            for x=((width/2^l)+1):(width/(2^(l-1)))
                xwave(y*2-1,x*2-1) = bitxor(wave(y*2-1,x*2-1),wave(y,x));
                xwave(y*2,x*2-1) = bitxor(wave(y*2,x*2-1),wave(y,x));
                xwave(y*2-1,x*2) = bitxor(wave(y*2-1,x*2),wave(y,x));
                xwave(y*2,x*2) = bitxor(wave(y*2,x*2),wave(y,x));
            end
        end
    end
%HORIZONTAL------------------------
%VERTICAL------------------------
    for l=2:level
        for y=((height/2^l)+1):(height/(2^(l-1)))
            for x=1:(width/2^l)
                xwave(y*2-1,x*2-1) = bitxor(wave(y*2-1,x*2-1),wave(y,x));
                xwave(y*2,x*2-1) = bitxor(wave(y*2,x*2-1),wave(y,x));
                xwave(y*2-1,x*2) = bitxor(wave(y*2-1,x*2),wave(y,x));
                xwave(y*2,x*2) = bitxor(wave(y*2,x*2),wave(y,x));
            end
        end
    end
%VERTICAL------------------------
%DIAGONAL------------------------
    for l=2:level
        for y=((height/2^l)+1):(height/(2^(l-1)))
            for x=((width/2^l)+1):(width/(2^(l-1)))
                xwave(y*2-1,x*2-1) = bitxor(wave(y*2-1,x*2-1),wave(y,x));
                xwave(y*2,x*2-1) = bitxor(wave(y*2,x*2-1),wave(y,x));
                xwave(y*2-1,x*2) = bitxor(wave(y*2-1,x*2),wave(y,x));
                xwave(y*2,x*2) = bitxor(wave(y*2,x*2),wave(y,x));
            end
        end
    end
%DIAGONAL------------------------
end
end

