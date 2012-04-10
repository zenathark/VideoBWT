function [ wave p ] = ippredict3x3( predicted, level )
%IPPREDICT3X3 Summary of this function goes here
%   Detailed explanation goes here
if level > 1
    xwave = predicted;
    height = size(predicted,1);
    width = size(predicted,2);
%HORIZONTAL------------------------
    for l=2:level
        for y=1:(height/2^l)
            for x=((width/2^l)+1)+1:(width/(2^(l-1)))-1
                if (predicted(y,x)+predicted(y,x+1)+predicted(y,x-1))>2
                    xwave(y*2-1,x*2-1) = 1;
                    xwave(y*2,x*2-1) = 1;
                    xwave(y*2-1,x*2) = 1;
                    xwave(y*2,x*2) = 1;
                end
            end
        end
    end
%HORIZONTAL------------------------
%VERTICAL------------------------
    for l=level:-1:2
        for y=((height/2^l)+1)+1:(height/(2^(l-1)))-2
            for x=1:(width/2^l)
                if (predicted(y,x)+predicted(y+1,x)+predicted(y-1,x))>2
                    xwave(y*2-1,x*2-1) = 1;
                    xwave(y*2,x*2-1) = 1;
                    xwave(y*2-1,x*2) = 1;
                    xwave(y*2,x*2) = 1;
                end
            end
        end
    end
%VERTICAL------------------------
%DIAGONAL------------------------
    for l=level:-1:2
        for y=((height/2^l)+1):(height/(2^(l-1)))
            for x=((width/2^l)+1):(width/(2^(l-1)))
                if (predicted(y,x)+predicted(y+1,x)+predicted(y-1,x)+predicted(y,x+1)+predicted(y,x-1))>3
                    xwave(y*2-1,x*2-1) = 1;
                    xwave(y*2,x*2-1) = 1;
                    xwave(y*2-1,x*2) = 1;
                    xwave(y*2,x*2) = 1;
                end
            end
        end
    end
%DIAGONAL------------------------
    p = xwave;
    wave = bitxor(predicted,xwave);
    th = height/(2^(level-1));
    tw = width/(2^(level-1));
    wave(1:th,1:tw) = predicted(1:th,1:tw);
end

end

