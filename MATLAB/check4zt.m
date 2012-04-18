function [ ztmtx ] = check4zt( wave, level )
%CHECK4ZT Summary of this function goes here
%   Detailed explanation goes here
[r c] = size(wave);
ztmtx = logical(zeros(r/2,c/2) == 1);
l = 1;
for y=r/(2^l):-1:r/(2^(l+1))-1
    for x=c/(2^l):-1:1
        if wave(y,x) == 0 && wave(2*y,2*x) == 0 && wave(2*y-1,2*x) == 0 && wave(2*y,2*x-1) == 0 && wave(2*y-1,2*x-1) == 0
            ztmtx(y,x) = 1;
        end
    end
end
for y=r/(2^l):-1:1
    for x=c/(2^l):-1:c/(2^(l+1))-1
        if wave(y,x) == 0 && wave(2*y,2*x) == 0 && wave(2*y-1,2*x) == 0 && wave(2*y,2*x-1) == 0 && wave(2*y-1,2*x-1) == 0
            ztmtx(y,x) = 1;
        end
    end
end
if level>2
    for l=2:level-1
        for y=r/(2^l):-1:r/(2^(l+1))-1
            for x=c/(2^l):-1:1
                if ztmtx(2*y,2*x) == 1 && ztmtx(2*y-1,2*x) == 1 && ztmtx(2*y,2*x-1) == 1 && ztmtx(2*y-1,2*x-1) == 1
                    ztmtx(y,x) = 1;
                end
            end
        end
        for y=r/(2^l):-1:1
            for x=c/(2^l):-1:c/(2^(l+1))-1
                if ztmtx(2*y,2*x) == 1 && ztmtx(2*y-1,2*x) == 1 && ztmtx(2*y,2*x-1) == 1 && ztmtx(2*y-1,2*x-1) == 1
                    ztmtx(y,x) = 1;
                end
            end
        end
    end
end
end

