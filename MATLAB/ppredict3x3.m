function [ predicted p ] = ppredict3x3( wave, level )
%PPREDICT3X3 Summary of this function goes here
%   Detailed explanation goes here
if level > 1
    xwave = wave;
    height = size(wave,1);
    width = size(wave,2);
%HORIZONTAL------------------------
    xwave(1:128,129:256) = 0;
    for l=2:level
        for y=2:(height/2^l)-1
            for x=((width/2^l)+1)+1:(width/(2^(l-1)))-1
                A=wave(y-1,x-1);
                B=wave(y-1,x);
                C=wave(y,x-1);
                d=wave(y*2-1,x*2-1);
                e=wave(y*2-1,x*2);
                f=wave(y*2-1,x*2+1);
                g=wave(y*2,x*2-1);
                j=wave(y*2+1,x*2-1);
                X = 2*(e + g) + (d + f + j) - 3*(B + C) - A;
                if X <= 4  && wave(y,x) == 1
                    rflag = 0;
                    while X<=4
                        if rflag == 0
                            X = X + 1;
                            xwave(y*2,x*2) = 1;
                           
                        end
                        if rflag == 1
                            X = X + 2;
                            xwave(y*2-1,x*2) = 1;
                        end
                        if rflag == 2
                            X = X + 2;
                            xwave(y*2,x*2-1) = 1;
                        end
                        if rflag == 3
                            X = X + 5;
                            xwave(y*2-1,x*2-1) = 1;
                        end
                    end
                elseif X<=4 && wave(y,x) == 0
                    xwave(y*2-1,x*2-1) = 0;
                    xwave(y*2,x*2-1) = 0;
                    xwave(y*2-1,x*2) = 0;
                    xwave(y*2,x*2) = 0;
                elseif X>=5 && wave(y,x) == 1
                    xwave(y*2-1,x*2-1) = 0;
                    xwave(y*2,x*2-1) = 0;
                    xwave(y*2-1,x*2) = 0;
                    xwave(y*2,x*2) = 0;
                else
                    xwave(y*2-1,x*2-1) = 0;
                    xwave(y*2,x*2-1) = 0;
                    xwave(y*2-1,x*2) = 0;
                    xwave(y*2,x*2) = 0;
                end
                %if (2*wave(y,x)+wave(y+1,x)+wave(y-1,x)+wave(y,x+1)+wave(y,x-1))>1
                %    xwave(y*2-1,x*2-1) = 1;
                %    xwave(y*2,x*2-1) = 1;
                %    xwave(y*2-1,x*2) = 1;
                %    xwave(y*2,x*2) = 1;
                %end
            end
        end
    end
%HORIZONTAL------------------------
%VERTICAL------------------------
    xwave(129:256,1:128) = 0;
    for l=2:level
        for y=((height/2^l)+1)+1:(height/(2^(l-1)))-2
            for x=2:(width/2^l)-1
                if (2*wave(y,x)+wave(y+1,x)+wave(y-1,x)+wave(y,x+1)+wave(y,x-1))>1
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
    xwave(129:256,129:256) = 0;
    for l=2:level
        for y=((height/2^l)+1):(height/(2^(l-1)))
            for x=((width/2^l)+1):(width/(2^(l-1)))
                if (wave(y,x)+wave(y+1,x)+wave(y-1,x)+wave(y,x+1)+wave(y,x-1))>1
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
    predicted = bitxor(wave,xwave);
    th = height/(2^(level-1));
    tw = width/(2^(level-1));
    predicted(1:th,1:tw) = wave(1:th,1:tw);
end

end

