function [ compressed_string ] = compressbwt( wave,level )
%COMPRESSBWT Summary of this function goes here
%   Detailed explanation goes here
[R C] = size(wave);
out = zeros(1,R*C);
bit_bucket=1;
zt = check4zt(wave,level);
for y=1:R/2^level
    for x=1:C/2^level
        out(bit_bucket) = wave(y,x);
        bit_bucket=bit_bucket+1;
    end 
end
for l=level:-1:1
    for y=1:R/2^l
        for x=(C/2^l)+1:C/2^(l-1)
            parent_x = ceil(x/2);
            parent_y = ceil(y/2);
            if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                out(bit_bucket) = wave(y,x);
                bit_bucket=bit_bucket+1;
                if wave(y,x) == 0 && y<=size(zt,1) && x<=size(zt,2)
                    if  zt(y,x) == 0
                        out(bit_bucket) = 0;
                        bit_bucket=bit_bucket+1;
                    else
                        if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                            out(bit_bucket) = 1;
                            bit_bucket=bit_bucket+1;
                        end
                    end
                end
            end
        end 
    end
    for y=(R/2^l)+1:R/2^(l-1)
        for x=1:C/2^(l-1)
            parent_x = ceil(x/2);
            parent_y = ceil(y/2);
            if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                out(bit_bucket) = wave(y,x);
                bit_bucket=bit_bucket+1;
                if wave(y,x) == 0 && y<=size(zt,1) && x<=size(zt,2)
                    if zt(y,x) == 0
                        out(bit_bucket) = 0;
                        bit_bucket=bit_bucket+1;
                    else
                        if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                            out(bit_bucket) = 1;
                            bit_bucket=bit_bucket+1;
                        end
                    end
                end
            end
        end 
    end
end
compressed_string = out(1:bit_bucket-1);
end

