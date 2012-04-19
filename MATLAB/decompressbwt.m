function [ out ] = decompressbwt( compressed_string, level, R, C )
%DECOMPRESSBWT Summary of this function goes here
%   Detailed explanation goes here
out = (zeros(R,C) == 1);
bit_bucket=1;
zt = (zeros(R/2,C/2) == 1);
for y=1:R/2^level
    for x=1:C/2^level
        out(y,x) = compressed_string(bit_bucket);
        bit_bucket=bit_bucket+1;
    end 
end
for l=level:-1:1
    for y=1:R/2^l
        for x=(C/2^l)+1:C/2^(l-1)
            parent_x = ceil(x/2);
            parent_y = ceil(y/2);
            if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                out(y,x) = compressed_string(bit_bucket);
                bit_bucket=bit_bucket+1;
                if (out(y,x) == 0 && y<=size(zt,1) && x<=size(zt,2))
                    nxt_bit = compressed_string(bit_bucket);
                    bit_bucket=bit_bucket+1;
                    if (nxt_bit == 1)
                        zt(y,x) = 1;
                    end
                end
            elseif zt(parent_y,parent_x) == 1 && y<=size(zt,1) && x<=size(zt,2)
                zt(y,x) = 1;
            end
        end 
    end
    for y=(R/2^l)+1:R/2^(l-1)
        for x=1:C/2^(l-1)
            parent_x = ceil(x/2);
            parent_y = ceil(y/2);
            if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                out(y,x) = compressed_string(bit_bucket);
                bit_bucket=bit_bucket+1;
                if (out(y,x) == 0 && y<=size(zt,1) && x<=size(zt,2))
                    nxt_bit = compressed_string(bit_bucket);
                    bit_bucket=bit_bucket+1;
                    if (nxt_bit == 1)
                        zt(y,x) = 1;
                    end
                end
            elseif zt(parent_y,parent_x) == 1 && y<=size(zt,1) && x<=size(zt,2)
                zt(y,x) = 1;
            end
        end 
    end
end

