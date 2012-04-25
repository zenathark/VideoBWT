function [ compressed_string rlez3] = compressbwt3z3(  wave,level )
%COMPRESSBWT3Z3 Summary of this function goes here
%   Detailed explanation goes here
[R C] = size(wave);
out = zeros(1,R*C);
bit_bucket=1;
zt = check4zt(wave,level);
rlez3_pass = zeros(1,256*256);
rlez3_c = 1;
rlez3_a = 0;
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
                if wave(y,x) == 0 && y<=size(zt,1) && x<=size(zt,2)
                    if  zt(y,x) == 0
                        if (rlez3_a > 1)
                            rlez3_pass(rlez3_c) = rlez3_a;
                            rlez3_a = 0;
                            rlez3_c = rlez3_c + 1;
                        end
                        out(bit_bucket) = wave(y,x);
                        bit_bucket=bit_bucket+1;
                        out(bit_bucket) = 1;
                        bit_bucket=bit_bucket+1;
                    else
                        if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                            if (zt(y+R/(2^l),x) == 1) && (zt(y+R/(2^l),x-C/(2^l))==1)
                                if (rlez3_a == 0)
                                    out(bit_bucket) = wave(y,x);
                                    bit_bucket=bit_bucket+1;
                                    out(bit_bucket) = 0;
                                    bit_bucket=bit_bucket+1;
                                    out(bit_bucket) = 1;
                                    bit_bucket=bit_bucket+1;
                                end
                                rlez3_a = rlez3_a + 1;
                            else
                                if (rlez3_a > 1)
                                    rlez3_pass(rlez3_c) = rlez3_a;
                                    rlez3_a = 0;
                                    rlez3_c = rlez3_c + 1;
                                end
                                out(bit_bucket) = wave(y,x);
                                bit_bucket=bit_bucket+1;
                                out(bit_bucket) = 0;
                                bit_bucket=bit_bucket+1;
                                out(bit_bucket) = 0;
                                bit_bucket=bit_bucket+1;
                            end
                        end
                    end
                else
                    if (rlez3_a > 1)
                        rlez3_pass(rlez3_c) = rlez3_a;
                        rlez3_a = 0;
                        rlez3_c = rlez3_c + 1;
                    end
                    out(bit_bucket) = wave(y,x);
                    bit_bucket=bit_bucket+1;
                end
            end
        end 
    end
    for y=(R/2^l)+1:R/2^(l-1)
        for x=1:C/2^(l-1)
            parent_x = ceil(x/2);
            parent_y = ceil(y/2);
            if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                if wave(y,x) == 0 && y<=size(zt,1) && x<=size(zt,2)
                    if zt(y,x) == 0
                        out(bit_bucket) = 0;
                        bit_bucket=bit_bucket+1;
                        out(bit_bucket) = 0;
                        bit_bucket=bit_bucket+1;
                    else
                        if (parent_x<=C/2^level && parent_y<=R/2^level) || zt(parent_y,parent_x) == 0
                            nb_y = y-R/2^l;
                            if x <= C/2^l
                                nb_x = x + C/2^l;
                                nb2_x = nb_x;
                            else
                                nb_x = x;
                                nb2_x = x - C/2^l;
                            end
                            if zt(nb_y,nb_x) == 0 ||  zt(y,nb2_x) == 0
                                out(bit_bucket) = 0;
                                bit_bucket=bit_bucket+1;
                                out(bit_bucket) = 1;
                                bit_bucket=bit_bucket+1;
                            end
                        end
                    end
                else
                    out(bit_bucket) = wave(y,x);
                    bit_bucket=bit_bucket+1;
                end
            end
        end 
    end
end
compressed_string = out(1:bit_bucket-1);
rlez3 = rlez3_pass(1:rlez3_c-1);
end

