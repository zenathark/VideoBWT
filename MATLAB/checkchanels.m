bpp = 0;
for i=1:8
    shk = zeros(512,512);
    shk(:,:) = lennash_5(i,:,:);
    if (sum(sum(shk))/(512*512) < 0.40)
        [c zt] = compressbwt3z3(shk,5);
        real_size = size(c,2) + sum(z3>255)*(16+2) + sum(z3>255)*(8+1);
        bpp = bpp + size(c,2) / (512*512);
    else
        bpp = bpp + 1;
    end
end