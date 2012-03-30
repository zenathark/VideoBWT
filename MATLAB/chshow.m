function [ ] = chshow( bch )
%CHSHOW Summary of this function goes here
%   Detailed explanation goes here
ch8(:,:) = bch(8,:,:);
ch7(:,:) = bch(7,:,:);
ch6(:,:) = bch(6,:,:);
ch5(:,:) = bch(5,:,:);
ch4(:,:) = bch(4,:,:);
ch3(:,:) = bch(3,:,:);
ch2(:,:) = bch(2,:,:);
ch1(:,:) = bch(1,:,:);
im = [ch8 ch7 ch6 ch5;
      ch4 ch3 ch2 ch1];
imshow(im);
end

