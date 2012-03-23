function [ wave ] = bwtg12d( img )
%BWTG12D Summary of this function goes here
%   Detailed explanation goes here
waverows = bwtg1(img);
wave = bwtg1(waverows');
end

