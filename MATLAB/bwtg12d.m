function [ wave ] = bwtg12d( img )
%BWTG12D Summary of this function goes here
%   Detailed explanation goes here
%get wavelet from rows
waverows = bwtg1(img);
%get wavelet from columns
wave = bwtg1(waverows');
wave = wave';
end

