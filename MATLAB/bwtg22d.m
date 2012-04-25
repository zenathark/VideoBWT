function [ wave ] = bwtg22d( img )
%BWTG22D Summary of this function goes here
%   Detailed explanation goes here
%get wavelet from rows
waverows = bwtg2(img);
%get wavelet from columns
wave = bwtg2(waverows');
wave = wave';
end

