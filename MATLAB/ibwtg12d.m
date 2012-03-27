function [ img ] = ibwtg12d( wave )
%IBWTG12D Summary of this function goes here
%   Detailed explanation goes here
%get wavelet from rows
wavecols = ibwtg1(wave');
%get wavelet from columns
img = ibwtg1(wavecols');
%wave = wave';

end

