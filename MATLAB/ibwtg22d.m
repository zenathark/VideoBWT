function [ img ] = ibwtg22d( wave )
%IBWTG22D Summary of this function goes here
%   Detailed explanation goes here
%get wavelet from rows
wavecols = ibwtg2(wave');
%get wavelet from columns
img = ibwtg2(wavecols');
%wave = wave';
end

