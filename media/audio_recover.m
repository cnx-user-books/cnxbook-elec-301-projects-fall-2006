%Authors: Weston Harper, Ramy-Claude James, John Passiak
%Date: Dec 19, 2006
%Adapted from matlab example for deconvwnr in matlab help for Image Processing Toolbox User's Guide
%Inputs:
%recording_in - audio signal (in the time domain) recorded in a room that will have the room's
%affect on it removed
%source_in - the known filter value of the room in the time domain
%noise_in - a time domain signal that is a example of the background noise
%of the room the recording was done in.
%original_in - either the original signal played and recorded to make
%   recording_in, or a signal that is statistically representative of the
%   type of signal that recording_in is
%Outputs:
%signal - a version of recording with the room affects removed
%Method:
%This function calculates the correlation of the ambient noise (noise_in)
%and the expected correlation of the input recording (from original_in) and
%uses these values in Wiener deconvolution to recover the audio signal.

function [signal] = audio_recover(recording_in,source_in,noise_in,original_in) 

recording=recording_in;
source=source_in; %determined room impulse response
noise=noise_in; %noise=wavread(noise_in);
original=original_in; %original=wavread(original_in);

% noise power
NP = abs(fftn(noise)).^2; 
% noise autocorrelation function, centered
NCORR = fftshift(real(ifftn(NP))); 

% original source power
originalP = abs(fftn(original)).^2;
% original autocorrelation function, centered
originalCORR = fftshift(real(ifftn(originalP)));

temp_signal=deconvwnr(recording,source,NCORR,originalCORR);%performs Wiener deconvolution
signal=temp_signal./(max(abs(temp_signal)));%normalizes the signal to values of -1 to 1 for easy output in matlab sound and wav commands