%Authors: Weston Harper, Ramy-Claude James, John Passiak
%Date: Dec 19, 2006
%Inputs:
%recording_in - output signal (in the time domain) that will be deconvolvedto find it's input.
%filter_in - the known filter value that will be deconvolved from the output signal
%Outputs:
%signal - the recovered input signal
%Method:
%this function performs deconvolution by division in the frequency domain
%to recover the input signal to a system.

function [signal] = normal_deconvolution(recording_in,filter_in) 

keylength=max([length(recording_in);length(filter_in)]);%finds which input is longer
recording_pad=fft([recording_in;zeros(keylength-length(recording_in),1)]);%zero pads recording if shorter than filter and takes the resulting FFT
filter_pad=fft([filter_in;zeros(keylength-length(filter_in),1)]);%zero pads the filter if shorter than the recording and takes the resulting FFT
signal=real(ifft(recording_pad./filter_pad));%divides the FFT and translates the resulting signal back to the time domain
