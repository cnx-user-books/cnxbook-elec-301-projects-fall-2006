function [f,magH] = vocalModel(vowelSig)
% This function takes in an isolated vowel signal and creates a model...
% ... of the vocal "system" from which it came

% Normalize speech signal max to 1
normalizedVowelSig = vowelSig / max(abs(vowelSig));

% Create an order 15 auto-regressive model of the speech signal...
% ... (essentially a model of the vocal tract of the person who...
% ... produced the signal) using the forward-backward approach. This...
% ... model is essentially an infinite-impulse response (IIR) filter.
model = ar(normalizedVowelSig,15);

% Determine the numerators and denominators of the transfer function...
% ... of the auto-regressive (AR) model.
[num,den] = tfdata(model,'v');

% Compute the frequency response (H) of the AR model.
[H,w] = freqz(num,den);

% Convert radians from frequency response into kHz for more intuitive graph.
f = w.*8000/(2000*pi);

% Compute magnitude of H as final y values of vocal model
magH = abs(H);

% Normalize magH
magH = magH/(sqrt(magH'*magH));