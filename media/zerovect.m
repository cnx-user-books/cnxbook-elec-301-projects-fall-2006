function Z = zerovect(vector)
% this function chops off the first 300 samples of a signal to remove
% the microphone turn on effect and centers the signal on zero, essentially
% removing any DC offset

Z = vector(300:end);

[A,B] = size(Z);

offset = mean(Z);

Z = Z - offset * ones(A,B);