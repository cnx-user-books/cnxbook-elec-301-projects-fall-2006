function graphFormant(f, magH, formant, magForm graphTitle)

% Graph frequency response of vocal model on a log scale for best view of...
% ... formants (peaks on graph).
 
figure;
semilogy(f, magH, 'r');
xlabel('Frequency(kHz)');
ylabel('Response');
title(graphTitle);

% Put circles at formants
hold on;
for g = 1:length(formant)
    stem(formant(g),magF(g))
end