function [strt,ct] = newenvdetect(signal, threshold)
% this function finds the index number of the beginning
% and end of each vowel envelope
% Inputs:
% signal - zero centered speech vector
% threshold - desired cutoff value (.2 is what we used)

% square signal to accentuate peaks
Z = signal .* signal;

% set highest peak equal to one
Z= Z/(max(abs(Z)));

% clip signal to .1, remove low energy components
for i = 1:length(Z)
    if Z(i) > .1
        Z(i) = .1;
    elseif Z(i) < -.1
        Z(i) = -.1;
    elseif abs(Z(i)) < .04
        Z(i) = 0;
    end
end

% generate time vector for plotting - fs = 8000Hz   
t=[0:1/8000:(length(Z)-1)/8000];

% Length 512 Boxcar Filter
S = ones(512,1);

% implement filter via fast convolution
env = fastconv(abs(Z),S);
env = env/max(abs(env)) * .1;

% generate slighty extended time vector due to convolution
t2 = [0:1/8000:(length(env)-1)/8000];


th = .1 * threshold;

% counter variables
s_index = 1;
c_index = 1;

pos_edge = 1; % are we searching for positive or negative zero crossing
loc = 1; % current sample number
strt = [ ]; % will hold indices of beginning of syllable in envelope
ct = [ ]; % will hold indices of end of syllable in envelope

while loc < length(env)
    
    if pos_edge == 1 % are we looking for a positive threshold crossing
        I = (loc-1)+ find(env(loc:end) > th , 1); % find crossing
        if isempty(I) % cancel search if no crossing found
          break
        end
        strt(s_index)= I; % keep track of crossing location
        s_index = s_index + 1;
        pos_edge = 0; % start looking for negative crossings
        loc = I + 100; % skip ahead a few samples to avoid noisy behavior
        
    elseif pos_edge == 0 % are we looking for a negative threshold crossing
        
        G = (loc-1) + find(env(loc:end) < th , 1);
        if isempty(G)
            break
        end
        ct(c_index)= G;
        c_index = c_index + 1;
        pos_edge = 1;
        loc = G + 100;
    end
    
end
