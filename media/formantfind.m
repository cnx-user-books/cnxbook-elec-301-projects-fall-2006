function [formantFreq, magForm] = formantfind(f, magH)
% Finds formants by finding peaks in vocal model by looking at...
% ... the derivative between points. Returns the formant frequencies...
%...and the magnitude

L = length(magH);

% Frequency step
delta = f(2) - f(1);

% Initialize first point at zero
deriv([1:2],1) = 0;

% Calculate derivative between each point
for k = 3:L
    deriv(k,1) = (magH(k) - magH(k-2))/delta;
end

% Redo first point, set equal to 2nd for smoothness
deriv([1:2],1) = deriv(3,1);

% Counter variables
index = 1;
f_index = 1;

% Constants
POS = 1;
NEG = -1;

% Sign of previous derivative - initialize as zero
prev = 0;

% Will store formant frequencies 
fmt_loc = [];


% While we still have values 
while index < L;

    % Look for perfect maxima
    if ((prev == POS) && (sign(deriv(index)) == 0))
        % Store formant frequencies and magnitude
        fmt_loc(f_index) = index;
        f_index = f_index + 1;
        % Approximate zero as a negative sign
        prev = NEG;
    
    % Look for sign changes - similar to maxima
    elseif ((prev == POS) && (sign(deriv(index)) == NEG))
        % Store formant frequencies and magnitude
        fmt_loc(f_index) = index;
        f_index = f_index + 1;
        prev = NEG;
        
    % In case of perfect minima
    elseif (sign(deriv(index)) == 0)
        prev = POS;
        
    else
        prev = sign(deriv(index));
        
    end
    
    index = index + 1;

end
    
% Make sure we have the maximum magnitude by looking at the 3 values on...
% ... both sides of each of our calculated formants, adjust the...
% ... corresponding frequencies
for h = 1:length(fmt_loc)
    [temp_mag(h),fmt_index(h)] = max(magH(fmt_loc(h)- 3 : fmt_loc(h) + 3));
    temp_f(h) = f(fmt_index(h)+ fmt_loc(h) - 4);
end

% Fill out final output vectors of calculated formants and corresponding...
% ... frequencies.  If the magnitude of a formant is less than the...
% ... magnitude of the next formant, discard the first formant
index = 1;
for t = 1: length(temp_f) - 1
    if temp_mag(t) > temp_mag(t+1)
        formantFreq(index)= temp_f(t);
        magForm(index) = temp_mag(t);
        index = index + 1;
    end
end







