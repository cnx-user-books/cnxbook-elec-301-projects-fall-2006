function main(mysterySignal, speakerDatabase)
% This is the main function that takes in the recorded mystery signal...
% ... from the person trying to "access the system".  The...
% ... speakerDatabase contains all the vowel models of people "allowed".

% This sets the threshhold for envelope detection
envThreshhold = .2;

error = '';

% Chop off first 300 samples (to eliminate microphone "turn-on" spike)...
% ... and zero the signal
mysterySignal = zerovect(mysterySignal);

% Normalize mysterySignal
mysterySignal = mysterySignal * ( 1 / sqrt(mysterySignal' * mysterySignal) );

% Run envelope detector to find beginning and end of each vowel in mysSig
[vowelStart, vowelEnd] = newenvdetect(mysterySignal, envThreshhold);

numberOfVowels = length(vowelStart);

% For each vowel in the mystery signal 
for k = 1 : numberOfVowels
    
    % Isolate kth vowel of mysterySignal
    mysVowelSig = mysterySignal( vowelStart(k) : vowelEnd(k) );
    
    % Create vocal model of kth vowel
    [freqkHz, mysVowelModelMag] = vocalModel(mysVowelSig);
    figure (k)
    subplot(511)
    semilogy(freqkHz, mysVowelModelMag);
    title('mysterySig')
    
    % Figure out which vowel (a, e, i, o...) the kth vowel is
    % Find formant frequencies of kth vowel
    [formantFreq, magForm] = formantfind(freqkHz, mysVowelModelMag);
    formantFreq

    % Database of known first and second formant values for vowel sounds
    % Note: we are using only select vowels appearing in the groups'...
    % ... names to emulate a security database  
    % format: vowel = [minFormant1 maxFormant1 minFormant2 maxFormant2]

     a = [ .670 .810 1.500 1.900 ]; % as in cat
    ah = [ .675 .830 .900 1.500]; % as in dog
    ay = [ .455 .570 1.875 2.500 ]; % as in pay
    ee = [ .200 .350 2.000 2.650 ]; % as in meet
    eh = [ .510 .675 1.620 1.960 ]; % as in bet
    ih = [ .350 .500 1.975 2.145 ]; % as in fish
    oh = [ .400 .600 .850 1.500 ]; % as in boat

    % Compare computed mystery formants with database of known values
    if ( formantFreq(1) > a(1) && formantFreq(1) < a(2)...
        && formantFreq(2) > a(3) && formantFreq(2) < a(4) )

        mysVowel(k) = 1; % 'a'

    elseif ( formantFreq(1) > ah(1) && formantFreq(1) < ah(2)...
        && formantFreq(2) > ah(3) && formantFreq(2) < ah(4) )

        mysVowel(k) = 2; % 'ah'

    elseif ( formantFreq(1) > ay(1) && formantFreq(1) < ay(2)...
        && formantFreq(2) > ay(3) && formantFreq(2) < ay(4) )

        mysVowel(k) = 3; % 'ay'

    elseif ( formantFreq(1) > ee(1) && formantFreq(1) < ee(2)...
        && formantFreq(2) > ee(3) && formantFreq(2) < ee(4) )

        mysVowel(k) = 4; % 'ee'

    elseif ( formantFreq(1) > eh(1) && formantFreq(1) < eh(2)...
        && formantFreq(2) > eh(3) && formantFreq(2) < eh(4) )

        mysVowel(k) = 5; % 'eh'

    elseif ( formantFreq(1) > ih(1) && formantFreq(1) < ih(2)...
        && formantFreq(2) > ih(3) && formantFreq(2) < ih(4) )

        mysVowel(k) = 6; % 'ih'

    elseif ( formantFreq(1) > oh(1) && formantFreq(1) < oh(2)...
        && formantFreq(2) > oh(3) && formantFreq(2) < oh(4) )

        mysVowel(k) = 7; % 'oh'
    else
        errmsg = sprintf('no vowel match in vowel number %d', k);
        disp(errmsg)
        continue

    end

    % Output to console which vowel was found
    disp(sprintf(' Found %d', mysVowel(k)));

    numberOfPeople = length(speakerDatabase);
    
    % For each person in speakerDatabase a.k.a. number of known speakers
    for l = 1 : numberOfPeople
        
        % Give each person a score for the current vowel and put it...
        % ... into a matrix with a row for each person and a column...
        % ... for each vowel
        matchMatrix(l, k) = speakerDatabase(l).vowel(mysVowel(k)).model'* mysVowelModelMag;
        figure(k)
        subplot(5,1,(l+1))
        semilogy(freqkHz , speakerDatabase(l).vowel(mysVowel(k)).model)
        title(sprintf('speaker number %d',l))
    end
    
end
            
% Add up the elements of each row in matchMatrix to get a total "score"...
% ... for each person
for m = 1 : size(matchMatrix, 1)
    scoreVector(m) = sum( matchMatrix(m, :));
end

% Get the highest score
[highestScore, highestScoreIndex] = max(scoreVector);

matchMatrix
scoreVector

% Set threshhold score that must be achieved for a positive match
thresholdMultiplier = length(nonzeros(matchMatrix(1,:)));
threshholdScore = .85 * thresholdMultiplier;


if (highestScore > threshholdScore)
    disp(['Hello ' speakerDatabase(highestScoreIndex).name]); 
    soundsc(speakerDatabase(highestScoreIndex).welcome,8000);

else
    disp('No Match');
end