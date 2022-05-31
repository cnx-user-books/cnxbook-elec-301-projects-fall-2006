function [formants, frmtModel ] = getFormantData(trial1)% , trial2, trial3, trial4, trial5)
% this function is used to get formant frequencies and vocal models for one vowel and one speaker
% code is correctly commented to take in one trial, can be expanded to five by uncommenting

%inputs: trial* = individual recordings of one vowel by one speaker

%outputs:   f1 = avg. value of first formant
%           f2 = avg. value of second formant
%           frmtModel = length 512 average model of formant graph

% Get Data for trial 1

%chop off first 300 samples and center
t1 = zerovect(trial1);
% get beginning and end of vowel
[strt , ct ] = envdect( t1 , .2);

%working with first and only vowel
vowelSig1 = t1(strt(1):ct(1));
% get formant model for vowel 
[freq , frmtModel] = vocalModel(vowelSig1);

[t1Freq , t1Mag] = formantfind(freq,frmtModel);

formants = t1Freq([1:2]);

% % Get Data for trial 2
% 
% %chop off first 300 samples and center
% t2 = zerovect(trial2);
% % get beginning and end of vowel
% [strt , ct ] = envdect( t2 , .2);
% 
% %working with first and only vowel
% vowelSig2 = t2(strt(1):ct(1));
% % get formant model for vowel 
% [freq , t2Model] = vocalModel(vowelSig2);
% 
% [t2Freq , t2Mag] = formantfind(freq, t2Model);
% 
% % Get Data for trial 3
% 
% %chop off first 300 samples and center
% t3 = zerovect(trial3);
% % get beginning and end of vowel
% [strt , ct ] = envdect( t3 , .2);
% 
% %working with first and only vowel
% vowelSig3 = t3(strt(1):ct(1));
% % get formant model for vowel 
% [freq , t3Model] = vocalModel(vowelSig3);
% 
% [t3Freq , t3Mag] = formantfind(freq, t3Model);
% 
% % Get Data for trial 4
% 
% %chop off first 300 samples and center
% t4 = zerovect(trial4);
% % get beginning and end of vowel
% [strt , ct ] = envdect( t4 , .2);
% 
% %working with first and only vowel
% vowelSig4 = t4(strt(1):ct(1));
% % get formant model for vowel 
% [freq , t4Model] = vocalModel(vowelSig4);
% 
% [t4Freq , t4Mag] = formantfind(freq, t4Model);
% 
% % Get Data for trial 5
% 
% %chop off first 300 samples and center
% t5 = zerovect(trial5);
% % get beginning and end of vowel
% [strt , ct ] = envdect( t5 , .2);
% 
% %working with first and only vowel
% vowelSig5 = t5(strt(1):ct(1));
% % get formant model for vowel 
% [freq , t5Model] = vocalModel(vowelSig5);
% 
% [t5Freq , t5Mag] = formantfind(freq, t5Model);

% f1 = (1/5) * (t1Freq(1) + t2Freq(1) + t3Freq(1) + t4Freq(1) + t5Freq(1));
% f2 = (1/5) * (t1Freq(2) + t2Freq(2) + t3Freq(2) + t4Freq(2) + t5Freq(2));
% formants = [f1,f2];
% frmtModel = (1/5) * (t1Model + t2Model + t3Model + t4Model + t5Model);