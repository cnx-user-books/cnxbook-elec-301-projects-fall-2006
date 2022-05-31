% databaseCreationScript
% We ran this after recording everybodies' vowel sounds to get the...
% ... models we needed for the speaker recognition.  Each vowel sound...
% ... is a field in a structure.  

% dhamodel = Damen Hattori's "a" model sound
knownSpeaker(1).name = 'Damen Hattori';
knownSpeaker(1).welcome = HelloDamen;
knownSpeaker(1).vowel(1).model   = dhamodel; % vowel 1 = a
knownSpeaker(1).vowel(2).model   = dhahmodel; % vowel 2 = ah
knownSpeaker(1).vowel(3).model   = dhaymodel; % vowel 3 = ay
knownSpeaker(1).vowel(4).model   = dheemodel; % vowel 4 = ee
knownSpeaker(1).vowel(5).model   = dhehmodel; % vowel 5 = eh
knownSpeaker(1).vowel(6).model   = dhihmodel; % vowel 6 = ih
knownSpeaker(1).vowel(7).model   = dhohmodel; % vowel 7 = oh
knownSpeaker(1).vowel(8).model   = dhohmodel; % vowel 8 = ue

knownSpeaker(2).name = 'Josh Long';
knownSpeaker(2).welcome = HelloJosh;
knownSpeaker(2).vowel(1).model   = jlamodel;
knownSpeaker(2).vowel(2).model  = jlahmodel;
knownSpeaker(2).vowel(3).model   = jlaymodel;
knownSpeaker(2).vowel(4).model   = jleemodel;
knownSpeaker(2).vowel(5).model   = jlehmodel;
knownSpeaker(2).vowel(6).model   = jlihmodel;
knownSpeaker(2).vowel(7).model   = jlohmodel;
knownSpeaker(2).vowel(8).model   = jluemodel;

knownSpeaker(3).name = 'Matt McDonell';
knownSpeaker(3).welcome = HelloMatt;
knownSpeaker(3).vowel(1).model   = mmamodel;
knownSpeaker(3).vowel(2).model  = mmahmodel;
knownSpeaker(3).vowel(3).model   = mmaymodel;
knownSpeaker(3).vowel(4).model   = mmeemodel;
knownSpeaker(3).vowel(5).model   = mmehmodel;
knownSpeaker(3).vowel(6).model   = mmihmodel;
knownSpeaker(3).vowel(7).model   = mmohmodel;
knownSpeaker(3).vowel(8).model   = mmuemodel;

knownSpeaker(4).name = 'Chris Pasich';
knownSpeaker(4).welcome = HelloChris;
knownSpeaker(4).vowel(1).model   = cpamodel;
knownSpeaker(4).vowel(2).model  = cpahmodel;
knownSpeaker(4).vowel(3).model   = cpaymodel;
knownSpeaker(4).vowel(4).model   = cpeemodel;
knownSpeaker(4).vowel(5).model   = cpehmodel;
knownSpeaker(4).vowel(6).model   = cpihmodel;
knownSpeaker(4).vowel(7).model   = cpohmodel;
knownSpeaker(4).vowel(8).model   = cpuemodel;

% Give the models names in the structure
for k = 1:4;
knownSpeaker(k).vowel(1).name = 'a' ;
knownSpeaker(k).vowel(2).name = 'ah' ;
knownSpeaker(k).vowel(3).name = 'ay' ;
knownSpeaker(k).vowel(4).name = 'ee' ;
knownSpeaker(k).vowel(5).name = 'eh' ;
knownSpeaker(k).vowel(6).name = 'ih' ;
knownSpeaker(k).vowel(7).name = 'oh' ;
knownSpeaker(k).vowel(8).name = 'ue' ;
end

% Create final database to use 
speakerDatabase = knownSpeaker(1 : 4);