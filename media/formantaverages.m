% this script calls getFormantData to get the formant frequencies and vocal model for all of the vowels for each speaker.
% it stores in the workplace the models used for matched filtering and vowel formant ranges for vowel recognition


% calculate values for ay
[ay(1,:),cpaymodel] = getFormantData(aycp1);%cp2,aycp3,aycp4,aycp5);
[ay(2,:),mmaymodel] = getFormantData(aymm1);%,aymm2,aymm3,aymm4,aymm5);
[ay(3,:),jlaymodel] = getFormantData(ayjl1);%,ayjl2,ayjl3,ayjl4,ayjl5);
[ay(4,:),dhaymodel] = getFormantData(aydh1);%,aydh2,aydh3,aydh4,aydh5);

ayAVG = .25 * sum(ay); % average value of vowel formants ( f1 & f2)
ayMAX = max(ay); % maximum value of vowel formants
ayMIN = min(ay); % minimum value of vowel formants

% calculate values for eh
[eh(1,:),cpehmodel] = getFormantData(ehcp1);%,ehcp2,ehcp3,ehcp4,ehcp5);
[eh(2,:),mmehmodel] = getFormantData(ehmm1);%,ehmm2,ehmm3,ehmm4,ehmm5);
[eh(3,:),jlehmodel] = getFormantData(ehjl1);%,ehjl2,ehjl3,ehjl4,ehjl5);
[eh(4,:),dhehmodel] = getFormantData(ehdh1);%,ehdh2,ehdh3,ehdh4,ehdh5);

ehAVG = .25 * sum(eh);
ehMAX = max(eh);
ehMIN = min(eh);

 %calculate values for a
[a(1,:),cpamodel] = getFormantData(acp1);%,acp2,acp3,acp4,acp5);
[a(2,:),mmamodel] = getFormantData(amm1);%,amm2,amm3,amm4,amm5);
[a(3,:),jlamodel] = getFormantData(ajl1);%,ajl2,ajl3,ajl4,ajl5);
[a(4,:),dhamodel] = getFormantData(adh1);%,adh2,adh3,adh4,adh5);

aAVG = .25 * sum(a);
aMAX = max(a);
aMIN = min(a);

 %calculate values for oh
[oh(1,:),cpohmodel] = getFormantData(ohcp1);%,ohcp2,ohcp3,ohcp4,ohcp5);
[oh(2,:),mmohmodel] = getFormantData(ohmm1);%,ohmm2,ohmm3,ohmm4,ohmm5);
[oh(3,:),jlohmodel] = getFormantData(ohjl1);%,ohjl2,ohjl3,ohjl4,ohjl5);
[oh(4,:),dhohmodel] = getFormantData(ohdh1);%,ohdh2,ohdh3,ohdh4,ohdh5);

ohAVG = .25 * sum(oh);
ohMAX = max(oh);
ohMIN = min(oh);

% calculate values for ee
[ee(1,:),cpeemodel] = getFormantData(eecp1);%,eecp2,eecp3,eecp4,eecp5);
[ee(2,:),mmeemodel] = getFormantData(eemm1);%,eemm2,eemm3,eemm4,eemm5);
[ee(3,:),jleemodel] = getFormantData(eejl1);%,eejl2,eejl3,eejl4,eejl5);
 [ee(4,:),dheemodel] = getFormantData(eedh1);%,eedh2,eedh3,eedh4,eedh5);

eeAVG = .25 * sum(ee);
eeMAX = max(ee);
eeMIN = min(ee);

% calculate values for ue
[ue(1,:),cpuemodel] = getFormantData(uecp1);%,uecp2,uecp3,uecp4,uecp5);
[ue(2,:),mmuemodel] = getFormantData(uemm1);%,uemm2,uemm3,uemm4,uemm5);
[ue(3,:),jluemodel] = getFormantData(uejl1);%,uejl2,uejl3,uejl4,uejl5);
 [ue(4,:),dhuemodel] = getFormantData(uedh1);%,uedh2,uedh3,uedh4,uedh5);

ueAVG = .25 * sum(ue);
ueMAX = max(ue);
ueMIN = min(ue);

% calculate values for ih
[ih(1,:),cpihmodel] = getFormantData(ihcp1);%,ihcp2,ihcp3,ihcp4,ihcp5);
[ih(2,:),mmihmodel] = getFormantData(ihmm1);%ihmm2,ihmm3,ihmm4,ihmm5);
[ih(3,:),jlihmodel] = getFormantData(ihjl1);%,ihjl2,ihjl3,ihjl4,ihjl5);
 [ih(4,:),dhihmodel] = getFormantData(ihdh1);%,ihdh2,ihdh3,ihdh4,ihdh5);

ihAVG = .25 * sum(ih);
ihMAX = max(ih);
ihMIN = min(ih);

% calculate values for ah
[ah(1,:),cpahmodel] = getFormantData(ahcp1);%,ahcp2,ahcp3,ahcp4,ahcp5);
[ah(2,:),mmahmodel] = getFormantData(ahmm1);%,ahmm2,ahmm3,ahmm4,ahmm5);
[ah(3,:),jlahmodel] = getFormantData(ahjl1);%,ahjl2,ahjl3,ahjl4,ahjl5);
[ah(4,:),dhahmodel] = getFormantData(ahdh1);%,ahdh2,ahdh3,ahdh4,ahdh5);

ahAVG = .25 * sum(ah);
ahMAX = max(ah);
ahMIN = min(ah);