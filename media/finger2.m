function [IDX1,IDX2,PVALUE] = finger2(filename)

% Input: MP3 file
% Output: Best matched music file index

%fprintf('Loading input file to memory...\n');
exten = filename(end-2:end);
if exten == '.au'
    [Y, FS, NBITS] = auread(filename);
else
    [Y, FS, NBITS] = mp3read(filename);
end

%fprintf('Sectioning and assigning spectral energies...\n');
SECS = makeframe(Y, FS);
%fprintf('Calculating frequency differences...\n');
fp = framedifference(SECS);

[IDX1, PMAX1, PVALUE] = matchprint2(fp);

PVALUE2 = PVALUE;
PVALUE2(IDX1) = -1;
[PMAX2, IDX2] = max(PVALUE2);

%fprintf('\nResults...\n');
load('fingernames.mat');
fprintf('Best match (score %g): %s\n', round(10*PMAX1)/10, eval(strcat('name_',num2str(IDX1))));
fprintf('Next match (score %g): %s\n', round(10*PMAX2)/10, eval(strcat('name_',num2str(IDX2))));

