function [IDX, PMAX, PVALUE] = matchprint2 (fp)

% Load database file
%fprintf('Loading audio fingerprint database...\n');
load('fingerprints_all.mat');

% Find longest fingerprint
%fprintf('Finding longest fingerprint...\n');
for ii = 1:fp_dbsize
    fp_len(ii) = length(eval(strcat('fp_',num2str(ii))));
end
maxlen = max(fp_len);

if maxlen < length(fp)
    maxlen = length(fp);
end

% Pre-determine strength ratio
ratio = length(fp) / maxlen;
if ratio > 1
    ratio = 1/ratio;
end

% Zero pad, change 0 to -1, normalize, rotate twice, fft
%fprintf('Processing input file...\n');
temp = zeros(24,maxlen);
temp(1:24, 1:size(fp,2)) = fp.*2-1;
temp = temp ./ matnorm(temp);
fp_fft = fft2(rot90(rot90(temp)));

% Zero pad, change 0 to -1, normalize, fft
% and match on the fly
%fprintf('Processing database and matching...\n');
for ii = 1:fp_dbsize
   % fprintf('%g',ii);
    temp = zeros(24, maxlen);
    temp(1:24, 1:size(eval(strcat('fp_',num2str(ii))),2)) = eval(strcat('fp_',num2str(ii))) .* 2-1;
    temp = temp ./ matnorm(temp);
    fprintf('.');
    PVALUE(ii) = max(max(ifft2(fp_fft .* fft2(temp))));
end

PVALUE = abs(PVALUE)/ratio;
[PMAX,IDX] = max(PVALUE);



