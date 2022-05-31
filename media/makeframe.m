function SECTIONS = makeframe (Y, FS)

% Input the raw signal array and sampling rate
% output: 25 critical frequency bands x number of timeframes
% output values are energy strength

Y = convertmono(Y);
len = length(Y);
window = round(0.37.*FS); % hanning window width
interval = round(0.0116.*FS); % 512
%interval = 1024;
len = len - window;
numintervals = len / interval;
kk=1;
hanw = hanning(window);


%SECTIONS = zeros(interval,numintervals);
for ii = 1:interval:len
%    if ii+interval > len
%        interval = len-ii+1;
%        disp('coming to end');
%        disp(interval);
%        disp(ii);
%        SECTIONS(ii:ii+interval-1,kk)=Y(ii:ii+interval-1);
%    else
       % SECTIONS(:,kk)=Y(ii:ii+interval-1);
      % (8:25pm)  SECTIONS(:,kk) = dividefreqs(Y(ii:ii+interval-1),FS).';
%      disp(ii)
%    fprintf('Grabbing interval %g to %g\n',ii,ii+window-1);
    SECTIONS(:,kk) = dividefreqs(Y(ii:ii+window-1).*hanw,FS).';
%    end

%    disp(size(SECTIONS));
    kk = kk+1;
end

