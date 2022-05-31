function SECENERGY = dividefreqs(X, FS)

bandsec = [
    1	100
    100	200
    200	300
    300	400
    400	510
    510	630
    630	770
    770	920
    920	1080
    1080	1270
    1270	1480
    1480	1720
    1720	2000
    2000	2320
    2320	2700
    2700	3150
    3150	3700
    3700	4400
    4400	5300
    5300	6400
    6400	7700
    7700	9500
    9500	12000
    12000	15500
    15500	20000];

FY = fft(X,44100);
FYP = FY(1:44100/2);
len = length(bandsec);

for ii=1:len
    SECENERGY(ii) = energy(FYP(bandsec(ii,1):bandsec(ii,2)));
end

SECENERGY = abs(SECENERGY);



