function w = hanning(n)
load('han.mat');
switch n
    case 2220
        w = h6000;
    case 2960
        w = h8000;
    case 4079
        w = h11025;
    case 5920
        w = h16000;
    case 8159
        w = h22050;
    case 11840
        w = h32000;
    case 16317
        w = h44100;
end

