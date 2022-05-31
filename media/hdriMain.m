% Taylor Johnson
% Sarah McGee
% Robert Ortman
% Tianhe Yang
% ELEC301
% Project - HDR Image Main Interface
% 2006-12-10

%hdriMain
%  This is the main interface to the outside world, and provides a common
%  gateway to create and then tonemap an image, given multiple LDR images.
%param algorithmCreate: string algorithm name to use for hdri creation
%param algorithmCreateParam: vector of input parameters to the creation algorithm
%param algorithmToneMap: string algorithm name to use for tone mapping
%param algorithmToneMapParam: vector of input parameters to the tone mapping algorithm
%param fileNPre: base file name prefix
%param fileNSuf: base file name suffix
%param startIndx: index to start at (3 instead of 1 for example)
%param N: number of input images
function result=hdriMain(algorithmCreate, algorithmCreateParam, algorithmToneMap, algorithmToneMapParam, fileNPre, fileNSuf, startIndx, N)
    bitrate = algorithmCreateParam(1);
    bitrateLDR = 8;

    %perform hdr image creation
    [imgRed imgGreen imgBlue] = hdriCreating(algorithmCreate, algorithmCreateParam, fileNPre, fileNSuf, startIndx, N);

    %perform tone mapping
    [imgRedMap imgGreenMap imgBlueMap] = hdriToneMapping(algorithmToneMap, algorithmToneMapParam, imgRed, imgGreen, imgBlue);

    %downsample to 8-bit
    imgRedMap = ceil(imgRedMap.*(2^bitrateLDR-1))-1;
    imgGreenMap = ceil(imgGreenMap.*(2^bitrateLDR-1))-1;
    imgBlueMap = ceil(imgBlueMap.*(2^bitrateLDR-1))-1;

    %convert back to 0.0 to 1.0 range
    imgRedMap = imgRedMap./(2^bitrateLDR-1);
    imgGreenMap = imgGreenMap./(2^bitrateLDR-1);
    imgBlueMap = imgBlueMap./(2^bitrateLDR-1);

    %concatenate tone-mapped arrays together into a single matrix (to use the image processing toolkit)
    cDisp = cat(3,imgRedMap,imgGreenMap,imgBlueMap);

    %display the image
    figure;
    image(cDisp);

    %write the new image to a file
    filename = strcat(fileNPre, 'c', '-', algorithmCreate);
    [n1 n2] = size(algorithmCreateParam);
    for (i=1:n2)
        filename = strcat(filename, 'p', num2str(i), '=', num2str(algorithmCreateParam(i)));
    end
    filename = strcat(filename, '-', algorithmToneMap);
    size(algorithmToneMapParam)
    [n1 n2] = size(algorithmToneMapParam);
    for (i=1:n2)
        filename = strcat(filename, 'p', num2str(i), '=', num2str(algorithmToneMapParam(i)));
    end
    filename = strcat('results/', filename, '.jpg');
    imwrite(cDisp, filename, 'Quality', 100, 'BitDepth', 8);

    %test tonal range
    %y = imread(filename);
    %[cL tL cH tH tLcL tHcH] = compareHistogram(y(:,:,1), y(:,:,1))
end

%compareHistogram
%Used to compare the histograms of two images to see which 
%contains more data
function [cL tL cH tH tLcL tHcH] = compareHistogram(imBefore, imAfter)
    bitrateL = 8;
    bitrateH = 12;
    hL = hist(double(imBefore),2^bitrateL);
    hH = hist(double(imBefore),2^bitrateH);

    [xhn yhn] = size(hL);
    c = 0;
    t = 0;
    for (x=1:xhn)
        not_zero = 0;
        for (y=1:yhn)
            if (hL(x,y) <= 0 && not_zero == 1)
                not_zero = 0;
            else
                not_zero = 1;
            end
        end

        if (not_zero == 0)
            c = c + 1;
        end
        t = t + 1;
    end
    cL = c;
    tL = t;

    [xhn yhn] = size(hH);
    c = 0;
    t = 0;
    for (x=1:xhn)
        not_zero = 0;
        for (y=1:yhn)
            if (hH(x,y) <= 0 && not_zero == 1)
                not_zero = 0;
            else
                not_zero = 1;
            end
        end

        if (not_zero == 0)
            c = c + 1;
        end
        t = t + 1;
    end
    cH = c;
    tH = t;

    tHcH = tH-cH;
    tLcL = tL-cL;
end
