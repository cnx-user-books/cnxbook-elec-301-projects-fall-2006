% Taylor Johnson
% Sarah McGee
% Robert Ortman
% Tianhe Yang
% ELEC301
% Project - HDR Image Creation Algorithms
% 2006-12-10

%hdriCreating Creates an HDR image
%param algorithm: string algorithm to use
%param fileNPre: base file name prefix
%param fileNSuf: base file name suffix
%param startIndx: index to start at
%param N: number of input images
function [mRed, mGreen, mBlue] = hdriCreating(algorithm, algorithmParam, fileNPre, fileNSuf, startIndx, N)
    bitrate = algorithmParam(1);

    %iterate over input images
    for (i=startIndx:N)
        %read the ith image: expects filenames like img001.jpg, img002.jpg, img003.jpg, ...
        strNum = strcat(floor(i/100) + '0', floor(i/10) + '0', mod(i,10) + '0');
        [imgData, imgMap] = imread(strcat(fileNPre, strNum, fileNSuf));

        [x y z] = size(imgData);

        reds(:,:,i-startIndx+1) = imgData(:,:,1);
        greens(:,:,i-startIndx+1) = imgData(:,:,2);
        blues(:,:,i-startIndx+1) = imgData(:,:,3);
    end

    switch lower(algorithm)
        case {'createaverage'}
            [imgRed, imgGreen, imgBlue] = createAverage(reds, greens, blues, bitrate);
        case {'createweighted'}
            [imgRed, imgGreen, imgBlue] = createWeighted(reds, greens, blues, bitrate);
        case {'createupsample'}
            [imgRed, imgGreen, imgBlue] = createRandomUpsample(reds, greens, blues, bitrate);
        case {'createpiece'}
            [imgRed, imgGreen, imgBlue] = createPiece(reds, greens, blues, bitrate);
        otherwise
            [imgRed, imgGreen, imgBlue] = createAverage(reds, greens, blues, bitrate);
    end

    %we'll return the composite matrix at this point
    mRed = imgRed;
    mGreen = imgGreen;
    mBlue = imgBlue;
end

%createAverage
%Divides each images' color data by the total number of images and adds them
%together, so each part of each image is given equal weight
function [mRed, mGreen, mBlue] = createAverage(reds, greens, blues, bitrate)
    [x y N] = size(reds);
    %set up output matrices
    pixels = x*y;
    imgRed = double(zeros(x,y));
    imgGreen = double(zeros(x,y));
    imgBlue = double(zeros(x,y));

    %iterate over input images
    for (i=1:N)
        %add luminance values from ith image to running total
        imgRed = imgRed + (double(reds(:,:,i)).*(2^bitrate) ./ (N));
        imgGreen = imgGreen + (double(greens(:,:,i)).*(2^bitrate) ./ (N));
        imgBlue = imgBlue + (double(blues(:,:,i)).*(2^bitrate) ./ (N));
    end

    %return the result
    mRed = imgRed;
    mGreen = imgGreen;
    mBlue = imgBlue;
end

%createWeighted
%Weights each images' color
function [mRed, mGreen, mBlue] = createWeighted(reds, greens, blues, bitrate)
    [x y N] = size(reds);
    %set up output matrices
    pixels = x*y;
    imgRed = (zeros(x,y));
    imgGreen = (zeros(x,y));
    imgBlue = (zeros(x,y));

    bitrateDefault = 8;

    %iterate over input images to find average and maximum luminance values
    for (i=1:N)
        [logAverage(i,:), lWhite(i,:)] = logAverageLuminance(reds(:,:,i), greens(:,:,i), blues(:,:,i));
    end

    lPercentage = logAverage ./ lWhite;

    lExponent = abs(lPercentage.^2 - lPercentage); %quadratic with roots at 0,1

    compositeLogAverage = sum(logAverage)./N;

    compositeMaxAverage = sum(lWhite)./N;

    compositePercentage = compositeLogAverage ./ compositeMaxAverage;

    thresholdMax = 0.9975;
    thresholdMin = 0.020;
    thresholdMiddleMin = 0.275;
    thresholdMiddleMax = 0.525;

    boundLower = 0.75;
    boundUpper = 1.25;

    excludeRed = 0;
    excludeGreen = 0;
    excludeBlue = 0;

    reds = double(reds);
    greens = double(greens);
    blues = double(blues);

    %iterate over input images
    for (i=1:N)
        %check for blownout high-lights
        if (lPercentage(i,1) >= thresholdMax)
            %'high-red'
            excludeRed = excludeRed + 1;
            reds(:,:,i) = reds(:,:,i) .^ 0.9;
        end
        if (lPercentage(i,2) >= thresholdMax)
            %'high-green'
            excludeGreen = excludeGreen + 1;
            greens(:,:,i) = greens(:,:,i) .^ 0.9;
        end
        if (lPercentage(i,3) >= thresholdMax)
            %'high-blue'
            excludeBlue = excludeBlue + 1;
            blues(:,:,i) = blues(:,:,i) .^ 0.9;
        end

        %check for detailless shadows
        if (lPercentage(i,1) <= thresholdMin)
            %'low-red'
            excludeRed = excludeRed + 1;
            reds(:,:,i) = reds(:,:,i) .^ 1.1;
        end
        if (lPercentage(i,2) <= thresholdMin)
            %'low-green'
            excludeGreen = excludeGreen + 1;
            greens(:,:,i) = greens(:,:,i) .^ 1.1;
        end
        if (lPercentage(i,3) <= thresholdMin)
            %'low-blue'
            excludeBlue = excludeBlue + 1;
            blues(:,:,i) = blues(:,:,i) .^ 1.1;
        end

        %enhance middle-tones
        if ((lPercentage(i,1) >= thresholdMiddleMin) && (lPercentage(i,1) <= thresholdMiddleMax))
            %'middle-red'
            reds(:,:,i) = reds(:,:,i) .* 1.1;
        end
        if ((lPercentage(i,2) >= thresholdMiddleMin) && (lPercentage(i,2) <= thresholdMiddleMax))
            %'middle-green'
            greens(:,:,i) = greens(:,:,i) .* 1.1;
        end
        if ((lPercentage(i,3) >= thresholdMiddleMin) && (lPercentage(i,3) <= thresholdMiddleMax))
            %'middle-blue'
            blues(:,:,i) = blues(:,:,i) .* 1.1;
        end

        imgRed = double(imgRed) + ((((double(reds(:,:,i))+(2^bitrateDefault*compositePercentage(1))^(lExponent(i,1)))).*(2^bitrate))+helperRandomIntRange((2^bitrate/2^bitrateDefault)*boundLower, (2^bitrate/bitrateDefault)*boundUpper))./N;
        imgGreen = double(imgGreen) + ((((double((greens(:,:,i)))+(2^bitrateDefault*compositePercentage(2))^(lExponent(i,2)))).*(2^bitrate))+helperRandomIntRange((2^bitrate/2^bitrateDefault)*boundLower, (2^bitrate/bitrateDefault)*boundUpper))./N;
        imgBlue = double(imgBlue) + ((((double((blues(:,:,i)))+(2^bitrateDefault*compositePercentage(3))^(lExponent(i,3)))).*(2^bitrate))+helperRandomIntRange((2^bitrate/2^bitrateDefault)*boundLower, (2^bitrate/bitrateDefault)*boundUpper))./N;
    end

    %return the result
    mRed = double(imgRed);
    mGreen = double(imgGreen);
    mBlue = double(imgBlue);
end

%createRandomUpsample
function [mRed, mGreen, mBlue] = createRandomUpsample(reds, greens, blues, bitrate)
    [xn yn N] = size(reds);
    %set up output matrices
    pixels = xn*yn;

    bitrateDefault = 8; %8 bits per channel default for most image types
    boundLower = 0.75;
    boundUpper = 1.25;

    imgRed = (zeros(xn,yn));
    imgGreen = (zeros(xn,yn));
    imgBlue = (zeros(xn,yn));

    %iterate over input images
    for (i=1:N)
        %iterate over every pixel
        imgRed = (imgRed) + ((double(reds(:,:,i)).*(2^bitrate)+helperRandomIntRange((2^bitrate/2^bitrateDefault)*boundLower, (2^bitrate/bitrateDefault)*boundUpper)))./N;
        imgGreen = (imgGreen) + ((double(greens(:,:,i)).*(2^bitrate)+helperRandomIntRange((2^bitrate/2^bitrateDefault)*boundLower, (2^bitrate/bitrateDefault)*boundUpper)))./N;
        imgBlue = (imgBlue) + ((double(blues(:,:,i)).*(2^bitrate)+helperRandomIntRange((2^bitrate/2^bitrateDefault)*boundLower, (2^bitrate/bitrateDefault)*boundUpper)))./N;
    end

    %return the result
    mRed = double(imgRed);
    mGreen = double(imgGreen);
    mBlue = double(imgBlue);
end

%createPiece
function [mRed, mGreen, mBlue] = createPiece(reds, greens, blues, bitrate)

end

%logAverageLuminance
function [logAverage, Lwhite] = logAverageLuminance(inputRed, inputGreen, inputBlue)
    matrix = cat(3,double(inputRed),double(inputGreen),double(inputBlue));

    for (index=1:1:3)
        Lw = matrix(:,:,index);
        [x y] = size(Lw);
        N = x*y; % number of pixels in matrix
        delta = 0.00001; % some small values for taking the log of a black pixel
        sum = 0;

        % finds the log-average luminance
        for (n = 1:x)
            for m = 1:y
                sum = sum + log(delta + Lw(n,m));
            end
        end
        logAverage(index) = exp(sum / N);

        % maximum luminance in scene to map to pure white
        Lwhite(index) = 0;
        for (n = 1:x)
            for (m = 1:y)
                if (matrix(n,m,index) > Lwhite(index))
                    Lwhite(index) = matrix(n,m,index);
                end
            end
        end
    end
end

%
function [rr] = helperRandomIntRange(min, max)
    rmax = (rand()*(max*(2/3)-1)) + min;
    rr = round(rmax-min);
end
