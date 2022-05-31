% Taylor Johnson
% Sarah McGee
% Robert Ortman
% Tianhe Yang
% ELEC301
% Project - HDR Image Tone Mapping Main Inteface and Operator Algorithms
% 2006-12-10

%hdriToneMapping
%Used to abstract the other algorithms away from the main interface
function [imgRed, imgGreen, imgBlue]=hdriToneMapping(algorithm, algorithmParam, red, green, blue)
    bitrateHDR = 32;
    bitrateLDR = 8;

    switch lower(algorithm)
        case {'localconvolution'}
            imgRed = hdriToneOperator_LocalConvolution(red, algorithmParam(1));
            imgGreen =  hdriToneOperator_LocalConvolution(green, algorithmParam(1));
            imgBlue = hdriToneOperator_LocalConvolution(blue, algorithmParam(1));
        case {'localaverage'}
            imgRed = hdriToneOperator_LocalAverage(red, algorithmParam(1), algorithmParam(2));
            imgGreen = hdriToneOperator_LocalAverage(green, algorithmParam(1), algorithmParam(2));
            imgBlue = hdriToneOperator_LocalAverage(blue, algorithmParam(1), algorithmParam(2));
        case {'localedge'}
            imgRed = hdriToneOperator_LocalEdge(red, algorithmParam(1), algorithmParam(2), algorithmParam(3), algorithmParam(4));
            imgGreen = hdriToneOperator_LocalEdge(green, algorithmParam(1), algorithmParam(2), algorithmParam(3), algorithmParam(4));
            imgBlue = hdriToneOperator_LocalEdge(blue, algorithmParam(1), algorithmParam(2), algorithmParam(3), algorithmParam(4));
        case {'localstochastic'}
            imgRed = hdriToneOperator_LocalStochastic(red, algorithmParam(1), algorithmParam(2), algorithmParam(3));
            imgGreen = hdriToneOperator_LocalStochastic(green, algorithmParam(1), algorithmParam(2), algorithmParam(3));
            imgBlue = hdriToneOperator_LocalStochastic(blue, algorithmParam(1), algorithmParam(2), algorithmParam(3));
        case {'wavelet'}
            imgRed = hdriToneOperator_LocalAverage(red, algorithmParam(1), algorithmParam(2));
            imgGreen = hdriToneOperator_LocalAverage(green, algorithmParam(1), algorithmParam(2));
            imgBlue = hdriToneOperator_LocalAverage(blue, algorithmParam(1), algorithmParam(2));
            imgRed = helperWaveletBlur(imgRed, algorithmParam(3), algorithmParam(4), algorithmParam(5));
            imgGreen = helperWaveletBlur(imgGreen, algorithmParam(3), algorithmParam(4), algorithmParam(5));
            imgBlue = helperWaveletBlur(imgBlue, algorithmParam(3), algorithmParam(4), algorithmParam(5));
        case {'localedgeconvolution'}
            imgRed = hdriToneOperator_LocalEdgeConvolution(red, algorithmParam(1), algorithmParam(2));
            imgGreen = hdriToneOperator_LocalEdgeConvolution(green, algorithmParam(1), algorithmParam(2));
            imgBlue = hdriToneOperator_LocalEdgeConvolution(blue, algorithmParam(1), algorithmParam(2));
        case {'global'}
            [imgRed, imgGreen, imgBlue, logAverage, Lwhite] = hdriToneOperator_Global(red, green, blue, algorithmParam(1));
        case {'sampling'}
            imgRed = hdriToneOperator_Sampling(red, algorithmParam(1), algorithmParam(2));
            imgGreen = hdriToneOperator_Sampling(green, algorithmParam(1), algorithmParam(2));
            imgBlue = hdriToneOperator_Sampling(blue, algorithmParam(1), algorithmParam(2));
        otherwise
            [imgRed, imgGreen, imgBlue, logAverage, Lwhite] = hdriToneOperator_Global(red, green, blue, algorithmParam(1));
    end
end

%hdriToneOperator_Global
%param inputMatrix: red, green, or blue input matrix
function [outRed, outGreen, outBlue, logAverage, Lwhite] = hdriToneOperator_Global(inputRed, inputGreen, inputBlue, brightness)
    matrix = cat(3,inputRed,inputGreen,inputBlue);

    for (index=1:1:3)
        Lw = matrix(:,:,index);
        [x y] = size(Lw);
        N = x*y;    % number of pixels in matrix
        a = 0.18*brightness;   %   key value, determined by user
        delta = 0.00001;    % some small values for taking the log of a black pixel
        sum = 0;

        % finds the log-average luminance
        for (n = 1:x)
            for m = 1:y
                sum = sum + log(delta + Lw(n,m));
            end
        end
        logAverage = exp(sum / N);

        % find scaled luminance
        matrixNew(:,:,index) = ones(x,y,1);
        for (n = 1:x)
            for (m = 1:y)
                matrixNew(n,m,index) = a*Lw(n,m)/logAverage;
            end
        end

        % maximum luminance in scene to map to pure white
        Lwhite = 0;
        for (n = 1:x)
            for (m = 1:y)
                if (matrixNew(n,m,index) > Lwhite)
                    Lwhite = matrixNew(n,m,index);
                end
            end
        end

        % tone mapping operator
        for (n = 1:x)
            for (m = 1:y)
                matrixNew(n,m,index) = (matrixNew(n,m,index)*(1 + matrixNew(n,m,index)/(Lwhite^2)))/(1 + matrixNew(n,m,index));
            end
        end

        % ensure no values are > 1.0 or < 0.0 due to rounding errors
        [xn yn zn] = size(matrixNew);
        for (i=1:xn)
            for (j=1:yn)
                if (matrixNew(i,j,index) > 1.0)
                    matrixNew(i,j,index) = 1.0;
                elseif (matrixNew(i,j,index) < 0.0)
                    matrixNew(i,j,index) = 0.0;
                end
            end
        end
    end

    outRed = matrixNew(:,:,1);
    outGreen = matrixNew(:,:,2);
    outBlue = matrixNew(:,:,3);
end

%hdriToneOperator_LocalConvolution
function [Lnew] = hdriToneOperator_LocalConvolution(Lw, brightness)
    [x y] = size(Lw);
    N = x*y;    % number of pixels in matrix
    a = 0.18*brightness;   %   key value, determined by user
    delta = 0.00001;    % some small values for taking the log of a black pixel
    sum = 0;
    Lnew = zeros(x,y);  %The final converted luminance (L_d)
    Lf = zeros(x,y);    %Var for the fft of L
    R1 = zeros(x,y);    %Temp var for FFT of Gaussian at i=1, s
    R2 = zeros(x,y);    %Temp var for FFT of Gaussian at i=2, s
    R_1 = zeros(x,y,8);  %Gaussian for i = 1
    R_2 = zeros(x,y,8); %Gaussian for i = 2
    V_1 = zeros(x,y,8); %Convolution result for i = 1
    V_2 = zeros(x,y,8); %Convolution result for i = 2
    S_m = zeros(x,y);  %Matrix of scales computed for each pixel
    scale = 1.6;
    alpha_1 = 1 / (2*sqrt(2));
    alpha_2 = alpha_1 * scale;
    phi = 8.0;
    eps = 0.05;


    % finds the log-average luminance
    for (n = 1:x)
        for m = 1:y
            sum = sum + log(delta + Lw(n,m));
        end
    end
    logAverage = exp(sum / N);

    % find scaled luminance
    L = ones(x,y);
    for (n = 1:x)
        for (m = 1:y)
            L(n,m) = a*Lw(n,m)/logAverage;
        end
    end


    %Method using convolution-based local adaption:

    for k = 1:8
        S(k) = round(scale^k);
    end

    %R_i(x,y,s):  Can be implemented as separate function
    for s = 1:8
        for n = 1:x
            for m = 1:y
                R_1(n,m,s) = (1 / (pi*(alpha_1*S(s))^2)) * exp(-1*((n^2+m^2) / ((alpha_1*S(s))^2)));
                R_2(n,m,s) = (1 / (pi*(alpha_2*S(s))^2)) * exp(-1*((n^2+m^2) / ((alpha_2*S(s))^2)));
            end
        end
    end

    for s = 1:8
        Lf = fft2(L);
        R1 = fft2(R_1(:,:,s));
        R2 = fft2(R_2(:,:,s));
        V_1(:,:,s) = ifft2(Lf .* R1);
        V_2(:,:,s) = ifft2(Lf .* R2);
    end

    for n = 1:x
        for m = 1:y
            s = 1;
            V = (V_1(n,m,s) - V_2(n,m,s)) / (2*phi*a/(S(s)^2) + V_1(n,m,s));
            while (s < 8)
                s = s + 1;
                V = (V_1(n,m,s) - V_2(n,m,s)) / (2*phi*a/(S(s)^2) + V_1(n,m,s));
                S_m(n,m) = s;
                if (abs(V) >= eps)  %If threshold of eps is exceeded take previous s value and set s to 10 to terminate while loop
                    S_m(n,m) = s - 1;
                    s = 10;
                end
            end
        end
    end

    for n = 1:x
        for m = 1:y
            Lnew(n,m) = L(n,m) / (1 + V_1(n,m,S_m(n,m)));
        end
    end

    % ensure no values are > 1.0 or < 0.0 due to rounding errors
    [xn yn] = size(Lnew);
    for (i=1:xn)
        for (j=1:yn)
            if (Lnew(i,j) > 1.0)
                Lnew(i,j) = 1.0;
            elseif (Lnew(i,j) < 0.0)
                Lnew = 0.0;
            end
        end
    end
end

%hdriToneOperator_LocalStochastic
%Iterates across the image with a window of specified size, randonly
%selecting pixels within that window to be averaged with respect to
%luminance.  This is similar to distributed ray-tracing.
function [Lnew] = hdriToneOperator_LocalStochastic(Lw, area, contrastLevel, truncatedLength)
    [x y] = size(Lw);
    localArea = area; % take the average of the values in a given area
    contrastLimit = contrastLevel;  % ignore vales too bright or too dark to avoid halos, value must be >1
    Lnew = ones(x,y); %  output image
    Lwhite = 0;
    logOrig = log(Lw)/log(contrastLimit);

    for n = 1:x
        for m = 1:y

            nei = helperFindNeighbors_Stochastic(Lw,n,m,localArea,truncatedLength);

            logNei = log(nei)/log(contrastLimit);

            logDiff = logNei - logOrig(n,m);

            weight = exp(-abs(logDiff).^(5*contrastLimit));

            weightedMean = sum(logDiff.*weight)/sum(weight);

            Lnew(n,m)=exp((weightedMean+logOrig(n,m)).*log(contrastLimit));

            if (Lnew(n,m) > Lwhite)
                Lwhite = Lnew(n,m);
            end

            % Lnew(n,m) = (Lnew(n,m)*(1 + Lnew(n,m)/(Lwhite^2)))/(1 + Lnew(n,m));

        end
    end

    [xn yn] = size(Lnew);
    for (i=1:xn)
        for (j=1:yn)
            Lnew(i,j) = Lnew(i,j)/Lwhite;
            if (Lnew(i,j) > 1.0)
                Lnew(i,j) = 1.0;
            elseif (Lnew(i,j) < 0.0)
                Lnew = 0.0;
            end
        end
    end
end

%hdriToneOperator_LocalEdge
%
function [Lnew] = hdriToneOperator_LocalEdge(Lw, area, contrastLevel, MaxDiff, MaxNumDiff)
    [x y] = size(Lw);
    localArea = area; % take the average of the values in a given area
    contrastLimit = contrastLevel;  % ignore vales too bright or too dark to avoid halos, value must be >1
    Lnew = ones(x,y); %  output image
    Lwhite = 0;
    logOrig = log(Lw)/log(contrastLimit);

    for n = 1:x
        for m = 1:y
            [difference,nei] = helperFindNeighbors_Edge(Lw,n,m,localArea, MaxDiff, MaxNumDiff);

            if(difference == 0)

                logNei = log(nei)/log(contrastLimit);

                logDiff = logNei - logOrig(n,m);

                weight = exp(-abs(logDiff).^(5*contrastLimit));

                weightedMean = sum(logDiff.*weight)/sum(weight);

                Lnew(n,m)=exp((weightedMean+logOrig(n,m)).*log(contrastLimit));

                if (Lnew(n,m) > Lwhite)
                    Lwhite = Lnew(n,m);
                end

            else
                Lnew(n,m) = Lw(n,m);
            end
        end
    end

    [xn yn] = size(Lnew);
    for (i=1:xn)
        for (j=1:yn)
            Lnew(i,j) = Lnew(i,j)/Lwhite;
            if (Lnew(i,j) > 1.0)
                Lnew(i,j) = 1.0;
            elseif (Lnew(i,j) < 0.0)
                Lnew = 0.0;
            end
        end
    end
end

%hdriToneOperator_LocalAverage
%
function [Lnew] = hdriToneOperator_LocalAverage(Lw, area, contrastLevel)
    [x y] = size(Lw);
    localArea = area; % take the average of the values in a given area
    contrastLimit = contrastLevel;  % ignore vales too bright or too dark to avoid halos, value must be >1
    Lnew = ones(x,y); %  output image
    Lwhite = 0;
    logOrig = log(Lw)/log(contrastLimit);

    for n = 1:x
        for m = 1:y
            nei = helperFindNeighbors_Average(Lw,n,m,localArea);

            logNei = log(nei)/log(contrastLimit);

            logDiff = logNei - logOrig(n,m);

            weight = exp(-abs(logDiff).^(5*contrastLimit));

            weightedMean = sum(logDiff.*weight)/sum(weight);

            Lnew(n,m)=exp((weightedMean+logOrig(n,m)).*log(contrastLimit));

            if (Lnew(n,m) > Lwhite)
                Lwhite = Lnew(n,m);
            end
        end
    end

    [xn yn] = size(Lnew);
    for (i=1:xn)
        for (j=1:yn)
            Lnew(i,j) = Lnew(i,j)/Lwhite;
            if (Lnew(i,j) > 1.0)
                Lnew(i,j) = 1.0;
            elseif (Lnew(i,j) < 0.0)
                Lnew = 0.0;
            end
        end
    end
end

%
function [Lnew] = hdriToneOperator_LocalEdgeConvolution(Lw, Car, Gw)
    [x y] = size(Lw);
    localArea = 1; % take the weighted average of the values in a given area
    Lnew = ones(x,y); %  output image
    Lwhite = 1;
    Gau = helperGaussian(Car, Gw);
    weight = sum(sum(Gau));
    convSiz = (2*(2*Car + 1)-1)^2;

    for n = 1:x
        for m = 1:y
            s = 0;

            [difference, nei] = helperFindNeighbors_EdgeConvolution(Lw,n,m,localArea,2.8,1);

            if (difference == 1)

                conved = helperConv2df(nei,Gau);

                convedSiz = size(conved);

                for e = 1:convedSiz(1)
                    for f = 1:convedSiz(2)
                        s = s + conved(e,f);
                    end
                end

                average = s/(convSiz+weight);

                Lnew(n,m) = average;

                if (Lnew(n,m) > Lwhite)
                    Lwhite = Lnew(n,m);
                end

            else
                Lnew(n,m) = Lw(n,m);
                if (Lnew(n,m) > Lwhite)
                    Lwhite = Lnew(n,m);
                end
            end
        end
    end

    [xn yn] = size(Lnew);
    for (i=1:xn)
        for (j=1:yn)
            Lnew(i,j) = Lnew(i,j)/(Lwhite);
            if (Lnew(i,j) > 1.0)
                Lnew(i,j) = 1.0;
            elseif (Lnew(i,j) < 0.0)
                Lnew(i,j) = 0.0;
            end
        end
    end
end

%Squishing function to original bitrate
%This gives a value between 0 and 1 in discrete steps of 2^bitrateLDR
%That is, it will map from a space 2^32 big into a space 2^8 big, assigning
%groups of values from the 2^32 space to signle values in the 2^8 space.
function [Lnew] = hdriToneOperator_Sampling(Lw, bitrateHDR, bitrateLDR)
    Lnew = (ceil(Lw./(2^bitrateHDR-1))-1)./(2^bitrateLDR-1);

    % ensure no values are > 1.0 or < 0.0 due to rounding errors
    [xn yn zn] = size(Lnew);
    for (i=1:xn)
        for (j=1:yn)
            if (Lnew(i,j) > 1.0)
                Lnew(i,j) = 1.0;
            elseif (Lnew(i,j) < 0.0)
                Lnew(i,j) = 0.0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isTooDifferent, neighbors] = helperFindNeighbors_Edge(img, row, col, area, MaxDiff, MaxNumDiff)
    [x,y] = size(img);
    neighbors = zeros((area-1)^2,1);
    isTooDifferent = 0;

    c = 1;
    numTooDifferent = 0;
    for n = row-area+1:row+area
        for m = col-area+1:col+area
            if (n>0 && n<x && m>0 && m<y)
                neighbors(c) = img(n,m);
                if (neighbors(c) < img(row, col)/MaxDiff || neighbors(c) > img(row, col)*MaxDiff)
                    if(numTooDifferent > MaxNumDiff)
                        isTooDifferent = 1;
                    end
                    numTooDifferent = numTooDifferent + 1;
                end
                c = c+1;
            end
        end
    end
end

%
function [neighbors] = helperFindNeighbors_Stochastic(img, row, col, area, truncatedLength)
    [x,y] = size(img);
    neighbors = [];
    index = 1;
    for n = 1:truncatedLength,
        r = helperRandomInt(2*area - 1);
        c = helperRandomInt(2*area - 1);
        if (row-area+1+r>0 && row-area+1+r<x && col-area+1+c>0 && col-area+1+c<y),
            neighbors(index) = img(row-area+1+r,col-area+1+c);
            index = index + 1;
        end
    end
    if (isempty(neighbors)),
        neighbors(1) = img(row,col);
    end
end

%
function [neighbors] = helperFindNeighbors_Average(img, row, col, area)
    [x,y] = size(img);

    c = 1;
    for n = row-area+1:row+area
        for m = col-area+1:col+area
            if (n>0 && n<x && m>0 && m<y)
            neighbors(c) = img(n,m);
            c = c+1;
            end
        end
    end
end

%
function [isTooDifferent, neighbors] = helperFindNeighbors_EdgeConvolution(img, row, col, Ar, MaxDiff, MaxNumDiff)
    [x,y] = size(img);
    neighbors = zeros((Ar-1)^2,1);
    isTooDifferent = 0;

    c = 1;
    numTooDifferent = 0;
    for n = row-Ar:row+Ar
        for m = col-Ar:col+Ar
            if (n>0 && n<x && m>0 && m<y)
                neighbors(c) = img(n,m);
                if (neighbors(c) < img(row, col)/MaxDiff || neighbors(c) > img(row, col)*MaxDiff)
                    if(numTooDifferent > MaxNumDiff)
                        isTooDifferent = 1;
                    end
                    numTooDifferent = numTooDifferent + 1;
                end
                c = c+1;
            end
        end
    end
end

%
function [rr] = helperRandomInt(range)
    r = (rand()*(range-1)) + 1;
    rr = round(r);
end

%helperWaveletBlur
% Wavelet implementation of an unsharp filter bank
% mappedInput - all values should be probably already be between 0 and 1
% filterBankLevel - choice between 1, 2, and 3
% threshold - the difference between mappedInput and the blurred image
% weightedDifference - by how much to make the dark pixels darker and the
% light pixels lighter
function [Lnew] = helperWaveletBlur(mappedInput, filterBankLevel, threshold, weightedDifference)
    %blah
    blur = helperWaveletFilterBank(mappedInput, filterBankLevel);
    diff = mappedInput - blur;
    [x y] = size(mappedInput);
    Lnew = ones(x,y);
    for m = 1:x,
        for n = 1:y,
            if (abs(diff(m,n)) > threshold),
                if ((mappedInput(m,n) < 0.5) && (mappedInput(m,n) > weightedDifference)),
                    Lnew(m,n) = mappedInput(m,n) - weightedDifference;
                elseif ((mappedInput(m,n) >= 0.5) && (mappedInput(m,n) < (1-weightedDifference))),
                    Lnew(m,n) = mappedInput(m,n) + weightedDifference;
                else
                    Lnew(m,n) = mappedInput(m,n);
                end
            else
                Lnew(m,n) = mappedInput(m,n);
            end
        end
    end
end

%helperWaveletFilterBank
%Implements levels 1, 2, and 3 discrete wavelet transform (DWT)
function [lowestSubband] = helperWaveletFilterBank(mappedInput, filterBankLevel)
    [x y] = size(mappedInput);
    lowestSubband = ones(x,y);
    switch(filterBankLevel)
        case 3,
            [cA,cH,cV,cD] = dwt2(mappedInput,'haar');
            [cA2,cH2,cV2,cD2] = dwt2(cA,'haar');
            [cA3,cH3,cV3,cD3] = dwt2(cA2,'haar');

            [x31 y31] = size(cH3);
            [x32 y32] = size(cV3);
            [x33 y33] = size(cD3);
            cH3 = zeros(x31,y31);
            cV3 = zeros(x32,y32);
            cD3 = zeros(x33,y33);
            cA2 = idwt2(cA3,cH3,cV3,cD3,'haar');

            [x21 y21] = size(cH2);
            [x22 y22] = size(cV2);
            [x23 y23] = size(cD2);
            cH2 = zeros(x21,y21);
            cV2 = zeros(x22,y22);
            cD2 = zeros(x23,y23);
            cA = idwt2(cA2,cH2,cV2,cD2,'haar');

            [x1 y1] = size(cH);
            [x2 y2] = size(cV);
            [x3 y3] = size(cD);
            cH = zeros(x1,y1);
            cV = zeros(x2,y2);
            cD = zeros(x3,y3);
            lowestSubband = idwt2(cA,cH,cV,cD,'haar');
        case 2,
            [cA,cH,cV,cD] = dwt2(mappedInput,'haar');
            [cA2,cH2,cV2,cD2] = dwt2(cA,'haar');

            [x21 y21] = size(cH2);
            [x22 y22] = size(cV2);
            [x23 y23] = size(cD2);
            cH2 = zeros(x21,y21);
            cV2 = zeros(x22,y22);
            cD2 = zeros(x23,y23);
            cA = idwt2(cA2,cH2,cV2,cD2,'haar');

            [x1 y1] = size(cH);
            [x2 y2] = size(cV);
            [x3 y3] = size(cD);
            cH = zeros(x1,y1);
            cV = zeros(x2,y2);
            cD = zeros(x3,y3);
            lowestSubband = idwt2(cA,cH,cV,cD,'haar');
        case 1,
            [cA,cH,cV,cD] = dwt2(mappedInput,'haar');
            [x1 y1] = size(cH);
            [x2 y2] = size(cV);
            [x3 y3] = size(cD);
            cH = zeros(x1,y1);
            cV = zeros(x2,y2);
            cD = zeros(x3,y3);
            lowestSubband = idwt2(cA,cH,cV,cD,'haar');
        otherwise
            [cA,cH,cV,cD] = dwt2(mappedInput,'haar');
            [x1 y1] = size(cH);
            [x2 y2] = size(cV);
            [x3 y3] = size(cD);
            cH = zeros(x1,y1);
            cV = zeros(x2,y2);
            cD = zeros(x3,y3);
            lowestSubband = idwt2(cA,cH,cV,cD,'haar');
    end
end

%
function [Gau] = helperGaussian(Ar, Gw)
    Gau = zeros(2*Ar+1,2*Ar+1);
    %[numRows, numCols] = size(img);
    guassWeight = Gw;
    siz = (2*Ar+1)^2;

    for n = 1:2*Ar+1
        for m = 1:2*Ar+1

            Gau(n,m) = exp(-guassWeight*sqrt((n-(Ar+1))^2+(m-(Ar+1))^2));

        end
    end
end

%
function [Z] = helperConv2df(x, y)
    sizex = size(x);
    sizey = size(y);
    rows = sizex(1) + sizey(1) - 1;
    cols = sizex(2) + sizey(2) - 1;
    Z = ifft2(fft2(x,rows,cols) .* fft2(y,rows,cols));
end
