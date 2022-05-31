% Taylor Johnson
% Sarah McGee
% Robert Ortman
% Tianhe Yang
% ELEC301
% Project - HDR Image Tests
% 2006-12-10

%Tests of combinations of creation and tone-mapping algorithms
function result=hdriTest()
    %input parameters: hdri create operator name, create operator parameters, tone mapping operator name,
    %                  tone mapping operator parameters, input file prefix, input file suffix,
    %                  input file start index, input file end index
    hdriMain('createAverage', [32], 'localconvolution', [3], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'localaverage', [2,5], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'localedge', [2,5,1.1,1], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'localstochastic', [2,5,10], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'localEdgeConvolution', [2,5], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'wavelet', [2,5,2,0.005,0.005], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'global', [2], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createAverage', [32], 'sampling', [32,8], 'images/memorial', 's.jpg', 1, 21);

    hdriMain('createUpsample', [32], 'localconvolution', [3], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'localaverage', [2,5], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'localedge', [2,5,1.1,1], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'localstochastic', [2,5,10], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'localEdgeConvolution', [2,5], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'wavelet', [2,5,2,0.005,0.005], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'global', [2], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createUpsample', [32], 'sampling', [32,8], 'images/memorial', 's.jpg', 1, 21);

    hdriMain('createweighted', [32], 'localconvolution', [3], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'localaverage', [2,5], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'localedge', [2,5,1.1,1], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'localstochastic', [2,5,10], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'localEdgeConvolution', [2,5], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'wavelet', [2,5,2,0.005,0.005], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'global', [2], 'images/memorial', 's.jpg', 1, 21);
    hdriMain('createweighted', [32], 'sampling', [32,8], 'images/memorial', 's.jpg', 1, 21);
end
