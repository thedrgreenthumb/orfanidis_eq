clear all;
close all;

%Check os
if(ispc)
    error('Windows based OS not supported.');
end;

%TODO: Fix it.
% The test programm call modification is required
exit(1)

%Include custom MATLAB Toolkit
addpath('./MDSPTK');

%Build configuration file
cfgFileName = 'orfanidis_eq.tstdat';
cfgFileMessage = '/* sample rate in Hz, number of bands, testing data vectors length */';
cfgSampleRateHz = 48000;
cfgNumberOfBands = 30; %Do not change it
cfgTestDataVectorsLength = cfgSampleRateHz/4;

cfgFId = fopen(cfgFileName,'w');
fprintf(cfgFId, '%s', cfgFileMessage);
fprintf(cfgFId, '%d, %d, %d\n', cfgSampleRateHz, cfgNumberOfBands, cfgTestDataVectorsLength);
fclose(cfgFId);

%Build gains configuration file
cfgGainsFileName = 'orfanidis_eq_gain.tstdat';
cfgGainsFileMessage = '/* gains for 1/3 octave eq, 30 values in dB */';
gainsdB = [0,0,0,0,3,0,5,0,0,16,0,0,0,0,0,0,0,0,0,0,0,-10,0,0,-5,0,0,-2,0,0];
csvwrite(cfgGainsFileName,gainsdB);

%Build and run test file
system('make');
system('./eq');

%Configure filter analizer from MDSPTK
fa = fanalyzer(cfgSampleRateHz);


%%% ---- eq1 class tests ----
fa.freqResp(csvread('butterworth_eq1.tstdat'),'log','eq1: butterworth');
fa.freqResp(csvread('chebyshev1_eq1.tstdat'),'log','eq1: chebyshev1');
fa.freqResp(csvread('chebyshev2_eq1.tstdat'),'log','eq1: chebyshev2');

%%% ---- eq2 class tests ----
fa.freqResp(csvread('butterworth_eq2.tstdat'),'log','eq2: butterworth');
fa.freqResp(csvread('chebyshev1_eq2.tstdat'),'log','eq2: chebyshev1');
fa.freqResp(csvread('chebyshev2_eq2.tstdat'),'log','eq2: chebyshev2');

%Clean up all
system('make clean');

