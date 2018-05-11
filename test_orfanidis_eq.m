clear all;
close all;

%Check os
if(ispc)
    error('Windows based OS is not supported.');
end

%Include custom MATLAB Toolkit
addpath('./MDSPTK');

%Build configuration file
sampleRateHz = 48000;
numberOfBands = 30; %Do not change it
testDataVectorsLength = sampleRateHz / 4;
gainsdB = '0 0 0 0 3 0 5 0 0 16 0 0 0 0 0 0 0 0 0 0 0 -10 0 0 -5 0 0 -2 0 0';

%Build and run test file
system('make');

run_cmd = sprintf('./eq -f %d -s %d -b %d -g %s', ...
    sampleRateHz, testDataVectorsLength, numberOfBands, gainsdB);
system(run_cmd);

%Configure filter analizer from MDSPTK
fa = fanalyzer(sampleRateHz);

% ---- eq class tests ----
fa.freqResp(csvread('butterworth.tstdat'),'log','eq2: butterworth');
fa.freqResp(csvread('chebyshev1.tstdat'),'log','eq2: chebyshev1');
fa.freqResp(csvread('chebyshev2.tstdat'),'log','eq2: chebyshev2');

%Clean up all
system('make clean');