clear all;
close all;

%Check os
if(ispc)
    error('Windows based OS not supported.');
end;

%Build configuration file
cfgFileName = 'orfanidis_eq.tstdat';
cfgFileMessage = '/* parameters for test bench: sample rate in Hz, number of bands, testing data vectors length */';
cfgSampleRateHz = 48000;
cfgNumberOfBands = 30;
cfgTestDataVectorsLength = cfgSampleRateHz/4;

cfgFId = fopen(cfgFileName,'w');
fprintf(cfgFId, '%s', cfgFileMessage);
fprintf(cfgFId, '%d, %d, %d\n', cfgSampleRateHz, cfgNumberOfBands, cfgTestDataVectorsLength);
fclose(cfgFId);

%Build and run test file
system('make');
system('./eq');

%Read data from files

% ------------ Butterworth ------------
%Data file template 'butterworth_x.tstdat'

%Read magnitudes
butterworth = zeros(cfgNumberOfBands,cfgTestDataVectorsLength + 1); %Hack: number of samples x + 1
for(n = 0:(cfgNumberOfBands - 1))
    fileName = strcat('butterworth_', num2str(n),'.tstdat');
    butterworth(n+1,:) = csvread(fileName);
end;

% ------------ Chebyshev type 1 ------------
%Data file template 'chebyshev1_x.tstdat'

%Read magnitudes
chebyshev1 = zeros(cfgNumberOfBands,cfgTestDataVectorsLength + 1); 
for(n = 0:(cfgNumberOfBands - 1))
    fileName = strcat('chebyshev1_', num2str(n),'.tstdat');
    chebyshev1(n+1,:) = csvread(fileName);
end;

% ------------ Chebyshev type 2 ------------
%Data file template 'chebyshev2_x.tstdat'

%Read magnitudes
chebyshev2 = zeros(cfgNumberOfBands,cfgTestDataVectorsLength + 1);
for(n = 0:(cfgNumberOfBands - 1))
    fileName = strcat('chebyshev2_', num2str(n),'.tstdat');
    chebyshev2(n+1,:) = csvread(fileName);
end;


%Plot all

%Calculate and plot it
N=length(cfgTestDataVectorsLength + 1); %Hack again
a=zeros(1,N);
a(1)=1;

% ------------ Butterworth ------------

mxsz = size(butterworth);
for n=1:mxsz(1)
    [H1(n,:), F]=freqz(butterworth(n,:), a, 65536, cfgSampleRateHz);
end;

figure1 = figure;
axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on', 'XMinorGrid','on');
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

plot(F,10*log(abs(H1)));

title('Butterworth');
xlabel('f->');
ylabel('Mag dB');

% ------------ Chebyshev 1 ------------

mxsz = size(chebyshev1);
for n=1:mxsz(1)
    [H1(n,:), F]=freqz(chebyshev1(n,:), a, 65536, cfgSampleRateHz);
end;

figure2 = figure;
axes1 = axes('Parent',figure2,'XScale','log','XMinorTick','on', 'XMinorGrid','on');
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

plot(F,10*log(abs(H1)));

title('Chebyshev Type 1');
xlabel('f->');
ylabel('Mag dB');

% ------------ Chebyshev 2 ------------

mxsz = size(chebyshev2);
for n=1:mxsz(1)
    [H1(n,:), F]=freqz(chebyshev2(n,:), a, 65536, cfgSampleRateHz);
end;

figure3 = figure;
axes1 = axes('Parent',figure3,'XScale','log','XMinorTick','on','XMinorGrid','on');
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

plot(F,10*log(abs(H1)));

title('Chebyshev Type 2');
xlabel('f->');
ylabel('Mag dB');


%Clean up all
system('make clean');

