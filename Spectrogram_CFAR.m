%% Produces Spectrogram (Doppler x Time) of Measured Data
clear all;
close all;

%% parameters
c = 299e6; %speed of light in m/s
cpi = 0.10; %coherent processing interval - seconds
fc = 2590e6; %Center frequency (connect VCO Vtune to +5)-- this depends on the cantenna used and Fc chosen
maxSpeed = 30; %maximum speed to display in m/s
%% Read WAV file
wavFile = 'Audi_A1_Driving_Away_45KPH.wav'; 
        % 'Audi_A1_Driving_Away_fast.wav'; % 'Audi_A1_Driving_Towards_Fast.wav';
[y,fs] = audioread(wavFile,'native'); %'native': Samples in the native format found in the file.
%% Derive parameters
N_block = fix(cpi*fs); %number of samples per pulse
FFT_size = N_block;
Overlap = round(FFT_size/4);
%Overlap = 8;
lambda = c/fc; %wavelength in metres
%ovsFrames = N_block*ovsDop; %number of oversampled samples
%% compute a Doppler window 
Win = hamming(N_block); 

%% Compute STFT
[S, t, speed] = myspectro(y,Win,Overlap,FFT_size, fs, lambda, cpi, maxSpeed);
%% Apply CFAR detection to STFT of y
[SLD, row_det, column_det, counter, CFAR_T, row_detection, col_detection] = CFAR(S);

%% plot the spectrogram
figure;
imagesc(t,speed,S); %display image with scaled colours
axis xy; 
axis tight; 
colormap(jet(256)); 
caxis(max(S(:)) + [-60 0]); % show 60 dB dynamic range
xlabel("Time(s)")
ylabel("Speed(m/s)")
colorbar;

%% plot the spectrogram 1
figure;
imagesc(t,speed,S); %display image with scaled colours
axis xy; 
axis tight; 
colormap(jet(256)); 
caxis(max(S(:)) + [-60 0]); % show 60 dB dynamic range
xlabel("Time(s)")
ylabel("Speed(m/s)")


hold on;
t1 = row_det;
t2 = column_det;
rowX = speed(t1);
columnX = t(t2);
plot(columnX, rowX,'kx', 'MarkerSize',8, 'LineWidth',2);
%plot(t1,t2,'kx', 'MarkerSize',8, 'LineWidth',2);
grid on;

colorbar;

%% plot the spectrogram 2
figure;
imagesc(t,speed,S); %display image with scaled colours
axis xy; 
axis tight; 
colormap(jet(256)); 
caxis(max(S(:)) + [-60 0]); % show 60 dB dynamic range
xlabel("Time(s)")
ylabel("Speed(m/s)")
title("Median velocity");

hold on;
t1 = row_detection;
t2 = col_detection;
rowX = speed(t1);
columnX = t(t2);
plot(columnX, rowX,'kx', 'MarkerSize',8, 'LineWidth',2);
%plot(t1,t2,'kx', 'MarkerSize',8, 'LineWidth',2);
grid on;

colorbar;


% %% Plot speed
% figure;
% t1 = row_detection;
% t2 = col_detection;
% rowX = speed(t1);
% columnX = t(t2);
% 
% % fit low order polynomial
% order = 3;
% p = polyfit(columnX,rowX,order);
% f = polyval(p,columnX);
% 
% plot(columnX, rowX,'kx', columnX, f, 'r-');
% grid on;
% axis([t(1)  t(end)  0 30]);
% xlabel("Time(s)")
% ylabel("Speed(m/s)")
% legend('raw speed estimate', 'smooth speed estimate');

%% Plot speed
figure;
t1 = row_detection;
t2 = col_detection;
rowX = speed(t1);
columnX = t(t2);

% fit low order polynomial
order = 3;
rowX_T = rowX.';
p = polyfit(columnX,rowX_T,order);
f = polyval(p,columnX);

plot(columnX, rowX,'kx', columnX, f, 'r-');
grid on;
axis([t(1)  t(end)  0 30]);
xlabel("Time(s)")
ylabel("Speed(m/s)")
legend('raw speed estimate', 'smooth speed estimate');
