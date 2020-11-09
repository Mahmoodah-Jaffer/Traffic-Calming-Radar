%% Produces Spectrogram (Doppler x Time) of Measured Data
%% parameters
c = 299e6; %speed of light in m/s
cpi = 0.10; %coherent processing interval - seconds
fc = 2590e6; %Center frequency (connect VCO Vtune to +5)-- this depends on the cantenna used and Fc chosen
maxSpeed = 30; %maximum speed to display in m/s
%% Read WAV file
wavFile = 'Audi_A1_Driving_Towards_Fast.wav';
[y,fs] = audioread(wavFile,'native'); %'native': Samples in the native format found in the file.
%% Derive parameters
N_block = fix(cpi*fs); %number of samples per pulse
FFT_size = N_block;
Overlap = round(FFT_size/24); % 4 8 12 24
%Overlap = 8;
lambda = c/fc; %wavelength in metres
%ovsFrames = N_block*ovsDop; %number of oversampled samples
%% compute a Doppler window 
Win = hamming(N_block); 

%% Compute STFT
[S, t, speed] = myspectro(y,Win,Overlap,FFT_size, fs, lambda, cpi, maxSpeed);
%[S, t, speed] = cantenna_dop_v3_yunus(wavFile);
%% Apply CFAR detection to STFT of y

[SLD, row_det, column_det, counter, CFAR_T] = CFAR(S);

%% plot the spectrogram
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
% rowX = speed(t1);
% columnX = t(t2);
%plot(rowX,columnX,'kx', 'MarkerSize',8, 'LineWidth',2);
plot(t1,t2,'kx', 'MarkerSize',8, 'LineWidth',2);
grid on;

colorbar;
