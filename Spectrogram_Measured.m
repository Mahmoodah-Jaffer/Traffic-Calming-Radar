%% Produces Spectrogram (Doppler x Time) of Measured Data
%% parameters
c = 299.8e6; %speed of light in m/s
cpi = 0.10; %coherent processing interval - seconds
fc = 2590e6; %Center frequency (connect VCO Vtune to +5)-- this depends on the cantenna used and Fc chosen
maxSpeed = 30; %maximum speed to display in m/s
%% Read WAV file
wavFile = '20130603 01 test 1 cw fahad to door.wav';
%'native': Samples in the native format found in the file.
[y,fs] = audioread(wavFile,'native');
%% Derive parameters
N_block = fix(cpi*fs); %number of samples per pulse
FFT_size = N_block;
Overlap = 8;
lambda = c/fc; %wavelength in metres
%ovsFrames = N_block*ovsDop; %number of oversampled samples
%% compute a Doppler window 
Win = hamming(N_block); 

%% Compute STFT
[S, t, speed] = myspectro(y,Win,Overlap,FFT_size, fs, lambda, cpi, maxSpeed);

%% plot the spectrogram
figure;
imagesc(t,speed,S'); %display image with scaled colours
axis xy; 
axis tight; 
colormap(jet(256)); 
caxis(max(S(:)) + [-60 0]); % show 60 dB dynamic range
xlabel("Time(s)")
ylabel("Speed(m/s)")
colorbar;

%% Function to compute STFT
function [STFT, t, speed] = myspectro(y,Win,Overlap,FFT_size, fs, lambda, cpi, maxSpeed)
    
    y = -y(:,2);%x must be a column vector and from 2nd channel 
     len_y = length(y);           % length of column vector y
%     len_w = length(Win);         % length of hamming window function
%     diff = len_y - len_w;
%     hop = len_w - Overlap;
    
    %nframes = 1+fix(diff/hop); %fix rounds down to nearest integer value

    nframes = 1+ fix(len_y / FFT_size * Overlap) - (Overlap) + 1; %number of overlapped frames

    %% compute axes parameters for the plot
    %Have not implemented oversampling
    doppler_freq = (0:FFT_size/2-1).' / (FFT_size) * fs; %Doppler frequency
    
    speed = (doppler_freq*lambda)/2; %velocity of wave at each Doppler frequency in m/s
    t = (1:nframes)*(cpi/Overlap); %time interval in seconds

    %limit the speed axis to a reasonable range
    speed = speed(speed <= maxSpeed);
    Speed_size = length(speed);
    STFT = zeros(Speed_size,nframes); %initialise STFT matrix
    
    % STFT (via time-localized FFT)
    for count = 0:nframes-1
        %windowing -- why does it need to be double precision for the
        %window to be applied?
        yIndx = ((1:FFT_size).' + (count)*floor(FFT_size/Overlap));
        yw = double(y(yIndx)).*Win;
        
        %FFT -- only positive values therefore no fftshift
        Y = fft(yw, FFT_size);

        %update of the stft matrix
        STFT(:, 1+count) = 20*log10(abs(Y(1:Speed_size)));%convert STFT to dB 
    end
    STFT = STFT.';
end
