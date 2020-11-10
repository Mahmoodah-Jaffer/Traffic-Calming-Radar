%% Function to compute STFT
function [STFT, t, speed] = myspectro(y,Win,Overlap,FFT_size, fs, lambda, cpi, maxSpeed)
    
     y = -y(:,2);%x must be a column vector and from 2nd channel 
     len_y = length(y);           % length of column vector y

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
    %STFT = STFT.';
end