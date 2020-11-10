function [SLD, row_det, column_det, counter, CFAR_T, row_detection2, col_detection2] = CFAR(S)
    PFA = 1e-5;
    signal = 10.^(S./20); %%convert from dB to linear unit
    SLD = abs(signal).^2; % Square Law Detector - produces power of signal
  
    [row, column] = size(SLD);
    %NumDataPoints = row*column;
    
    %% Arrays to store detection positions
    row_det = []; %detection positions in row
    column_det = []; %detection positions in column
    
    row_detection(column).vector = []; 
    col_detection = nan(1,column);
    
    %% Parameters
    window = 12; %window size 32
    N = window*2; %Number of Reference Cells
    guard_cells = 6; %Guard cells - 3/4  2
    
    CFAR_T =  []; % intialise array for each CUT
    counter = 1;
    counter2 = 1;
    
    %% Determine CFAR for each CUT
    %% Apply algorithm along each column in spectrogram
    for c = 1:column
        counter2 = 1;
        power = SLD(1:row,c); %Single column
        
        %No false alarm outside the window and guard cell region
        region = window + guard_cells +1;
        for r = region:row-region
            
            CUT = power(r); %power of CUT
            %CA-CFAR
            lag_window = power(r-window-guard_cells:r-1-guard_cells);
            lead_window = power(r+1+guard_cells:r+window+guard_cells);
            
            %Calculate CA-CFAR Threshold
            %1. Calculate the interference statistic
            g_CA = (sum(lag_window)+ sum(lead_window))./N;
            %2. Calculate CFAR constant
            alpha_CA = N*(PFA^(-1/(N))-1);
            %CA-CFAR Threshold
            threshold = g_CA*alpha_CA;
            
            CFAR_T = [CFAR_T;threshold];
            
            if (threshold < CUT)
                row_det = [row_det; r];
                column_det = [column_det;c];
                
                row_detection(c).vector(counter2) = r;
                col_detection(c) = c;
                
                counter = counter + 1;
                counter2 = counter2 + 1;
            end 

        end
  
    end
    
    x = 1;
    
    for c = 1:column
        row_detection1(c) = round(median(row_detection(c).vector));
    end
    
    row_detection2 = row_detection1(~isnan(row_detection1))';
    col_detection2 = col_detection(~isnan(col_detection))';
    
end