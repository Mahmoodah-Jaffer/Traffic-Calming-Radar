function [SLD, row_det, column_det, counter, CFAR_T] = CFAR(S)
    PFA = 1e-3;
    %signal = S;
    signal = 10.^(S./20);
    SLD = abs(signal).^2; % Square Law Detector
    %SLD = SLD1(:);
  
    [row, column] = size(SLD);
    %NumDataPoints = row*column;
    
    %% Arrays to store detection positions
    row_det = []; %detection positions in row
    column_det = []; %detection positions in column
    
    %% Parameters
    window = 12; %window size 32
    N = window*2; %Number of Reference Cells
    guard_cells = 6; %Guard cells - 3/4  2
    
    CFAR_T =  []; % intialise array for each CUT
    counter = 0;
    
    %% Determine CFAR for each CUT
    %% Apply algorithm along each column in spectrogram
    for c = 1:column
        
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
%             test = sum(lag_window)+ sum(lead_window);
%             tester = sum(test);
%             g_CA = (tester)./N;
            g_CA = (sum(lag_window)+ sum(lead_window))./N;
            %2. Calculate CFAR constant
            alpha_CA = N*(PFA^(-1/(N))-1);
            %CA-CFAR Threshold
            threshold = g_CA*alpha_CA;
            
            CFAR_T = [CFAR_T;threshold];
            
            if (threshold < CUT)
                row_det = [row_det; r];
                column_det = [column_det;c];
                counter = counter + 1;
            end 

        end
  
    end
    
%% Plot r(CUT) and threshold - 100000*10e-3 = 100 false alarms
% figure
% plot(test(1:NumDataPoints), SLD(1:NumDataPoints),'color','r')
% title('CUT and Threshold')
% 
% hold on
% 
% plot(test(1:NumDataPoints),CFAR_T(1:NumDataPoints),  'color','g')
% 
% hold off
% grid on;
% legend('Data', 'Threshold');
    
end