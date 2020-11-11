%% Cell Averaging CFAR Detection Algorithm
clear all;  
close all;

PFA = 1e-3; % Abdul Gaffar 1e-3 = 1/1000             
 
row = 1000;
column = 100;

a = normrnd(0,1,row,column); %Chosen by Gaussian PDF
b = normrnd(0,1,row,column); %Chosen by Gaussian PDF
x = (a + (1i*b))/sqrt(2); %a - I, b - Q --> NOISE SO all detections are false alarms??

SLD = abs(x).^2; % Square Law Detector
  
    %[row, column] = size(SLD);
    NumDataPoints = row*column;
    
    %% Arrays to store detection positions
    row_det = []; %detection postiions in row
    column_det = []; %detection postiions in column
    
    %% Parameters
    window = 32; %window size
    N = window*2; %Number of Reference Cells
    guard_cells = 2; %Guard cells - 3/4
    
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
    PFA_sim = counter/NumDataPoints;
    Error = abs(((PFA_sim-PFA)/PFA)*100); % Abdul Gaffar abs()
% NumDataPoints = rows*columns;
% test = 1:NumDataPoints;%100 000 data points
% 
% 
% a = normrnd(0,1,rows,columns); %Chosen by Gaussian PDF
% b = normrnd(0,1,rows,columns); %Chosen by Gaussian PDF
% x = (a + (1i*b))/sqrt(2); %a - I, b - Q --> NOISE SO all detections are false alarms??
% 
% SLD = abs(x).^2; % Square Law Detector
% %SLD = SLD1(:);
% i_power = mean(SLD); %Interference power of Square Law Detector
% 
% window = 16; %window size
% N = window*2; %Number of Reference Cells
% guard_cells = 3; %Guard cells - 3/4
% 
% len_x = NumDataPoints;
% CFAR_T = zeros(rows,columns); % intialise matrix for each CUT
% counter = 0; %keep track of number of false alarms
% 
% %% Determine CFAR for each CUT
% for CUT = 1:numel(SLD)
%     %No false alarm outside the window and guard cell region
%     region = window + guard_cells +1;
%     if (CUT < region || CUT > (len_x - region))
%         CFAR_T(CUT) = 0;
%         continue;
%     end
%     
%     %CA-CFAR
%     lag_window = SLD(CUT-window-guard_cells:CUT-1-guard_cells);
%     lead_window = SLD(CUT+1+guard_cells:CUT+window+guard_cells);
%      
%      
%     %Calculate CA-CFAR Threshold
%     %1. Calculate the interference statistic
%     g_CA = (sum(lag_window)+ sum(lead_window))./N;
%     %2. Calculate CFAR constant
%     alpha_CA = N*(PFA^(-1/(N))-1);
%     %CA-CFAR Threshold
%     threshold = g_CA*alpha_CA;
%     
%     %CFAR Threshold
%     %threshold = -log(PFA)*i_power; %page 631 of textbook
%     
%     CFAR_T(CUT) = threshold;
%     
%     if (CFAR_T(CUT)<SLD(CUT))
%         counter = counter + 1;
%     end 
%  
% end


% for r = 1:rows   
%     for c = 1:columns
%         
%         %No false alarm outside the window and guard cell region
%         region = window + guard_cells +1;
%         if (r < region || r > (rows - region) || c < region || c > (columns - region) )
%             CFAR_T(r,c) = 0;
%             continue;
%         end
%         
%         %CUT = r*c;
%         %CA-CFAR
% %           lag_window = SLD(CUT-window-guard_cells:CUT-1-guard_cells);
% %           lead_window = SLD(CUT+1+guard_cells:CUT+window+guard_cells);
%         lag_window = SLD(r-window-guard_cells:r-1-guard_cells,c-window-guard_cells:c-1-guard_cells );
%         %lag_window = SLD(r-window-guard_cells,c-window-guard_cells:r-1-guard_cells,c-1-guard_cells );
%         lead_window = SLD(r+1+guard_cells:r+window+guard_cells, c+1+guard_cells:c+window+guard_cells);
%          
%        %Calculate CA-CFAR Threshold
%        %1. Calculate the interference statistic
%        test = sum(lag_window)+ sum(lead_window);
% %        tester = sum(test);
% %        g_CA = (tester)./N;
%        
%        %2. Calculate CFAR constant
%         alpha_CA = N*(PFA^(-1/(N))-1);
%         %CA-CFAR Threshold
%         threshold = g_CA*alpha_CA;
%        
%         CFAR_T(r,c) = threshold;
%         why = SLD(r,c);
%         
%         if (CFAR_T(r,c)<SLD(r,c))
%             counter = counter + 1;
%         end           
%     end
% end
% 
%% Plot r(CUT) and threshold - 100000*10e-3 = 100 false alarms
% PFA_sim = counter/NumDataPoints;
% Error = abs(((PFA_sim-PFA)/PFA)*100); % Abdul Gaffar abs()
% 
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