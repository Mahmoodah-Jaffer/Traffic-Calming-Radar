%% Cell Averaging CFAR Detection Algorithm
clear all;  
close all;

PFA = 1e-3; % Abdul Gaffar 1e-3 = 1/1000
             
NumDataPoints = 1000000; 
test = 1:NumDataPoints;%100 000 data points

a = normrnd(0,1,1,NumDataPoints); %Chosen by Gaussian PDF
b = normrnd(0,1,1,NumDataPoints); %Chosen by Gaussian PDF
x = (a + (1i*b))/sqrt(2); %a - I, b - Q --> NOISE SO all detections are false alarms??

SLD = abs(x).^2; % Square Law Detector
i_power = mean(SLD); %Interference power of Square Law Detector

window = 32; %window size
N = window*2; %Number of Reference Cells
guard_cells = 2; %Guard cells - 3/4

len_x = length(SLD);
CFAR_T = zeros(len_x,1); % intialise matrix for each CUT
counter = 0; %keep track of number of false alarms

%% Determine CFAR for each CUT
for CUT = 1:len_x
    %No false alarm outside the window and guard cell region
    region = window + guard_cells +1;
    if (CUT < region || CUT > (len_x - region))
        CFAR_T(CUT) = 0;
        continue;
    end
    
    %CA-CFAR
    lag_window = SLD(CUT-window-guard_cells:CUT-1-guard_cells);
    lead_window = SLD(CUT+1+guard_cells:CUT+window+guard_cells);
     
     
    %Calculate CA-CFAR Threshold
    %1. Calculate the interference statistic
    g_CA = (sum(lag_window)+ sum(lead_window))./N;
    %2. Calculate CFAR constant
    alpha_CA = N*(PFA^(-1/(N))-1);
    %CA-CFAR Threshold
    threshold = g_CA*alpha_CA;
    
    %CFAR Threshold
    %threshold = -log(PFA)*i_power; %page 631 of textbook
    
    CFAR_T(CUT) = threshold;
    
    if (CFAR_T(CUT)<SLD(CUT))
        counter = counter + 1;
    end 
 
end


%% Plot r(CUT) and threshold - 100000*10e-3 = 100 false alarms
PFA_sim = counter/NumDataPoints;
Error = abs(((PFA_sim-PFA)/PFA)*100); % Abdul Gaffar abs()

figure
plot(test(1:10000), SLD(1:10000),'color','r')
title('CUT and Threshold')

hold on

plot(test(1:10000),CFAR_T(1:10000),  'color','g')

hold off
grid on;
legend('Data', 'Threshold');




