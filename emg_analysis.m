%{
Author: Kalpajyoti Hazarika
Date: '01-Jul-2023'

Desciption: This file analyse the electromyography data from visual
reactive key Pressing task

%}

%%
% data extraction from the nev and nsx files

clear; close all; clc;
addpath('D:\my drive\OneDrive - Indian Institute of Science\lab works\lab works\emg_analysis\blackrock analysis codes');

openNEV
openNSx('uV');  
openNSx('uV'); 
%%
ns4 = NS4.Data{3};
ns6 = NS6.Data{2};
data_ns4= ns4(1:16,:)/4;      % 10kHz sampling
data_ns6 = ns6(1:16,:)/4;      % raw data
%%
% channel data
ch_01_uf = data_ns4(12,:);   
ch_02_uf = data_ns4(13,:);   
ch_03_uf = data_ns4(14,:);   
ch_04_uf = data_ns4(15,:);   
ch_05_uf = data_ns4(16,:);

% parameters
Fs = 10000;             % Sampling frequency (10 kHz)                  
T  = 1/Fs;              % Sampling period (0.0001s or 0.1 ms)      
L  = size(data_ns4,2);  % Length of signal
t  = (0:L-1)*T;         % in seconds or t = (0:1:(L-1))./10000; 

%%
date = char(datetime(strtok(NEV.MetaTags.DateTime,' '),'Format','yyMMdd'));
fname_br = sprintf('Analyzed_EMG_%s_%s',date,NEV.MetaTags.Filename); 

% Save file
% Identify folder where the analyzed monkeylogic file is saved, and save the analyzed blackrock file
filepath_br = uigetdir('D:\lab works\save_emg_data\','Choose folder to store file');
% cd(filepath_br)

save(fullfile(filepath_br,fname_br),"fname_br","filepath_br",'Fs','t',"NS6","NS4","NEV",'ch_01_uf','ch_02_uf','ch_03_uf','ch_04_uf', ...
'ch_05_uf','-v7.3'); 

% save(fullfile(filepath_br,fname_br),"fname_br","filepath_br",'Fs','t',"NS6","NS4","NEV",'ch_01_uf','ch_02_uf','ch_03_uf','ch_04_uf', ...
% 'ch_05_uf','ch_06_uf','ch_07_uf','ch_08_uf','ch_09_uf','ch_10_uf','ch_11_uf','ch_12_uf','ch_13_uf','ch_14_uf','ch_15_uf','ch_16_uf','-v7.3'); 

addpath(genpath('D:\lab works\save_emg_data\'))
savepath
%%
% Even Marker
br_ttl = NEV.Data.SerialDigitalIO.TimeStampSec(1:end);
%%
% key_pressing time and report time for stimuli changes based on blackrock ttl information
% 

nb_trial = 30;

diff_br_ttl = diff(br_ttl);

response_time = diff_br_ttl(1:4:end)';

percieved_time = diff_br_ttl(3:4:end)'/2;

%%
% visualation of emg signal
plot(ch_01_uf)



% for ii = 1:nb_trial
%     for jj = 1:length(br_ttl)
%     response_time(ii) = br_ttl(jj+1) - br_ttl(jj);
%     end
% end


%%
end_time = zeros(nb_trial,1);

for i = 1:nb_trial
    end_time(i) = eval(strcat('Trial',num2str(i),'.BehavioralCodes.CodeTimes(8)'));
end


%%
data= double(data_ns4);
% addpath('J:\adiLab\kalpa\satya_code')
addpath('D:\trial_15\removeLineNoise_SpectrumEstimation');
% emg_data = linenoise_remove_emg(data,[2,550]);
data = data(12:16,:);

%{
In the data set the emg signal is recorded 3.3741 second later than the
experiment began. So, I have added 33741 zeros before the signal starts.
%}

diff_time = br_ttl(2) - 10.9888;

% event marker
br_ttl = NEV.Data.SerialDigitalIO.TimeStampSec(1:end);
br_ttlNew = (br_ttl - diff_time)*10000;

% Filter the signal using bandpass filter [10 500]Hz

% Bandpass Filter Design Parameters
order = 4;              % Filter order
fs = 10000;              % Sampling frequency (Hz)
fpass = [10 500];       % Passband frequencies (Hz)
Wpass = fpass / (fs/2); % Normalized passband frequencies

% Butterworth Bandpass Filter Design
[b, a] = butter(order, Wpass, 'bandpass');

% Filter Signal

filtered_signal = zeros(5,length(data));

for ii = 1:5
    input_signal = data(ii,:);
    filtered_signal(ii,:) = filter(b, a, input_signal);
    
    % Generate time vector
    t = (0:length(input_signal)-1) / fs;
    
    % Plot input signal
    figure(ii)
    subplot(2,2,1);
    plot(t, input_signal);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Input Signal');
    
    % Plot filtered signal
    subplot(2, 2, 2);
    plot(t, filtered_signal(ii,:));
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Filtered Signal');
    
    % Adjust plot layout
    sgtitle('IIR Bandpass Filtered Signals');
    
    % frequency response
    
    % Calculate frequency response of original signal
    original_spectrum = abs(fft(input_signal));
    original_spectrum = original_spectrum(1:length(original_spectrum)/2+1);
    frequencies = (0:length(original_spectrum)-1) * fs / length(original_spectrum);
    
    % Calculate frequency response of filtered signal
    filtered_spectrum = abs(fft(filtered_signal(ii,:)));
    filtered_spectrum = filtered_spectrum(1:length(filtered_spectrum)/2+1);
    
    % Plot frequency response    
    subplot(2,2,3)
    plot(frequencies, original_spectrum, 'b', 'LineWidth', 2);
   
    subplot(2,2,4)
    plot(frequencies, filtered_spectrum, 'r', 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Frequency Response');
    grid on
 
end
    
% segment the data based on TTL
emg_data = cell(120,50000);
seg_data = cell(120,50000);

diff_sc_rt = diff(br_ttlNew);

for_idx = mean(diff_sc_rt(1:2:end));
back_idx = mean(diff_sc_rt(2:2:end));


for jj = 1:60
    % segment the original raw data
   seg_data{jj} = data(:,round(br_ttlNew(jj*2))-round(for_idx):round(br_ttlNew(jj*2))+round(back_idx));
   emg_extract = filtered_signal(:,round(br_ttlNew(jj*2))-round(for_idx):round(br_ttlNew(jj*2))+round(back_idx));
   emg_data{jj} = removeLineNoise_SpectrumEstimation(emg_extract,10000,'LF = 50, NH = 5, HW = 4, M = 2048');
end


% paramters to plot
Fs= 10000;
ts = 1/Fs;
L = length(emg_data{1}(1,:));
t = (0:(L-1))*ts;

for i = 1:5
    subplot(5,1,i)
    plot(t,emg_data{1}(i,:),'r');
    hold on;
    plot(t,seg_data{1}(i,:),'b')
    xlabel('Sec')
    ylabel('\muV')
    axis 'tight'
    legend('filtered','Original')
end

for ii = 1:5

    original_spectrum = abs(fft(seg_data{1}(ii,:)));
    original_spectrum = original_spectrum(1:length(original_spectrum)/2+1);
    frequencies = (0:length(original_spectrum)-1) * fs / length(original_spectrum);
    
    % Calculate frequency response of filtered signal
    filtered_spectrum = abs(fft(emg_data{1}(ii,:)));
    filtered_spectrum = filtered_spectrum(1:length(filtered_spectrum)/2+1);
    
    % Plot frequency response    
    subplot(5,1,ii)
    plot(frequencies, original_spectrum, 'b', 'LineWidth', 2);
    hold on;
    
    plot(frequencies, filtered_spectrum, 'r', 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Frequency Response');
    legend('Original','Filtered')    
    grid on

end

%%
% Load EMG data (replace 'emg_data.mat' with your actual data file)

% Set threshold as mean ± 3 standard deviations

threshold_p = mean(emg_data{1}(1,:)) + 3 * std(emg_data{1}(1,:));
threshold_n = mean(emg_data{1}(1,:)) - 3 * std(emg_data{1}(1,:));

% Find spike indices that cross the threshold
spikeIdx_p = find(emg_data{1}(1,:) > threshold_p);
spikeIdx_n = find(emg_data{1}(1,:) < threshold_n);

%%
% Sample data (replace with your actual EMG signal and time vector)
% emg_signal = [your EMG signal];
% time_vector = [your time vector];

% Step 1: Calculate mean and standard deviation
mean_emg = mean(emg_data{1}(1,:));
std_emg = std(emg_data{1}(1,:));

% Step 2: Define time bins for the PSTH (you can adjust the bin size)
bin_size = 0.1; % in seconds
time_bins = min(t):bin_size:max(t);
time_vector = t;
emg_signal = emg_data{1}(1,:);
% Step 3: Divide the EMG signal into segments corresponding to time bins
emg_segments = cell(length(time_bins) - 1, 1);
for i = 1:length(time_bins)-1
    idx = time_vector >= time_bins(i) & time_vector < time_bins(i+1);
    emg_segments{i} = emg_signal(idx);
end

% Step 4: Calculate mean and standard deviation for each time bin
mean_psth = zeros(size(time_bins) - 1);
sd_psth = zeros(size(time_bins) - 1);
for i = 1:length(emg_segments)
    mean_psth(i) = mean(emg_segments{i});
    sd_psth(i) = std(emg_segments{i});
end

% Step 5: Plot the PSTH with mean +/- 2 SD
figure;
plot(time_bins(1:end-1), mean_psth, 'b');
hold on;
plot(time_bins(1:end-1), mean_psth + 3 * sd_psth, 'r--');
plot(time_bins(1:end-1), mean_psth - 3 * sd_psth, 'r--');
xlabel('Time (s)');
ylabel('EMG Amplitude');
title('Peri-Stimulus Time Histogram (PSTH)');
legend('Mean', 'Mean + 2 SD', 'Mean - 2 SD');
grid on;


%%

% Create spike times vector based on sampling frequency (replace 'fs' with your sampling frequency)
spike_times = spike_indices / fs;

% Create spike labels (optional)
spike_labels = ones(size(spike_times));

% Plot raster plot
figure;
scatter(spike_times, spike_labels, 'k.');
xlabel('Time (s)');
ylabel('Spike');
title('Raster Plot of EMG Spikes');







%%
% extracting the event marker

files = dir(filepath_br);
for k = 1:numel(files)
    if contains(files(k).name,"Analyzed_ML",'IgnoreCase',true) 
       load(files(k).name)
    end
end

%%
A = cell(size(channels,1),1); 
% M = zeros(size(cond_rt,1),50000); 
M = 120;
A(:,1) = {M};

% unfiltered data, all trials, aligned on target onset & movement onset 
all_trials_uf_gocue = A;    all_trials_uf_mov = A;   

% filtered data, all trials, aligned on target onset & movement onset 
gocue_wise_all_ch_filt = A; mov_wise_all_ch_filt = A;    

%%


