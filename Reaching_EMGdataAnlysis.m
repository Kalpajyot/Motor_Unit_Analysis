%
% Author : Kalpajyoti Hazarika
% Date: '08-Jul-2023'
% This code analyse emg data
% Preprocessing and motor unit extraction
%%
clear; close all; clc;
addpath('D:\my drive\OneDrive - Indian Institute of Science\lab works\lab works\Reaching_LibetTask\emgData');

openNEV
openNSx('uV');  
openNSx('uV'); 

%%
ns4 = NS4.Data;
ns6 = NS6.Data;
data_ns4= ns4(1:16,:)/4;      % 10kHz sampling
data_ns6 = ns6(1:16,:)/4;      % raw data
%%
% segregate the data into in and out resposne field trial from the
% behavioural data
addpath('J:\onedrive\OneDrive - Indian Institute of Science\lab works\lab works\Reaching_LibetTask\behavioral data')
load('230705_sadhvika_reach_task_fix.mat')

nb_trial = TrialRecord.CurrentTrialNumber;
trial_cond = zeros(nb_trial,1);

for i = 1:nb_trial
    trial_cond(i) = eval(strcat('Trial',num2str(i),'.Condition'));
end

%%
% Even Marker from blackrock 
% Scene :
% Dial start and ends

br_ttl = NEV.Data.SerialDigitalIO.TimeStampSec(1:end);
% duration of scene in each trial
diff_ttl = diff(br_ttl);
time_start_end = diff_ttl(1:2:end);

% monkey logic event marker timing

dial_start =zeros(nb_trial,1);
movement_end = zeros(nb_trial,1);
go_cueTime = cell(1,nb_trial);


for i = 1:nb_trial
    dial_start(i) = round(eval(strcat('Trial',num2str(i),'.BehavioralCodes.CodeTimes(5)')));
    movement_end(i) = round(eval(strcat('Trial',num2str(i),'.BehavioralCodes.CodeTimes(6)')));    
    go_cueTime{i} = (dial_start(i):movement_end(i));
end

% Scene duration in each trial 
scene_duration  = (movement_end - dial_start)/1000;

% compare black rock event marker with the moknkey logic event marker
diff_ml_br = time_start_end(1:20)'- scene_duration(1:20);


%{
It is important to check the duration of each trial. 
In my first analysis of proactive task, I found that after 11 trial( incorrect
trial), there is a mismatch in the duration with the actual scene duration.
So, I will proceed with the 10 trial
%}

%%
% channel data
ch_01_uf = data_ns4(12,:);   
ch_02_uf = data_ns4(13,:);   
ch_03_uf = data_ns4(14,:);   
ch_04_uf = data_ns4(15,:);   
ch_05_uf = data_ns4(16,:);

emg_data = [ch_01_uf;ch_02_uf;ch_03_uf;ch_04_uf;ch_05_uf];


% parameters
Fs = 10000;             % Sampling frequency (10 kHz)                  
T  = 1/Fs;              % Sampling period (0.0001s or 0.1 ms)      
L  = size(data_ns4,2);  % Length of signal
t  = (0:L-1)*T;         % in seconds or t = (0:1:(L-1))./10000; 

% extract first 10 trial data from each channel

stimulus_time = [br_ttl(1:2:19)];

data1 = cell(10,1);
% data is taken 1.5 sec before the onset of the stimulus and the 6 second
% after the stimulus onset


for i = 1:10
    data1{i} = emg_data(:,round(stimulus_time(i)*10000)-19999:round(stimulus_time(i)*10000)+60000);
end

% In this section last trial is removed as there are not any sample points
% after the movement. If all the trials have enough data then keep only 
% st1 = [br_ttl(22:2:end)]

st1 = [br_ttl(22:2:end)];

data12 = cell(length(st1),1);

for j = 1:length(st1)
    data12{j} = emg_data(:,round(st1(j)*10000)-19999:round(st1(j)*10000)+60000);
end

%  in this data the 11th trial(as it is incomplete trial) and last trial are removed.
final_data = [data1;data12];

time = (-2:T:6-T);


for i = 1:length(final_data)
    subplot(5,6,i);
    % channel 1
    plot(time,final_data{i}(1,:));
    title(sprintf('Trial %1.0f',i));
    ylim([-500,500])
    
end

%%
% Trial no to remove ( Bad trial)
% Reason :
% Urge time is more than the reaction time

% trial_toRemove = [2,11,21,27];
% trial_cond(trial_toRemove) = [];

%%
% Frequency spectrum

% Filter the signal using bandpass filter [10 500]Hz

% Bandpass Filter Design Parameters
order = 4;              % Filter order
fs = 10000;              % Sampling frequency (Hz)
fpass = [10 550];       % Passband frequencies (Hz)
Wpass = fpass / (fs/2); % Normalized passband frequencies

% Butterworth Bandpass Filter Design
[b, a] = butter(order, Wpass, 'bandpass');

% Filter Signal

filtered_signal = cell(length(final_data),1);

for jj =1:length(final_data)
    for ii = 1:5

        input_signal = final_data{jj}(ii,:);
        filtered_signal{jj}(ii,:) = filter(b, a, input_signal);
        
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
        plot(t, filtered_signal{jj}(ii,:));
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Filtered Signal');

        % Adjust plot layout
        sgtitle('Signals Time and Frequency Domain');
        
        % frequency response        
        % Calculate frequency response of original signal
        original_spectrum = abs(fft(input_signal));
        original_spectrum = original_spectrum(1:length(original_spectrum)/2+1);
        frequencies = (0:length(original_spectrum)-1) * fs / length(original_spectrum);

        % Calculate frequency response of filtered signal
        filtered_spectrum = abs(fft(filtered_signal{jj}(ii,:)));
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
     
end

   %%
% Removing the line noise
% filterd_data
Noise_filtData = cell(length(filtered_signal),1);

for ii = 1:length(filtered_signal)
    Noise_filtData{ii} = removeLineNoise_SpectrumEstimation(filtered_signal{ii},10000,'LF = 50, NH = 12, HW = 16');
end

%%
% Create a sample cell array
cellArray = Noise_filtData;


% Convert cell data to a three-dimensional matrix
numCells = numel(cellArray);
[row, col] = size(cellArray{1});

matrix = zeros(row, col, numCells);
for i = 1:numCells
    matrix(:, :, i) = cellArray{i};
end

% Display the matrix size
disp('Matrix size:');
disp(size(matrix));



%%

% Visualization of the signal after applying line noise removal method
chose_trial = 29;

for i = 1:5
    noiseFilt_spectrum = abs(fft(final_data{chose_trial}(i,:)));
    noiseFilt_spectrum = noiseFilt_spectrum(1:length(noiseFilt_spectrum)/2+1);
    frequencies = (0:length(noiseFilt_spectrum)-1) * fs / length(noiseFilt_spectrum);
    
    % Calculate frequency response of filtered signal
    filtered_spectrum = abs(fft(Noise_filtData{chose_trial}(i,:)));
    filtered_spectrum = filtered_spectrum(1:length(filtered_spectrum)/2+1);
    t = (0:length(final_data{chose_trial})-1) / fs;
    
    % Plot frequency response
    figure(i)
    subplot(221)
    plot(t,final_data{chose_trial},'k');
    xlabel('Time(sec)');
    ylabel('\muV');
    title('Raw')

    subplot(222)
    plot(t,Noise_filtData{chose_trial},'r');
    xlabel('Time(sec)');
    ylabel('\muV');
    title('Filtered')

    subplot(223)
    plot(frequencies, noiseFilt_spectrum, 'b', 'LineWidth', 1.2);
    title('Originnal')
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    
    subplot(224)
    plot(frequencies, filtered_spectrum, 'r', 'LineWidth', 1.2);
    title('Line Noise Removed');
    grid on
end



%%
% trial_toRemove = [2,11,21,27,30];
trial_toRemove = [11,22];
% trial_cond(trial_toRemove) = [];
temp_trialCond = trial_cond;

temp_trialCond(trial_toRemove) = [];

in_rf = find(temp_trialCond == 1);
out_rf = find(temp_trialCond ==2);

%%

in_spikeTime = cell(length(in_rf),1);
out_spikeTime = cell(length(out_rf),1);
% trial to remove (bad trials)
% trial_toDelete = [2,20,26];
trial_toDelete = 21;

Noise_filtData(trial_toDelete) =[]; 


%%
ch_1 = 1;
ch_2 = 2;

% For in response field spike time
data_inRF = cell(length(in_rf),1);
data_outRF = cell(length(out_rf),1);

for jj = 1:length(in_rf)   
    
    trial_signal1 = Noise_filtData{in_rf(jj)}(ch_1,:);
    trial_signal2 = Noise_filtData{in_rf(jj)}(ch_2,:);
    SD_data = trial_signal1 - trial_signal2;
    % SD_12 = Noise_filtData{in_rf(jj)}(ch_1,1:14999);
    SD_12 = SD_data(1,1:20000);

    threshold_p = mean(SD_12) + 3 * std(SD_12);
    threshold_n = mean(SD_12) - 3 * std(SD_12);

    % peak location and amplitude 

    [pks_p,locs_p] = findpeaks(SD_data,'MinPeakHeight',threshold_p,'MinPeakDistance',1);  
    [pks_n,locs_n] = findpeaks(-SD_data,'MinPeakHeight',-threshold_n,'MinPeakDistance',1);       
    idx = sort([locs_p,locs_n]);   
    in_spikeTime{jj} = t(idx);
    idx = sort([locs_p,locs_n]);

    temp = SD_data;
    temp(~ismember(1:numel(temp), idx)) = 0;    
    data_inRF{jj} = temp;

end


%  for out response field spike time
for jj = 1:length(out_rf)

    trial_signal1 = Noise_filtData{out_rf(jj)}(ch_1,:);
    trial_signal2 = Noise_filtData{out_rf(jj)}(ch_2,:);
    SD_data = trial_signal1 - trial_signal2;
    % SD_12 = Noise_filtData{in_rf(jj)}(ch_1,1:14999);
    SD_12 = SD_data(1,1:20000);

    threshold_p = mean(SD_12) + 3 * std(SD_12);
    threshold_n = mean(SD_12) - 3 * std(SD_12);
    % peak location and amplitude 

    [pks_p,locs_p] = findpeaks(SD_data,'MinPeakHeight',threshold_p,'MinPeakDistance',1);  
    [pks_n,locs_n] = findpeaks(-SD_data,'MinPeakHeight',-threshold_n,'MinPeakDistance',1);       
    idx = sort([locs_p,locs_n]);   
    out_spikeTime{jj} = t(idx); 

end


% psth plot
% in response
[H1,timeVals1] = psthplot(in_spikeTime,1,[0,7],1);
% out rf
[H2,timeVals2] = psthplot(out_spikeTime,1,[0,7],1);

% smoothing 
sigma = 2; % standard deviation of Gaussian filter
win_size = ceil(2*sigma)*2+2; % window size of filter

% win_size = 50; it means 5 ms
% Apply Gaussian filter
gauss_filter = gausswin(win_size, sigma); % create Gaussian filter
smoothed_H1 = filtfilt(gauss_filter, 1, H1)/1000; % apply filter
smoothed_H2 = filtfilt(gauss_filter, 1, H2)/1000; % apply filter


figure(2)
% subplot(211)
plot(timeVals1,smoothed_H1,'r');
xlabel('time (s)');
ylabel('spikes/s');
xlim([0,7])
% title("In resposne field")

% subplot(212)
hold on;
plot(timeVals2,smoothed_H2,'k');
% xline(2.5,'--b',{'Stimulus Onset'},'linewidth',2,'fontweight','bold')
xlim([0,7])
% title("Out- repsonse field")
xlabel('time (s)');
ylabel('spikes/s');

xline(2,'--b',{'Stimulus Onset'},'linewidth',2,'fontweight','bold')
legend('In response','Out response','Stimulus onset')


%%

% subplot(2,1,1)
% rasterplot(in_spikeTime);
% title('Raster Plot')
% xlim([1,7])
% xlabel('Sec')
% ylabel('Trial')
% title('In response field')
% 
% subplot(2,1,2)
% rasterplot(out_spikeTime);
% title('Out response field')
% xlim([1,6])
% xlabel('Sec')
% ylabel('Trial')
%%

% extract the urge time
urge_time = sadhvikaProactiveanalysisS2.CorrectUrgeTime;
urge_time = urge_time(~isnan(urge_time));
reaction_time = sadhvikaProactiveanalysisS2.correctReactionTime;
reaction_time = reaction_time(~isnan(reaction_time));
diff_rt_ut = reaction_time - urge_time;


% extract data from urge time( urge time is taken as cente)
% urge_time(26) = [];
% reaction_time(26) = [];
data_urge = cell(length(urge_time),1);

% time values for plotting
t_dataUrge = (1:30000);
 
% allign with the urge time

for jj = 1:length(data_urge)
    data_urge{jj} = Noise_filtData{jj}(:,round(20000+urge_time(jj)*10000-14999):round(20000+urge_time(jj)*10000+15000));
    subplot(6,5,jj)   
    plot(t_dataUrge,data_urge{jj}');
    xline(15000,'--r','linewidth',1.2)
    xline(15000+diff_rt_ut(jj)*10000,'--k','linewidth',1.2')
    ylabel('\muV');
    xlabel('time(second)');  
end

%%

% BIPOLAR SIGNAL generation

% only parts of the signal which crosses the  mean +/3 sd of the signal can
% be taken for analysis

% this is three dimensional matrix ( channel by time samples by trial)
all_data = matrix;
% bipolar signal
all_data = diff(all_data,1);

% remove the bad trials from the all_data matrix
all_data(:,:,trial_toDelete) =[];

data_extract = cell(size(all_data,3),1);

chan_name = ["SD:1-2","SD:2-3","SD:3-4","SD:4-5"];


% in response field
data_in = cell(length(in_rf),1);

[spike_3,locs_3,t_3,data_3] = MeanSd_signal(all_data,in_rf,urge_time,3,fs);

[spike_5,locs_5,t_5,data_5] = MeanSd_signal(all_data,in_rf,urge_time,5,fs);

[spike_7,locs_7,t_7,data_7] = MeanSd_signal(all_data,in_rf,urge_time,7,fs);


% out response field 
data_out = cell(length(out_rf),1);

[spikeOut_3,locsOut_3,tOut_3,dataOut_3] = MeanSd_signal(all_data,out_rf,urge_time,3,fs);

[spikeOut_5,locsOut_5,tOut_5,dataOut_5] = MeanSd_signal(all_data,out_rf,urge_time,5,fs);

[spikeOut_7,locsOut_7,tOut_7,dataOut_7] = MeanSd_signal(all_data,out_rf,urge_time,7,fs);

%%
% small, medium and big small motor unit
% choose channels
ch_nb = 4;
trial_nb = 6;
% subplot(3,1,1)

for i = 1:4
    % ch_nb = 4;
    % locs_35 = intersect(locs_3{14}{i},locs_5{14}{i});
    locs_35 = setdiff(locs_3{trial_nb}{i},locs_5{trial_nb}{i});
    locs_57 = setdiff(locs_5{trial_nb}{i},locs_7{trial_nb}{i});
    plot(t_3,data_3{trial_nb}(i,:));
    hold on;
    % plot(t_3(locs_3{trial_nb}{i}),data_3{trial_nb}(i,locs_3{trial_nb}{i}),'^k','MarkerSize',3,'MarkerFaceColor','g');
    % hold on;
    plot(t_3(locs_35),data_3{trial_nb}(i,locs_35),'^','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r');
    hold on;
    plot(t_3(locs_57),data_3{trial_nb}(i,locs_57),'^','MarkerSize',3,'MarkerFaceColor','g','MarkerEdgeColor','g');
    hold on;
    plot(t_7(locs_7{trial_nb}{i}),data_3{trial_nb}(i,locs_7{trial_nb}{i}),'^','MarkerSize',3,'MarkerFaceColor','b','MarkerEdgeColor','b');
    % title('3 SD');
    xline(0,'--k',{'Urge Time'},'linewidth',1.2);
    xlabel('Time(second)');
    xlim([-1,1])
    ylabel('\muV');
     
end

legend('channle data','between 3 and 5(Small)','between 5 and 7(Medium) ','7 and above(Large)','location','southeast','Fontsize',12,'fontweight','bold');

%% 
% create a 3 dimensional matrix 
%  spike times of in response field
numCells_in = numel(spike_3);
[row, col] = size(data_3{1});
sp_3_in = zeros(row, col, numCells_in);
sp_5_in = zeros(row, col, numCells_in);
sp_7_in =zeros(row, col, numCells_in);

for ii = 1:numCells_in
    sp_3_in(:, :, ii) = spike_3{ii};
    sp_5_in(:,:,ii) = spike_5{ii};
    sp_7_in(:,:,ii) = spike_7{ii};

end

% spike times of out response field
numCells_out = numel(spikeOut_3);
[row, col] = size(dataOut_3{1});
sp_3_out = zeros(row, col, numCells_out);
sp_5_out = zeros(row, col, numCells_out);
sp_7_out = zeros(row, col, numCells_out);

for ii = 1:numCells_out
    sp_3_out(:, :, ii) = spikeOut_3{ii};
    sp_5_out(:,:,ii) = spikeOut_5{ii};
    sp_7_out(:,:,ii) = spikeOut_7{ii};

end

% average across the trials (in response field spike times)
H_3_in= mean(sp_3_in,3);
H_5_in = mean(sp_5_in,3);
H_7_in = mean(sp_7_in,3);

% average across the trials (out response field spike times)
H_3_out= mean(sp_3_out,3);
H_5_out = mean(sp_5_out,3);
H_7_out = mean(sp_7_out,3);


% Plotting the psth
% smoothing 
sigma = 2; % standard deviation of Gaussian filter
% win_size = ceil(2*sigma)*2+2; % window size of filter

% it means 5 ms
win_size = 100;

% channel no
chn_id = 1;

% Apply Gaussian filter
gauss_filter = gausswin(win_size, sigma); % create Gaussian filter
psth_35_in = filtfilt(gauss_filter, 1, (H_3_in(chn_id,:)-H_5_in(chn_id,:)))/1000; % apply filter
psth_35_out = filtfilt(gauss_filter, 1,(H_3_out(chn_id,:)-H_5_out(chn_id,:)))/1000; % apply filter

psth_57_in = filtfilt(gauss_filter, 1, (H_5_in(chn_id,:)-H_7_in(chn_id,:)))/1000; % apply filter
psth_57_out = filtfilt(gauss_filter, 1,(H_5_out(chn_id,:)-H_7_out(chn_id,:)))/1000; % apply filter


psth_7_in = filtfilt(gauss_filter, 1, H_7_in(chn_id,:))/1000; % apply filter
psth_7_out = filtfilt(gauss_filter, 1,H_7_out(chn_id,:))/1000; % apply filter


clf;
ts = 1/fs;
% time for psth plot
t_psth = (-1:ts:1.5-ts);
subplot(3,1,1)
plot(t_psth,psth_35_in,'-r');
hold on;
plot(t_psth,psth_35_out,'-k')
xline(0,'--k',{'Urge time'},'linewidth',1.2,'fontsize',10,'fontweight','bold')
xlabel('Second');
ylabel('spikes/second');
title('3 minus 5 (small)')
legend('In resposne field','Out resposne field','location','best','fontsize',10,'fontweight','bold')
subplot(3,1,2)
plot(t_psth,psth_57_in,'-r');
hold on;
plot(t_psth,psth_57_out,'-k');
xline(0,'--k',{'Urge time'},'linewidth',1.2,'fontsize',10,'fontweight','bold')
xlabel('Second');
ylabel('spikes/second')
title('5 minus 7 (medium)')
legend('In resposne field','Out resposne field','location','best','fontsize',10,'fontweight','bold')
subplot(3,1,3)
plot(t_psth,psth_7_in,'-r');
hold on;
plot(t_psth,psth_7_out,'-k');
xline(0,'--k',{'Urge time'},'linewidth',1.2,'fontsize',10,'fontweight','bold')
xlabel('Second');
ylabel('spikes/second');
title('7 and above (large)')

legend('In resposne field','Out resposne field','location','best','fontsize',10,'fontweight','bold')


%%

% subplot(3,1,2)
% for i = 1:4
%     % ch_nb = 4;
%     locs_35 = setdiff(locs_3{trial_nb}{i},locs_5{trial_nb}{i});
%     plot(t_3,data_3{trial_nb}(i,:));
%     hold on;
%     plot(t_3(locs_35),data_3{trial_nb}(i,locs_35),'ob','MarkerSize',3);
%     hold on;
%     title('3 SD minus 5 SD');
%     xline(0,'--r',{'Urge Time'},'linewidth',1.2);
%     xlabel('Time(second)')
% 
% end
% 
% subplot(3,1,3)
% for i = 1:4
%     % ch_nb = 4;
%     locs_57 = setdiff(locs_5{trial_nb}{i},locs_7{trial_nb}{i});
%     plot(t_3,data_3{trial_nb}(i,:));
%     hold on;
%     plot(t_3(locs_57),data_3{trial_nb}(i,locs_57),'ob','MarkerSize',3);
%     hold on;
%     title('5 SD minus 7 SD');
%     xline(0,'--r',{'Urge Time'},'linewidth',1.2);
%     xlabel('Time(second)')
% 
% end
%%
% Create a sample cell array

% data_in35 = cell(length(in_rf),1);
% 
% for i = 1:length(in_rf)
%     data_in35{i} = data_in3{i}(ch);
% end
% 
% 
% % Convert cell data to a three-dimensional matrix
% numCells = numel(data_in3);
% [row, col] = size(data_in3{1});
% 
% matrix_3 = zeros(row, col, numCells);
% for i = 1:numCells
%     matrix_3(:, :, i) = data_in3{i};
% end
% 
% matrix_5 = zeros(row, col, numCells);
% for i = 1:numCells
%     matrix_5(:, :, i) = data_in5{i};
% end
% 
% 
% matrix_7 = zeros(row, col, numCells);
% for i = 1:numCells
%     matrix_7(:, :, i) = data_in7{i};
% end
% 
% % Display the matrix size
% disp('Matrix size:');
% disp(size(matrix_5));

% in and out response field data extraction
% for ii = 1:length(cond_trial)
    % clf;
%     figure(ii)
%     for jj = 1:4
% 
%         befor_stim = all_data(jj,1:15000,rf_type(ii));
% 
%         threshold_p = mean(befor_stim) + 7 * std(befor_stim);
%         threshold_n = mean(befor_stim) - 7 * std(befor_stim);
% 
%         % peak location and amplitude
%         data_analysis = all_data(jj,round(20000+urge_time(ii)*10000-9999):round(20000+urge_time(ii)*10000+15000),rf_type(ii));
%         data_temp = data_analysis;
%         ts = 1/fs;
%         t_data = (-1:ts:1.5-ts);
%         % 
%         % [pks_p,locs_p] = findpeaks(data_temp,'MinPeakHeight',threshold_p,'MinPeakDistance',100);  
%         % [pks_n,locs_n] = findpeaks(-data_temp,'MinPeakHeight',-threshold_n,'MinPeakDistance',100);
%         % 
%         % idx = sort([locs_p,locs_n]);
%         pos_peak = find(data_temp>threshold_p);
%         neg_peak = find(data_temp<threshold_n);
%         data_temp([pos_peak,neg_peak]) = 1;
%         data_temp(~ismember(1:numel(data_temp),[pos_peak,neg_peak])) = 0;
% 
%         data_in{ii}(jj,:) = data_temp;
%         valuesToMark = sort([pos_peak,neg_peak]);
%         plot(t_data,data_analysis);
%         hold on;
%         plot(t_data(valuesToMark), data_analysis(valuesToMark), 'ro', 'MarkerSize', 1.2, 'MarkerFaceColor', 'r');
%         xline(t_data(10000),'--b',{'Urge Time'},'linewidth',1.2)
%         % xline(t_data(10000 + round(diff_rt_ut(rf_type(ii))*10000)),'--r',{'Reaction Time'},'linewidth',1.2)
%         title(sprintf('Trial %1.0f',rf_type(ii)));
%         % % Mark specific values
%         % valuesToMark = sort([pos_peak,neg_peak])./10000;
%         % for i = 1:length(valuesToMark)
%         %     value = valuesToMark(i);
%         %     % index = find(t_data == value);
%         %     if ~isempty(index)
%         %         plot(t_data(valuesToMark), y(index), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
%         %         % text(x(index), y(index), num2str(value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%         %     end
%         % end
% 
%         % idx = sort([locs_p,locs_n]);
%         % figure(1)
%         % % Create the plot
%         % subplot(2,2,jj)
%         % plot(1:length(data_urge{ii}),data_temp);        
%         % % Add markers at points 3, 5 and 10
%         % hold on;
%         % scatter(t(idx)*10000, data_temp(1,idx), 'filled', 'MarkerFaceColor', 'r');
%         % xline(10000,'--r',{'Urge Time'},'linewidth',1.2)
%         % xline(10000+diff_rt_ut(ii)*10000,'--k',{'Reaction Time'},'linewidth',1.2')
%         % ylabel('\muV');
%         % xlabel('time(second)');
%         % title(sprintf('Bipolar %s',chan_name(jj)));
%         % 
%         % % temp = all_data(jj,:,ii);
%         % % temp(~ismember(1:numel(temp), idx)) = 0;    
%         % % data_extract{ii} = temp;
%         % figure(2)
%         % plot(1:length(data_urge{ii}),data_temp);
%         % hold on;
% 
% 
%     end
% 
% end

%
% scatter(urge_time,reaction_time)
% 
% % Convert cell data to a three-dimensional matrix
% numCells1 = numel(data_in);
% [row, col] = size(data_in{1});
% 
% in_emg= zeros(row, col, numCells1);
% for i = 1:numCells1
%     in_emg(:, :, i) = data_in{i};
% end
% 
% numCells2 = numel(data_out);
% [row, col] = size(data_out{1});
% 
% out_emg= zeros(row, col, numCells2);
% for i = 1:numCells2
%     out_emg(:, :, i) = data_out{i};
% end
% 
% % Display the matrix size
% disp('Matrix size:');
% disp(size(in_emg));
% disp(size(out_emg));
% 
% subplot(211)
% for i = 1:length(in_rf)
%     plot(i*in_emg(1,1:10:end,i),'|k');
%     ylim([0.5,length(in_rf)])
%     hold on;
% end
% 
% subplot(212)
% for i = 1:length(out_rf)
%     plot(i*out_emg(1,1:10:end,i),'|r');
%     ylim([0.5,length(out_rf)])
%     hold on;
% end
% 
% %%
% % average across trials
% H1 = mean(in_emg,3);
% % average across trials
% H2 = mean(out_emg,3);
% 
% plot(H1(1,:),'r');
% hold on;
% plot(H2(1,:),'k')
% 
% 
% % smoothing 
% sigma = 3; % standard deviation of Gaussian filter
% % win_size = ceil(2*sigma)*2+2; % window size of filter
% 
% % it means 5 ms
% win_size = 20; 
% 
% % Apply Gaussian filter
% gauss_filter = gausswin(win_size, sigma); % create Gaussian filter
% smoothed_H1 = filtfilt(gauss_filter, 1, H1(1,:))/length(in_rf); % apply filter
% smoothed_H2 = filtfilt(gauss_filter, 1, H2(1,:))/length(out_rf); % apply filter
% 
% ts = 1/fs;
% % time for psth plot
% t_psth = (-1:ts:1.5-ts);
% 
% plot(t_psth,smoothed_H1,'-r');
% hold on;
% plot(t_psth,smoothed_H2,'-k')
% xlabel('Second');
% ylabel('spikes/second')

%%
% in rf
% [H1,timeVals1] = psthplot(data_in,1,[0,7],1);
% % out rf
% [H2,timeVals2] = psthplot(data_out,1,[0,7],1);
% 
% H1 = mean(in_emg,3);
% 
% figure(2)
% % subplot(211)
% plot(timeVals1,smoothed_H1,'r');
% xlabel('time (s)');
% ylabel('spikes/s');
% % xlim([1,3])
% % title("In resposne field")
% 
% % subplot(212)
% hold on;
% plot(timeVals2,smoothed_H2,'k');
% % xline(2.5,'--b',{'Stimulus Onset'},'linewidth',2,'fontweight','bold')
% % xlim([1,3])
% % title("Out- repsonse field")
% xlabel('time (s)');
% ylabel('spikes/s');
% 

%{
**** remove all those trials whose difference between the urge time and reaction
time is less than 400 ms.

%}

% for i = 1:length(urge_time)
%     % 400 milliseconds before urge time and 985 millisecond after urge time
%     data_urge{i} = Noise_filtData{i}(:,round(urge_time(i)*10000) - 4000):round(urge_time(i)*10000)+round(max(diff_rt_ut)*10000);
% 
% end
%%


