%
% File Description 

%{
 Author: Kalpajyoti Hazarika
 Date : '05-Jul-2023'
 This code will analyse the urge time and response time from the proactive
 reaching task

%}
%%
clear,clc
addpath('J:\onedrive\OneDrive - Indian Institute of Science\lab works\lab works\Reaching_LibetTask\behavioral data')
load('230705_sadhvika_reach_task_fix.mat')
% load('230705_sadhvika3_Proactive_reach_task_fix.mat')

% number of trial
nb_trial = TrialRecord.CurrentTrialNumber;
trial_dialPosition = zeros(nb_trial,1);

% extracting the trial with no error
Trial_id = zeros(nb_trial,1);

for i = 1: nb_trial
    Trial_errInfo = eval(strcat('Trial',num2str(i),'.TrialError'));
    if Trial_errInfo == 0
        Trial_id(i) = i;
    end
    trial_dialPosition(i) = eval(strcat('Trial',num2str(i),'.UserVars.dial_position'));
end

Trial_id = Trial_id(Trial_id ~=0);


%%
% Calculate the velcity Profile

nb_correctTrial = length(Trial_id);

x = cell(1,nb_correctTrial);
y = cell(1,nb_correctTrial);

for i =1:length(Trial_id)
    x{i} = eval(strcat('Trial',num2str(Trial_id(i)),'.AnalogData.Touch(:,1)'));
    y{i} = eval(strcat('Trial',num2str(Trial_id(i)),'.AnalogData.Touch(:,2)'));
end



x_gocue_onwards = cell(1,nb_correctTrial);
y_gocue_onwards = cell(1,nb_correctTrial);

% stop_pos = [3000	3750	2250	4250	3500	3500	3500	5000	1500	1000	1500	1000	3500	2000	4000	4000	5000	4250	2000	3750	2000	2250	1500	2250	2000	5000	3500	3000	1000	3000].*(60/1000);

dial_start =zeros(nb_correctTrial,1);
movement_end = zeros(nb_correctTrial,1);
go_cueTime = cell(1,nb_correctTrial);


for i = 1:nb_correctTrial
    dial_start(i) = round(eval(strcat('Trial',num2str(Trial_id(i)),'.BehavioralCodes.CodeTimes(5)')));
    movement_end(i) = round(eval(strcat('Trial',num2str(Trial_id(i)),'.BehavioralCodes.CodeTimes(6)')));
    x_gocue_onwards{i} = x{i}(dial_start(i):movement_end(i));
    y_gocue_onwards{i} = y{i}(dial_start(i):movement_end(i));
    lenX = length(x_gocue_onwards{i});
    go_cueTime{i} = [dial_start(i):movement_end(i)];
end


%%
% Visualization of the velocity profile

diff_factor = 17;
plot(1:diff_factor:length(x_gocue_onwards{1}),x_gocue_onwards{1}(1:diff_factor:end))
hold on;
plot(1:diff_factor:length(y_gocue_onwards{1}),y_gocue_onwards{1}(1:diff_factor:end))

gocue_x = x_gocue_onwards;
gocue_y = y_gocue_onwards;

sm_vel = cell(nb_correctTrial,1);


rt_vel_method = zeros(nb_correctTrial,1);

for i = 1:nb_correctTrial
    
    fx = diff(gocue_x{i}(1:diff_factor:end));
    fy = diff(gocue_y{i}(1:diff_factor:end));
    vel = hypot(fx,fy);

    subplot(6,5,i)
    plot(vel,'k')
    hold on;
    ff_x = filtfilt(ones(5,1)/5,1,gocue_x{i}(1:diff_factor:end));
    ff_y = filtfilt(ones(5,1)/5,1,gocue_y{i}(1:diff_factor:end));
    sm_vel{i} = (hypot(diff(ff_x),diff(ff_y))'./diff_factor)*1000;
    t_new = go_cueTime{i}(1:diff_factor:end);

    plot(sm_vel{i},'r');

    max_vel = max(sm_vel{i});
    idx_max_vel = find(sm_vel{i} == max(sm_vel{i}), 1, 'first');

    [TF_peaks,peak_vel] = islocalmax(sm_vel{i},'MinProminence',20,'MinSeparation',50);    % for multiple peaks
    idx_peaks = find(TF_peaks ~= 0);

    if isempty(idx_peaks)
       idx_peaks = idx_max_vel;
    end

    idx_vel_10 = find(sm_vel{i} > 10, 1, 'first');

    vel_backwards = flip(sm_vel{i}(1:idx_vel_10));
    
    idx_vel_0 = find(vel_backwards == 0, 1, 'first');

    idx_vel_onset = idx_vel_10 - idx_vel_0 + 1;

    if isempty(idx_vel_onset)
       idx_vel_onset = 1;
    end   
    

    rt_vel_method(i) = t_new(idx_vel_onset)-go_cueTime{i}(1);

end


%%
% I tried this way also 
% rt_vel_method = zeros(nb_correctTrial,1);
% 
% 
% for i = 1:nb_correctTrial
% 
%     % Direct calculattion of x and y coordiantes using Diff function
%     xx = diff(x_gocue_onwards{i});
%     yy = diff(y_gocue_onwards{i});  
% 
%     % Calculate the length of the variable without zero elements
%     lengthWithoutZero_x = length(xx(xx ~= 0));
%     lengthWithoutZero_y = length(yy(yy ~= 0));
% 
% %     max_xy = max(lengthWithoutZero_x,lengthWithoutZero_y);
% 
%     if lengthWithoutZero_x>lengthWithoutZero_y        
%         id_x = find(xx~=0);
%         xx_n = x_gocue_onwards{i}(id_x);
%         yy_n = y_gocue_onwards{i}(id_x);
%         go_cue_t = go_cueTime{i}(id_x);        
%     else
%         id_y = find(yy~=0);        
%         yy_n = y_gocue_onwards{i}(id_y);
%         xx_n = x_gocue_onwards{i}(id_y);
%         go_cue_t = go_cueTime{i}(id_y);
% 
%     end
% 
%     subplot(6,5,i);
%     scatter(xx_n,yy_n);
% 
%     vel_xy = hypot(diff(xx_n),diff(yy_n));
% 
%     ff_x = filtfilt(ones(5,1)/5,1,xx_n);
%     ff_y = filtfilt(ones(5,1)/5,1,yy_n);
%     vel_sm = ((diff(ff_x).^2 + diff(ff_y).^2)./diff(go_cue_t)')*1000;
%     t_new = (go_cue_t(2:end));
% 
%     % 
%     % subplot(6,5,i);    
%     % plot(vel_xy(1:end),'-*r');
%     % hold on;
%     % plot(vel_sm(1:end),'-ob','MarkerFaceColor','g','MarkerEdgeColor','b');
%     % legend('raw','filtered')
% 
% %     xk = x_gocue_onwards{i}(find(xx~=0));
% %     yk = y_gocue_onwards{i}(find(yy~=0));
% %     if length(xk)> length(yk)
% %         velocity = sqrt(diff(xk(1:length(yk))).^2 + diff(yk.^2));
% %         
% %     else
% %         velocity = sqrt(diff(xk).^2 + diff(yk(1:length(xk))).^2);
% %     end
% %     subplot(5,4,i);
% %     plot(velocity,'-*');
% 
% 
%     max_vel = max(vel_sm);
%     idx_max_vel = find(vel_sm == max(vel_sm), 1, 'first');
%     [TF_peaks,peak_vel] = islocalmax(vel_sm,'MinProminence',20,'MinSeparation',50);    % for multiple peaks
%     idx_peaks = find(TF_peaks ~= 0);
% 
%     if isempty(idx_peaks)
%        idx_peaks = idx_max_vel;
%     end
% 
%     idx_vel_10 = find(vel_sm > 10, 1, 'first');
% 
%     % Searching backwards
%     vel_backwards = flip(vel_sm(1:idx_vel_10));
%     idx_vel_0 = find(vel_backwards == 0, 1, 'first');
% 
%     idx_vel_onset = idx_vel_10 - idx_vel_0 + 1;
% 
%     if isempty(idx_vel_onset)
%        idx_vel_onset = 1;
%     end   
% 
%     rt_vel_method(i) = t_new(idx_vel_onset)-go_cue_t(1);    
% 
% end
% 

% Sample stepwise x and y coordinates with random step sizes

% rt_vel_method = zeros(nb_correctTrial,1);
% 
% for i = 1:nb_correctTrial
%     x = x_gocue_onwards{i}; % Replace with your own stepwise x-coordinate data
%     y = y_gocue_onwards{i}; % Replace with your own stepwise y-coordinate data
% 
%     % Identify the indices of the first steps
%     x_diff = diff(x);
%     y_diff = diff(y);
%     first_step_indices = find(x_diff ~= 0 | y_diff ~= 0) + 1;
% 
%     % Extract the values at the first steps
%     x_first_step = x(first_step_indices);
%     y_first_step = y(first_step_indices);
% 
%     % Display the values at the first steps
% %     disp('First Step Values:');
% %     disp('X-coordinate:'); disp(x_first_step);
% %     disp('Y-coordinate:'); disp(y_first_step);
% 
% %     plot(x_first_step,y_first_step,'or','MarkerSize',10,'MarkerFaceColor','b')
% %     xlabel('x coordinate');
% %     ylabel('y coordinate')
% 
%     velocity = hypot(diff(x_first_step),diff(y_first_step))./diff(go_cueTime{i})*1000;
% 
%     % smoothed version
%     ff_x1 = filtfilt(ones(5,1)/5,1,x_first_step);
%     ff_y1 = filtfilt(ones(5,1)/5,1,y_first_step);
%     vel_sm1 = (hypot(diff(ff_x1),diff(ff_y1))./diff(go_cueTime{i}))*1000;
%     t_new = (go_cueTime{i}(2:end));
% 
%     % subplot(6,5,i)
%     % plot(velocity,'*-r');
%     % hold on;
%     % plot(vel_sm1,'*-k')
% 
%     max_vel = max(vel_sm1);
%     idx_max_vel = find(vel_sm1 == max(vel_sm1), 1, 'first');
%     [TF_peaks,peak_vel] = islocalmax(vel_sm1,'MinProminence',20,'MinSeparation',50);    % for multiple peaks
%     idx_peaks = find(TF_peaks ~= 0);
% 
%     if isempty(idx_peaks)
%        idx_peaks = idx_max_vel;
%     end
% 
%     idx_vel_10 = find(vel_sm1 > 10, 1, 'first');
% 
%     % Searching backwards
%     vel_backwards = flip(vel_sm1(1:idx_vel_10));
%     idx_vel_0 = find(vel_backwards == 0, 1, 'first');
%     idx_vel_onset = idx_vel_10 - idx_vel_0 + 1;
%     if isempty(idx_vel_onset)
%        idx_vel_onset = 1;
%     end   
% 
%     % rt_vel_method(i) = t_new(idx_vel_onset);
%     rt_vel_method(i) = t_new(idx_vel_onset) - go_cueTime{i}(1);
% 
% end



%%
urge_pos = [4.25	6.25	11.5	5.75	8.25	12.75	7.25	4.5	5.75	6.75	6.25	12.75	8.25	5.25	10.25	8.5	11.25	8.5	8.25	4.25	6.25	12.25	7.25	5.75	12.5	9.5	12.75	7.25	3.25];

dial_position_seq(11) = [];
dial_start = dial_position_seq/5;

diff_ut_rt = zeros(1,length(urge_pos));

for i = 1:length(dial_start)
    if dial_start(i)< urge_pos(i)
        diff_ut_rt(i) =  (urge_pos(i) - dial_start(i))*0.5;
    else
        diff_ut_rt(i) = ((12-dial_start(i))+urge_pos(i))*0.5;
    end

end
%%





















