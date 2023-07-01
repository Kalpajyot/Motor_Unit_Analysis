nb_trial = TrialRecord.CurrentTrialNumber;

x = cell(1,nb_trial);
y = cell(1,nb_trial);

for i =1:nb_trial
    x{i} = eval(strcat('Trial',num2str(i),'.AnalogData.Touch(:,1)'));
    y{i} = eval(strcat('Trial',num2str(i),'.AnalogData.Touch(:,2)'));
end



x_gocue_onwards = cell(1,nb_trial);
y_gocue_onwards = cell(1,nb_trial);

stop_pos = [3000	3750	2250	4250	3500	3500	3500	5000	1500	1000	1500	1000	3500	2000	4000	4000	5000	4250	2000	3750	2000	2250	1500	2250	2000	5000	3500	3000	1000	3000].*(60/1000);

bhv_code_dial =zeros(nb_trial,1);

for i = 1:nb_trial
    bhv_code_dial(i) = round(eval(strcat('Trial',num2str(i),'.BehavioralCodes.CodeTimes(6)')));
    x_gocue_onwards{i} = x{i}(bhv_code_dial(i):end);
    y_gocue_onwards{i} = y{i}(bhv_code_dial(i):end);
end

dialAppear_time = bhv_code_dial + stop_pos(1:20)';




% x_diff = diff(x_gocue_onwards{2}(1:17:end));
% y_diff = diff(y_gocue_onwards{2}(1:17:end));


for i = 1:nb_trial
    xx = diff(x_gocue_onwards{i});
    yy = diff(y_gocue_onwards{i});  
    
    % Calculate the length of the variable without zero elements
    lengthWithoutZero_x = length(xx(xx ~= 0));
    lengthWithoutZero_y = length(yy(yy ~= 0));

%     max_xy = max(lengthWithoutZero_x,lengthWithoutZero_y);
    
    if lengthWithoutZero_x>lengthWithoutZero_y        
        id_x = find(xx~=0);
        xx_n = x_gocue_onwards{i}(id_x);
        yy_n = y_gocue_onwards{i}(id_x);
    else
        id_y = find(yy~=0);        
        yy_n = y_gocue_onwards{i}(id_y);
        xx_n = x_gocue_onwards{i}(id_y);
    end
    
    subplot(5,4,i);
%     scatter(xx_n,yy_n);
    
    vel_xy = hypot(diff(xx_n),diff(yy_n));
     
    ff_x = filtfilt(ones(4,1)/4,1,xx_n);
    ff_y = filtfilt(ones(4,1)/4,1,yy_n);
    vel_sm = diff(ff_x).^2 + diff(ff_y).^2;
    subplot(5,4,i);    
    plot(vel_xy(1:2:end),'*r');
    hold on;
    plot(vel_sm(1:2:end),'-ob','MarkerFaceColor','g','MarkerEdgeColor','b');
    legend('raw','filtered')

%     xk = x_gocue_onwards{i}(find(xx~=0));
%     yk = y_gocue_onwards{i}(find(yy~=0));
%     if length(xk)> length(yk)
%         velocity = sqrt(diff(xk(1:length(yk))).^2 + diff(yk.^2));
%         
%     else
%         velocity = sqrt(diff(xk).^2 + diff(yk(1:length(xk))).^2);
%     end
%     subplot(5,4,i);
%     plot(velocity,'-*');
    
end



% Sample stepwise x and y coordinates with random step sizes
for i = 1:nb_trial
    x = x_gocue_onwards{i}; % Replace with your own stepwise x-coordinate data
    y = y_gocue_onwards{i}; % Replace with your own stepwise y-coordinate data

    % Identify the indices of the first steps
    x_diff = diff(x);
    y_diff = diff(y);
    first_step_indices = find(x_diff ~= 0 | y_diff ~= 0) + 1;

    % Extract the values at the first steps
    x_first_step = x(first_step_indices);
    y_first_step = y(first_step_indices);

    % Display the values at the first steps
%     disp('First Step Values:');
%     disp('X-coordinate:'); disp(x_first_step);
%     disp('Y-coordinate:'); disp(y_first_step);

%     plot(x_first_step,y_first_step,'or','MarkerSize',10,'MarkerFaceColor','b')
%     xlabel('x coordinate');
%     ylabel('y coordinate')

    velocity = sqrt(diff(x_first_step).^2 + diff(y_first_step).^2);
    subplot(5,4,i)
    plot(velocity,'*-r');
end








