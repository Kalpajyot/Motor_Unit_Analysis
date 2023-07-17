function [outputData,valuesToMark,t_data,data_analysis] = MeanSd_signal(input_data,rf_type,urge_time,sd_value,fs)

    outputData = cell(length(rf_type),1);
    data_analysis = cell(length(rf_type),1);
    valuesToMark = cell(length(rf_type),1);

    for ii = 1:length(outputData)
        % clf;
        % figure(ii)
        for jj = 1:4
            
            befor_stim = input_data(jj,1:15000,rf_type(ii));
    
            threshold_p = mean(befor_stim) + sd_value * std(befor_stim);
            threshold_n = mean(befor_stim) - sd_value * std(befor_stim);
        
            % peak location and amplitude
            data_analysis{ii}(jj,:) = input_data(jj,round(20000+urge_time(rf_type(ii))*10000-9999):round(20000+urge_time(rf_type(ii))*10000+15000),rf_type(ii));
            data_temp = data_analysis{ii}(jj,:);
            ts = 1/fs;
            t_data = (-1:ts:1.5-ts);
            % 
            % [pks_p,locs_p] = findpeaks(data_temp,'MinPeakHeight',threshold_p,'MinPeakDistance',100);  
            % [pks_n,locs_n] = findpeaks(-data_temp,'MinPeakHeight',-threshold_n,'MinPeakDistance',100);
            % 
            % idx = sort([locs_p,locs_n]);
            pos_peak = find(data_temp>threshold_p);
            neg_peak = find(data_temp<threshold_n);
            data_temp([pos_peak,neg_peak]) = 1;
            data_temp(~ismember(1:numel(data_temp),[pos_peak,neg_peak])) = 0;

            outputData{ii}(jj,:) = data_temp;
            % all trials plus all channels
            valuesToMark{ii}{jj} = sort([pos_peak,neg_peak]);
            % subplot(2,2,jj)
            % plot(t_data,data_analysis);
            % hold on;
            % plot(t_data(valuesToMark), data_analysis(valuesToMark), 'ro', 'MarkerSize', 1.2, 'MarkerFaceColor', 'r');
            % xline(t_data(10000),'--b',{'Urge Time'},'linewidth',1.2)
            % % xline(t_data(10000 + round(diff_rt_ut(rf_type(ii))*10000)),'--r',{'Reaction Time'},'linewidth',1.2)
            % title(sprintf('Trial %1.0f',rf_type(ii)));
            % % Mark specific values
            % valuesToMark = sort([pos_peak,neg_peak])./10000;
            % for i = 1:length(valuesToMark)
            %     value = valuesToMark(i);
            %     % index = find(t_data == value);
            %     if ~isempty(index)
            %         plot(t_data(valuesToMark), y(index), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            %         % text(x(index), y(index), num2str(value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            %     end
            % end
            
            % idx = sort([locs_p,locs_n]);
            % figure(1)
            % % Create the plot
            % subplot(2,2,jj)
            % plot(1:length(data_urge{ii}),data_temp);        
            % % Add markers at points 3, 5 and 10
            % hold on;
            % scatter(t(idx)*10000, data_temp(1,idx), 'filled', 'MarkerFaceColor', 'r');
            % xline(10000,'--r',{'Urge Time'},'linewidth',1.2)
            % xline(10000+diff_rt_ut(ii)*10000,'--k',{'Reaction Time'},'linewidth',1.2')
            % ylabel('\muV');
            % xlabel('time(second)');
            % title(sprintf('Bipolar %s',chan_name(jj)));
            % 
            % % temp = all_data(jj,:,ii);
            % % temp(~ismember(1:numel(temp), idx)) = 0;    
            % % data_extract{ii} = temp;
            % figure(2)
            % plot(1:length(data_urge{ii}),data_temp);
            % hold on;
    
    
        end
    
    end
end



