function LFPpop_sim_temp(exp)
% Only to extract and plot spectrograms and coherence over time.

%% load, extract, and save relevant variables

% path = 'D:\output_data';
% cd(path)

% fprintf(['Time:  ' num2str(clock) '\n']);
% fnames = dir('experiments*.mat');
% cnt=1;
% for i = 1:length(fnames)
%     fprintf(['****   Loading file ' num2str(fnames(i).name) '   ****' '\n'])
%     load(fnames(i).name);
%     for sess = 1:length(experiments.sessions)
%         %         if ~isempty(experiments.sessions(sess).lfps(1).stats)
%         if ~isempty(experiments.sessions(sess).lfps(1).trials)
%             %% Read areas
%             if experiments.sessions(sess).monk_id == 53
%                 for nlfps = 1:length(experiments.sessions(1).lfps)
%                     indx_MST(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'MST');
%                     indx_PFC(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'PFC');
%                     indx_PPC(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'PPC');
%                 end
%             else
%                 for nlfps = 1:length(experiments.sessions(1).lfps)
%                     indx_MST(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'MST');
%                     indx_PPC(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'PPC');
%                 end
%             end
%             
%             %% extract per area
%             if experiments.sessions(sess).monk_id == 53
%                 exp(cnt).area.MST.lfps.stats = [experiments.sessions(sess).lfps(indx_MST).stats];
%                 exp(cnt).area.PFC.lfps.stats = [experiments.sessions(sess).lfps(indx_PFC).stats];
%                 exp(cnt).area.PPC.lfps.stats = [experiments.sessions(sess).lfps(indx_PPC).stats];
%                 % %                 exp(cnt).area.MST.signal = experiments.sessions(sess).lfps(indx_MST); % extract lfp signal and band passed signal
%                 % %                 exp(cnt).area.PPC.signal = experiments.sessions(sess).lfps(indx_PPC); % extract lfp signal and band passed signal
%                 % %                 exp(cnt).area.PFC.signal = experiments.sessions(sess).lfps(indx_PFC); % extract lfp signal and band passed signal
%                 exp(cnt).behavior = [experiments.sessions(sess).behaviours];
%             else
%                 exp(cnt).area.MST.lfps.stats = [experiments.sessions(sess).lfps(indx_MST).stats]; 
%                 exp(cnt).area.PPC.lfps.stats = [experiments.sessions(sess).lfps(indx_PPC).stats];
%                 % %                 exp(cnt).area.MST.signal = []; exp(cnt).area.MST.signal = experiments.sessions(sess).lfps(indx_MST); % extract lfp signal and band passed signal
%                 % %                 exp(cnt).area.PPC.signal = experiments.sessions(sess).lfps(indx_PPC); % extract lfp signal and band passed signal
%                 exp(cnt).behavior = [experiments.sessions(sess).behaviours];
%             end
%             exp(cnt).pop = experiments.sessions(sess).populations.lfps.stats;
%             exp(cnt).monk_id = experiments.sessions(sess).monk_id; % extract lfp signal and band passed signal
%         end
%         cnt=cnt+1;
%     end
%     disp('Clearing experiments... . . .')
%     clear experiments
% end
% 
% % save
% disp('Saving... . . .')
% fprintf(['Time:  ' num2str(clock) '\n']);
% save('exp_out_lfp_stats_pop_2020_08_20_coherence_align_stop','exp', '-v7.3');
% % load train
% % sound(y,Fs)
% disp('           Saved! ')
% fprintf(['Time:  ' num2str(clock) '\n']);
% disp('Extracting... . . .')

%% Choose what to analyze and save
save_behavior = true;
extract_lfp_raw = false; % raw and per trial lfps
save_lfp_raw = false; % raw and per trial lfps
do_PSD = false;  % extract power spectral densities
save_spectro = false; % save spectrogram file?
save_spectro_per_trial = false;
save_spectro_per_trial_align_stop = false;
avg_monks = false; % average for all monkeys?
do_cohero = true; % extract coherograms
doCSD = false; % Perform CSD analysis for MST recordings?
do_ERP = false; % extract ERPs (evoked LFPs)
save_lfp_band_pass = false; % extract band passed lfp signal only
do_band_passed = false; % extract and analyse band passed signal
do_band_passed_vs_accuracy = false;

monk = []; all_monks=[];
%% Save behavior
if save_behavior
    monks = unique([exp.monk_id]);
    for ii = 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            monk(ii).behavior(j) = p_monk(j).behavior;
        end
    end
end
%% PSD
%% avg per session
if do_PSD
    trialtype = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype); events = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype.all.events);
    for i = 1:length(exp) % num of sessions
        nareas = numel(fieldnames(exp(i).area)); areas = fieldnames(exp(i).area);
        for a = 1:length(areas) % num of areas
            for type = 1:length(trialtype) % num of trial types
                if ~isempty(exp(i).area.(areas{a}).lfps.stats)
                    nconds = length(exp(i).area.(areas{a}).lfps.stats(1).trialtype.(trialtype{type})); clear cond
                    for cond = 1:nconds
                        clear pow pow_eye
                        for ch = 1:length(exp(i).area.(areas{a}).lfps.stats) % extract per channel
                            if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
                                pow_eye(ch,:) = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).spectrum.psd;
                                freq_eye = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).spectrum.freq;
                            else
                                pow(ch,:) = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).spectrum.psd;
                            end
                        end
                        % store
                        freq = exp(i).area.(areas{a}).lfps.stats(1).trialtype.all.spectrum.freq;
                        if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
                            psden_eye(i).area.(areas{a}).(trialtype{type})(cond).ch = pow_eye;
                            psden_eye(i).area.(areas{a}).(trialtype{type})(cond).mu = nanmean(pow_eye);
                            psden_eye(i).area.(areas{a}).(trialtype{type})(cond).sem = nanstd(pow_eye)/sqrt(length(ch));
                            psden_eye(i).area.(areas{a}).(trialtype{type})(cond).max = max(max(nanmean(pow_eye)));
                        else
                            psden(i).area.(areas{a}).(trialtype{type})(cond).ch = pow;
                            psden(i).area.(areas{a}).(trialtype{type})(cond).mu = nanmean(pow);
                            psden(i).area.(areas{a}).(trialtype{type})(cond).sem = nanstd(pow)/sqrt(length(ch));
                            psden(i).area.(areas{a}).(trialtype{type})(cond).max = max(max(nanmean(pow)));
                        end
                    end
                end
            end
        end
        psden(i).monk_id = exp(i).monk_id;
    end
    
    %% normalize by max
    % get max values
    % cnt=1;
    % for i = 1:length(exp) % num of sessions
    %     nareas = numel(fieldnames(exp(i).area)); areas = fieldnames(exp(i).area);
    %     for a = 1:length(areas) % num of areas
    %         for type = 1:length(trialtype) % num of trial types
    %             if ~isempty(exp(i).area.(areas{a}).lfps.stats)
    %                 nconds = length(exp(i).area.(areas{a}).lfps.stats(1).trialtype.(trialtype{type})); clear cond
    %                 for cond = 1:nconds
    %                     if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
    %                         max_psd(cnt) = psden_eye(i).area.(areas{a}).(trialtype{type})(cond).max;
    %                     else
    %                         max_psd(cnt) = psden(i).area.(areas{a}).(trialtype{type})(cond).max;
    %                     end
    %                     cnt = cnt+1;
    %                 end
    %             end
    %         end
    %     end
    % end
    % %maxPSDval = max(max_psd(max_psd<1000));
    % maxPSDval = 1; % take the max between 0 and 50 Hz
    %
    % % normalize for max of max (between all conditions)
    % for i = 1:length(exp) % num of sessions
    %     nareas = numel(fieldnames(exp(i).area)); areas = fieldnames(exp(i).area);
    %     for a = 1:length(areas) % num of areas
    %         for type = 1:length(trialtype) % num of trial types
    %             if ~isempty(exp(i).area.(areas{a}).lfps.stats)
    %                 nconds = length(exp(i).area.(areas{a}).lfps.stats(1).trialtype.(trialtype{type})); clear cond
    %                 for cond = 1:nconds
    %                     if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
    %                         psden_eye(i).area.(areas{a}).(trialtype{type})(cond).mu_norm = psden_eye(i).area.(areas{a}).(trialtype{type})(cond).mu; %/maxPSDval;
    %                     else
    %                         psden(i).area.(areas{a}).(trialtype{type})(cond).mu_norm = psden(i).area.(areas{a}).(trialtype{type})(cond).mu; %/maxPSDval;
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
    
    %% avg between sessions for each monkey
    monks = unique([psden.monk_id]);
    for i = 1:length(monks)
        m = [psden.monk_id] == monks(i); p_monk = psden(m); p_monk_eye = psden_eye(m);
        nareas = numel(fieldnames(p_monk(1).area)); areas = fieldnames(p_monk(1).area);
        for a = 1:length(areas)
            trialtype = [fieldnames(p_monk(1).area.(areas{a})); fieldnames(p_monk_eye(1).area.(areas{a}))];
            for type = 1:length(trialtype)
                if type == 6 | type == 7;  % eyes
                    nconds = 1; clear cond
                else
                    nconds = length(p_monk(1).area.(areas{a}).(trialtype{type})); clear cond
                end
                for cond = 1:nconds
                    clear psd_eye_all psd_all
                    for k = 1:length(p_monk)
                        if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
                            psd_eye_all(k,:) = p_monk_eye(k).area.(areas{a}).(trialtype{type})(cond).mu;
                        else
                            psd_all(k,:) = p_monk(k).area.(areas{a}).(trialtype{type})(cond).mu;
                        end
                    end
                    if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
                        monk(i).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = nanmean(psd_eye_all); % average per trial type
                        monk(i).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = nanstd(psd_eye_all);
                    else
                        if size(psd_all,1)<2
                            monk(i).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = psd_all;
                            monk(i).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = psd_all;
                        else
                            monk(i).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = nanmean(psd_all);  % average per trial type
                            monk(i).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = nanstd(psd_all);
                        end
                    end
                end
            end
        end
        monk(i).pw.freq = exp(i).area.(areas{a}).lfps.stats(1).trialtype.all.spectrum.freq;
        monk(i).pw.freq_eye = exp(i).area.(areas{a}).lfps.stats(1).trialtype.eyesfree.spectrum.freq;
        monk(i).pw.monk_id = unique([psden(m).monk_id]);
    end
    
    %% avg across monkeys
    if avg_monks
        clear all_m all_m_eye
        for i=1:length(monks), area_m{i} = fieldnames(monk(i).pw.area); end
        areas = intersect(string(area_m{1}), string(area_m{2}));
        for a = 1:length(areas)
            for type = 1:length(trialtype)
                nconds = length(monk(i).pw.area.(areas{a}).(trialtype{type})); clear cond
                for cond = 1:nconds
                    for i = 1:length(monk)
                        if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
                            all_m_eye(i,:) = monk(i).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess;
                        else
                            all_m(i,:) = monk(i).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess;
                        end
                        
                        if strcmp((trialtype{type}), 'eyesfree') | strcmp((trialtype{type}), 'eyesfixed')
                            all_monks.pw.area.(areas{a}).(trialtype{type})(cond).mu = nanmean(all_m_eye);
                            all_monks.pw.area.(areas{a}).(trialtype{type})(cond).sem = nanstd(all_m_eye);
                        else
                            all_monks.pw.area.(areas{a}).(trialtype{type})(cond).mu = nanmean(all_m);
                            all_monks.pw.area.(areas{a}).(trialtype{type})(cond).sem = nanstd(all_m);
                        end
                    end
                end
            end
        end
        %%
        all_monks.pw.freq = monk(1).pw.freq;
        all_monks.pw.freq_eye = monk(1).pw.freq_eye;
    end
end
%% spectrogram
% average across channels
if save_spectro
    trialtype = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype); events = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype.all.events);
    for i = 1:length(exp) % num of sessions
        nareas = numel(fieldnames(exp(i).area)); areas = fieldnames(exp(i).area);
        for a = 1:length(areas) % num of areas
            for type = 1:3 % length(trialtype) % num of trial types
                if ~isempty(exp(i).area.(areas{a}).lfps.stats)
                    nconds = length(exp(i).area.(areas{a}).lfps.stats(1).trialtype.(trialtype{type})); clear cond
                    for cond = 1:nconds
                        for ev = 1:length(events)
                            clear spec
                            for ch = 1:length(exp(i).area.(areas{a}).lfps.stats) % extract per channel
                                spec(ch,:,:) = real(exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).p_spectrogram');
                            end
                            % store and normalize by max
                            p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu = squeeze(nanmean(spec,1))/max(max(squeeze(nanmean(spec,1))));
                            p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu_no_norm = squeeze(nanmean(spec,1));
                            p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).std = squeeze(nanstd(spec,1)/sqrt(length(ch)));
                            p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).freq_spectrogram;
                            p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).ts_spectrogram;
                        end
                    end
                end
            end
        end
        p_spec(i).monk_id = exp(i).monk_id;
    end
    
    % plot -- sanity check MST
    % monkey_num = 1; area = 'PPC' ; ev = 'move'
    % freq = p_spec(monkey_num).area.(area).all.events.(ev).freq;
    % ts =  p_spec(monkey_num).area.(area).all.events.(ev).ts;
    % p_spectro = p_spec(monkey_num).area.(area).all.events.(ev).pow_mu;
    % % plot
    % figure('Name',['Monk ' num2str(monkey_num) ' ' area]);
    % imagesc(ts-1,freq,p_spectro, [0 160]); axis xy; colorbar;
    % set(gca,'xlim',[-0.25 0.25], 'ylim',[4 80], 'FontSize', 22)
    % xlabel('time (s)'); ylabel('frequency (Hz)')
    
    %% avg between sessions for each monkey
    monks = unique([p_spec.monk_id]);
    for i = 1:length(monks)
        m = [p_spec.monk_id] == monks(i); p_monk = p_spec(m);
        nareas = numel(fieldnames(p_monk(1).area)); areas = fieldnames(p_monk(1).area);
        for a = 1:length(areas)
            trialtype = fieldnames(p_monk(1).area.(areas{a}));
            for type = 1:length(trialtype)
                nconds = length(p_monk(1).area.(areas{a}).(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ev = 1:length(events)
                        clear spec spec_no_norm spec_f spec_ts
                        for k = 1:length(p_monk)
                            spec(k,:,:) = p_monk(k).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu;
                            spec_no_norm(k,:,:) = p_monk(k).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu_no_norm;
                            spec_f(k,:,:) = p_monk(k).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq;
                            spec_ts(k,:,:) = p_monk(k).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts;
                        end
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).mu_sess = squeeze(nanmean(spec,1));  % average per trial type
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).mu_sess_no_norm = squeeze(nanmean(spec_no_norm,1));  % average per trial type
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).std_sess = squeeze(nanstd(spec,1));
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq_sess = squeeze(spec_f(1,:));
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts_sess = squeeze(spec_ts(1,:));
                    end
                end
                
            end
        end
    end
    
    %% avg across monkeys
    if avg_monks
        clear all_m
        for i=[1 3] %1:length(monk), area_m{i} = fieldnames(monk(i).spec.area); end  % if size(area_m{i},1) == 1,  end
            areas = intersect(string(area_m{1}), string(area_m{3})); % hardcoded for now, only look at Quigley(44) and Schro(51) % change to 2 if less than 3 monkeys.
            for a = 1:length(areas)
                for type = 1:length(trialtype)
                    nconds = length(monk(3).spec.area.(areas{a}).(trialtype{type})); clear cond
                    for cond = 1:nconds
                        for ev = 1:length(events)
                            for i = [1 3] %1:length(monk)
                                all_m(i,:,:) = monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).mu_sess;
                            end
                            all_m(2,:,:) = []; % hardcoded, change!
                            all_monks.spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).mu = squeeze(nanmean(all_m,1));
                        end
                    end
                end
            end
        end
    end
end
%% Gather spectrograms for each trial per electrode per monkey aligned to target onset
if save_spectro_per_trial
    monks = unique([exp.monk_id]); clear p_monk
    for i = 1:length(monks)
        m = [exp.monk_id] == monks(i); p_monk = exp(m);
        for sess = 1:length(p_monk) % num of sessions
            trialtype = fieldnames(p_monk(sess).pop.trialtype);
            for type = 1:length(trialtype) % num of trial types
                for cond = 2%1:length(p_monk(sess).pop.trialtype.(trialtype{type}));
                    areas = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area);
                    for a=1:length(areas);
                        for ch = 1:length(p_monk(sess).pop.trialtype.reward(1).area.(areas{a}).ch)
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).ch(ch).trl = p_monk(sess).pop.trialtype.reward(cond).area.(areas{a}).ch(ch).trl;
                        end
                    end
                end
            end
        end
    end
    
    %% Gather spectrograms for all electrodes each trial per monkey aligned to target onset
    monks = unique([exp.monk_id]); clear p_monk
    for i = 1:length(monks)
        m = [exp.monk_id] == monks(i); p_monk = exp(m);
        for sess = 1:length(p_monk) % num of sessions
            trialtype = fieldnames(p_monk(sess).pop.trialtype);
            for type = 1:length(trialtype) % num of trial types
                for cond = 1:length(p_monk(sess).pop.trialtype.(trialtype{type}))
                    areas = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area);
                    for a=1:length(areas);
                        for trl = 1:length(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl)
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).pw_trl(trl).ts = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl(trl).ts_spectrogram;
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).pw_trl(trl).freq = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl(trl).freq_spectrogram;
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).pw_trl(trl).spectro = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl(trl).spectrogram;
                        end
                    end
                end
            end
        end
    end
end

%% Gather spectrograms for each trial per electrode per monkey aligned to movement stop
if save_spectro_per_trial_align_stop
    monks = unique([exp.monk_id]); clear p_monk
    for i = 1:length(monks)
        m = [exp.monk_id] == monks(i); p_monk = exp(m);
        for sess = 1:length(p_monk) % num of sessions
            trialtype = fieldnames(p_monk(sess).pop.trialtype);
            for type = 1:length(trialtype) % num of trial types
                for cond = 2%1:length(p_monk(sess).pop.trialtype.(trialtype{type}));
                    areas = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area);
                    for a=1:length(areas)
                        for ch = 1:length(p_monk(sess).pop.trialtype.reward(1).area.(areas{a}).ch)
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).ch(ch).trl = p_monk(sess).pop.trialtype.reward(cond).area.(areas{a}).ch(ch).trl;
                        end
                    end
                end
            end
        end
    end
    
    %% Gather spectrograms for all electrodes each trial per monkey aligned to target onset
    monks = unique([exp.monk_id]); clear p_monk
    for i = 1:length(monks)
        m = [exp.monk_id] == monks(i); p_monk = exp(m);
        for sess = 1:length(p_monk) % num of sessions
            trialtype = fieldnames(p_monk(sess).pop.trialtype);
            for type = 1:length(trialtype) % num of trial types
                for cond = 1:length(p_monk(sess).pop.trialtype.(trialtype{type}))
                    areas = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area);
                    for a=1:length(areas)
                        for trl = 1:length(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl)
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).pw_trl(trl).ts = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl(trl).ts_spectrogram_stop;
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).pw_trl(trl).freq = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl(trl).freq_spectrogram_stop;
                            monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).pw_trl(trl).spectro = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).pop_trl(trl).spectrogram_stop;
                        end
                    end
                end
            end
        end
    end
end

%% coherogram
% average for each session
if do_cohero
    monks = unique([exp.monk_id]);
    for i = 1:length(monks) % [1 3]
        m = [exp.monk_id] == monks(i); p_monk = exp(m);
        coh_areas = fieldnames(p_monk(i).pop.trialtype.reward(1).events.stop);
        trialtype = fieldnames(p_monk(i).pop.trialtype);
        for sess = 1:length(p_monk)
            for type = 1:length(trialtype)
                nconds = length(p_monk(sess).pop.trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    events = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).events);
                    for ev = 1:length(events)
                        clear coh_ar coh_ev coh_phi coh_ts coh_freq
                        for coh_ar = 1:length(coh_areas)
                            coh_ev(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher;
                            coh_phi(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coherPhi;
                            coh_ts(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher_ts;
                            coh_freq(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher_freq;
                        end
                        % plot -- sanity
                        %             figure('Name', 'PPC-->MST');
                        %             imagesc(coh_ts(2,:)-1, coh_freq(2,:),squeeze(coh_ev(2,:,:)));
                        %             axis xy; set(gca,'xlim',[-0.25 0.25], 'ylim', [4 50],'FontSize', 22); colorbar;
                        %
                        
                        
                        % coherogram mean
                        for coh_ar = 1:length(coh_areas)
                            monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher = squeeze(coh_ev(coh_ar,:,:));
                            monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher_phi = squeeze(coh_phi(coh_ar,:,:));
                            monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher_ts = squeeze(coh_ts(coh_ar,:,:));
                            monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).coher_freq = squeeze(coh_freq(coh_ar,:,:));
                        end

                    end
                end
            end
        end
    end
end

%% Compute CSD for each session
if doCSD
    monks = unique([exp.monk_id]);
    for ii = [1 3]; %1:length(monks)
        areas = fieldnames(p_monk(ii).area); trialtype = fieldnames(p_monk(ii).pop.trialtype); events = fieldnames(p_monk(ii).area.MST.lfps.stats(1).trialtype.all.events);
        m = [exp.monk_id] == monks(ii); sess = exp(m);
        for s = 1:length(sess)
            for type = 1:length(trialtype)
                nconds = length(p_monk(ii).area.MST.lfps.stats(1).trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ev = 1:length(events)
                        t = sess(s).area.MST.lfps.stats(1).trialtype.(trialtype{type})(cond).events.(events{ev}).time; %t_win = t(t >= -0.5 & t <= 0.5);
                        clear wave
                        for ch = 1:length(p_monk(s).area.MST.lfps.stats)
                            wave(ch,:) = sess(s).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu; %(t >= -0.5 & t <= 0.5);
                        end
                        % compute csd
                        monk(ii).csd_sess(s).area.MST.trialtype.(trialtype{type})(cond).events.(events{ev}).csd = computecsd(wave, t, 0.1); % 0.1 electrode distance
                        monk(ii).csd_sess(s).area.MST.trialtype.(trialtype{type})(cond).events.(events{ev}).time = t;
                    end
                end
            end
        end
    end
end
%% Gather pop analyses for monks with MST rec
if do_ERP
    monks = unique([exp.monk_id]);
    for ii = [1 3]; % 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            monk(ii).pop(j) = p_monk(j).pop;
            for type = 1:length(trialtype)
                nconds = length(p_monk(j).area.PPC.lfps.stats(1).trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ch = 1:length(p_monk(j).area.MST.lfps.stats)
                        % cont
                        %                     monk(ii).cont.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous; %gather continuous
                        for ev = 1:length(events)
                            % ERPs
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).lfp_align = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).lfp_align;
                            if extract_lfp_raw
                                monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                            end
                        end
                    end
                    for ch = 1:length(p_monk(j).area.PPC.lfps.stats)
                        % cont
                        %                     monk(ii).cont.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous;
                        for ev = 1:length(events)
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            if extract_lfp_raw
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                            end
                        end
                    end
                end
            end
        end
    end
end
%% Gather pop analyses for monks without MST rec
if do_ERP
    monks = unique([exp.monk_id]);
    for ii = 2  %Bruno 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            monk(ii).pop(j) = p_monk(j).pop;
            for type = 1:length(trialtype)
                nconds = length(p_monk(j).area.PPC.lfps.stats(1).trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ch = 1:length(p_monk(j).area.PPC.lfps.stats)
                        % cont
                        %                     monk(ii).cont.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous;
                        for ev = 1:length(events)
                            % ERPs
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            if extract_lfp_raw
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                            end
                        end
                    end
                end
            end
        end
    end
end
%% Gather erp's for PFC on Schro for now
if do_ERP
    monks = unique([exp.monk_id]);
    for ii = 3  %3 Schro  1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            monk(ii).pop(j) = p_monk(j).pop;
            for type = 1:length(trialtype)
                nconds = length(p_monk(j).area.PFC.lfps.stats(1).trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ch = 1:length(p_monk(j).area.PFC.lfps.stats)
                        % cont
                        %                     monk(ii).cont.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous;
                        for ev = 1:length(events)
                            % ERPs
                            monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            if extract_lfp_raw
                                monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Get ERPs separated by move on time MST
    for ii = 1:length(monks)
        monks = unique([exp.monk_id]);
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            % select trials
            correct = p_monk(j);
        end
    end
end

if save_lfp_band_pass
    monks = unique([exp.monk_id]); clear p_monk
    for ii = 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for s = 1:length(p_monk)
            areas = fieldnames(p_monk(s).area);
            for a = 1:length(areas)
                for trl = 1:length(p_monk)
                    monk(ii).sess(s).area.(areas{a}).lfp = p_monk(s).area.(areas{a}).signal;
                end
            end
        end
    end
end

%% Analyses on band passed trials
if do_band_passed
    disp('                 Do band passed . . .     ')
    fprintf(['Time:  ' num2str(clock) '\n']);
    monks = unique([exp.monk_id]);
    for m = 1 %1:length(monks)
        ii = [exp.monk_id] == monks(m); p_monk = exp(ii);
        nsess = 1:length(p_monk);
        for s = 1:length(nsess)
            corr = p_monk(s).behavior.stats.trialtype.reward(strcmp({p_monk(s).behavior.stats.trialtype.reward.val},'rewarded')).trlindx;
            incorr = p_monk(s).behavior.stats.trialtype.reward(strcmp({p_monk(s).behavior.stats.trialtype.reward.val},'unrewarded')).trlindx;
            behv_corr = p_monk(s).behavior.trials(corr);
            behv_incorr = p_monk(s).behavior.trials(incorr);
            areas = fieldnames(p_monk(s).area); nareas = length(areas);
            for ar = 1:nareas
                for nlfp = 1:length(p_monk(s).area.(areas{ar}).signal)
                    if size(p_monk(s).area.(areas{ar}).signal,2) ~= 0
                        trl_corr = p_monk(s).area.(areas{ar}).signal(nlfp).trials(corr);
                        trl_incorr = p_monk(s).area.(areas{ar}).signal(nlfp).trials(incorr);
                        % extract t_stop corr
                        clear trl_end_corr trl_end_incorr
                        for j = 1:length(behv_corr), trl_end_corr(j) = behv_corr(j).events.t_end; end
                        for j = 1:length(behv_incorr), trl_end_incorr(j) = behv_incorr(j).events.t_end; end
                        
                        %% corr theta
                        % align to t_stop and extract around stop time
                        for trl = 1:length(behv_corr)
                            this_ts = behv_corr(trl).continuous.ts-trl_end_corr(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr(trl).lfp_theta(this_ts>-1 & this_ts<1); theta_corr(trl,:) = this_lfp_corr(1:333,:);  %% Choose band
                                % monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr(trl).lfp = real(theta_corr(trl,:));
                                this_ninety = prctile(real(theta_corr(trl,:)),95); theta_corr_indx_95(trl,:) = real(theta_corr(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr(trl).tspk = ts(trl,theta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.ts_95_corr] = Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr, [-1:0.02:1], 0.05);
                        %% corr beta
                        % align to t_stop and extract around stop time
                        for trl = 1:length(behv_corr)
                            this_ts = behv_corr(trl).continuous.ts-trl_end_corr(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_corr(trl,:) = this_lfp_corr(1:333,:);  %% Choose band
                                % monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr(trl).lfp = real(beta_corr(trl,:));
                                this_ninety = prctile(real(beta_corr(trl,:)),95); beta_corr_indx_95(trl,:) = real(beta_corr(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr(trl).tspk = ts(trl,beta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.ts_95_corr] = Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr, [-1:0.02:1], 0.05);
                        %% incorr theta
                        for trl = 1:length(behv_incorr)
                            this_ts = behv_incorr(trl).continuous.ts-trl_end_incorr(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_incorr = trl_incorr(trl).lfp_theta(this_ts>-1 & this_ts<1); theta_incorr(trl,:) = this_lfp_incorr(1:333,:);  %% Choose band
                                % monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr(trl).lfp = real(theta_incorr(trl,:));
                                this_ninety = prctile(real(theta_incorr(trl,:)),95); theta_incorr_indx_95(trl,:) = real(theta_incorr(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr(trl).tspk = ts(trl,theta_incorr_indx_95(trl,:))';
                            else
                                if trl ==1
                                    monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_corr(trl).tspk)));
                                else
                                    monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr(trl-1).tspk)));
                                end
                            end
                            monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr(trl).ts = ts;
                            % take a psth of the 95th prcnt
                            [monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.rate_95_incorr,  monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.ts_95_incorr] = Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).trl_incorr, [-1:0.02:1], 0.05);
                            %% incorr beta
                            for trl = 1:length(behv_incorr)
                                this_ts = behv_incorr(trl).continuous.ts-trl_end_incorr(trl);
                                if this_ts(end)>0.9
                                    this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                    this_lfp_incorr = trl_incorr(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_incorr(trl,:) = this_lfp_incorr(1:333,:);  %% Choose band
                                    % monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr(trl).lfp = real(beta_incorr(trl,:));
                                    this_ninety = prctile(real(beta_incorr(trl,:)),95); beta_incorr_indx_90(trl,:) = real(beta_incorr(trl,:))>this_ninety;
                                    monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr(trl).tspk = ts(trl,beta_incorr_indx_90(trl,:))';
                                else
                                    if trl ==1
                                        monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_corr(trl).tspk)));
                                    else
                                        monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr(trl-1).tspk)));
                                    end
                                end
                            end
                            monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr(trl).ts = ts;
                            %% take a psth of the 95th prcnt
                            [monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.rate_95_incorr,  monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.ts_95_incorr] = Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).trl_incorr, [-1:0.02:1], 0.05);
                        end
                    end
                    disp(['>>>>>       Monkey ' num2str(m) ' area ' (areas{ar}) ' channel ' num2str(nlfp)  ' done       <<<<<' ])
                    fprintf(['Time:  ' num2str(clock) '\n']);
                end
                % save per monkey per session to have smaller files and save time
                if ~isempty(nlfp)
                    dataToSave = monk(m).sess(s);
                    disp(['    Saving monkey ' num2str(m)])
                    fprintf(['Time:  ' num2str(clock) '\n']);
                    save(['monk ' num2str(m) ' band_passed area ' (areas{ar}) ' sess ' num2str(s) ], 'dataToSave', '-v7.3');
                    disp(['    Done! Clearing memory... .  .  ' (m)])
                    clear dataToSave
                    monk(m).sess(s) = [];
                end
            end
        end
        
    end
end

if do_band_passed_vs_accuracy
    disp('                 Do band passed vs accuracy . . .     ')
    fprintf(['Time:  ' num2str(clock) '\n']);
    monks = unique([exp.monk_id]);
    for m = 1:length(monks)
        ii = [exp.monk_id] == monks(m); p_monk = exp(ii);
        nsess = 1:length(p_monk);
        for s = 1:length(nsess)
            corr = p_monk(s).behavior.stats.trialtype.reward(strcmp({p_monk(s).behavior.stats.trialtype.reward.val},'rewarded')).trlindx;
            incorr = p_monk(s).behavior.stats.trialtype.reward(strcmp({p_monk(s).behavior.stats.trialtype.reward.val},'unrewarded')).trlindx;
            behv_corr = p_monk(s).behavior.trials(corr);
            % get monkey and target position only for correct trials
            monk_abs_x = p_monk(s).behavior.stats.pos_abs.x_monk(corr); ntrls_correct = length(monk_abs_x);
            monk_abs_y = p_monk(s).behavior.stats.pos_abs.y_monk(corr);
            r_targ = p_monk(s).behavior.stats.pos_final.r_targ(corr);
            theta_targ = p_monk(s).behavior.stats.pos_final.theta_targ(corr)*pi/180;
            x_targ = r_targ.*sin(theta_targ); % convert to cartesian
            y_targ = r_targ.*cos(theta_targ);
            % get end position for monkey
            for indx = 1:ntrls_correct, x_monk(indx) = monk_abs_x{indx}(end-130); y_monk(indx) = monk_abs_y{indx}(end-130) + 37.5; end
            % take out distance between monk position and target position in each trial
            clear monk2targ
            for indx = 1:ntrls_correct
                monk2targ(indx) = sqrt((x_monk(indx) - x_targ(indx))^2 + (y_monk(indx) - y_targ(indx))^2);
            end
            % select close to target and far away from target
            %             low_third_indx = monk2targ < 21.66;
            %             middle_third_indx = monk2targ > 21.66 & monk2targ < 43.32;
            %             upper_third_indx = monk2targ > 43.32;
            
            low_third_indx = monk2targ < 32.5;
            middle_third_indx = monk2targ > 32.5;
            upper_third_indx = monk2targ > 32.5;
            
            % plot histogram
            %             histogram(monk2targ,100);
            %             title('End point to Target')
            areas = fieldnames(p_monk(s).area); nareas = length(areas);
            for ar = 1:nareas
                for nlfp = 1:length(p_monk(s).area.(areas{ar}).signal)
                    if size(p_monk(s).area.(areas{ar}).signal,2) ~= 0
                        
                        trls_corr = p_monk(s).area.(areas{ar}).signal(nlfp).trials(corr);
                        trl_corr_low = trls_corr(low_third_indx); behv_low = behv_corr(low_third_indx);
                        trl_corr_middle = trls_corr(middle_third_indx); behv_middle = behv_corr(middle_third_indx);
                        trl_corr_upper = trls_corr(upper_third_indx); behv_upper = behv_corr(upper_third_indx);
                        % extract t_stop corr
                        clear trl_end_corr_low trl_end_corr_middle trl_end_corr_upper
                        for j = 1:length(trl_corr_low), trl_end_corr_low(j) = behv_low(j).events.t_end; end
                        for j = 1:length(trl_corr_middle), trl_end_corr_middle(j) = behv_middle(j).events.t_end; end
                        for j = 1:length(trl_corr_upper), trl_end_corr_upper(j) = behv_upper(j).events.t_end; end
                        
                        %% corr theta low
                        % align to t_stop and extract around stop time
                        clear theta_corr_indx_95 ts
                        for trl = 1:length(trl_corr_low)
                            this_ts = behv_low(trl).continuous.ts-trl_end_corr_low(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr_low(trl).lfp_theta(this_ts>-1 & this_ts<1); theta_corr_low(trl,:) = this_lfp_corr(1:333,:);  %%%%%%%%%%%%%%% Choose band
                                this_ninety = prctile(real(theta_corr_low(trl,:)),95); theta_corr_indx_95(trl,:) = real(theta_corr_low(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).low.trl_corr(trl).tspk = ts(trl,theta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).low.trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).low.trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).low.trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.low.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.low.ts_95_corr] = ...
                            Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).low.trl_corr, [-1:0.02:1], 0.05);
                        clear theta_corr_indx_95 ts
                        %% corr theta middle
                        % align to t_stop and extract around stop time
                        for trl = 1:length(trl_corr_middle)
                            this_ts = behv_middle(trl).continuous.ts-trl_end_corr_middle(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr_middle(trl).lfp_theta(this_ts>-1 & this_ts<1); theta_corr_middle(trl,:) = this_lfp_corr(1:333,:); %%%%%%%%%%%%%%% Choose band
                                this_ninety = prctile(real(theta_corr_middle(trl,:)),95); theta_corr_indx_95(trl,:) = real(theta_corr_middle(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).middle.trl_corr(trl).tspk = ts(trl,theta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).middle.trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).middle.trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).middle.trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.middle.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.middle.ts_95_corr] = ...
                            Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).middle.trl_corr, [-1:0.02:1], 0.05);
                        clear theta_corr_indx_95 ts
                        %% corr theta upper
                        % align to t_stop and extract around stop time
                        for trl = 1:length(trl_corr_upper)
                            this_ts = behv_upper(trl).continuous.ts-trl_end_corr_upper(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr_upper(trl).lfp_theta(this_ts>-1 & this_ts<1); theta_corr_upper(trl,:) = this_lfp_corr(1:333,:);  %%%%%%%%%%%%%%% Choose band
                                this_ninety = prctile(real(theta_corr_upper(trl,:)),95); theta_corr_indx_95(trl,:) = real(theta_corr_upper(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).upper.trl_corr(trl).tspk = ts(trl,theta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).upper.trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).upper.trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).upper.trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.upper.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).psth.upper.ts_95_corr] = ...
                            Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).upper.trl_corr, [-1:0.02:1], 0.05);
                        clear theta_corr_indx_95 ts
                        %% corr beta low
                        % align to t_stop and extract around stop time
                        for trl = 1:length(trl_corr_low)
                            this_ts = behv_low(trl).continuous.ts-trl_end_corr_low(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr_low(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_corr_low(trl,:) = this_lfp_corr(1:333,:);  %%%%%%%%%%%%%%% Choose band
                                this_ninety = prctile(real(beta_corr_low(trl,:)),95); beta_corr_indx_95(trl,:) = real(beta_corr_low(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).low.trl_corr(trl).tspk = ts(trl,beta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).low.trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).low.trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).low.trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.low.rate_95_corr_low, monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.low.ts_95_corr] = ...
                            Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.theta.lfp(nlfp).low.trl_corr, [-1:0.02:1], 0.05);
                        clear beta_corr_indx_95 ts
                        %% corr beta middle
                        % align to t_stop and extract around stop time
                        for trl = 1:length(trl_corr_middle)
                            this_ts = behv_middle(trl).continuous.ts-trl_end_corr_middle(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr_middle(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_corr_middle(trl,:) = this_lfp_corr(1:333,:); %%%%%%%%%%%%%%% Choose band
                                this_ninety = prctile(real(beta_corr_middle(trl,:)),95); beta_corr_indx_95(trl,:) = real(beta_corr_middle(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).middle.trl_corr(trl).tspk = ts(trl,beta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).middle.trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).middle.trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).middle.trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.middle.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.middle.ts_95_corr] = ...
                            Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).middle.trl_corr, [-1:0.02:1], 0.05);
                        clear beta_corr_indx_95 ts
                        %% corr beta  upper
                        % align to t_stop and extract around stop time
                        for trl = 1:length(trl_corr_upper)
                            this_ts = behv_upper(trl).continuous.ts-trl_end_corr_upper(trl);
                            if this_ts(end)>0.9
                                this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                                this_lfp_corr = trl_corr_upper(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_corr_upper(trl,:) = this_lfp_corr(1:333,:);  %%%%%%%%%%%%%%% Choose band
                                this_ninety = prctile(real(beta_corr_upper(trl,:)),95); beta_corr_indx_95(trl,:) = real(beta_corr_upper(trl,:))>this_ninety;
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).upper.trl_corr(trl).tspk = ts(trl,beta_corr_indx_95(trl,:))';
                            else
                                monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).upper.trl_corr(trl).tspk = single(NaN(size(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).upper.trl_corr(trl-1).tspk)));
                            end
                        end
                        monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).upper.trl_corr(trl).ts = ts;
                        % take a psth of the 95th prcnt
                        [monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.upper.rate_95_corr, monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).psth.upper.ts_95_corr] = ...
                            Spiketimes2Rate(monk(m).sess(s).area.(areas{ar}).band_passed.beta.lfp(nlfp).upper.trl_corr, [-1:0.02:1], 0.05);
                        
                        disp(['>>>>>       Monkey ' num2str(m) ' area ' (areas{ar}) ' channel ' num2str(nlfp)  ' done       <<<<<' ])
                        fprintf(['Time:  ' num2str(clock) '\n']);
                    end
                end
                % save per monkey per session to have smaller files and save time
                if ~isempty(nlfp)
                    dataToSave = monk(m).sess(s);
                    disp(['    Saving monkey ' num2str(m)])
                    % fprintf(['Time:  ' num2str(clock) '\n']);
                    save(['monk ' num2str(m) ' band_passedVSaccuracyHalf area ' (areas{ar}) ' sess ' num2str(s) ], 'dataToSave', '-v7.3');
                    disp(['    Done! Clearing memory... .  .  '])
                    clear dataToSave
                    monk(m).sess(s) = [];
                end
            end
        end
    end 
end

%% Save
    disp('                 Done, saving . . .     ')
    fprintf(['Time:  ' num2str(clock) '\n']);
    save('lfp_pop_sim_2020_08_20_coherence_align_stop', 'monk', 'all_monks', '-v7.3');
    load train
    sound(y,Fs)
    disp('                     Saved!')
    fprintf(['Time:  ' num2str(clock) '\n']);
