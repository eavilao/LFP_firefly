function LFPpop_sim_all(exp)

% Function to  load, extract, and save relevant variables from
% experiments.m
% If running for the first time, you don't need (exp) just hit run
%% Choose what to analyze and save
extract_exp_out = false; % load experiments.m file and extract. If saved 'exp_out', make it false.
save_exp_out = false; % save mat file without raw lfp signal
save_pop = false; % if this is true it will only extract pop

extract_lfp_raw = false; % raw and per trial lfps
save_lfp_raw = false; % raw and per trial lfps
do_PSD = true;  % extract power spectral densities
save_spectro = false;
save_spectro_per_trial = false;
save_spectro_per_trial_align_stop = false;
avg_monks = true; % average for all monkeys?
do_cohero = false; % extract coherograms
do_cohero_band_passed = false; % extract coherograms per band
doCSD = false; % Perform CSD analysis for MST recordings?
do_ERP = true; % extract ERPs (evoked LFPs)
save_band_pass_analysis = false; % extract band passed lfp signal only (used only for plotting)
do_band_passed_pop = false;  % needs pop
do_phases = false; % needs pop

name_output_exp_out_file = 'exp_out_psd_erp_all_2022_04_15';
name_output_file = 'lfp_psd_erp_all_2022_04_15';

%% Extract
if extract_exp_out
    %      path = 'D:\Output';
    %      cd(path)
    fprintf(['Time:  ' num2str(clock) '\n']); 
    fnames = dir('experiments*.mat');
    cnt=1;
    for i = 1:length(fnames)
        fprintf(['****   Loading file ' num2str(fnames(i).name) '   ****' '\n'])
        load(fnames(i).name);
        for sess = 1:length(experiments.sessions)
            if ~isempty(experiments.sessions(sess).lfps(1).trials)
                %% Read areas
                for nlfps = 1:length(experiments.sessions(sess).lfps)
                    indx_MST(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'MST');
                    indx_PFC(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'PFC');
                    indx_PPC(nlfps) = strcmp(experiments.sessions(sess).lfps(nlfps).brain_area, 'PPC');
                end
                %% extract per area
                if save_pop
                    exp(cnt).behavior = [experiments.sessions(sess).behaviours];
                    exp(cnt).pop = experiments.sessions(sess).populations.lfps.stats;
                    exp(cnt).monk_id = experiments.sessions(sess).monk_id;
                else
                    if save_band_pass_analysis
                        mst_ch = find(indx_MST); ppc_ch = find(indx_PPC); pfc_ch = find(indx_PFC);
                        if experiments.sessions(sess).monk_id == 44 || experiments.sessions(sess).monk_id == 53 ...
                                for ch = 1:sum(indx_MST), exp(cnt).area.MST.band_passed(ch) = experiments.sessions(sess).lfps(mst_ch(ch)).stats.band_passed;end
                        end % extract band passed signal
                        if experiments.sessions(sess).monk_id == 53 || experiments.sessions(sess).monk_id == 71 ...
                                for ch = 1:sum(indx_PFC), exp(cnt).area.PFC.band_passed(ch) = experiments.sessions(sess).lfps(pfc_ch(ch)).stats.band_passed;end
                        end
                        for ch = 1:sum(indx_PPC), exp(cnt).area.PPC.band_passed(ch) = experiments.sessions(sess).lfps(ppc_ch(ch)).stats.band_passed;end
                    else
                        
                        mst_ch = find(indx_MST); ppc_ch = find(indx_PPC); pfc_ch = find(indx_PFC);
                        if experiments.sessions(sess).monk_id == 44 || experiments.sessions(sess).monk_id == 53 ...
                                for ch = 1:sum(indx_MST), exp(cnt).area.MST.lfps(ch).stats.trialtype.all = experiments.sessions(sess).lfps(mst_ch(ch)).stats.trialtype.all;end
                        end % extract band passed signal
                        if experiments.sessions(sess).monk_id == 53 || experiments.sessions(sess).monk_id == 71 ...
                                for ch = 1:sum(indx_PFC), exp(cnt).area.PFC.lfps(ch).stats.trialtype.all = experiments.sessions(sess).lfps(pfc_ch(ch)).stats.trialtype.all;end
                        end
                        for ch = 1:sum(indx_PPC), exp(cnt).area.PPC.lfps(ch).stats.trialtype.all = experiments.sessions(sess).lfps(ppc_ch(ch)).stats.trialtype.all;end
                        
                        %
                        %                         if experiments.sessions(sess).monk_id == 44 || experiments.sessions(sess).monk_id == 53 ...
                        %                                 exp(cnt).area.MST.lfps.stats = [experiments.sessions(sess).lfps(indx_MST).stats];
                        %                         end
                        %                         if experiments.sessions(sess).monk_id == 53 || experiments.sessions(sess).monk_id == 71 ...
                        %                                 exp(cnt).area.PFC.lfps.stats = [experiments.sessions(sess).lfps(indx_PFC).stats];
                        %                         end
                        %                         exp(cnt).area.PPC.lfps.stats = [experiments.sessions(sess).lfps(indx_PPC).stats];
                    end
                    
                    exp(cnt).behavior = [experiments.sessions(sess).behaviours];
                    exp(cnt).monk_id = experiments.sessions(sess).monk_id;
                end
            end
            cnt=cnt+1;
        end
        disp('Clearing experiments... . . .')
        clear experiments
    end
    
    if save_exp_out
        % save
        disp('Saving... . . .')
        fprintf(['Time:  ' num2str(clock) '\n']);
        save(name_output_exp_out_file,'exp', '-v7.3');
        % load train
        % sound(y,Fs)
        disp('           Saved! ')
        fprintf(['Time:  ' num2str(clock) '\n']);
        disp('Extracting... . . .')
    end
end

monk = []; all_monks=[];
%% Save behavior

monks = unique([exp.monk_id]);
for ii = 1:length(monks)
    m = [exp.monk_id] == monks(ii); p_monk = exp(m);
    for j = 1:length(p_monk)
        monk(ii).behavior(j) = p_monk(j).behavior;
    end
end

%% PSD
%% avg per session
if do_PSD
    trial_type = fieldnames(exp(1).area.PPC.lfps(1).stats.trialtype); events = fieldnames(exp(1).area.PPC.lfps(1).stats.trialtype.all.events);
    for i = 1:length(exp) % num of sessions
        nareas = numel(fieldnames(exp(i).area)); areas = fieldnames(exp(i).area);
        for a = 1:length(areas) % num of areas
            for type = 1:length(trial_type) % num of trial types
                if ~isempty(exp(i).area.(areas{a}).lfps)
                    nconds = length(exp(i).area.(areas{a}).lfps(1).stats.trialtype.(trial_type{type})); clear cond
                    for cond = 1:nconds
                        clear pow pow_eye
                        for ch = 1:length(exp(i).area.(areas{a}).lfps) % extract per channel
                            if strcmp((trial_type{type}), 'eyesfree') | strcmp((trial_type{type}), 'eyesfixed')
                                pow_eye(ch,:) = exp(i).area.(areas{a}).lfps(ch).stats.trialtype.(trial_type{type})(cond).spectrum.psd;
                                freq_eye = exp(i).area.(areas{a}).lfps(ch).stats.trialtype.(trial_type{type})(cond).spectrum.freq;
                            else
                                pow(ch,:) = exp(i).area.(areas{a}).lfps(ch).stats.trialtype.(trial_type{type})(cond).spectrum.psd;
                            end
                        end
                        % store
                        freq = exp(i).area.(areas{a}).lfps(1).stats.trialtype.(trial_type{type})(cond).spectrum.freq;
                        if strcmp((trial_type{type}), 'eyesfree') | strcmp((trial_type{type}), 'eyesfixed')
                            psden_eye(i).area.(areas{a}).(trial_type{type})(cond).ch = pow_eye;
                            psden_eye(i).area.(areas{a}).(trial_type{type})(cond).mu = nanmean(pow_eye);
                            psden_eye(i).area.(areas{a}).(trial_type{type})(cond).sem = nanstd(pow_eye)/sqrt(length(ch));
                            psden_eye(i).area.(areas{a}).(trial_type{type})(cond).max = max(max(nanmean(pow_eye)));
                        else
                            psden(i).area.(areas{a}).(trial_type{type})(cond).ch = pow;
                            psden(i).area.(areas{a}).(trial_type{type})(cond).mu = nanmean(pow);
                            psden(i).area.(areas{a}).(trial_type{type})(cond).sem = nanstd(pow)/sqrt(length(ch));
                            psden(i).area.(areas{a}).(trial_type{type})(cond).max = max(max(nanmean(pow)));
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
        m = [psden.monk_id] == monks(i); p_monk = psden(m); % p_monk_eye = psden_eye(m);
        nareas = numel(fieldnames(p_monk(1).area)); areas = fieldnames(p_monk(1).area);
        for a = 1:length(areas)
            % trialtype = [fieldnames(p_monk(1).area.(areas{a})); fieldnames(p_monk_eye(1).area.(areas{a}))];
            trialtype = fieldnames(p_monk(1).area.(areas{a}));
            for type = 1:length(trialtype)
                if type == 6 | type == 7  % eyes
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
        monk(i).pw.freq = exp(i).area.(areas{a}).lfps(1).stats.trialtype.all.spectrum.freq;
        % monk(i).pw.freq_eye = exp(i).area.(areas{a}).lfps.stats(1).trialtype.eyesfree.spectrum.freq;
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
        %all_monks.pw.freq_eye = monk(1).pw.freq_eye;
    end
end
%% spectrogram
% average across channels
if save_spectro
    trialtype = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype); events = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype.all.events);
    for i = 1:length(exp) % num of sessions
        nareas = numel(fieldnames(exp(i).area)); areas = fieldnames(exp(i).area);
        for a = 1:length(areas) % num of areas
            for type = 2% 1:3 % length(trialtype) % num of trial types
                if ~isempty(exp(i).area.(areas{a}).lfps.stats)
                    nconds = length(exp(i).area.(areas{a}).lfps.stats(1).trialtype.(trialtype{type})); clear cond
                    for cond = 1:nconds
                        for ev = 1:length(events)
                            clear spec
                            if ev == 4 && cond == 1
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu = NaN;
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu_no_norm = NaN;
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).std = NaN;
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq = NaN;
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts = NaN;
                            else
                                for ch = 1:length(exp(i).area.(areas{a}).lfps.stats) % extract per channel
                                    spec(ch,:,:) = real(exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).all_freq.p_spectrogram');
                                end
                                
                                % store and normalize by max
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu = squeeze(nanmean(spec,1))/max(max(squeeze(nanmean(spec,1))));
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).pow_mu_no_norm = squeeze(nanmean(spec,1));
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).std = squeeze(nanstd(spec,1)/sqrt(length(ch)));
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).all_freq.freq_spectrogram;
                                p_spec(i).area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts = exp(i).area.(areas{a}).lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).all_freq.ts_spectrogram;
                            end
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
            for type = 1 %1:length(trialtype)
                nconds = 2 %length(p_monk(1).area.(areas{a}).(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ev = [2 3 4] %1:length(events)
                        clear spec spec_no_norm spec_f spec_ts
                        for k = 1:length(p_monk)
                            ts = round(p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).ts);
                            spec(k,:,:) = p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).pow_mu(:,ts>=-1.56 & ts<=1.56);
                            spec_no_norm(k,:,:) = p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).pow_mu_no_norm(:,ts>-1.51 & ts<1.51);
                            %spec_f(k,:,:) = p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).freq;
                            %spec_ts(k,:,:) = p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).ts(ts>-1.51 & ts<1.51);
                        end
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).mu_sess = squeeze(nanmean(spec,1));  % average per trial type
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).mu_sess_no_norm = squeeze(nanmean(spec_no_norm,1));  % average per trial type
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).std_sess = squeeze(nanstd(spec,1));
                        % monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq_sess = squeeze(spec_f(1,:));
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).freq_sess = p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).freq;
                        % monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts_sess = squeeze(spec_ts(1,:));
                        monk(i).spec.area.(areas{a}).(trialtype{type})(cond).events.(events{ev}).ts_sess = p_monk(k).area.(areas{a}).reward(cond).events.(events{ev}).ts(ts>-1.51 & ts<1.51)
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
                for type = 2 %1:length(trialtype)
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
    %     monks = unique([exp.monk_id]); clear p_monk
    %     for i = 1:length(monks)
    %         m = [exp.monk_id] == monks(i); p_monk = exp(m);
    %         for sess = 1:length(p_monk) % num of sessions
    %             trialtype = fieldnames(p_monk(sess).pop.trialtype);
    %             for type = 1:length(trialtype) % num of trial typesM
    %                 for cond = 2%1:length(p_monk(sess).pop.trialtype.(trialtype{type}));
    %                     areas = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area);
    %                     for a=1:length(areas)
    %                         events = fieldnames(exp(1).pop.trialtype.reward(2).area.PPC.events);
    %                         for ev = 1:length(events)
    %                             for chn = 1:length(p_monk(sess).pop.trialtype.reward(cond).area.(areas{a}).events.(events{ev}).ch)
    %                                 monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).ch(chn).trl = p_monk(sess).pop.trialtype.reward(cond).area.(areas{a}).events.(events{ev}).ch(chn).trl;
    %                             end
    %                         end
    %                     end
    %                 end
    %             end
    %         end
    %     end
    
    %% Gather spectrograms for all electrodes each trial per monkey aligned to event onset
    monks = unique([exp.monk_id]); clear p_monk
    for i = 1:length(monks)
        m = [exp.monk_id] == monks(i); p_monk = exp(m);
        for sess = 1:length(p_monk) % num of sessions
            trialtype = fieldnames(p_monk(sess).pop.trialtype);
            for type = 1:length(trialtype) % num of trial types
                for cond = 1:length(p_monk(sess).pop.trialtype.(trialtype{type}))
                    areas = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area);
                    for a=1:length(areas)
                        events = fieldnames(exp(1).pop.trialtype.reward(2).area.PPC.events);
                        for ev = 1:length(events)
                            if strcmp((events{ev}),'reward') && cond==1 || strcmp((events{ev}),'reward') && cond==3 || strcmp((events{ev}),'reward') && cond==4
                                monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl.ts = NaN;
                                monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl.freq = NaN;
                                monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl.spectro = NaN;
                            else
                                for trl = 1:length(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pop_trl)
                                    % Z normalization and store
                                    monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(trl).ts = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pop_trl(trl).ts_spectrogram;
                                    monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(trl).freq = p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pop_trl(trl).freq_spectrogram;
                                    p_spectro = (p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pop_trl(trl).spectrogram - nanmean(nanmean(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pop_trl(trl).spectrogram)))./...
                                        nanstd(nanstd(p_monk(sess).pop.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pop_trl(trl).spectrogram));
                                    monk(i).sess(sess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(trl).spectro = p_spectro;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% average across monkeys
    if avg_monks
        areas = {'PPC', 'PFC', 'MST'};
        for a = 1:length(areas)
            for type = 1:length(trialtype)
                clear cond
                for cond = 1:length(monk(i).sess(sess).trialtype.(trialtype{type}))
                    for ev = 2:4 %1:length(events)
                        clear spectro_trl_freq monk_sess_mu monk_mu
                        if strcmp((areas{a}),'PFC')
                            for i = [3 4] % change if adding one more monkey
                                for nsess = 1:length(monk(i).sess)  
                                    for trl = 1:length(monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl)
                                        ts = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).ts;
                                        for freq = 1:length(monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq)
                                            pw_trl =  monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(trl).spectro(ts>-1.51 & ts<1.514,freq)'; %extract 1.5s around event
                                            %if ev == 3; spectro_trl_freq(trl,freq,:) = pw_trl(1:62);end
                                            spectro_trl_freq(trl,freq,:) = pw_trl;
                                        end
                                    end
                                    
                                    % average session per monkey
                                    monk_sess_mu(nsess,:,:) = squeeze(nanmean(spectro_trl_freq));
                                end
                                monk_mu(i,:,:) = squeeze(nanmean(monk_sess_mu)); % store avg per monkey
                                %% sanity plot per monkey
%                                 if ~isnan(monk_mu(:,1,1))
%                                      J = customcolormap_preset('black_teal_white');
%                                     f = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq;
%                                     figure('Name',['monk ' num2str(i) ]'); hold on; colormap(J);
%                                     ts_plot = ts(ts>-1.51 & ts<1.514);
%                                     if ev == 2, ts_plot = ts_plot-0.3; x_lim = [-0.5 0.5]; elseif ev == 3,  x_lim = [-1.5 1.5]; else x_lim = [-1 1]; end
%                                     imagesc(ts_plot,f,squeeze(monk_mu(i,:,:)), [-0.2 2.4]); axis xy; colorbar;
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'ylim',[4 50], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('frequency (Hz)'); axis square
%                                     title(['cond ' num2str(cond)]); if ev == 2, vline([-0.3 0], 'w'); else vline(0, 'w'); end
%                                     % theta
%                                     if cond == 1, colorline = [0 0 0]; else colorline = [0 1 0]; end
%                                     figure; hold on
%                                     shadedErrorBar(ts_plot,nanmean(squeeze(monk_mu(i,f>=3.5 & f<12,:))), ...
%                                         nanstd(squeeze(monk_mu(i,f>=3.5 & f<12,:)))/sqrt(size(nanmean(squeeze(monk_mu(i,f>=3.5 & f<12,:))),1)),'lineprops',{'Color', colorline})
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('Normalized power'); axis square
%                                     title('PPC theta'); if ev == 2, vline([-0.3 0], 'k'); else vline(0, '--k'); end
%                                     % beta
%                                     figure; hold on
%                                     %plot(ts_s,p_spectro_corr(freq>=11 & freq<=21,:)','Color',[0 0.5 0]); alpha(0.1)
%                                     shadedErrorBar(ts_plot,nanmean(squeeze(monk_mu(i,f>=11 & f<21,:))), ...
%                                         nanstd(squeeze(monk_mu(i,f>=11 & f<21,:)))/sqrt(size(nanmean(squeeze(monk_mu(i,f>=11 & f<21,:))),1)),'lineprops',{'Color', colorline})
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('Normalized power'); axis square
%                                     title('PPC beta'); if ev == 2, vline([-0.3 0], 'k'); else vline(0, '--k'); end
%                                 end
                            end
                            % average for all monkeys
                            monk_sess_mu(1,:,:) = []; monk_sess_mu(2,:,:) = []; % hardcoded because the first and second monks have no MST rec
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_mu = squeeze(nanmean(monk_sess_mu));
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_ts = ts(ts>-1.51 & ts<1.514);
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_freq = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq;
                        elseif strcmp((areas{a}),'MST')
                            for i = [1 3] % change if adding one more monkey
                                for nsess = 1:length(monk(i).sess)
                                    ts = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).ts;
                                    for trl = 1:length(monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl)
                                        for freq = 1:length(monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq)
                                            pw_trl =  monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(trl).spectro(ts>-1.51 & ts<1.514,freq)'; %extract 1.5s around event
                                            %if ev == 3; spectro_trl_freq(trl,freq,:) = pw_trl(1:62);end
                                            spectro_trl_freq(trl,freq,:) = pw_trl;
                                        end
                                    end
                                    % average session per monkey
                                    monk_sess_mu(nsess,:,:) = squeeze(nanmean(spectro_trl_freq));
                                end
                                monk_mu(i,:,:) = squeeze(nanmean(monk_sess_mu)); % store avg per monkey
                                %% sanity plot per monkey
%                                 if ~isnan(monk_mu(:,1,1))
%                                     J = customcolormap_preset('black_teal_white');
%                                     f = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq;
%                                     figure('Name',['monk ' num2str(i) ]'); hold on; colormap(J);
%                                     ts_plot = ts(ts>-1.51 & ts<1.514);
%                                     if ev == 2, ts_plot = ts_plot-0.3; x_lim = [-0.5 0.5]; elseif ev == 3,  x_lim = [-1.5 1.5]; else x_lim = [-1 1]; end
%                                     imagesc(ts_plot,f,squeeze(monk_mu(i,:,:)), [-0.2 2.4]); axis xy; colorbar;
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'ylim',[4 50], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('frequency (Hz)'); axis square
%                                     title(['cond ' num2str(cond)]); if ev == 2, vline([-0.3 0], 'w'); else vline(0, 'w'); end
%                                     % theta
%                                     if cond == 1, colorline = [0 0 0]; else colorline = [0 1 0]; end
%                                     figure; hold on
%                                     shadedErrorBar(ts_plot,nanmean(squeeze(monk_mu(i,f>=3.5 & f<12,:))), ...
%                                         nanstd(squeeze(monk_mu(i,f>=3.5 & f<12,:)))/sqrt(size(nanmean(squeeze(monk_mu(i,f>=3.5 & f<12,:))),1)),'lineprops',{'Color', colorline})
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('Normalized power'); axis square
%                                     title('PPC theta'); if ev == 2, vline([-0.3 0], 'k'); else vline(0, '--k'); end
%                                     % beta
%                                     figure; hold on
%                                     %plot(ts_s,p_spectro_corr(freq>=11 & freq<=21,:)','Color',[0 0.5 0]); alpha(0.1)
%                                     shadedErrorBar(ts_plot,nanmean(squeeze(monk_mu(i,f>=11 & f<21,:))), ...
%                                         nanstd(squeeze(monk_mu(i,f>=11 & f<21,:)))/sqrt(size(nanmean(squeeze(monk_mu(i,f>=11 & f<21,:))),1)),'lineprops',{'Color', colorline})
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('Normalized power'); axis square
%                                     title('PPC beta'); if ev == 2, vline([-0.3 0], 'k'); else vline(0, '--k'); end
%                                 end
                            end
                            % average for all monkeys
                            monk_sess_mu(2,:,:) = []; % hardcoded because the second monk has no MST rec
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_mu = squeeze(nanmean(monk_sess_mu));
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_ts = ts(ts>-1.51 & ts<1.514);
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_freq = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq;
                        else
                            for i = 1:length(monk) % change if adding one more monkey
                                for nsess = 1:length(monk(i).sess)
                                    ts = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).ts;
                                    for trl = 1:length(monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl)
                                        for freq = 1:length(monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq)
                                            pw_trl =  monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(trl).spectro(ts>-1.51 & ts<1.514,freq)'; %extract 1.5s around event
                                            %if ev == 3; spectro_trl_freq(trl,freq,:) = pw_trl(1:62);end
                                            spectro_trl_freq(trl,freq,:) = pw_trl;
                                        end
                                    end
                                    % average session per monkey
                                    monk_sess_mu(nsess,:,:) = squeeze(nanmean(spectro_trl_freq));
                                end
                                monk_mu(i,:,:) = squeeze(nanmean(monk_sess_mu)); % store avg per monkey
                                %% sanity plot per monkey
%                                 if ~isnan(monk_mu(:,1,1))
%                                    J = customcolormap_preset('black_teal_white');
%                                     f = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq;
%                                     figure('Name',['monk ' num2str(i) ]'); hold on; colormap(J);
%                                     ts_plot = ts(ts>-1.51 & ts<1.514);
%                                     if ev == 2, ts_plot = ts_plot-0.3; x_lim = [-0.5 0.5]; elseif ev == 3,  x_lim = [-1.5 1.5]; else x_lim = [-1 1]; end
%                                     imagesc(ts_plot,f,squeeze(monk_mu(i,:,:)), [-0.2 2.4]); axis xy; colorbar;
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'ylim',[4 50], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('frequency (Hz)'); axis square
%                                     title(['cond ' num2str(cond)]); if ev == 2, vline([-0.3 0], 'w'); else vline(0, 'w'); end
%                                     % theta
%                                     if cond == 1, colorline = [0 0 0]; else colorline = [0 1 0]; end
%                                     figure; hold on
%                                     shadedErrorBar(ts_plot,nanmean(squeeze(monk_mu(i,f>=3.5 & f<12,:))), ...
%                                         nanstd(squeeze(monk_mu(i,f>=3.5 & f<12,:)))/sqrt(size(nanmean(squeeze(monk_mu(i,f>=3.5 & f<12,:))),1)),'lineprops',{'Color', colorline})
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('Normalized power'); axis square
%                                     title('PPC theta'); if ev == 2, vline([-0.3 0], 'k'); else vline(0, '--k'); end
%                                     % beta
%                                     figure; hold on
%                                     %plot(ts_s,p_spectro_corr(freq>=11 & freq<=21,:)','Color',[0 0.5 0]); alpha(0.1)
%                                     shadedErrorBar(ts_plot,nanmean(squeeze(monk_mu(i,f>=11 & f<21,:))), ...
%                                         nanstd(squeeze(monk_mu(i,f>=11 & f<21,:)))/sqrt(size(nanmean(squeeze(monk_mu(i,f>=11 & f<21,:))),1)),'lineprops',{'Color', colorline})
%                                     set(gca,'xlim',[x_lim(1) x_lim(2)], 'FontSize', 22)
%                                     xlabel([(events{ev}) ' time (s)']); ylabel('Normalized power'); axis square
%                                     title('PPC beta'); if ev == 2, vline([-0.3 0], 'k'); else vline(0, '--k'); end
%                                 end
                            end
                            % average for all monkeys
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_mu = squeeze(nanmean(monk_sess_mu));
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_ts = ts;
                            all_monks.trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).all_monks_freq = monk(i).sess(nsess).trialtype.(trialtype{type})(cond).area.(areas{a}).events.(events{ev}).pw_trl(1).freq;
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
        coh_areas = fieldnames(p_monk(i).pop.trialtype.reward(2).events.stop.coherogram);
        trialtype = fieldnames(p_monk(i).pop.trialtype);
        for sess = 1:length(p_monk)
            for type = 2 %1:length(trialtype)
                nconds = length(p_monk(sess).pop.trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    events = fieldnames(p_monk(sess).pop.trialtype.(trialtype{type})(cond).events);
                    for ev = 1:length(events)
                        clear coh_ar coh_ev coh_phi coh_ts coh_freq
                        for coh_ar = 1:length(coh_areas)
                            coh_ev(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).coherogram.(coh_areas{coh_ar}).coher;
                            coh_phi(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).coherogram.(coh_areas{coh_ar}).coherPhi;
                            coh_ts(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).coherogram.(coh_areas{coh_ar}).coher_ts;
                            coh_freq(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).coherogram.(coh_areas{coh_ar}).coher_freq;
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

if do_cohero_band_passed
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
                        bands = fieldnames(p_monk(1).pop.trialtype.reward(1).events.stop.(coh_areas{1})); % per band
                        for b = 1:length(bands)
                            clear coh_ar coh_ev coh_phi coh_ts coh_freq
                            for coh_ar = 1:length(coh_areas)
                                coh_ev(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher;
                                coh_phi(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coherPhi;
                                coh_ts(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher_ts;
                                coh_freq(coh_ar,:,:) = p_monk(sess).pop.trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher_freq;
                            end
                            % plot -- sanity
                            %             figure('Name', 'PPC-->MST');
                            %             imagesc(coh_ts(2,:)-1, coh_freq(2,:),squeeze(coh_ev(2,:,:)));
                            %             axis xy; set(gca,'xlim',[-0.25 0.25], 'ylim', [4 50],'FontSize', 22); colorbar;
                            %
                            
                            
                            % coherogram mean
                            for coh_ar = 1:length(coh_areas)
                                monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher = squeeze(coh_ev(coh_ar,:,:));
                                monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher_phi = squeeze(coh_phi(coh_ar,:,:));
                                monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher_ts = squeeze(coh_ts(coh_ar,:,:));
                                monk(i).coher.sess(sess).trialtype.(trialtype{type})(cond).events.(events{ev}).(coh_areas{coh_ar}).(bands{b}).coher_freq = squeeze(coh_freq(coh_ar,:,:));
                            end
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
    for ii = [1 3] %1:length(monks)
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
    for ii = [1 3]% 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            for type = 1:length(trialtype)
                nconds = length(p_monk(j).area.PPC.lfps(1).stats.trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    events = fieldnames(p_monk(j).area.MST.lfps(1).stats.trialtype.all.events);
                    for ch = 1:length(p_monk(j).area.MST.lfps)
                        
                        % cont
                        %                     monk(ii).cont.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous; %gather continuous
                        for ev = 1:length(events)
                            % ERPs
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.MST.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.MST.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.MST.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            % monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).lfp_align = p_monk(j).area.MST.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).all_freq.lfp_align;
                            if extract_lfp_raw
                                monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.MST.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                % monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.MST.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.MST.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).beta;
                            end
                        end
                    end
                    for ch = 1:length(p_monk(j).area.PPC.lfps)
                        % cont
                        %                     monk(ii).cont.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous;
                        for ev = 1:length(events)
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            if extract_lfp_raw
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                % monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).beta;
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
    for ii = [2 4]  %Bruno Vik 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            for type = 1:length(trialtype)
                nconds = length(p_monk(j).area.PPC.lfps(1).stats.trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ch = 1:length(p_monk(j).area.PPC.lfps)
                        % cont
                        %                     monk(ii).cont.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous;
                        for ev = 1:length(events)
                            % ERPs
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            if extract_lfp_raw
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).beta;
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
    for ii = [3 4] %3 Schro Vik  1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for j = 1:length(p_monk)
            for type = 1:length(trialtype)
                nconds = length(p_monk(j).area.PFC.lfps(1).stats.trialtype.(trialtype{type})); clear cond
                for cond = 1:nconds
                    for ch = 1:length(p_monk(j).area.PFC.lfps)
                        % cont
                        %                     monk(ii).cont.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).continuous = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype{type})(cond).continuous;
                        for ev = 1:length(events)
                            % ERPs
                            monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_mu = p_monk(j).area.PFC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_mu;
                            monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_sem = p_monk(j).area.PFC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).potential_sem;
                            monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).erp_time = p_monk(j).area.PFC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).time;
                            if extract_lfp_raw
                                monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).theta = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).theta;
                                monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).alpha = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).alpha;
                                monk(ii).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype{type})(cond).events.(events{ev}).beta = p_monk(j).area.PPC.lfps(ch).stats.trialtype.(trialtype{type})(cond).events.(events{ev}).beta;
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

%% Analyses on band passed trials

if save_band_pass_analysis
    monks = unique([exp.monk_id]); clear p_monk
    for ii = 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for s = 1:length(p_monk)
            areas = fieldnames(p_monk(s).area);
            for a = 1:length(areas)
                for trl = 1:length(p_monk)
                    monk(ii).sess(s).area.(areas{a}).lfp = p_monk(s).area.(areas{a}).band_passed;
                end
            end
        end
    end
end

%% only extract data for raster
if do_band_passed_pop
    monks = unique([exp.monk_id]); clear p_monk
    for ii = 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for s = 1:length(p_monk)
            areas = fieldnames(p_monk(s).pop.area);
            for a = 1:length(areas)
                for trl = 1:length(p_monk)
                    monk(ii).sess(s).pop = p_monk(s).pop;
                end
            end
        end
    end
end

if do_phases
    monks = unique([exp.monk_id]); clear p_monk
    for ii = 1:length(monks)
        m = [exp.monk_id] == monks(ii); p_monk = exp(m);
        for s = 1:length(p_monk)
            monk(ii).sess(s).pop = p_monk(s).pop;
        end
    end
end

%% Save
% This file is needed to plot in plotLFPpop_sim.m

%%
disp('                 Done, now saving . . .     ')
fprintf(['Time:  ' num2str(clock) '\n']);
save(name_output_file, 'monk','all_monks', '-v7.3')
% load train
% sound(y,Fs)
disp('                     Saved!')
fprintf(['Time:  ' num2str(clock) '\n']);
