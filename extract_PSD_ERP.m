function extract_PSD_ERP
% Function to  load, extract, and save PSD and ERP from
% experiments.m
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
            
            exp(sess).behavior = [experiments.sessions(sess).behaviours];
            exp(sess).monk_id = experiments.sessions(sess).monk_id;
            exp(sess).area.MST.lfps.stats = [experiments.sessions(sess).lfps(indx_MST).stats];
            exp(sess).area.PFC.lfps.stats = [experiments.sessions(sess).lfps(indx_PFC).stats];
            exp(sess).area.PPC.lfps.stats = [experiments.sessions(sess).lfps(indx_PPC).stats];
            
            % extract PSD
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
                                monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = nanmean(psd_eye_all); % average per trial type
                                monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = nanstd(psd_eye_all);
                            else
                                if size(psd_all,1)<2
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = psd_all;
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = psd_all;
                                else
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = nanmean(psd_all);  % average per trial type
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = nanstd(psd_all);
                                end
                            end
                        end
                    end
                end
                monk(cnt).pw.freq = exp(i).area.(areas{a}).lfps.stats(1).trialtype.all.spectrum.freq;
                % monk(i).pw.freq_eye = exp(i).area.(areas{a}).lfps.stats(1).trialtype.eyesfree.spectrum.freq;
                monk(cnt).pw.monk_id = unique([psden(m).monk_id]);
            end
            
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
                                monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = nanmean(psd_eye_all); % average per trial type
                                monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = nanstd(psd_eye_all);
                            else
                                if size(psd_all,1)<2
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = psd_all;
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = psd_all;
                                else
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).mu_sess = nanmean(psd_all);  % average per trial type
                                    monk(cnt).pw.area.(areas{a}).(trialtype{type})(cond).std_sess = nanstd(psd_all);
                                end
                            end
                        end
                    end
                end
                monk(cnt).pw.freq = exp(i).area.(areas{a}).lfps.stats(1).trialtype.all.spectrum.freq;
                % monk(i).pw.freq_eye = exp(i).area.(areas{a}).lfps.stats(1).trialtype.eyesfree.spectrum.freq;
                monk(cnt).pw.monk_id = unique([psden(m).monk_id]);
            end
            
            %% Gather pop analyses for monks with MST rec
            monks = unique([exp.monk_id]);
            for ii = 1:length(monks)
                m = [exp.monk_id] == monks(ii); p_monk = exp(m);
                % trialtype = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype); events = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype.all.events);
                trialtype = 'reward'; events = fieldnames(exp(1).area.PPC.lfps.stats(1).trialtype.all.events);
                for j = 1:length(p_monk)
                    for type = 1:length(trialtype)
                        nconds = length(p_monk(j).area.PPC.lfps.stats(1).trialtype.reward); clear cond
                        for cond = 1:nconds
                            %% MST
                            if ~isempty(p_monk(j).area.MST.lfps.stats)
                                for ch = 1:length(p_monk(j).area.MST.lfps.stats)
                                    for ev = 1:length(events)
                                        % ERPs
                                        monk(cnt).erp.MST.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_mu = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).potential_mu;
                                        monk(cnt).erp.MST.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_sem = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).potential_sem;
                                        monk(cnt).erp.MST.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_time = p_monk(j).area.MST.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).time;
                                    end
                                end
                            end
                            %% PPC
                            if ~isempty(p_monk(j).area.PPC.lfps.stats)
                                for ch = 1:length(p_monk(j).area.PPC.lfps.stats)
                                    for ev = 1:length(events)
                                        monk(cnt).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_mu = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).potential_mu;
                                        monk(cnt).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_sem = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).potential_sem;
                                        monk(cnt).erp.PPC.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_time = p_monk(j).area.PPC.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).time;
                                    end
                                end
                            end
                            %% PFC 
                            if ~isempty(p_monk(j).area.PFC.lfps.stats)
                                for ch = 1:length(p_monk(j).area.PFC.lfps.stats)
                                    for ev = 1:length(events)
                                        monk(cnt).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_mu = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).potential_mu;
                                        monk(cnt).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_sem = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).potential_sem;
                                        monk(cnt).erp.PFC.sess(j).lfps(ch).trialtype.(trialtype)(cond).events.(events{ev}).erp_time = p_monk(j).area.PFC.lfps.stats(ch).trialtype.(trialtype)(cond).events.(events{ev}).time;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    disp('Saving individual monk files... . . .')
    save(['monk_' num2str(cnt)], 'monk','-v7.3')
    disp('Clearing experiments... . . .')
    clearvars -except fnames monk cnt
    cnt = cnt+1;
end
disp('               All Done, now saving . . .     ')
fprintf(['Time:  ' num2str(clock) '\n']);
save('PSD_ERP', 'monk', '-v7.3')
% load train
% sound(y,Fs)
disp('                     Saved!')
fprintf(['Time:  ' num2str(clock) '\n']);
end