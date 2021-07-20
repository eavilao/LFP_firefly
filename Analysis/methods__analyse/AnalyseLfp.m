function stats = AnalyseLfp(trials_lfps,stationary_lfps,mobile_lfps,eyesfixed_lfps,eyesfree_lfps,trials_behv,behv_stats,prs)

stats = [];
%% load analysis params
dt = prs.dt; % sampling resolution (s)
temporal_binwidth = prs.temporal_binwidth;
corr_lag = prs.corr_lag;
duration_zeropad = prs.duration_zeropad;
nbootstraps = prs.nbootstraps;
peaktimewindow = prs.peaktimewindow;
minpeakprominence = prs.minpeakprominence.neural;
mintrialsforstats = prs.mintrialsforstats;
event_potential = prs.event_potential;
ntrls = length(trials_lfps);
fixateduration = prs.fixateduration;
around_event = prs.around_event;
compute_psd = prs.compute_psd;
analyse_theta = prs.analyse_theta;
analyse_alpha = prs.analyse_alpha;
analyse_beta = prs.analyse_beta;
extract_band_pass = prs.extract_band_passed;
analyse_phase = prs.analyse_phase;
analyse_phase_within_area;


%% load cases
trialtypes = fields(behv_stats.trialtype);
events = cell2mat({trials_behv.events});
continuous = cell2mat({trials_behv.continuous});

spectralparams.tapers = prs.spectrum_tapers;
spectralparams.trialave = prs.spectrum_trialave;
spectralparams.Fs = 1/dt;
%% event-aligned, trial-averaged LFP
if event_potential
    gettuning = prs.tuning_events;
    for i=1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).events = stats.trialtype.all.events;
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                trials_lfps_temp = trials_lfps(trlindx);
                
                %% aligned to movement onset
                if any(strcmp(gettuning,'move'))
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_move], 'lfp');
                    lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.move)';
                    stats.trialtype.(trialtypes{i})(j).events.move.potential_mu = nanmean(lfps_temp2);
                    stats.trialtype.(trialtypes{i})(j).events.move.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                    stats.trialtype.(trialtypes{i})(j).events.move.time = prs.ts.move;
                    
                    if analyse_theta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_move], 'lfp_theta');
                        lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.move)';
                        stats.trialtype.(trialtypes{i})(j).events.move.theta.potential_mu = nanmean(lfps_temp2);
                        stats.trialtype.(trialtypes{i})(j).events.move.theta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                        stats.trialtype.(trialtypes{i})(j).events.move.theta.time = prs.ts.move;
                    end
                    if analyse_beta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_move], 'lfp_beta');
                        lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.move)';
                        stats.trialtype.(trialtypes{i})(j).events.move.beta.potential_mu = nanmean(lfps_temp2);
                        stats.trialtype.(trialtypes{i})(j).events.move.beta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                        stats.trialtype.(trialtypes{i})(j).events.move.beta.time = prs.ts.move;
                    end
                    %% compute spectrogram
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_move], 'lfp');
                    stats.trialtype.(trialtypes{i})(j).events.move.all_freq.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_move])-0.5 & ts < 2,:);
                    stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_lfp_align = ts(ts > -max([events_temp.t_move])-0.5 & ts < 2);
                    [stats.trialtype.(trialtypes{i})(j).events.move.all_freq.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.move.all_freq.freq_spectrogram] = ...
                        mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.move.all_freq.lfp_align,prs.spectrogram_movingwin,spectralparams);
                    %                                         figure; imagesc(stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_spectrogram-abs(stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_lfp_align(1)),...
                    %                                             stats.trialtype.(trialtypes{i})(j).events.move.all_freq.freq_spectrogram,real(stats.trialtype.(trialtypes{i})(j).events.move.p_spectrogram')); axis xy;
                    %                                         set(gca,'xlim',[-1.5 1.5], 'ylim',[1 80])
                    stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_spectrogram = stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_spectrogram-...
                        abs(stats.trialtype.(trialtypes{i})(j).events.move.all_freq.ts_lfp_align(1));
                    
                    %% align band passed signal to event theta
                    if analyse_theta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_move], 'lfp_theta');
                        stats.trialtype.(trialtypes{i})(j).events.move.theta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_move])-0.5 & ts < 2,:);
                        stats.trialtype.(trialtypes{i})(j).events.move.theta.ts_lfp_align = ts(ts > -max([events_temp.t_move])-0.5 & ts < 2);
                    end
                    
                    %% align band passed signal to event beta
                    if analyse_beta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_move], 'lfp_beta');
                        stats.trialtype.(trialtypes{i})(j).events.move.beta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_move])-0.5 & ts < 2,:);
                        stats.trialtype.(trialtypes{i})(j).events.move.beta.ts_lfp_align = ts(ts > -max([events_temp.t_move])-0.5 & ts < 2);
                    end
                    %%
                    if extract_band_pass && strcmp(trialtypes{i},'reward')
                        if j ==1   % incorrect
                            % for every trial take 95th percentile value
                            for ntrl = 1:sum(trlindx)
                                stats.band_passed.move.err.ts = stats.trialtype.(trialtypes{i})(j).events.move.beta.ts_lfp_align;
                                % theta
                                this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.move.theta.lfp_align(:,ntrl));
                                this_95_theta = prctile(this_theta,95);
                                % stats.band_passed.move.err.trl(ntrl).theta.lfp_move = this_theta;
                                stats.band_passed.move.err.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                                
                                % beta
                                this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.move.beta.lfp_align(:,ntrl));
                                this_95_beta = prctile(this_beta,95);
                                % stats.band_passed.move.err.trl(ntrl).beta.lfp_move = this_beta;
                                stats.band_passed.move.err.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                            end
                        else % correct
                            for ntrl = 1:sum(trlindx)
                                stats.band_passed.move.corr.ts = stats.trialtype.(trialtypes{i})(j).events.move.beta.ts_lfp_align;
                                % theta
                                this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.move.theta.lfp_align(:,ntrl));
                                this_95_theta = prctile(this_theta,95);
                                % stats.band_passed.move.corr.trl(ntrl).theta.lfp_move = this_theta;
                                stats.band_passed.move.corr.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                                
                                % beta
                                this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.move.beta.lfp_align(:,ntrl));
                                this_95_beta = prctile(this_beta,95);
                                % stats.band_passed.move.corr.trl(ntrl).beta.lfp_move = this_beta;
                                stats.band_passed.move.corr.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                            end
                        end
                    end
                    
                    %%
                    if analyse_phase && strcmp(trialtypes{i},'reward')
                        freq = [prs.lfp_freqmin:prs.lfp_freqmax];
                        this_lfp = fillmissing(stats.trialtype.(trialtypes{i})(j).events.move.all_freq.lfp_align,'constant',0);
                        for f_indx = 1:length(freq)
                            if f_indx ~= length(freq)
                                % extract analytic form for each freq
                                [b,a] = butter(prs.lfp_filtorder,[freq(f_indx) freq(f_indx+1)]/(prs.fs_lfp/2));
                                lfp_filt = filtfilt(b,a,this_lfp);
                                lfp_angle = angle(hilbert(lfp_filt)); % extract phase for that freq
                                stats.trialtype.(trialtypes{i})(j).events.move.all_freq.lfp_angle(f_indx,:,:) = lfp_angle;
                                % Compute phase clustering for all trials
                                % in one timepoint: abs(mean(exp(1i*angles_at_one_time_point_across_trials)))-->(Lachaux et al 1999)
                                if analyse_phase_within_area
                                    for tmp = 1:size(lfp_angle,1)
                                        stats.trialtype.(trialtypes{i})(j).events.move.all_freq.plv(f_indx,tmp) = abs( mean( exp( 1i * lfp_angle(tmp,:) )));
                                    end
                                end
                            end
                        end
                    end
                    
                end
                %% aligned to target onset
                if any(strcmp(gettuning,'target'))
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_targ],'lfp');
                    lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.target)';
                    stats.trialtype.(trialtypes{i})(j).events.target.potential_mu = nanmean(lfps_temp2);
                    stats.trialtype.(trialtypes{i})(j).events.target.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                    stats.trialtype.(trialtypes{i})(j).events.target.time = prs.ts.target;
                    
                    if analyse_theta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_targ],'lfp_theta');
                        lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.target)';
                        stats.trialtype.(trialtypes{i})(j).events.target.theta.potential_mu = nanmean(lfps_temp2);
                        stats.trialtype.(trialtypes{i})(j).events.target.theta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                        stats.trialtype.(trialtypes{i})(j).events.target.theta.time = prs.ts.target;
                    end
                    if analyse_beta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_targ],'lfp_beta');
                        lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.target)';
                        stats.trialtype.(trialtypes{i})(j).events.target.beta.potential_mu = nanmean(lfps_temp2);
                        stats.trialtype.(trialtypes{i})(j).events.target.beta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                        stats.trialtype.(trialtypes{i})(j).events.target.beta.time = prs.ts.target;
                    end
                    %% compute spectrogram
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_targ],'lfp');
                    stats.trialtype.(trialtypes{i})(j).events.target.all_freq.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_targ])-0.5 & ts < 2,:);
                    stats.trialtype.(trialtypes{i})(j).events.target.all_freq.ts_lfp_align = ts(ts > -max([events_temp.t_targ])-0.5 & ts < 2);
                    [stats.trialtype.(trialtypes{i})(j).events.target.all_freq.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.target.all_freq.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.target.all_freq.freq_spectrogram] = ...
                        mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.target.all_freq.lfp_align,prs.spectrogram_movingwin,spectralparams);
                    %                     figure; imagesc((stats.trialtype.(trialtypes{i})(j).events.target.ts_spectrogram),stats.trialtype.(trialtypes{i})(j).events.target.freq_spectrogram,real(stats.trialtype.(trialtypes{i})(j).events.target.p_spectrogram'), [0 0.8e-04]); axis xy;
                    %                     set(gca,'xlim',[-0.5 0.5]) %'ylim',[4 80])
                    
                    stats.trialtype.(trialtypes{i})(j).events.target.all_freq.ts_spectrogram = stats.trialtype.(trialtypes{i})(j).events.target.all_freq.ts_spectrogram -...
                        abs(stats.trialtype.(trialtypes{i})(j).events.target.all_freq.ts_lfp_align(1));
                    %% align band passed to event lfp theta
                    if analyse_theta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_targ], 'lfp_theta');
                        stats.trialtype.(trialtypes{i})(j).events.target.theta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_targ])-0.5 & ts < 2,:);
                        stats.trialtype.(trialtypes{i})(j).events.target.theta.ts_lfp_align = ts(ts > -max([events_temp.t_targ])-0.5 & ts < 2);
                    end
                    %% align band passed to event lfp beta
                    if analyse_beta
                        [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_targ], 'lfp_beta');
                        stats.trialtype.(trialtypes{i})(j).events.target.beta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_targ])-0.5 & ts < 2,:);
                        stats.trialtype.(trialtypes{i})(j).events.target.beta.ts_lfp_align = ts(ts > -max([events_temp.t_targ])-0.5 & ts < 2);
                    end
                    %% band passed
                    if extract_band_pass && strcmp(trialtypes{i},'reward')
                        if j ==1   % incorrect
                            % for every trial take 95th percentile value
                            for ntrl = 1:sum(trlindx)
                                stats.band_passed.target.err.ts = stats.trialtype.(trialtypes{i})(j).events.target.beta.ts_lfp_align;
                                % theta
                                this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.target.theta.lfp_align(:,ntrl));
                                this_95_theta = prctile(this_theta,95);
                                % stats.band_passed.target.err.trl(ntrl).theta.lfp_target = this_theta;
                                stats.band_passed.target.err.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                                
                                % beta
                                this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.target.beta.lfp_align(:,ntrl));
                                this_95_beta = prctile(this_beta,95);
                                % stats.band_passed.target.err.trl(ntrl).beta.lfp_target = this_beta;
                                stats.band_passed.target.err.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                            end
                        else % correct
                            for ntrl = 1:sum(trlindx)
                                stats.band_passed.target.corr.ts = stats.trialtype.(trialtypes{i})(j).events.target.beta.ts_lfp_align;
                                % theta
                                this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.target.theta.lfp_align(:,ntrl));
                                this_95_theta = prctile(this_theta,95);
                                % stats.band_passed.target.corr.trl(ntrl).theta.lfp_target = this_theta;
                                stats.band_passed.target.corr.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                                
                                % beta
                                this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.target.beta.lfp_align(:,ntrl));
                                this_95_beta = prctile(this_beta,95);
                                % stats.band_passed.target.corr.trl(ntrl).beta.lfp_target = this_beta;
                                stats.band_passed.target.corr.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                            end
                        end
                    end
                    
                    %% phase
                    if analyse_phase && strcmp(trialtypes{i},'reward')
                        freq = [prs.lfp_freqmin:prs.lfp_freqmax];
                        this_lfp = fillmissing(stats.trialtype.(trialtypes{i})(j).events.target.all_freq.lfp_align,'constant',0);
                        for f_indx = 1:length(freq)
                            if f_indx ~= length(freq)
                                % extract analytic form for each freq
                                [b,a] = butter(prs.lfp_filtorder,[freq(f_indx) freq(f_indx+1)]/(prs.fs_lfp/2));
                                lfp_filt = filtfilt(b,a,this_lfp);
                                lfp_angle = angle(hilbert(lfp_filt)); % extract phase for that freq
                                stats.trialtype.(trialtypes{i})(j).events.target.all_freq.lfp_angle(f_indx,:,:) = lfp_angle;
                                % Compute phase clustering for all trials
                                % in one timepoint: abs(mean(exp(1i*angles_at_one_time_point_across_trials)))-->(Lachaux et al 1999)
                                if analyse_phase_within_area
                                    for tmp = 1:size(lfp_angle,1)
                                        stats.trialtype.(trialtypes{i})(j).events.target.all_freq.plv(f_indx,tmp) = abs( mean( exp( 1i * lfp_angle(tmp,:) )));
                                    end
                                end
                            end
                        end
                    end
                end
            end
            %% aligned to movement stop
            if any(strcmp(gettuning,'stop'))
                [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_stop], 'lfp');
                lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.stop)';
                stats.trialtype.(trialtypes{i})(j).events.stop.potential_mu = nanmean(lfps_temp2);
                stats.trialtype.(trialtypes{i})(j).events.stop.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                stats.trialtype.(trialtypes{i})(j).events.stop.time = prs.ts.stop;
                
                if analyse_theta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_stop],'lfp_theta');
                    lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.stop)';
                    stats.trialtype.(trialtypes{i})(j).events.stop.theta.potential_mu = nanmean(lfps_temp2);
                    stats.trialtype.(trialtypes{i})(j).events.stop.theta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                    stats.trialtype.(trialtypes{i})(j).events.stop.theta.time = prs.ts.stop;
                end
                if analyse_beta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_stop],'lfp_beta');
                    lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.stop)';
                    stats.trialtype.(trialtypes{i})(j).events.stop.beta.potential_mu = nanmean(lfps_temp2);
                    stats.trialtype.(trialtypes{i})(j).events.stop.beta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                    stats.trialtype.(trialtypes{i})(j).events.stop.beta.time = prs.ts.stop;
                end
                %% compute spectrogram lfp not band passed
                [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_stop], 'lfp');
                stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_stop])-0.5 & ts < 2,:);
                stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_lfp_align = ts(ts > -max([events_temp.t_stop])-0.5 & ts < 2);
                [stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.freq_spectrogram] = ...
                    mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.lfp_align,prs.spectrogram_movingwin,spectralparams);
                %                                         figure; imagesc(stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_spectrogram-abs(stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_lfp_align(1)),...
                %                                             stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.freq_spectrogram,real(stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.p_spectrogram')); axis xy;
                %                                         set(gca,'xlim',[-1.5 1.5])
                stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_spectrogram = stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_spectrogram-...
                    abs(stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.ts_lfp_align(1));
                
                %% compute spectrogram lfp theta
                if analyse_theta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_stop], 'lfp_theta');
                    stats.trialtype.(trialtypes{i})(j).events.stop.theta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_stop])-0.5 & ts < 2,:);
                    stats.trialtype.(trialtypes{i})(j).events.stop.theta.ts_lfp_align = ts(ts > -max([events_temp.t_stop])-0.5 & ts < 2);
                    % %                     [stats.trialtype.(trialtypes{i})(j).events.stop.theta.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.theta.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.theta.freq_spectrogram] = ...
                    % %                         mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.stop.theta.lfp_align,prs.spectrogram_movingwin,spectralparams);
                end
                %% compute spectrogram lfp beta
                if analyse_beta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_stop], 'lfp_beta');
                    stats.trialtype.(trialtypes{i})(j).events.stop.beta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_stop])-0.5 & ts < 2,:);
                    stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_lfp_align = ts(ts > -max([events_temp.t_stop])-0.5 & ts < 2);
                    [stats.trialtype.(trialtypes{i})(j).events.stop.beta.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.beta.freq_spectrogram] = ...
                        mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.stop.beta.lfp_align,prs.spectrogram_movingwin,spectralparams);
                    %                     figure; imagesc(stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_spectrogram-1,...
                    %                         stats.trialtype.(trialtypes{i})(j).events.stop.beta.freq_spectrogram,real(stats.trialtype.(trialtypes{i})(j).events.stop.beta.p_spectrogram')); axis xy;
                    %                     set(gca,'xlim',[-0.75 0.75]) %'ylim',[4 80])
                end
                %% band passed
                if extract_band_pass && strcmp(trialtypes{i},'reward')
                    if j == 1   % incorrect
                        % for every trial take 95th percentile value
                        for ntrl = 1:sum(trlindx)
                            stats.band_passed.stop.err.ts = stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_lfp_align;
                            % theta
                            this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.stop.theta.lfp_align(:,ntrl));
                            this_95_theta = prctile(this_theta,95);
                            % stats.band_passed.stop.err.trl(ntrl).theta.lfp_stop = this_theta;
                            stats.band_passed.stop.err.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                            
                            % beta
                            this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.stop.beta.lfp_align(:,ntrl));
                            this_95_beta = prctile(this_beta,95);
                            % stats.band_passed.stop.err.trl(ntrl).beta.lfp_stop = this_beta;
                            stats.band_passed.stop.err.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                        end
                    else % correct
                        for ntrl = 1:sum(trlindx)
                            stats.band_passed.stop.corr.ts = stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_lfp_align;
                            % theta
                            this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.stop.theta.lfp_align(:,ntrl));
                            this_95_theta = prctile(this_theta,95);
                            % stats.band_passed.stop.corr.trl(ntrl).theta.lfp_stop = this_theta;
                            stats.band_passed.stop.corr.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                            
                            % beta
                            this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.stop.beta.lfp_align(:,ntrl));
                            this_95_beta = prctile(this_beta,95);
                            % stats.band_passed.stop.corr.trl(ntrl).beta.lfp_stop = this_beta;
                            stats.band_passed.stop.corr.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                        end
                    end
                end
                
                %% phase
                if analyse_phase && strcmp(trialtypes{i},'reward')
                    freq = [prs.lfp_freqmin:prs.lfp_freqmax];
                    this_lfp = fillmissing(stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.lfp_align,'constant',0);
                    for f_indx = 1:length(freq)
                        if f_indx ~= length(freq)
                            % extract analytic form for each freq
                            [b,a] = butter(prs.lfp_filtorder,[freq(f_indx) freq(f_indx+1)]/(prs.fs_lfp/2));
                            lfp_filt = filtfilt(b,a,this_lfp);
                            lfp_angle = angle(hilbert(lfp_filt)); % extract phase for that freq
                            stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.lfp_angle(f_indx,:,:) = lfp_angle;
                            % Compute phase clustering for all trials
                            % in one timepoint: abs(mean(exp(1i*angles_at_one_time_point_across_trials)))-->(Lachaux et al 1999)
                            if analyse_phase_within_area
                                for tmp = 1:size(lfp_angle,1)
                                    stats.trialtype.(trialtypes{i})(j).events.stop.all_freq.plv(f_indx,tmp) = abs( mean( exp( 1i * lfp_angle(tmp,:) )));
                                end
                            end
                        end
                    end
                end
            end
            %% aligned to reward
            if any(strcmp(gettuning,'reward'))
                [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_rew], 'lfp');
                lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.reward)';
                stats.trialtype.(trialtypes{i})(j).events.reward.potential_mu = nanmean(lfps_temp2);
                stats.trialtype.(trialtypes{i})(j).events.reward.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                stats.trialtype.(trialtypes{i})(j).events.reward.time = prs.ts.reward;
                if analyse_theta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_rew],'lfp_theta');
                    lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.reward)';
                    stats.trialtype.(trialtypes{i})(j).events.reward.theta.potential_mu = nanmean(lfps_temp2);
                    stats.trialtype.(trialtypes{i})(j).events.reward.theta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                    stats.trialtype.(trialtypes{i})(j).events.reward.theta.time = prs.ts.reward;
                end
                if analyse_beta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_rew],'lfp_beta');
                    lfps_temp2 = interp1(ts,(trials_lfps_temp2),prs.ts.reward)';
                    stats.trialtype.(trialtypes{i})(j).events.reward.beta.potential_mu = nanmean(lfps_temp2);
                    stats.trialtype.(trialtypes{i})(j).events.reward.beta.potential_sem = nanstd(lfps_temp2)/sqrt(size(lfps_temp2,1));
                    stats.trialtype.(trialtypes{i})(j).events.reward.beta.time = prs.ts.reward;
                end
                %% compute spectrogram lfp not band passed
                [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_rew], 'lfp');
                stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_rew])-0.5 & ts < 2,:);
                stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_lfp_align = ts(ts > -max([events_temp.t_rew])-0.5 & ts < 2);
                [stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.freq_spectrogram] = ...
                    mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.lfp_align,prs.spectrogram_movingwin,spectralparams);
                %                     figure; imagesc(stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_spectrogram,...
                %                         stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.freq_spectrogram,real(stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.p_spectrogram')); axis xy;
                %                     set(gca,'xlim',[-0.75 0.75])
                if ~isempty(stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_lfp_align)
                    stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_spectrogram = stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_spectrogram - ...
                        abs(stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.ts_lfp_align(1));
                end
                %% compute spectrogram lfp theta
                if analyse_theta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_rew], 'lfp_theta');
                    stats.trialtype.(trialtypes{i})(j).events.reward.theta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_rew])-0.5 & ts < 2,:);
                    stats.trialtype.(trialtypes{i})(j).events.reward.theta.ts_lfp_align = ts(ts > -max([events_temp.t_rew])-0.5 & ts < 2);
                    % %                     [stats.trialtype.(trialtypes{i})(j).events.reward.theta.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.reward.theta.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.reward.theta.freq_spectrogram] = ...
                    % %                         mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.reward.theta.lfp_align,prs.spectrogram_movingwin,spectralparams);
                    %
                end
                %% compute spectrogram lfp beta
                if analyse_beta
                    [trials_lfps_temp2,ts] = ShiftLfps(trials_lfps_temp,continuous_temp,[events_temp.t_rew], 'lfp_beta');
                    stats.trialtype.(trialtypes{i})(j).events.reward.beta.lfp_align = trials_lfps_temp2(ts > -max([events_temp.t_stop])-0.5 & ts < 2,:);
                    stats.trialtype.(trialtypes{i})(j).events.reward.beta.ts_lfp_align = ts(ts > -max([events_temp.t_stop])-0.5 & ts < 2);
                    [stats.trialtype.(trialtypes{i})(j).events.reward.beta.p_spectrogram, stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_spectrogram, stats.trialtype.(trialtypes{i})(j).events.reward.beta.freq_spectrogram] = ...
                        mtspecgramc(stats.trialtype.(trialtypes{i})(j).events.reward.beta.lfp_align,prs.spectrogram_movingwin,spectralparams);
                    %                     figure; imagesc(stats.trialtype.(trialtypes{i})(j).events.stop.beta.ts_spectrogram-1,...
                    %                         stats.trialtype.(trialtypes{i})(j).events.stop.beta.freq_spectrogram,real(stats.trialtype.(trialtypes{i})(j).events.stop.beta.p_spectrogram')); axis xy;
                    %                     set(gca,'xlim',[-0.75 0.75]) %'ylim',[4 80])
                end
                
                %% band passed
                if extract_band_pass && strcmp(trialtypes{i},'reward')
                    if j ==1   % incorrect
                        % for every trial take 95th percentile value
                        for ntrl = 1:sum(trlindx)
                            stats.band_passed.reward.err.ts = stats.trialtype.(trialtypes{i})(j).events.reward.beta.ts_lfp_align;
                            % theta
                            this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.reward.theta.lfp_align(:,ntrl));
                            this_95_theta = prctile(this_theta,95);
                            % stats.band_passed.reward.err.trl(ntrl).theta.lfp_reward = this_theta;
                            stats.band_passed.reward.err.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                            
                            % beta
                            this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.reward.beta.lfp_align(:,ntrl));
                            this_95_beta = prctile(this_beta,95);
                            % stats.band_passed.reward.err.trl(ntrl).beta.lfp_reward = this_beta;
                            stats.band_passed.reward.err.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                        end
                    else % correct
                        for ntrl = 1:sum(trlindx)
                            stats.band_passed.reward.corr.ts = stats.trialtype.(trialtypes{i})(j).events.reward.beta.ts_lfp_align;
                            % theta
                            this_theta = abs(stats.trialtype.(trialtypes{i})(j).events.reward.theta.lfp_align(:,ntrl));
                            this_95_theta = prctile(this_theta,95);
                            % stats.band_passed.reward.corr.trl(ntrl).theta.lfp_reward = this_theta;
                            stats.band_passed.reward.corr.trl(ntrl).theta_95_indx = this_theta > this_95_theta;
                            
                            % beta
                            this_beta = abs(stats.trialtype.(trialtypes{i})(j).events.reward.beta.lfp_align(:,ntrl));
                            this_95_beta = prctile(this_beta,95);
                            % stats.band_passed.reward.corr.trl(ntrl).beta.lfp_reward = this_beta;
                            stats.band_passed.reward.corr.trl(ntrl).beta_95_indx = this_beta > this_95_beta;
                        end
                    end
                end
                
                %% phase
                if analyse_phase && strcmp(trialtypes{i},'reward')
                    if j == 2
                        freq = [prs.lfp_freqmin:prs.lfp_freqmax];
                        this_lfp = fillmissing(stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.lfp_align,'constant',0);
                        for f_indx = 1:length(freq)
                            if f_indx ~= length(freq)
                                % extract analytic form for each freq
                                [b,a] = butter(prs.lfp_filtorder,[freq(f_indx) freq(f_indx+1)]/(prs.fs_lfp/2));
                                lfp_filt = filtfilt(b,a,this_lfp);
                                lfp_angle = angle(hilbert(lfp_filt)); % extract phase for that freq
                                stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.lfp_angle(f_indx,:,:) = lfp_angle;
                                % Compute phase clustering for all trials
                                % in one timepoint: abs(mean(exp(1i*angles_at_one_time_point_across_trials)))-->(Lachaux et al 1999)
                                if analyse_phase_within_area
                                    for tmp = 1:size(lfp_angle,1)
                                        stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.plv(f_indx,tmp) = abs( mean( exp( 1i * lfp_angle(tmp,:) )));
                                    end
                                end
                            end
                        end
                    else
                        stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.plv = NaN;
                        stats.trialtype.(trialtypes{i})(j).events.reward.all_freq.lfp_angle = NaN;
                    end
                end
            end
            
        end
    end
end

%% power spectral density
if compute_psd
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    % during trials
    for i=1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).spectrum = stats.trialtype.all.spectrum;
            else
                sMarkers = [];
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                trials_lfps_temp = trials_lfps(trlindx);
                %%
                lfp_concat = cell2mat({trials_lfps_temp.lfp}'); % concatenate trials
                triallen = cellfun(@(x) length(x), {trials_lfps_temp.lfp});
                sMarkers(:,1) = cumsum([1 triallen(1:end-1)]); sMarkers(:,2) = cumsum(triallen); % demarcate trial onset and end
                [stats.trialtype.(trialtypes{i})(j).spectrum.psd , stats.trialtype.(trialtypes{i})(j).spectrum.freq] = ...
                    mtspectrumc_unequal_length_trials(lfp_concat, prs.spectrum_movingwin , spectralparams, sMarkers); % needs http://chronux.org/
            end
        end
    end
    
    %% stationary period
    stationary_lfps_temp = []; sMarkers = [];
    for i=1:length(stationary_lfps)
        if ~isempty(stationary_lfps(i).lfp) % gather available inter-trials
            stationary_lfps_temp(end+1).lfp = stationary_lfps(i).lfp;
        end
    end
    lfp_concat = cell2mat({stationary_lfps_temp.lfp}); % concatenate trials
    triallen = cellfun(@(x) length(x), {stationary_lfps_temp.lfp});
    sMarkers(:,1) = cumsum([1 triallen(1:end-1)]); sMarkers(:,2) = cumsum(triallen); % demarcate trial onset and end
    [stats.trialtype.stationary.spectrum.psd , stats.trialtype.stationary.spectrum.freq] = ...
        mtspectrumc_unequal_length_trials(lfp_concat(:), [1 1] , spectralparams, sMarkers); % needs http://chronux.org/
    
    %% mobile period
    mobile_lfps_temp = []; sMarkers = [];
    trlindx = behv_stats.trialtype.all.trlindx; mobile_lfps = mobile_lfps(trlindx);
    for i=1:length(mobile_lfps)
        if ~isempty(mobile_lfps(i).lfp) % gather available inter-trials
            mobile_lfps_temp(end+1).lfp = mobile_lfps(i).lfp;
        end
    end
    lfp_concat = cell2mat({mobile_lfps_temp.lfp}); % concatenate trials
    triallen = cellfun(@(x) length(x), {mobile_lfps_temp.lfp});
    sMarkers(:,1) = cumsum([1 triallen(1:end-1)]); sMarkers(:,2) = cumsum(triallen); % demarcate trial onset and end
    [stats.trialtype.mobile.spectrum.psd , stats.trialtype.mobile.spectrum.freq] = ...
        mtspectrumc_unequal_length_trials(lfp_concat(:), [1 1] , spectralparams, sMarkers); % needs http://chronux.org/
    
    %     %% eyes-fixed period
    %     eyesfixed_lfps_temp = []; sMarkers = [];
    %     for i=1:length(eyesfixed_lfps)
    %         if ~isempty(eyesfixed_lfps(i).lfp)
    %             eyesfixed_lfps_temp(end+1).lfp = eyesfixed_lfps(i).lfp;
    %         end
    %     end
    %     lfp_concat = cell2mat({eyesfixed_lfps_temp.lfp}); % concatenate trials
    %     triallen = cellfun(@(x) length(x), {eyesfixed_lfps_temp.lfp});
    %     sMarkers(:,1) = cumsum([1 triallen(1:end-1)]); sMarkers(:,2) = cumsum(triallen); % demarcate trial onset and end
    %     [stats.trialtype.eyesfixed.spectrum.psd , stats.trialtype.eyesfixed.spectrum.freq] = ...
    %         mtspectrumc_unequal_length_trials(lfp_concat(:), [fixateduration fixateduration] , spectralparams, sMarkers); % needs http://chronux.org/
    %
    %     %% eyes-free period
    %     eyesfree_lfps_temp = []; sMarkers = [];
    %     for i=1:length(eyesfree_lfps)
    %         if ~isempty(eyesfree_lfps(i).lfp)
    %             eyesfree_lfps_temp(end+1).lfp = eyesfree_lfps(i).lfp;
    %         end
    %     end
    %     lfp_concat = cell2mat({eyesfree_lfps_temp.lfp}); % concatenate trials
    %     triallen = cellfun(@(x) length(x), {eyesfree_lfps_temp.lfp});
    %     sMarkers(:,1) = cumsum([1 triallen(1:end-1)]); sMarkers(:,2) = cumsum(triallen); % demarcate trial onset and end
    %     [stats.trialtype.eyesfree.spectrum.psd , stats.trialtype.eyesfree.spectrum.freq] = ...
    %         mtspectrumc_unequal_length_trials(lfp_concat(:), [fixateduration fixateduration] , spectralparams, sMarkers); % needs http://chronux.org/
end

%% theta LFP
trials_theta(ntrls) = struct();
if analyse_theta
    for i=1:ntrls
        trials_theta(i).lfp = trials_lfps(i).lfp_theta(:); % read as column vector
        theta_freq = [(1/dt)/(2*pi)*diff(unwrap(angle(trials_theta(i).lfp))) ; nan];
        theta_freq(theta_freq<prs.lfp_theta(1) | theta_freq>prs.lfp_theta(2)) = nan;
        trials_theta(i).freq = theta_freq;
    end
    
    
    for i=1%:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        for j=1:nconds
            trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
            events_temp = events(trlindx);
            continuous_temp = continuous(trlindx);
            trials_theta_temp = trials_theta(trlindx);
            %% define time windows for computing tuning
            timewindow_move = [[events_temp.t_move]' [events_temp.t_stop]']; % when the subject is moving
            %% linear velocity, v
            stats.trialtype.(trialtypes{i})(j).continuous.v.thetafreq = ...
                ComputeTuning({continuous_temp.v},{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method);
            %% angular velocity, w
            stats.trialtype.(trialtypes{i})(j).continuous.w.thetafreq = ...
                ComputeTuning({continuous_temp.w},{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method);
            %% vw
            stats.trialtype.(trialtypes{i})(j).continuous.vw.thetafreq = ...
                ComputeTuning2D({continuous_temp.v},{continuous_temp.w},{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,prs.tuning,prs.tuning_method);
            %% horizontal eye velocity
            %             heye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false); % average both eyes (if available)
            %             heyevel = cellfun(@(x) [0 ; diff(x)'/dt],heye,'UniformOutput',false);
            %             stats.trialtype.(trialtypes{i})(j).continuous.heyevel.thetafreq = ...
            %                 ComputeTuning(heyevel,{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.heye_vel);
            %% vertical velocity
            %             veye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false); % average both eyes (if available)
            %             veyevel = cellfun(@(x) [0 ; diff(x)'/dt],veye,'UniformOutput',false);
            %             stats.trialtype.(trialtypes{i})(j).continuous.veyevel.thetafreq = ...
            %                 ComputeTuning(veyevel,{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.veye_vel);
        end
    end
end
%% alpha LFP
trials_alpha(ntrls) = struct();
if analyse_alpha
    for i=1:ntrls
        trials_alpha(i).lfp = trials_lfps(i).lfp_alpha(:); % read as column vector
        alpha_freq = [(1/dt)/(2*pi)*diff(unwrap(angle(trials_alpha(i).lfp))) ; nan];
        alpha_freq(alpha_freq<prs.lfp_alpha(1) | alpha_freq>prs.lfp_alpha(2)) = nan;
        trials_alpha(i).freq = alpha_freq;
    end
    for i=1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        for j=1:nconds
            trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
            events_temp = events(trlindx);
            continuous_temp = continuous(trlindx);
            trials_alpha_temp = trials_alpha(trlindx);
            %% define time windows for computing tuning
            timewindow_move = [[events_temp.t_move]' [events_temp.t_stop]']; % when the subject is moving
            %% linear velocity, v
            stats.trialtype.(trialtypes{i})(j).continuous.v.alphafreq = ...
                ComputeTuning({continuous_temp.v},{continuous_temp.ts},{trials_alpha_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method);
            %% angular velocity, w
            stats.trialtype.(trialtypes{i})(j).continuous.w.alphafreq = ...
                ComputeTuning({continuous_temp.w},{continuous_temp.ts},{trials_alpha_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method);
            %% vw
            stats.trialtype.(trialtypes{i})(j).continuous.vw.alphafreq = ...
                ComputeTuning2D({continuous_temp.v},{continuous_temp.w},{continuous_temp.ts},{trials_alpha_temp.freq},timewindow_move,prs.tuning,prs.tuning_method);
            %% horizontal eye velocity
            %             heye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false); % average both eyes (if available)
            %             heyevel = cellfun(@(x) [0 ; diff(x)'/dt],heye,'UniformOutput',false);
            %             stats.trialtype.(trialtypes{i})(j).continuous.heyevel.thetafreq = ...
            %                 ComputeTuning(heyevel,{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.heye_vel);
            %% vertical velocity
            %             veye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false); % average both eyes (if available)
            %             veyevel = cellfun(@(x) [0 ; diff(x)'/dt],veye,'UniformOutput',false);
            %             stats.trialtype.(trialtypes{i})(j).continuous.veyevel.thetafreq = ...
            %                 ComputeTuning(veyevel,{continuous_temp.ts},{trials_theta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.veye_vel);
        end
    end
end

%% beta LFP
trials_beta(ntrls) = struct();
if analyse_beta
    for i=1:ntrls
        trials_beta(i).lfp = trials_lfps(i).lfp_beta(:); % read as column vector
        beta_freq = [(1/dt)/(2*pi)*diff(unwrap(angle(trials_beta(i).lfp))) ; nan];
        beta_freq(beta_freq<prs.lfp_beta(1) | beta_freq>prs.lfp_beta(2)) = nan;
        trials_beta(i).freq = beta_freq;
    end
    for i=1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        for j=1:nconds
            trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
            events_temp = events(trlindx);
            continuous_temp = continuous(trlindx);
            trials_beta_temp = trials_beta(trlindx);
            %% define time windows for computing tuning
            timewindow_move = [[events_temp.t_move]' [events_temp.t_stop]']; % when the subject is moving
            %% linear velocity, v
            stats.trialtype.(trialtypes{i})(j).continuous.v.betafreq = ...
                ComputeTuning({continuous_temp.v},{continuous_temp.ts},{trials_beta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method);
            %% angular velocity, w
            stats.trialtype.(trialtypes{i})(j).continuous.w.betafreq = ...
                ComputeTuning({continuous_temp.w},{continuous_temp.ts},{trials_beta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method);
            %% vw
            stats.trialtype.(trialtypes{i})(j).continuous.vw.betafreq = ...
                ComputeTuning2D({continuous_temp.v},{continuous_temp.w},{continuous_temp.ts},{trials_beta_temp.freq},timewindow_move,prs.tuning,prs.tuning_method);
            %% horizontal eye velocity
            %             heye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false); % average both eyes (if available)
            %             heyevel = cellfun(@(x) [0 ; diff(x)'/dt],heye,'UniformOutput',false);
            %             stats.trialtype.(trialtypes{i})(j).continuous.heyevel.betafreq = ...
            %                 ComputeTuning(heyevel,{continuous_temp.ts},{trials_beta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.heye_vel);
            %% vertical velocity
            %             veye = cellfun(@(x,y) nanmean([x(:)' ; y(:)']),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false); % average both eyes (if available)
            %             veyevel = cellfun(@(x) [0 ; diff(x)'/dt],veye,'UniformOutput',false);
            %             stats.trialtype.(trialtypes{i})(j).continuous.veyevel.betafreq = ...
            %                 ComputeTuning(veyevel,{continuous_temp.ts},{trials_beta_temp.freq},timewindow_move,duration_zeropad,corr_lag,nbootstraps,prs.tuning,prs.tuning_method,prs.binrange.veye_vel);
        end
    end
end
