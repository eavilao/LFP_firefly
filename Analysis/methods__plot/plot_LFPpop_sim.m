function plot_LFPpop_sim(monk, plot_type)

% plot from LFPpop_sim
% names of variables
% monk: contains data per monkey
%
%% Plot directory

% 'trial_times'
% 'behavior_steering'
% 'move_on_target_on_distribution'
% 'erp_single_sess_all'
% 'erp_single_sess_reward_densities'
% 'erp_all_MST' (and max amp and time)
% 'erp_all_PPC' (and max amp and time)
% 'erp_reward_densities_MST' (and max amp and time)
% 'erp_before_after_move_MST'
% 'erp_all_PPC' (and max amp and time)
% 'erp_reward_densities_PPC' (and max amp and time)
% 'erp_move_before_after_PPC'
% 'PSD_all'
% 'PSD_all_together'
% 'PSD_reward'
% 'PSD_move'
% 'PSD_move_together'
% 'PSD_move_together_modulation'
% 'PSD_eye'
% 'PSD_reward'
% 'spectrogram_move'
% 'spectrogram_target'
% 'spectrogram_stop'
% 'spectrogram_reward_density'
% 'spectrogram_reward_density_all_theta'
% 'spectrogram_reward_density_all_theta_diff'
% 'spectrogram_reward_density_all_beta'
% 'spectrogram_reward_density_all_beta_diff'
% 'spectrogram_reward_density_diff'
% 'spectrogram_trl'
% 'spectrogram_trl_session_align_target'
% 'spectrogram_trl_session_align_stop'
% 'coherogram_move'
% 'coherogram_target'
% 'coherogram_stop'
% 'coherogram_reward_PPC_MST'
% 'coherogram_target_stop_PPC_MST'
% 'coherence_dist'
% 'coherence_dist_all'
% 'coherence_dist_reward_density'
% 'sim_coherence_between_areas'
% 'speed_dependent_LFP_activity_MST'
% 'speed_dependent_LFP_activity_PPC'
% 'trial_band_passed'
% 'pop_psth_band_passed'
% 'psth_peak_pop'
% 'band_passed_vs_accuracy'
% 'distance_vs_reward'

% print('-painters', '-depsc2', 'raster')

%% Plots per monkey

switch plot_type
    
    case 'trial_times'
        for m = 1:length(monk)
            for sess = 1:length(monk(m).behavior)
                trl_indx = monk(m).behavior(sess).stats.trialtype.reward(2).trlindx; corr = find(trl_indx);
                for ii = 1:length(corr)
                    trl_t(m,sess,ii) = monk(m).behavior(sess).trials(corr(ii)).events.t_rew - monk(m).behavior(sess).trials(corr(ii)).events.t_beg_correction;
                    trl_rew(m,sess,ii) = monk(m).behavior(sess).trials(corr(ii)).events.t_rew;
                end
                max_trl_rew = max(squeeze(trl_rew(m,sess,:)))
                figure; hold on;
                title(['Monk ' num2str(m) ' sess ' num2str(sess)])
                histogram(squeeze(trl_rew(m,sess,:)),100)
                xlim([0.5 5]); ylim([0 80])
            end
        end
        
        
    case 'behavior_steering'
        %% behavioural data
        for m = 1:length(monk)
            for sess = 1:length(monk(m).behavior)
                correct = monk(m).behavior(sess).stats.trialtype.reward(strcmp({monk(m).behavior(sess).stats.trialtype.reward.val},'rewarded')).trlindx;
                incorrect = monk(m).behavior(sess).stats.trialtype.reward(strcmp({monk(m).behavior(sess).stats.trialtype.reward.val},'unrewarded')).trlindx;
                crazy = ~(correct | incorrect); ntrls = sum(~crazy);
                behv_correct = monk(m).behavior(sess).trials(correct); ntrls_correct = length(behv_correct);
                behv_incorrect = monk(m).behavior(sess).trials(incorrect); ntrls_incorrect = length(behv_incorrect);
                
                % velocity
                for j = 1:length(behv_correct)
                    this_ts = behv_correct(j).continuous.ts;
                    this_ts_move = behv_correct(j).continuous.ts-behv_correct(j).events.t_move;
                    this_ts_stop = behv_correct(j).continuous.ts-behv_correct(j).events.t_stop;
                    this_ts_reward = behv_correct(j).continuous.ts-behv_correct(j).events.t_rew;
                    ts_win = this_ts(this_ts>-0.51 & this_ts<0.51); ts_win_targ = this_ts(this_ts>-0.81 & this_ts<0.81);
                    
                    move_v = behv_correct(j).continuous.v(this_ts_move>-0.51 & this_ts_move<0.51); move_w = behv_correct(j).continuous.w(this_ts_move>-0.51 & this_ts_move<0.51);
                    targ_v = behv_correct(j).continuous.v(this_ts>-0.81 & this_ts<0.81); targ_w = behv_correct(j).continuous.w(this_ts>-0.81 & this_ts<0.81);
                    stop_v = behv_correct(j).continuous.v(this_ts_stop>-0.51 & this_ts_stop<0.51); stop_w = behv_correct(j).continuous.w(this_ts_stop>-0.51 & this_ts_stop<0.51);
                    rew_v = behv_correct(j).continuous.v(this_ts_reward>-0.51 & this_ts_reward<0.51); rew_w = behv_correct(j).continuous.w(this_ts_reward>-0.51 & this_ts_reward<0.51);
                    
                    % store
                    %                     v_this_move(j,:) = move_v(1:169);
                    %                     v_this_targ(j,:) = targ_v(1:269); % so dirrrty...yuuuk. Extra sample, somewhere.
                    %                     v_this_stop(j,:) = stop_v(1:169);
                    %                     v_this_rew(j,:) = rew_v(1:169);
                    %
                    %                     w_this_move(j,:) = abs(move_w(1:169));
                    %                     w_this_targ(j,:) = abs(targ_w(1:269));
                    %                     w_this_stop(j,:) = abs(stop_w(1:169));
                    %                     w_this_rew(j,:) = abs(rew_w(1:169));
                    
                    %% velocity distribution per session
                    t_v_targ(j) = behv_correct(j).continuous.v(this_ts>-0.001 & this_ts<0.005);
                    t_w_targ(j) = behv_correct(j).continuous.w(this_ts>-0.001 & this_ts<0.005);
                end
                % mean for a session
                                ts = ts_win(1:169); % ts_win(1:169);
                                ts_win_targ = ts_win_targ(1:269);
                
                                v_sess_move(sess,:) = nanmean(v_this_move);
                                v_sess_targ(sess,:) = nanmean(v_this_targ);
                                v_sess_stop(sess,:) = nanmean(v_this_stop);
                                v_sess_reward(sess,:) = nanmean(v_this_rew);
                
                                w_sess_move(sess,:) = nanmean(w_this_move);
                                w_sess_targ(sess,:) = nanmean(w_this_targ);
                                w_sess_stop(sess,:) = nanmean(w_this_stop);
                                w_sess_reward(sess,:) = nanmean(w_this_rew);
                
                % plot per distribution of velocties when target is on per
                % session
                figure; hold on
                histogram(t_v_targ,35, 'DisplayStyle', 'bar');
                set(gca,'TickDir', 'out', 'FontSize', 22); box off;
                xlabel('Linear velocity cm/s');
                
                figure; hold on
                histogram(t_w_targ,35, 'DisplayStyle', 'bar');
                set(gca,'TickDir', 'out', 'FontSize', 22); box off;
                xlabel('Angular velocity cm/s');
            end
        end
        
        %% plot move
        figure; hold on;
        [ax,~,~] = plotyy(ts, mean(v_sess_move),ts, mean(w_sess_move))
        set(ax(1),'YLim',[0 200]);
        set(ax(2),'YLim',[0 45]);
        set(gca,'xlim',[-0.5 0.5],'xTick',[],'TickDir', 'out', 'FontSize', 22); box off;
        vline(0,'k'); axis square;
        
        %% plot target
        figure; hold on;
        [ax,~,~] = plotyy(ts_win_targ, mean(v_sess_targ),ts_win_targ, mean(w_sess_targ));
        set(ax(1),'YLim',[0 200]);
        set(ax(2),'YLim',[0 45]);
        set(gca,'xlim',[-0.2 0.8],'xTick',[],'TickDir', 'out', 'FontSize', 22); box off;
        vline(0.3,'k'); vline(0,'k'); axis square;
        
        figure; plot(ts_win_targ, mean(v_sess_targ))
        set(gca,'xlim',[0.05 0.55],'xTick',[],'TickDir', 'out', 'FontSize', 22); box off; axis square; ylim([0 200]); vline(0.3,'k'); vline(0,'k');
        figure;plot(ts_win_targ, mean(w_sess_targ))
        set(gca,'xlim',[0.05 0.55],'xTick',[],'TickDir', 'out', 'FontSize', 22); box off; axis square; ylim([0 45]); vline(0.3,'k'); vline(0,'k');
        %% Velocity distribution at target onset per session
        
        
        %% plot stop
        figure; hold on;
        [ax,~,~] = plotyy(ts, mean(v_sess_stop),ts, mean(w_sess_stop))
        set(ax(1),'YLim',[0 200]);
        set(ax(2),'YLim',[0 45]);
        set(gca,'xlim',[-0.5 0.5],'xTick',[],'TickDir', 'out', 'FontSize', 22); box off;
        vline(0,'k'); axis square;
        
        %% plot reward
        figure; hold on;
        [ax,~,~] = plotyy(ts, mean(v_sess_reward),ts, mean(w_sess_reward))
        set(ax(1),'YLim',[0 200]);
        set(ax(2),'YLim',[0 45]);
        set(gca,'xlim',[-0.5 0.5],'xTick',[],'TickDir', 'out', 'FontSize', 22); box off;
        vline(0,'k'); axis square;
        
    case 'move_on_target_on_distribution'
        for m = 1:length(monk)
            for sess = 1:length(monk(m).behavior)
                correct = monk(m).behavior(sess).stats.trialtype.reward(strcmp({monk(m).behavior(sess).stats.trialtype.reward.val},'rewarded')).trlindx;
                incorrect = monk(m).behavior(sess).stats.trialtype.reward(strcmp({monk(m).behavior(sess).stats.trialtype.reward.val},'unrewarded')).trlindx;
                crazy = ~(correct | incorrect); ntrls = sum(~crazy);
                behv_correct = monk(m).behavior(sess).trials(correct); ntrls_correct = length(behv_correct);
                behv_incorrect = monk(m).behavior(sess).trials(incorrect); ntrls_incorrect = length(behv_incorrect);
                
                for j = 1:length(behv_correct)
                    target_t(j) = behv_correct(j).events.t_targ;
                    move_t(j) = behv_correct(j).events.t_move;
                end
                
                figure; hold on
                histogram(move_t,30, 'DisplayStyle', 'stairs');
                set(gca,'xlim',[-2 2],'TickDir', 'out', 'FontSize', 22); box off;
                vline(mean(move_t),'--r'); vline(median(move_t), '--g')
                xlabel('Time from target onset');
                
                figure; hold on;
                move_after = move_t>0;
                move_before = move_t<0;
                bar(1,sum(move_after));
                bar(2,sum(move_before));
                set(gca,'xTick', [],'TickDir', 'out', 'FontSize', 22); box off;
                legend({'after' 'before'}, 'box', 'off');
            end
        end
        
    case 'erp_single_sess_all'
        type = 'all'
        ts_move_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.move.erp_time;
        ts_target_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.target.erp_time;
        ts_stop_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.stop.erp_time;
        ts_reward_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.reward.erp_time;
        ts_move_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.move.erp_time;
        ts_target_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.target.erp_time;
        ts_stop_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.stop.erp_time;
        ts_reward_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.reward.erp_time;
        
        % Mean all channels per session MST
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for ch = 1:length(monk(1).cont.MST.sess(sess).lfps)
                    monk(m).sess(sess).erp_move_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.move.erp_mu*1000; % move
                    monk(m).sess(sess).erp_target_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.target.erp_mu*1000; % move
                    monk(m).sess(sess).erp_stop_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.stop.erp_mu*1000; % move
                    monk(m).sess(sess).erp_reward_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.reward.erp_mu*1000; % move
                end
            end
        end
        % Mean all channels per session PPC
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for ch = 1:length(monk(1).cont.PPC.sess(sess).lfps)
                    monk(m).sess(sess).erp_move_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.move.erp_mu; % move
                    monk(m).sess(sess).erp_target_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.target.erp_mu; % move
                    monk(m).sess(sess).erp_stop_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.stop.erp_mu; % move
                    monk(m).sess(sess).erp_reward_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.reward.erp_mu; % move
                end
            end
        end
        
        % average across channels MST
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                monk(m).sess(sess).erp_move_mu_ch_MST = nanmean(monk(m).sess(sess).erp_move_MST);
                monk(m).sess(sess).erp_target_mu_ch_MST = nanmean(monk(m).sess(sess).erp_target_MST);
                monk(m).sess(sess).erp_stop_mu_ch_MST = nanmean(monk(m).sess(sess).erp_stop_MST);
                monk(m).sess(sess).erp_reward_mu_ch_MST = nanmean(monk(m).sess(sess).erp_reward_MST);
            end
        end
        
        % average across channels PPC
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                monk(m).sess(sess).erp_move_mu_ch_PPC = nanmean(monk(m).sess(sess).erp_move_PPC);
                monk(m).sess(sess).erp_target_mu_ch_PPC = nanmean(monk(m).sess(sess).erp_target_PPC);
                monk(m).sess(sess).erp_stop_mu_ch_PPC = nanmean(monk(m).sess(sess).erp_stop_PPC);
                monk(m).sess(sess).erp_reward_mu_ch_PPC = nanmean(monk(m).sess(sess).erp_reward_PPC);
            end
        end
        
        %% plot one session
        figure; % move
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_move_MST, smooth(monk(m).sess(s).erp_move_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_move_PPC, smooth(monk(m).sess(s).erp_move_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 1],'ylim',[-30 30],'yTick',[-30 0 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP target monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        % target
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_target_MST, smooth(monk(m).sess(s).erp_target_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_target_PPC, smooth(monk(m).sess(s).erp_target_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 1],'ylim',[-30 30],'yTick',[-30 0 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP target monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        % stop
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_stop_MST, smooth(monk(m).sess(s).erp_stop_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_stop_PPC, smooth(monk(m).sess(s).erp_stop_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 0.5],'ylim',[-30 30],'yTick',[-30 0 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP target monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        % reward
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_reward_MST, smooth(monk(m).sess(s).erp_reward_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_reward_PPC, smooth(monk(m).sess(s).erp_reward_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 0.5],'ylim',[-30 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP target monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
    case 'erp_single_sess_reward_densities'
        type = 'reward'
        ncond = length(monk(1).erp.MST.sess(1).lfps(1).trialtype.(type));
        ts_move_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time;
        ts_target_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time;
        ts_stop_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time;
        ts_reward_MST = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time;
        ts_move_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time;
        ts_target_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time;
        ts_stop_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time;
        ts_reward_PPC = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time;
        
        % Mean all channels per session MST
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for cond = 1:ncond
                    for ch = 1:length(monk(1).cont.MST.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_target_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_stop_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_reward_MST(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu*1000; % move
                    end
                end
            end
        end
        % Mean all channels per session PPC
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for cond = 1:ncond
                    for ch = 1:length(monk(1).cont.PPC.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_target_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_stop_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_reward_PPC(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu; % move
                    end
                end
            end
        end
        
        % average across channels MST
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                for cond = 1:ncond
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch_MST = nanmean(monk(m).sess(sess).(type)(cond).erp_move_MST);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch_MST = nanmean(monk(m).sess(sess).(type)(cond).erp_target_MST);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch_MST = nanmean(monk(m).sess(sess).(type)(cond).erp_stop_MST);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch_MST = nanmean(monk(m).sess(sess).(type)(cond).erp_reward_MST);
                end
            end
        end
        
        % average across channels PPC
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                for cond = 1:ncond
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch_PPC = nanmean(monk(m).sess(sess).(type)(cond).erp_move_PPC);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch_PPC = nanmean(monk(m).sess(sess).(type)(cond).erp_target_PPC);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch_PPC = nanmean(monk(m).sess(sess).(type)(cond).erp_stop_PPC);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch_PPC = nanmean(monk(m).sess(sess).(type)(cond).erp_reward_PPC);
                end
            end
        end
        
        %% plot one session
        figure; % move
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_move_MST, smooth(monk(m).sess(s).(type)(cond).erp_move_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_move_PPC, smooth(monk(m).sess(s).(type)(cond).erp_move_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 1],'ylim',[-30 30],'yTick',[-30 0 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP target monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        % target
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_target_MST, smooth(monk(m).sess(s).(type)(cond).erp_target_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_target_PPC, smooth(monk(m).sess(s).(type)(cond).erp_target_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 1],'ylim',[-30 30],'yTick',[-30 0 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP target monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        % stop
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_stop_MST, smooth(monk(m).sess(s).(type)(cond).erp_stop_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_stop_PPC, smooth(monk(m).sess(s).(type)(cond).erp_stop_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 0.5],'ylim',[-30 30],'yTick',[-30 0 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP stop monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        % reward
        for m = 1:length(monk)
            for s = 1:length(monk(m).sess)
                figure; hold on;
                plot(ts_reward_MST, smooth(monk(m).sess(s).(type)(cond).erp_reward_mu_ch_MST),'-b', 'LineWidth', 1.5);
                plot(ts_reward_PPC, smooth(monk(m).sess(s).(type)(cond).erp_reward_mu_ch_PPC),'-c', 'LineWidth', 1.5);
                set(gca,'xlim',[-1 0.5],'ylim',[-30 30], 'TickDir', 'out', 'FontSize', 22); box off;
                ylabel('µV'); xlabel('Time (s)')
                title(['ERP reward monk ' num2str(m)]); legend({'MST', 'PPC'},'box', 'off');axis square
            end
        end
        
        
        
    case 'erp_all_MST'
        type = 'all'
        ts_move = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = [1 3]; %1:length(monk)
            for sess = 1:length(monk(m).pop)
                for ncond = 1
                    for ch = 1:length(monk(1).cont.MST.sess(sess).lfps)
                        monk(m).sess(sess).erp_move(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.move.erp_mu*1000; % move
                        monk(m).sess(sess).erp_target(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.target.erp_mu*1000; % move
                        monk(m).sess(sess).erp_stop(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.stop.erp_mu*1000; % move
                        monk(m).sess(sess).erp_reward(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type).events.reward.erp_mu*1000; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = [1 3] % 1:length(monk)
            for sess = 1:length(monk(m).sess)
                monk(m).sess(sess).erp_move_mu_ch = nanmean(monk(m).sess(sess).erp_move);
                monk(m).sess(sess).erp_target_mu_ch = nanmean(monk(m).sess(sess).erp_target);
                monk(m).sess(sess).erp_stop_mu_ch = nanmean(monk(m).sess(sess).erp_stop);
                monk(m).sess(sess).erp_reward_mu_ch = nanmean(monk(m).sess(sess).erp_reward);
            end
        end
        
        % average across sessions
        for m = [1 3] % 1:length(monk)
            clear th_v th_w bet_v bet_w
            for sess = 1:length(monk(m).sess)
                erp_move_sess(sess,:) =  monk(m).sess(sess).erp_move_mu_ch; ...
                    [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                erp_target_sess(sess,:) = monk(m).sess(sess).erp_target_mu_ch; ...
                    [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                erp_stop_sess(sess,:) = monk(m).sess(sess).erp_stop_mu_ch; ...
                    [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                erp_reward_sess(sess,:) = monk(m).sess(sess).erp_reward_mu_ch; ...
                    [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
            end
            % mean
            monk(m).erp.move = nanmean(erp_move_sess);     monk(m).erp.move_std = nanmean(nanstd(erp_move_sess));
            monk(m).erp.target = nanmean(erp_target_sess); monk(m).erp.target_std = nanmean(nanstd(erp_target_sess));
            monk(m).erp.stop = nanmean(erp_stop_sess);     monk(m).erp.stop_std = nanmean(nanstd(erp_stop_sess));
            monk(m).erp.reward = nanmean(erp_reward_sess); monk(m).erp.reward_std = nanmean(nanstd(erp_reward_sess));
            
            monk(m).erp.move_max = nanmean(max_move_sess); monk(m).erp.move_max_amp_std = nanstd(max_move_sess);
            monk(m).erp.target_max = nanmean(max_target_sess);  monk(m).erp.target_max_amp_std = nanstd(max_target_sess);
            monk(m).erp.stop_max = nanmean(max_stop_sess); monk(m).erp.stop_max_amp_std = nanstd(max_stop_sess);
            monk(m).erp.reward_max = nanmean(max_reward_sess); monk(m).erp.reward_max_amp_std = nanstd(max_reward_sess);
            
            monk(m).erp.move_max_time = nanmean(max_move_time); monk(m).erp.move_max_std = nanstd(max_move_time);
            monk(m).erp.target_max_time = nanmean(max_target_time);  monk(m).erp.target_max_std = nanstd(max_target_time);
            monk(m).erp.stop_max_time = nanmean(max_stop_time); monk(m).erp.stop_max_std = nanstd(max_stop_time);
            monk(m).erp.reward_max_time = nanmean(max_reward_time); monk(m).erp.reward_max_std = nanstd(max_reward_time);
        end
        
        
        
        %% plot
        %  move
        figure;
        for m = [1 3]; %1:length(monk)
            plot(ts_move, smooth(monk(m).erp.move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP move'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        %% target
        figure;
        for m = [1 3]; %1:length(monk)
            plot(ts_target, smooth(monk(m).erp.target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP target'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        % avg for all monks
        for m = [1 3];
            erp_targ(m,:) = monk(m).erp.target; % erp_targ(m,:) = abs(monk(m).erp.target); %
        end
        erp_targ = erp_targ([1 3],:);
        plot(ts_target-0.3, smooth(nanmean(erp_targ)),'r', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP target all');
        
        %% stop
        figure;
        for m = [1 3]; %1:length(monk)
            plot(ts_stop, smooth(monk(m).erp.stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP stop'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        % avg for all monks
        for m = [1 3];
            erp_stop(m,:) = monk(m).erp.stop; % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        erp_stop = erp_stop([1 3],:);
        plot(ts_stop, smooth(nanmean(erp_stop)),'r', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP stop all');
        
        %% reward
        figure;
        for m = [1 3]; %1:length(monk)
            plot(ts_reward, smooth(monk(m).erp.reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP reward'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        % avg for all monks
        for m = [1 3];
            erp_reward(m,:) = monk(m).erp.reward; % erp_reward(m,:) = abs(monk(m).erp.target); %
        end
        erp_reward = erp_reward([1 3],:);
        plot(ts_reward, smooth(nanmean(erp_reward)),'r', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward all');
        
        % plot max amp
        figure; hold on;
        errorbar(1.2,monk(1).erp.move_max, monk(1).erp.move_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(1).erp.target_max, monk(1).erp.target_max_amp_std,'ob', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(1).erp.stop_max, monk(1).erp.stop_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(1).erp.reward_max, monk(1).erp.reward_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1.2,monk(2).erp.move_max, monk(2).erp.move_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(2).erp.target_max, monk(2).erp.target_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(2).erp.stop_max, monk(2).erp.stop_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(2).erp.reward_max, monk(2).erp.reward_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[0 40], 'yTick',[0 20 40], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Max ERP (µV)'); title('MST'); axis square;
        
        
        % plot max time
        figure; hold on;
        errorbar(1,monk(1).erp.move_max_time, monk(1).erp.move_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2,monk(1).erp.target_max_time, monk(1).erp.target_max_std,'ob', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3,monk(1).erp.stop_max_time, monk(1).erp.stop_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4,monk(1).erp.reward_max_time, monk(1).erp.reward_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1,monk(2).erp.move_max_time, monk(2).erp.move_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2,monk(2).erp.target_max_time, monk(2).erp.target_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3,monk(2).erp.stop_max_time, monk(2).erp.stop_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4,monk(2).erp.reward_max_time, monk(2).erp.reward_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[-0.6 0.8], 'yTick',[-0.6 0 0.6], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Time from event (s)'); hline(0, '--k'); title('MST'); axis square;
        
    case 'erp_reward_densities_MST'
        type = 'reward'
        ncond = length(monk(1).erp.MST.sess(1).lfps(1).trialtype.(type));
        ts_move = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = [1 3] % 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for cond = 1:ncond
                    for ch = 1:length(monk(1).cont.MST.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_target(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_stop(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_reward(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu*1000; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = [1 3] %1:length(monk)
            for sess = 1:length(monk(m).sess)
                for cond = 1:ncond
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_move);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_target);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_stop);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_reward);
                end
            end
        end
        
        % average across sessions
        for m = [1 3] %1:length(monk)
            clear th_v th_w bet_v bet_w
            for cond = 1:ncond
                for sess = 1:length(monk(m).sess)
                    
                    erp_move_sess(sess,:) =  monk(m).sess(sess).(type)(cond).erp_move_mu_ch; ...
                        [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).(type)(cond).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                    erp_target_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_target_mu_ch; ...
                        [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).(type)(cond).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                    erp_stop_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_stop_mu_ch; ...
                        [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).(type)(cond).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                    erp_reward_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_reward_mu_ch; ...
                        [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).(type)(cond).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
                    
                end
                
                % mean
                monk(m).erp.(type)(cond).move = nanmean(erp_move_sess);     monk(m).erp.(type)(cond).move_std = nanmean(nanstd(erp_move_sess));
                monk(m).erp.(type)(cond).target = nanmean(erp_target_sess); monk(m).erp.(type)(cond).target_std = nanmean(nanstd(erp_target_sess));
                monk(m).erp.(type)(cond).stop = nanmean(erp_stop_sess);     monk(m).erp.(type)(cond).stop_std = nanmean(nanstd(erp_stop_sess));
                monk(m).erp.(type)(cond).reward = nanmean(erp_reward_sess); monk(m).erp.(type)(cond).reward_std = nanmean(nanstd(erp_reward_sess));
                
                monk(m).erp.(type)(cond).move_max = nanmean(max_move_sess);      monk(m).erp.(type)(cond).move_max_amp_std = nanstd(max_move_sess);
                monk(m).erp.(type)(cond).target_max = nanmean(max_target_sess);  monk(m).erp.(type)(cond).target_max_amp_std = nanstd(max_target_sess);
                monk(m).erp.(type)(cond).stop_max = nanmean(max_stop_sess);      monk(m).erp.(type)(cond).stop_max_amp_std = nanstd(max_stop_sess);
                monk(m).erp.(type)(cond).reward_max = nanmean(max_reward_sess);  monk(m).erp.(type)(cond).reward_max_amp_std = nanstd(max_reward_sess);
                
                monk(m).erp.(type)(cond).move_max_time = nanmean(max_move_time);      monk(m).erp.(type)(cond).move_max_std = nanstd(max_move_time);
                monk(m).erp.(type)(cond).target_max_time = nanmean(max_target_time);  monk(m).erp.(type)(cond).target_max_std = nanstd(max_target_time);
                monk(m).erp.(type)(cond).stop_max_time = nanmean(max_stop_time);      monk(m).erp.(type)(cond).stop_max_std = nanstd(max_stop_time);
                monk(m).erp.(type)(cond).reward_max_time = nanmean(max_reward_time);  monk(m).erp.(type)(cond).reward_max_std = nanstd(max_reward_time);
            end
        end
        
        
        %% plot
        %  move
        figure;
        for m = [1 3] %1:length(monk)
            for cond = 1:ncond
                plot(ts_move, smooth(monk(m).erp.(type)(cond).move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP move'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = [1 3];
            erp_move(m,:) = smooth(monk(m).erp.(type)(2).move); % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        erp_move = erp_move([1 3],:);
        figure; hold on;
        plot(ts_move,erp_move,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_move, smooth(nanmean(erp_move)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP move all'); axis square;
        %% target
        figure;
        for m = [1 3] %1:length(monk)
            for cond = 1:ncond
                plot(ts_target, smooth(monk(m).erp.(type)(cond).target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP target'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = [1 3];
            erp_targ(m,:) = monk(m).erp.(type)(2).target; % erp_targ(m,:) = abs(monk(m).erp.target); %
        end
        erp_targ = erp_targ([1 3],:);
        figure; hold on;
        plot(ts_target-0.3,erp_targ,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_target-0.3, smooth(nanmean(erp_targ)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP target all'); axis square;
        %% stop
        figure;
        for m = [1 3] %1:length(monk)
            for cond = 1:ncond
                plot(ts_stop, smooth(monk(m).erp.(type)(cond).stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP stop'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = [1 3];
            erp_stop(m,:) = monk(m).erp.(type)(2).stop; % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        erp_stop = erp_stop([1 3],:);
        figure; hold on;
        plot(ts_stop,erp_stop,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_stop, smooth(nanmean(erp_stop)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP stop all'); axis square;
        %% reward
        figure;
        for m = [1 3] %1:length(monk)
            for cond = 1:ncond
                plot(ts_reward, smooth(monk(m).erp.(type)(cond).reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP reward'); legend({'Monk Q rew', 'Monk S rew'}, 'box', 'off');
            end
        end
        
        % avg for all monks
        for m = [1 3];
            erp_reward(m,:) = monk(m).erp.(type)(2).reward; % erp_reward(m,:) = abs(monk(m).erp.target); %
        end
        erp_reward = erp_reward([1 3],:);
        figure; hold on;
        plot(ts_reward,erp_reward,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_reward, smooth(nanmean(erp_reward)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward all'); axis square;
        
        %% plot max amp
        cond = 2;
        figure; hold on;
        errorbar(1.2,monk(1).erp.(type)(cond).move_max, monk(1).erp.(type)(cond).move_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(1).erp.(type)(cond).target_max, monk(1).erp.(type)(cond).target_max_amp_std,'ob', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(1).erp.(type)(cond).stop_max, monk(1).erp.(type)(cond).stop_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(1).erp.(type)(cond).reward_max, monk(1).erp.(type)(cond).reward_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1.2,monk(2).erp.(type)(cond).move_max, monk(2).erp.(type)(cond).move_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(2).erp.(type)(cond).target_max, monk(2).erp.(type)(cond).target_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(2).erp.(type)(cond).stop_max, monk(2).erp.(type)(cond).stop_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(2).erp.(type)(cond).reward_max, monk(2).erp.(type)(cond).reward_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[0 40], 'yTick',[0 20 40], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Max ERP (µV)'); title('MST'); axis square;
        
        
        % plot max time
        figure; hold on;
        errorbar(1,monk(1).erp.(type)(cond).move_max_time, monk(1).erp.(type)(cond).move_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2,monk(1).erp.(type)(cond).target_max_time, monk(1).erp.(type)(cond).target_max_std,'ob', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3,monk(1).erp.(type)(cond).stop_max_time, monk(1).erp.(type)(cond).stop_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4,monk(1).erp.(type)(cond).reward_max_time, monk(1).erp.(type)(cond).reward_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1,monk(2).erp.(type)(cond).move_max_time, monk(2).erp.(type)(cond).move_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2,monk(2).erp.(type)(cond).target_max_time, monk(2).erp.(type)(cond).target_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3,monk(2).erp.(type)(cond).stop_max_time, monk(2).erp.(type)(cond).stop_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4,monk(2).erp.(type)(cond).reward_max_time, monk(2).erp.(type)(cond).reward_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[-0.6 0.8], 'yTick',[-0.6 0 0.6], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Time from event (s)'); hline(0, '--k'); title('MST'); axis square;
        
    case 'erp_before_after_move_MST'
        
        type = 'reward'
        
        ncond = length(monk(1).erp.MST.sess(1).lfps(1).trialtype.(type));
        ts_move = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = [1 3] % 1:length(monk)
            for sess = 1:length(monk(m).erp.MST.sess)
                for cond = [3 4]
                    for ch = 1:length(monk(m).erp.MST.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_target(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_stop(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu*1000; % move
                        monk(m).sess(sess).(type)(cond).erp_reward(ch,:) = monk(m).erp.MST.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu*1000; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = [1 3] %1:length(monk)
            for sess = 1:length(monk(m).erp.MST.sess)
                for cond = [3 4]
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_move);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_target);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_stop);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_reward);
                end
            end
        end
        
        % average across sessions
        for m = [1 3] %1:length(monk)
            clear th_v th_w bet_v bet_w
            for cond = [3 4]
                for sess = 1:length(monk(m).erp.MST.sess)
                    
                    erp_move_sess(sess,:) =  monk(m).sess(sess).(type)(cond).erp_move_mu_ch; ...
                        [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).(type)(cond).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                    erp_target_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_target_mu_ch; ...
                        [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).(type)(cond).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                    erp_stop_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_stop_mu_ch; ...
                        [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).(type)(cond).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                    erp_reward_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_reward_mu_ch; ...
                        [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).(type)(cond).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
                    
                end
                
                % mean
                monk(m).erp.(type)(cond).move = nanmean(erp_move_sess);     monk(m).erp.(type)(cond).move_std = nanmean(nanstd(erp_move_sess));
                monk(m).erp.(type)(cond).target = nanmean(erp_target_sess); monk(m).erp.(type)(cond).target_std = nanmean(nanstd(erp_target_sess));
                monk(m).erp.(type)(cond).stop = nanmean(erp_stop_sess);     monk(m).erp.(type)(cond).stop_std = nanmean(nanstd(erp_stop_sess));
                monk(m).erp.(type)(cond).reward = nanmean(erp_reward_sess); monk(m).erp.(type)(cond).reward_std = nanmean(nanstd(erp_reward_sess));
                
                monk(m).erp.(type)(cond).move_max = nanmean(max_move_sess);      monk(m).erp.(type)(cond).move_max_amp_std = nanstd(max_move_sess);
                monk(m).erp.(type)(cond).target_max = nanmean(max_target_sess);  monk(m).erp.(type)(cond).target_max_amp_std = nanstd(max_target_sess);
                monk(m).erp.(type)(cond).stop_max = nanmean(max_stop_sess);      monk(m).erp.(type)(cond).stop_max_amp_std = nanstd(max_stop_sess);
                monk(m).erp.(type)(cond).reward_max = nanmean(max_reward_sess);  monk(m).erp.(type)(cond).reward_max_amp_std = nanstd(max_reward_sess);
                
                monk(m).erp.(type)(cond).move_max_time = nanmean(max_move_time);      monk(m).erp.(type)(cond).move_max_std = nanstd(max_move_time);
                monk(m).erp.(type)(cond).target_max_time = nanmean(max_target_time);  monk(m).erp.(type)(cond).target_max_std = nanstd(max_target_time);
                monk(m).erp.(type)(cond).stop_max_time = nanmean(max_stop_time);      monk(m).erp.(type)(cond).stop_max_std = nanstd(max_stop_time);
                monk(m).erp.(type)(cond).reward_max_time = nanmean(max_reward_time);  monk(m).erp.(type)(cond).reward_max_std = nanstd(max_reward_time);
            end
        end
        
        
        %% plot
        %  move
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_move, smooth(monk(m).erp.(type)(cond).move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP move'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = [1 3];
            erp_move_before(m,:) = smooth(monk(m).erp.(type)(3).move);
            erp_move_after(m,:) = smooth(monk(m).erp.(type)(4).move);
        end
        erp_move_before = erp_move_before([1 3],:); erp_move_after = erp_move_after([1 3],:);
        figure; hold on;
        plot(ts_move,erp_move_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_move,erp_move_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_move, smooth(nanmean(erp_move_before)),'k', 'LineWidth', 2); hold on; plot(ts_move, smooth(nanmean(erp_move_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP move before/after'); axis square;
        %% target
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_target, smooth(monk(m).erp.(type)(cond).target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP target'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = [1 3];
            erp_targ_before(m,:) = monk(m).erp.(type)(3).target;
            erp_targ_after(m,:) = monk(m).erp.(type)(4).target;
        end
        erp_targ_before = erp_targ_before([1 3],:); erp_targ_after = erp_targ_after([1 3],:);
        figure; hold on;
        plot(ts_target-0.3,erp_targ_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_target-0.3,erp_targ_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_target-0.3, smooth(nanmean(erp_targ_before)),'k', 'LineWidth', 2); hold on; plot(ts_target-0.3, smooth(nanmean(erp_targ_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off; vline(-0.3,'k');
        title('ERP target before/after'); axis square;
        %% stop
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_stop, smooth(monk(m).erp.(type)(cond).stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP stop'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = [1 3];
            erp_stop_before(m,:) = smooth(monk(m).erp.(type)(3).stop);
            erp_stop_after(m,:) = smooth(monk(m).erp.(type)(4).stop);
        end
        erp_stop_before = erp_stop_before([1 3],:); erp_stop_after = erp_stop_after([1 3],:);
        figure; hold on;
        plot(ts_stop,erp_stop_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_stop,erp_stop_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_stop, smooth(nanmean(erp_stop_before)),'k', 'LineWidth', 2); hold on; plot(ts_stop, smooth(nanmean(erp_stop_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP stop before/after'); axis square;
        %% reward
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_reward, smooth(monk(m).erp.(type)(cond).reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP reward'); legend({'Monk Q rew', 'Monk S rew'}, 'box', 'off');
        %             end
        %         end
        
        % avg for all monks
        for m = [1 3];
            erp_reward_before(m,:) = smooth(monk(m).erp.(type)(3).reward);
            erp_reward_after(m,:) = smooth(monk(m).erp.(type)(4).reward);
        end
        erp_reward_before = erp_reward_before([1 3],:); erp_reward_after = erp_reward_after([1 3],:);
        figure; hold on;
        plot(ts_reward,erp_reward_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_reward,erp_reward_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_reward, smooth(nanmean(erp_reward_before)),'k', 'LineWidth', 2); hold on; plot(ts_reward, smooth(nanmean(erp_reward_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward before/after'); axis square;
        
        
    case 'erp_all_PPC'
        type = 'all'
        ts_move = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(1).erp.MST.sess(1).lfps(1).trialtype.(type).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        % Mean all channels per session
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for ch = 1:length(monk(1).cont.PPC.sess(sess).lfps)
                    monk(m).sess(sess).erp_move(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.move.erp_mu; % move
                    monk(m).sess(sess).erp_target(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.target.erp_mu; % move
                    monk(m).sess(sess).erp_stop(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.stop.erp_mu; % move
                    monk(m).sess(sess).erp_reward(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type).events.reward.erp_mu; % move
                end
            end
        end
        
        
        % average across channels
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                monk(m).sess(sess).erp_move_mu_ch = nanmean(monk(m).sess(sess).erp_move);
                monk(m).sess(sess).erp_target_mu_ch = nanmean(monk(m).sess(sess).erp_target);
                monk(m).sess(sess).erp_stop_mu_ch = nanmean(monk(m).sess(sess).erp_stop);
                monk(m).sess(sess).erp_reward_mu_ch = nanmean(monk(m).sess(sess).erp_reward);
            end
        end
        
        % average across sessions
        for m = 1:length(monk)
            clear th_v th_w bet_v bet_w
            for sess = 1:length(monk(m).sess)
                erp_move_sess(sess,:) =  monk(m).sess(sess).erp_move_mu_ch; ...
                    [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                erp_target_sess(sess,:) = monk(m).sess(sess).erp_target_mu_ch; ...
                    [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                erp_stop_sess(sess,:) = monk(m).sess(sess).erp_stop_mu_ch; ...
                    [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                erp_reward_sess(sess,:) = monk(m).sess(sess).erp_reward_mu_ch; ...
                    [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
            end
            % mean
            monk(m).erp.move = nanmean(erp_move_sess);     monk(m).erp.move_std = nanmean(nanstd(erp_move_sess));
            monk(m).erp.target = nanmean(erp_target_sess); monk(m).erp.target_std = nanmean(nanstd(erp_target_sess));
            monk(m).erp.stop = nanmean(erp_stop_sess);     monk(m).erp.stop_std = nanmean(nanstd(erp_stop_sess));
            monk(m).erp.reward = nanmean(erp_reward_sess); monk(m).erp.reward_std = nanmean(nanstd(erp_reward_sess));
            
            monk(m).erp.move_max = nanmean(max_move_sess); monk(m).erp.move_max_amp_std = nanstd(max_move_sess);
            monk(m).erp.target_max = nanmean(max_target_sess);  monk(m).erp.target_max_amp_std = nanstd(max_target_sess);
            monk(m).erp.stop_max = nanmean(max_stop_sess); monk(m).erp.stop_max_amp_std = nanstd(max_stop_sess);
            monk(m).erp.reward_max = nanmean(max_reward_sess); monk(m).erp.reward_max_amp_std = nanstd(max_reward_sess);
            
            monk(m).erp.move_max_time = nanmean(max_move_time); monk(m).erp.move_max_std = nanstd(max_move_time);
            monk(m).erp.target_max_time = nanmean(max_target_time);  monk(m).erp.target_max_std = nanstd(max_target_time);
            monk(m).erp.stop_max_time = nanmean(max_stop_time); monk(m).erp.stop_max_std = nanstd(max_stop_time);
            monk(m).erp.reward_max_time = nanmean(max_reward_time); monk(m).erp.reward_max_std = nanstd(max_reward_time);
        end
        
        
        %  move
        figure;
        for m = 1:length(monk)
            plot(ts_move, smooth(monk(m).erp.move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP move'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        %% target
        figure;
        for m = 1:length(monk)
            plot(ts_target, smooth(monk(m).erp.target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP target'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        % avg for all monks
        for m = [1 3];
            erp_targ(m,:) = monk(m).erp.target; % erp_targ(m,:) = abs(monk(m).erp.target); %
        end
        erp_targ = erp_targ([1 3],:);
        plot(ts_target-0.3, smooth(nanmean(erp_targ)),'r', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP target all');
        %% stop
        figure;
        for m = 1:length(monk)
            plot(ts_stop, smooth(monk(m).erp.stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP stop'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        % avg for all monks
        for m = [1 3];
            erp_stop(m,:) = monk(m).erp.stop; % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        erp_stop = erp_stop([1 3],:);
        plot(ts_stop, smooth(nanmean(erp_stop)),'r', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP stop all');
        %% reward
        figure;
        for m = 1:length(monk)
            plot(ts_reward, smooth(monk(m).erp.reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
            set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
            title('ERP reward'); legend({'Monk Q', 'Monk S'}, 'box', 'off');
        end
        % avg for all monks
        for m = [1 3];
            erp_reward(m,:) = monk(m).erp.reward; % erp_reward(m,:) = abs(monk(m).erp.target); %
        end
        erp_reward = erp_reward([1 3],:);
        plot(ts_reward, smooth(nanmean(erp_reward)),'r', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward all');
        %% plot max amp
        figure; hold on;
        errorbar(1.2,monk(1).erp.move_max, monk(1).erp.move_max_amp_std, 'oc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(1).erp.target_max, monk(1).erp.target_max_amp_std,'oc', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(1).erp.stop_max, monk(1).erp.stop_max_amp_std, 'oc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(1).erp.reward_max, monk(1).erp.reward_max_amp_std, 'oc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1.2,monk(2).erp.move_max, monk(2).erp.move_max_amp_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(2).erp.target_max, monk(2).erp.target_max_amp_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(2).erp.stop_max, monk(2).erp.stop_max_amp_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(2).erp.reward_max, monk(2).erp.reward_max_amp_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[0 40], 'yTick',[0 20 40], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Max ERP (µV)'); title('PPC'); axis square;
        
        % plot max time
        figure; hold on;
        errorbar(1.2,monk(1).erp.move_max_time, monk(1).erp.move_max_std, 'oc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(1).erp.target_max_time, monk(1).erp.target_max_std,'oc', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(1).erp.stop_max_time, monk(1).erp.stop_max_std, 'oc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(1).erp.reward_max_time, monk(1).erp.reward_max_std, 'oc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1.2,monk(2).erp.move_max_time, monk(2).erp.move_max_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(2).erp.target_max_time, monk(2).erp.target_max_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(2).erp.stop_max_time, monk(2).erp.stop_max_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(2).erp.reward_max_time, monk(2).erp.reward_max_std, 'sc','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[-0.6 0.8], 'yTick',[-0.6 0 0.6], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Time from event (s)'); hline(0,'--k'); title('PPC'); axis square;
        
    case 'erp_reward_densities_PPC'
        type = 'reward'
        ncond = length(monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type));
        ts_move = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for cond = 1:ncond
                    for ch = 1:length(monk(m).cont.PPC.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_target(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_stop(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_reward(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                for cond = 1:ncond
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_move);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_target);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_stop);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_reward);
                end
            end
        end
        
        % average across sessions
        for m = 1:length(monk)
            clear th_v th_w bet_v bet_w
            for cond = 1:ncond
                for sess = 1:length(monk(m).sess)
                    
                    erp_move_sess(sess,:) =  monk(m).sess(sess).(type)(cond).erp_move_mu_ch; ...
                        [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).(type)(cond).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                    erp_target_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_target_mu_ch; ...
                        [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).(type)(cond).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                    erp_stop_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_stop_mu_ch; ...
                        [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).(type)(cond).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                    erp_reward_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_reward_mu_ch; ...
                        [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).(type)(cond).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
                    
                end
                
                % mean
                monk(m).erp.(type)(cond).move = nanmean(erp_move_sess);     monk(m).erp.(type)(cond).move_std = nanmean(nanstd(erp_move_sess));
                monk(m).erp.(type)(cond).target = nanmean(erp_target_sess); monk(m).erp.(type)(cond).target_std = nanmean(nanstd(erp_target_sess));
                monk(m).erp.(type)(cond).stop = nanmean(erp_stop_sess);     monk(m).erp.(type)(cond).stop_std = nanmean(nanstd(erp_stop_sess));
                monk(m).erp.(type)(cond).reward = nanmean(erp_reward_sess); monk(m).erp.(type)(cond).reward_std = nanmean(nanstd(erp_reward_sess));
                
                monk(m).erp.(type)(cond).move_max = nanmean(max_move_sess);      monk(m).erp.(type)(cond).move_max_amp_std = nanstd(max_move_sess);
                monk(m).erp.(type)(cond).target_max = nanmean(max_target_sess);  monk(m).erp.(type)(cond).target_max_amp_std = nanstd(max_target_sess);
                monk(m).erp.(type)(cond).stop_max = nanmean(max_stop_sess);      monk(m).erp.(type)(cond).stop_max_amp_std = nanstd(max_stop_sess);
                monk(m).erp.(type)(cond).reward_max = nanmean(max_reward_sess);  monk(m).erp.(type)(cond).reward_max_amp_std = nanstd(max_reward_sess);
                
                monk(m).erp.(type)(cond).move_max_time = nanmean(max_move_time);      monk(m).erp.(type)(cond).move_max_std = nanstd(max_move_time);
                monk(m).erp.(type)(cond).target_max_time = nanmean(max_target_time);  monk(m).erp.(type)(cond).target_max_std = nanstd(max_target_time);
                monk(m).erp.(type)(cond).stop_max_time = nanmean(max_stop_time);      monk(m).erp.(type)(cond).stop_max_std = nanstd(max_stop_time);
                monk(m).erp.(type)(cond).reward_max_time = nanmean(max_reward_time);  monk(m).erp.(type)(cond).reward_max_std = nanstd(max_reward_time);
            end
        end
        
        
        %% plot
        %  move
        figure;
        for m = 1:length(monk)
            for cond = 1:ncond
                plot(ts_move, smooth(monk(m).erp.(type)(cond).move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP move'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 1:length(monk)
            erp_move(m,:) = smooth(monk(m).erp.(type)(2).move); % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        figure; hold on;
        plot(ts_move,erp_move,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_move, smooth(nanmean(erp_move)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP move all'); axis square;
        %% target
        figure;
        for m = 1:length(monk)
            for cond = 1:ncond
                plot(ts_target, smooth(monk(m).erp.(type)(cond).target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP target'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 1:length(monk)
            erp_targ(m,:) = monk(m).erp.(type)(2).target; % erp_targ(m,:) = abs(monk(m).erp.target); %
        end
        figure; hold on;
        plot(ts_target-0.3,erp_targ,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_target-0.3, smooth(nanmean(erp_targ)),'k', 'LineWidth', 2);
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP target all'); axis square;
        %% stop
        figure;
        for m = 1:length(monk)
            for cond = 1:ncond
                plot(ts_stop, smooth(monk(m).erp.(type)(cond).stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP stop'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 1:length(monk)
            erp_stop(m,:) = monk(m).erp.(type)(2).stop; % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        figure; hold on;
        plot(ts_stop,erp_stop,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_stop, smooth(nanmean(erp_stop)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off; set(gca,'YTick', [0 10]);
        title('ERP stop all'); axis square;
        %% reward
        figure;
        for m = 1:length(monk)
            for cond = 1:ncond
                plot(ts_reward, smooth(monk(m).erp.(type)(cond).reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP reward'); legend({'Monk Q rew', 'Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 1:length(monk)
            erp_reward(m,:) = monk(m).erp.(type)(2).reward; % erp_reward(m,:) = abs(monk(m).erp.target); %
        end
        figure; hold on;
        plot(ts_reward,erp_reward,'Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_reward, smooth(nanmean(erp_reward)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward all'); axis square;
        %% plot max amp
        cond = 2;
        figure; hold on;
        errorbar(1.2,monk(1).erp.(type)(cond).move_max, monk(1).erp.(type)(cond).move_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(1).erp.(type)(cond).target_max, monk(1).erp.(type)(cond).target_max_amp_std,'ob', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(1).erp.(type)(cond).stop_max, monk(1).erp.(type)(cond).stop_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(1).erp.(type)(cond).reward_max, monk(1).erp.(type)(cond).reward_max_amp_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1.2,monk(2).erp.(type)(cond).move_max, monk(2).erp.(type)(cond).move_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2.2,monk(2).erp.(type)(cond).target_max, monk(2).erp.(type)(cond).target_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3.2,monk(2).erp.(type)(cond).stop_max, monk(2).erp.(type)(cond).stop_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4.2,monk(2).erp.(type)(cond).reward_max, monk(2).erp.(type)(cond).reward_max_amp_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[0 40], 'yTick',[0 20 40], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Max ERP (µV)'); title('MST'); axis square;
        
        
        % plot max time
        figure; hold on;
        errorbar(1,monk(1).erp.(type)(cond).move_max_time, monk(1).erp.(type)(cond).move_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2,monk(1).erp.(type)(cond).target_max_time, monk(1).erp.(type)(cond).target_max_std,'ob', 'LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3,monk(1).erp.(type)(cond).stop_max_time, monk(1).erp.(type)(cond).stop_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4,monk(1).erp.(type)(cond).reward_max_time, monk(1).erp.(type)(cond).reward_max_std, 'ob','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        
        errorbar(1,monk(2).erp.(type)(cond).move_max_time, monk(2).erp.(type)(cond).move_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(2,monk(2).erp.(type)(cond).target_max_time, monk(2).erp.(type)(cond).target_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(3,monk(2).erp.(type)(cond).stop_max_time, monk(2).erp.(type)(cond).stop_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        errorbar(4,monk(2).erp.(type)(cond).reward_max_time, monk(2).erp.(type)(cond).reward_max_std, 'sb','LineWidth',1,'MarkerSize',14, 'Capsize',0);
        set(gca,'xlim',[0.5 4.5],'xTick',[1 2 3 4],'xTickLabel',[{'move' 'target' 'stop' 'reward'}],'yLim',[-0.6 0.8], 'yTick',[-0.6 0 0.6], 'TickDir', 'out', 'FontSize', 22); box off;
        ylabel('Time from event (s)'); hline(0, '--k'); title('MST'); axis square;
        
        
    case 'erp_before_after_move_PPC'
        
        type = 'reward'
        
        ncond = length(monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type));
        ts_move = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = 1:length(monk)
            for sess = 1:length(monk(m).erp.PPC.sess)
                for cond = [3 4]
                    for ch = 1:length(monk(m).erp.PPC.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_target(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_stop(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_reward(ch,:) = monk(m).erp.PPC.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = 1:length(monk)
            for sess = 1:length(monk(m).erp.PPC.sess)
                for cond = [3 4]
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_move);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_target);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_stop);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_reward);
                end
            end
        end
        
        % average across sessions
        for m = 1:length(monk)
            clear th_v th_w bet_v bet_w
            for cond = [3 4]
                for sess = 1:length(monk(m).erp.PPC.sess)
                    
                    erp_move_sess(sess,:) =  monk(m).sess(sess).(type)(cond).erp_move_mu_ch; ...
                        [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).(type)(cond).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                    erp_target_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_target_mu_ch; ...
                        [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).(type)(cond).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                    erp_stop_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_stop_mu_ch; ...
                        [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).(type)(cond).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                    erp_reward_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_reward_mu_ch; ...
                        [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).(type)(cond).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
                    
                end
                
                % mean
                monk(m).erp.(type)(cond).move = nanmean(erp_move_sess);     monk(m).erp.(type)(cond).move_std = nanmean(nanstd(erp_move_sess));
                monk(m).erp.(type)(cond).target = nanmean(erp_target_sess); monk(m).erp.(type)(cond).target_std = nanmean(nanstd(erp_target_sess));
                monk(m).erp.(type)(cond).stop = nanmean(erp_stop_sess);     monk(m).erp.(type)(cond).stop_std = nanmean(nanstd(erp_stop_sess));
                monk(m).erp.(type)(cond).reward = nanmean(erp_reward_sess); monk(m).erp.(type)(cond).reward_std = nanmean(nanstd(erp_reward_sess));
                
                monk(m).erp.(type)(cond).move_max = nanmean(max_move_sess);      monk(m).erp.(type)(cond).move_max_amp_std = nanstd(max_move_sess);
                monk(m).erp.(type)(cond).target_max = nanmean(max_target_sess);  monk(m).erp.(type)(cond).target_max_amp_std = nanstd(max_target_sess);
                monk(m).erp.(type)(cond).stop_max = nanmean(max_stop_sess);      monk(m).erp.(type)(cond).stop_max_amp_std = nanstd(max_stop_sess);
                monk(m).erp.(type)(cond).reward_max = nanmean(max_reward_sess);  monk(m).erp.(type)(cond).reward_max_amp_std = nanstd(max_reward_sess);
                
                monk(m).erp.(type)(cond).move_max_time = nanmean(max_move_time);      monk(m).erp.(type)(cond).move_max_std = nanstd(max_move_time);
                monk(m).erp.(type)(cond).target_max_time = nanmean(max_target_time);  monk(m).erp.(type)(cond).target_max_std = nanstd(max_target_time);
                monk(m).erp.(type)(cond).stop_max_time = nanmean(max_stop_time);      monk(m).erp.(type)(cond).stop_max_std = nanstd(max_stop_time);
                monk(m).erp.(type)(cond).reward_max_time = nanmean(max_reward_time);  monk(m).erp.(type)(cond).reward_max_std = nanstd(max_reward_time);
            end
        end
        
        
        %% plot
        %  move
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_move, smooth(monk(m).erp.(type)(cond).move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP move'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = 1:length(monk)
            erp_move_before(m,:) = smooth(monk(m).erp.(type)(3).move);
            erp_move_after(m,:) = smooth(monk(m).erp.(type)(4).move);
        end
        figure; hold on;
        plot(ts_move,erp_move_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_move,erp_move_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_move, smooth(nanmean(erp_move_before)),'k', 'LineWidth', 2); hold on; plot(ts_move, smooth(nanmean(erp_move_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP move before/after'); axis square;
        %% target
        %         figure;
        %         for m = [1 3] % 1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_target, smooth(monk(m).erp.(type)(cond).target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP target'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = 1:length(monk)
            erp_targ_before(m,:) = monk(m).erp.(type)(3).target;
            erp_targ_after(m,:) = monk(m).erp.(type)(4).target;
        end
        figure; hold on;
        plot(ts_target-0.3,erp_targ_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_target-0.3,erp_targ_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_target-0.3, smooth(nanmean(erp_targ_before)),'k', 'LineWidth', 2); hold on; plot(ts_target-0.3, smooth(nanmean(erp_targ_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off; vline(-0.3,'k');
        title('ERP target before/after'); axis square;
        %% stop
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_stop, smooth(monk(m).erp.(type)(cond).stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP stop'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = 1:length(monk)
            erp_stop_before(m,:) = smooth(monk(m).erp.(type)(3).stop);
            erp_stop_after(m,:) = smooth(monk(m).erp.(type)(4).stop);
        end
        figure; hold on;
        plot(ts_stop,erp_stop_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_stop,erp_stop_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_stop, smooth(nanmean(erp_stop_before)),'k', 'LineWidth', 2); hold on; plot(ts_stop, smooth(nanmean(erp_stop_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP stop before/after'); axis square;
        %% reward
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_reward, smooth(monk(m).erp.(type)(cond).reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP reward'); legend({'Monk Q rew', 'Monk S rew'}, 'box', 'off');
        %             end
        %         end
        
        % avg for all monks
        for m = 1:length(monk)
            erp_reward_before(m,:) = smooth(monk(m).erp.(type)(3).reward);
            erp_reward_after(m,:) = smooth(monk(m).erp.(type)(4).reward);
        end
        figure; hold on;
        plot(ts_reward,erp_reward_before,'Color',[0.5 0.5 0.5], 'LineWidth', 2); plot(ts_reward,erp_reward_after,'--','Color',[0.5 0.5 0.5], 'LineWidth', 2);
        plot(ts_reward, smooth(nanmean(erp_reward_before)),'k', 'LineWidth', 2); hold on; plot(ts_reward, smooth(nanmean(erp_reward_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward before/after'); axis square;
        
    case 'erp_reward_densities_PFC'
        type = 'reward' % Only Schro has recordings in PFC for now (Apr 2020)
        ncond = length(monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type));
        ts_move = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = 3 % 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for cond = 1:ncond
                    for ch = 1:length(monk(m).cont.PPC.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_target(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_stop(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_reward(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = 3 % 1:length(monk)
            for sess = 1:length(monk(m).sess)
                for cond = 1:ncond
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_move);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_target);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_stop);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_reward);
                end
            end
        end
        
        % average across sessions
        for m = 3 %1:length(monk)
            clear th_v th_w bet_v bet_w
            for cond = 1:ncond
                for sess = 1:length(monk(m).sess)
                    
                    erp_move_sess(sess,:) =  monk(m).sess(sess).(type)(cond).erp_move_mu_ch; ...
                        [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).(type)(cond).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                    erp_target_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_target_mu_ch; ...
                        [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).(type)(cond).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                    erp_stop_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_stop_mu_ch; ...
                        [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).(type)(cond).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                    erp_reward_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_reward_mu_ch; ...
                        [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).(type)(cond).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
                    
                end
                
                % mean
                monk(m).erp.(type)(cond).move = nanmean(erp_move_sess);     monk(m).erp.(type)(cond).move_std = nanmean(nanstd(erp_move_sess));
                monk(m).erp.(type)(cond).target = nanmean(erp_target_sess); monk(m).erp.(type)(cond).target_std = nanmean(nanstd(erp_target_sess));
                monk(m).erp.(type)(cond).stop = nanmean(erp_stop_sess);     monk(m).erp.(type)(cond).stop_std = nanmean(nanstd(erp_stop_sess));
                monk(m).erp.(type)(cond).reward = nanmean(erp_reward_sess); monk(m).erp.(type)(cond).reward_std = nanmean(nanstd(erp_reward_sess));
                
                monk(m).erp.(type)(cond).move_max = nanmean(max_move_sess);      monk(m).erp.(type)(cond).move_max_amp_std = nanstd(max_move_sess);
                monk(m).erp.(type)(cond).target_max = nanmean(max_target_sess);  monk(m).erp.(type)(cond).target_max_amp_std = nanstd(max_target_sess);
                monk(m).erp.(type)(cond).stop_max = nanmean(max_stop_sess);      monk(m).erp.(type)(cond).stop_max_amp_std = nanstd(max_stop_sess);
                monk(m).erp.(type)(cond).reward_max = nanmean(max_reward_sess);  monk(m).erp.(type)(cond).reward_max_amp_std = nanstd(max_reward_sess);
                
                monk(m).erp.(type)(cond).move_max_time = nanmean(max_move_time);      monk(m).erp.(type)(cond).move_max_std = nanstd(max_move_time);
                monk(m).erp.(type)(cond).target_max_time = nanmean(max_target_time);  monk(m).erp.(type)(cond).target_max_std = nanstd(max_target_time);
                monk(m).erp.(type)(cond).stop_max_time = nanmean(max_stop_time);      monk(m).erp.(type)(cond).stop_max_std = nanstd(max_stop_time);
                monk(m).erp.(type)(cond).reward_max_time = nanmean(max_reward_time);  monk(m).erp.(type)(cond).reward_max_std = nanstd(max_reward_time);
            end
        end
        
        
        %% plot
        %  move
        figure;
        for m = 3 %1:length(monk)
            for cond = 1:ncond
                plot(ts_move, smooth(monk(m).erp.(type)(cond).move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP move'); legend({'Monk unrew', 'Monk  rew', 'Monk  unrew','Monk  rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 3;
            erp_move(m,:) = smooth(monk(m).erp.(type)(2).move); % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_move(m,:) = monk(m).erp.(type)(2).move; % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        figure; hold on;
        for ii = 1:length(erp_move_sess(:,1)),plot(ts_move,smooth(erp_move_sess(ii,:)),'Color',[0.5 0.5 0.5], 'LineWidth', 2);end
        plot(ts_move, smooth(nanmean(erp_move)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off; set(gca,'YTick', [0 10]);
        title('ERP move all'); axis square;
        %% target
        figure;
        for m = 3 %1:length(monk)
            for cond = 1:ncond
                plot(ts_target, smooth(monk(m).erp.(type)(cond).target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP target'); legend({'Monk unrew', 'Monk rew', 'Monk unrew','Monk rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_targ(m,:) = monk(m).erp.(type)(2).target; % erp_targ(m,:) = abs(monk(m).erp.target); %
        end
        erp_targ = erp_targ(3,:);
        figure; hold on;
        for ii = 1:length(erp_target_sess(:,1)),plot(ts_target-0.3,smooth(erp_target_sess(ii,:)),'Color',[0.5 0.5 0.5], 'LineWidth', 2);end
        plot(ts_target-0.3, smooth(erp_targ),'k', 'LineWidth', 2);
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP target all'); axis square;
        %% stop
        figure;
        for m = 3 %1:length(monk)
            for cond = 1:ncond
                plot(ts_stop, smooth(monk(m).erp.(type)(cond).stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP stop'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_stop(m,:) = monk(m).erp.(type)(2).stop; % erp_stop(m,:) = abs(monk(m).erp.stop); %
        end
        figure; hold on;
        for ii = 1:length(erp_stop_sess(:,1)),plot(ts_stop,smooth(erp_stop_sess(ii,:)),'Color',[0.5 0.5 0.5], 'LineWidth', 2);end
        plot(ts_stop, smooth(nanmean(erp_stop)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off; set(gca,'YTick', [0 10]);
        title('ERP stop all'); axis square;
        %% reward
        figure;
        for m = 3 %1:length(monk)
            for cond = 1:ncond
                plot(ts_reward, smooth(monk(m).erp.(type)(cond).reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
                set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
                title('ERP reward'); legend({'Monk Q rew', 'Monk S rew'}, 'box', 'off');
            end
        end
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_reward(m,:) = monk(m).erp.(type)(2).reward; % erp_reward(m,:) = abs(monk(m).erp.target); %
        end
        figure; hold on;
        for ii = 1:length(erp_reward_sess(:,1)),plot(ts_reward,smooth(erp_reward_sess(ii,:)),'Color',[0.5 0.5 0.5], 'LineWidth', 2);end
        plot(ts_reward, smooth(nanmean(erp_reward)),'k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward all'); axis square;
        
    case 'erp_before_after_move_PFC'
        
        type = 'reward'
        
        ncond = length(monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type));
        ts_move = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.move.erp_time; ts_move_win = ts_move(ts_move>=-0.5 & ts_move<=0.5);
        ts_target = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.target.erp_time; ts_target_win = ts_target(ts_target>=-0.5 & ts_target<=0.5);
        ts_stop = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.stop.erp_time; ts_stop_win = ts_stop(ts_stop>=-0.5 & ts_stop<=0.5);
        ts_reward = monk(3).erp.PFC.sess(1).lfps(1).trialtype.(type)(ncond).events.reward.erp_time; ts_reward_win = ts_reward(ts_reward>=-0.5 & ts_reward<=0.5);
        
        % Mean all channels per session
        for m = 3 % 1:length(monk)
            for sess = 1:length(monk(m).erp.PPC.sess)
                for cond = [3 4]
                    for ch = 1:length(monk(m).erp.PPC.sess(sess).lfps)
                        monk(m).sess(sess).(type)(cond).erp_move(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.move.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_target(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.target.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_stop(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.stop.erp_mu; % move
                        monk(m).sess(sess).(type)(cond).erp_reward(ch,:) = monk(m).erp.PFC.sess(sess).lfps(ch).trialtype.(type)(cond).events.reward.erp_mu; % move
                    end
                end
            end
        end
        
        % average across channels
        for m = 3 %1:length(monk)
            for sess = 1:length(monk(m).erp.PPC.sess)
                for cond = [3 4]
                    monk(m).sess(sess).(type)(cond).erp_move_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_move);
                    monk(m).sess(sess).(type)(cond).erp_target_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_target);
                    monk(m).sess(sess).(type)(cond).erp_stop_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_stop);
                    monk(m).sess(sess).(type)(cond).erp_reward_mu_ch = nanmean(monk(m).sess(sess).(type)(cond).erp_reward);
                end
            end
        end
        
        % average across sessions
        for m = 3 % 1:length(monk)
            clear th_v th_w bet_v bet_w
            for cond = [3 4]
                for sess = 1:length(monk(m).erp.PFC.sess)
                    
                    erp_move_sess(sess,:) =  monk(m).sess(sess).(type)(cond).erp_move_mu_ch; ...
                        [max_move_sess(sess,:),indx_move] = max(abs(monk(m).sess(sess).(type)(cond).erp_move_mu_ch(1,ts_move>=-0.5 & ts_move<=0.5))); max_move_time(sess,:) = ts_move_win(indx_move);
                    erp_target_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_target_mu_ch; ...
                        [max_target_sess(sess,:),indx_target] = max(abs(monk(m).sess(sess).(type)(cond).erp_target_mu_ch(1,ts_target>=-0.5 & ts_target<=0.5))); max_target_time(sess,:) = ts_target_win(indx_target);
                    erp_stop_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_stop_mu_ch; ...
                        [max_stop_sess(sess,:),indx_stop] = max(abs(monk(m).sess(sess).(type)(cond).erp_stop_mu_ch(1,ts_stop>=-0.5 & ts_stop<=0.5))); max_stop_time(sess,:) = ts_stop_win(indx_stop);
                    erp_reward_sess(sess,:) = monk(m).sess(sess).(type)(cond).erp_reward_mu_ch; ...
                        [max_reward_sess(sess,:),indx_reward] = max(abs(monk(m).sess(sess).(type)(cond).erp_reward_mu_ch(1,ts_reward>=-0.5 & ts_reward<=0.5))); max_reward_time(sess,:) = ts_reward_win(indx_reward);
                    
                end
                
                % mean
                monk(m).erp.(type)(cond).move = nanmean(erp_move_sess);     monk(m).erp.(type)(cond).move_std = nanmean(nanstd(erp_move_sess));
                monk(m).erp.(type)(cond).target = nanmean(erp_target_sess); monk(m).erp.(type)(cond).target_std = nanmean(nanstd(erp_target_sess));
                monk(m).erp.(type)(cond).stop = nanmean(erp_stop_sess);     monk(m).erp.(type)(cond).stop_std = nanmean(nanstd(erp_stop_sess));
                monk(m).erp.(type)(cond).reward = nanmean(erp_reward_sess); monk(m).erp.(type)(cond).reward_std = nanmean(nanstd(erp_reward_sess));
                
                monk(m).erp.(type)(cond).move_max = nanmean(max_move_sess);      monk(m).erp.(type)(cond).move_max_amp_std = nanstd(max_move_sess);
                monk(m).erp.(type)(cond).target_max = nanmean(max_target_sess);  monk(m).erp.(type)(cond).target_max_amp_std = nanstd(max_target_sess);
                monk(m).erp.(type)(cond).stop_max = nanmean(max_stop_sess);      monk(m).erp.(type)(cond).stop_max_amp_std = nanstd(max_stop_sess);
                monk(m).erp.(type)(cond).reward_max = nanmean(max_reward_sess);  monk(m).erp.(type)(cond).reward_max_amp_std = nanstd(max_reward_sess);
                
                monk(m).erp.(type)(cond).move_max_time = nanmean(max_move_time);      monk(m).erp.(type)(cond).move_max_std = nanstd(max_move_time);
                monk(m).erp.(type)(cond).target_max_time = nanmean(max_target_time);  monk(m).erp.(type)(cond).target_max_std = nanstd(max_target_time);
                monk(m).erp.(type)(cond).stop_max_time = nanmean(max_stop_time);      monk(m).erp.(type)(cond).stop_max_std = nanstd(max_stop_time);
                monk(m).erp.(type)(cond).reward_max_time = nanmean(max_reward_time);  monk(m).erp.(type)(cond).reward_max_std = nanstd(max_reward_time);
            end
        end
        
        
        %% plot
        %  move
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_move, smooth(monk(m).erp.(type)(cond).move),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP move'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = 3 % 1:length(monk)
            erp_move_before(m,:) = smooth(monk(m).erp.(type)(3).move);
            erp_move_after(m,:) = smooth(monk(m).erp.(type)(4).move);
        end
        figure; hold on;
        plot(ts_move, smooth(nanmean(erp_move_before)),'k', 'LineWidth', 2); hold on; plot(ts_move, smooth(nanmean(erp_move_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP move before/after'); axis square;
        %% target
        %         figure;
        %         for m = [1 3] % 1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_target, smooth(monk(m).erp.(type)(cond).target),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP target'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_targ_before(m,:) = monk(m).erp.(type)(3).target;
            erp_targ_after(m,:) = monk(m).erp.(type)(4).target;
        end
        figure; hold on;
        plot(ts_target-0.3, smooth(nanmean(erp_targ_before)),'k', 'LineWidth', 2); hold on; plot(ts_target-0.3, smooth(nanmean(erp_targ_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off; vline(-0.3,'k');
        title('ERP target before/after'); axis square;
        %% stop
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_stop, smooth(monk(m).erp.(type)(cond).stop),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP stop'); legend({'Monk Q unrew', 'Monk Q rew', 'Monk S unrew','Monk S rew'}, 'box', 'off');
        %             end
        %         end
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_stop_before(m,:) = smooth(monk(m).erp.(type)(3).stop);
            erp_stop_after(m,:) = smooth(monk(m).erp.(type)(4).stop);
        end
        figure; hold on;
        plot(ts_stop, smooth(nanmean(erp_stop_before)),'k', 'LineWidth', 2); hold on; plot(ts_stop, smooth(nanmean(erp_stop_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP stop before/after'); axis square;
        %% reward
        %         figure;
        %         for m = [1 3] %1:length(monk)
        %             for cond = [3 4]
        %                 plot(ts_reward, smooth(monk(m).erp.(type)(cond).reward),'color',[1 3 2]== m, 'LineWidth', 2); hold on;
        %                 set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        %                 title('ERP reward'); legend({'Monk Q rew', 'Monk S rew'}, 'box', 'off');
        %             end
        %         end
        
        % avg for all monks
        for m = 3 %1:length(monk)
            erp_reward_before(m,:) = smooth(monk(m).erp.(type)(3).reward);
            erp_reward_after(m,:) = smooth(monk(m).erp.(type)(4).reward);
        end
        figure; hold on;
        plot(ts_reward, smooth(nanmean(erp_reward_before)),'k', 'LineWidth', 2); hold on; plot(ts_reward, smooth(nanmean(erp_reward_after)),'--k', 'LineWidth', 2); hold on;
        set(gca,'xlim',[-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); box off;
        title('ERP reward before/after'); axis square;
        
    case 'PSD_all'
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).pw.area);
            f = monk(1).pw.freq; % frequency
            for narea = 1:length(areas)
                psd = monk(nmonk).pw.area.(areas{narea}).all.mu_sess; psd_sem = monk(nmonk).pw.area.(areas{narea}).all.std_sess;
                
                figure; %plot(f,psd,'LineWidth',2,'Color','k');
                shadedErrorBar(f,psd,psd_sem, 'lineprops','k');
                %set(s.mainLine,'LineWidth', 1);
                %set(s.edge,'LineStyle', 'none');
                %set(s.patch, 'FaceAlpha', 0.1);
                xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)');
                set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
                title(['all' areas(narea) 'Monk' nmonk])
            end
        end
        
    case 'PSD_all_together'
        f = monk(1).pw.freq; % frequency
        for nmonk = 3 % 1:length(monk)  [1 3]
            psd_mst(nmonk,:) = monk(nmonk).pw.area.MST.all.mu_sess;
            psd_mst_sem(nmonk,:) = monk(nmonk).pw.area.MST.all.std_sess;
            psd_ppc(nmonk,:) = monk(nmonk).pw.area.PPC.all.mu_sess;
            psd_ppc_sem(nmonk,:) = monk(nmonk).pw.area.PPC.all.std_sess;
            psd_pfc(nmonk,:) = monk(nmonk).pw.area.PFC.all.mu_sess;
            psd_pfc_sem(nmonk,:) = monk(nmonk).pw.area.PFC.all.std_sess;
        end
        
        figure; %  MST
        shadedErrorBar(f,mean(psd_mst),mean(psd_mst_sem), 'lineprops','k');
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 0.08*10e-4])
        set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        title('MST')
        
        % Separate for each monkey
        figure; %  MST
        plot(f,psd_mst,'k');
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 0.08*10e-4])
        set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        title('MST')
        
        figure; %  PPC
        shadedErrorBar(f,mean(psd_ppc),mean(psd_ppc_sem), 'lineprops','k');
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 60])
        set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        title('PPC')
        
        %% PFC
        max_pfc = max(max(psd_pfc));
        figure;
        shadedErrorBar(f,mean(psd_pfc)./max_pfc,mean(psd_pfc_sem)./max_pfc, 'lineprops','k');
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 0.01])
        set(gca,'TickDir', 'out','yTick',[0 0.005 0.01], 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        title('PFC')
        
        
    case 'PSD_move'
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).pw.area);
            f = monk(1).pw.freq; % frequency
            for narea = 1:length(areas)
                psd1 = monk(nmonk).pw.area.(areas{narea}).stationary.mu_sess; psd1_sem = monk(nmonk).pw.area.(areas{narea}).stationary.std_sess; % stationary
                psd2 = monk(nmonk).pw.area.(areas{narea}).mobile.mu_sess; psd2_sem = monk(nmonk).pw.area.(areas{narea}).mobile.std_sess; % move
                %plot
                figure; subplot(1,2,1); hold on;
                shadedErrorBar(f,psd1,psd1_sem, 'lineprops','r');
                shadedErrorBar(f,psd2,psd2_sem, 'lineprops','k');
                xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('PSD (dB)'); %ylabel('Power spectral density (\muV^2/Hz)');
                set(gca,'TickDir', 'out', 'FontSize', 22); box off; %set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
                % ratio
                subplot(1,2,2); hold on; plot(psd1./psd2, 'k'); plot(nanmean(psd1./psd2), 'LineWidth',2, 'Color','k');
                axis([1 50 0 6]); hline(1,'k'); xlabel('Frequency (Hz)'); ylabel('PSD');
                box off; set(gca,'TickDir', 'out', 'FontSize', 22); title('motion');
                title(['Move vs stationary' areas(narea) 'monk ' nmonk])
                
                figure; hold on;
                % theta
                l_lim = 4 ; h_lim = 7.9;
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(1,nanmean(psd1_band),std(psd1_band),'or','MarkerFaceColor', 'r','LineWidth',1,'CapSize',0);
                errorbar(2,nanmean(psd2_band),std(psd2_band),'ok','MarkerFaceColor', 'k','LineWidth',1,'CapSize',0);
                % alpha
                l_lim = 8 ; h_lim = 11.9; % beta band
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(3,nanmean(psd1_band),std(psd1_band),'or','MarkerFaceColor', 'r','LineWidth',1,'CapSize',0);
                errorbar(4,nanmean(psd2_band),std(psd2_band),'ok','MarkerFaceColor', 'k','LineWidth',1,'CapSize',0);
                % Beta
                l_lim = 12 ; h_lim = 19; % beta band
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(5,nanmean(psd1_band),std(psd1_band),'or','MarkerFaceColor', 'r','LineWidth',1,'CapSize',0);
                errorbar(6,nanmean(psd2_band),std(psd2_band),'ok','MarkerFaceColor', 'k','LineWidth',1,'CapSize',0); ylabel('PSD');
                set(gca,'xlim', [0 6],'xTick',[1 3 5 ],'xTickLabel',[{'theta'},{'alpha'}, {'beta'}], 'TickDir', 'out', 'FontSize',22);
                title(['Monk ' areas(narea) num2str(nmonk)]); axis square;
                
            end
        end
        
    case 'PSD_move_together'
        f = monk(1).pw.freq; % frequency
        for nmonk = 3;  %1:length(monk) [1 3]
            %             psd1_mst(nmonk,:) = monk(nmonk).pw.area.MST.stationary.mu_sess;
            %             psd1_mst_sem(nmonk,:) = monk(nmonk).pw.area.MST.stationary.std_sess;
            %             psd2_mst(nmonk,:) = monk(nmonk).pw.area.MST.mobile.mu_sess;
            %             psd2_mst_sem(nmonk,:) = monk(nmonk).pw.area.MST.mobile.std_sess;
            %
            %             psd1_ppc(nmonk,:) = monk(nmonk).pw.area.PPC.stationary.mu_sess;
            %             psd1_ppc_sem(nmonk,:) = monk(nmonk).pw.area.PPC.stationary.std_sess;
            %             psd2_ppc(nmonk,:) = monk(nmonk).pw.area.PPC.mobile.mu_sess;
            %             psd2_ppc_sem(nmonk,:) = monk(nmonk).pw.area.PPC.mobile.std_sess;
            
            psd1_pfc(nmonk,:) = monk(nmonk).pw.area.PFC.stationary.mu_sess;
            psd1_pfc_sem(nmonk,:) = monk(nmonk).pw.area.PFC.stationary.std_sess;
            psd2_pfc(nmonk,:) = monk(nmonk).pw.area.PFC.mobile.mu_sess;
            psd2_pfc_sem(nmonk,:) = monk(nmonk).pw.area.PFC.mobile.std_sess;
            
        end
        
        %% MST
        %         figure;
        %         shadedErrorBar(f,mean(psd_mst),mean(psd_mst_sem), 'lineprops','k');
        %         xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 0.08*10e-4])
        %         set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        %         title('MST')
        
        figure; hold on;
        shadedErrorBar(f,nanmean(psd1_mst([1 3],:)), nanmean(psd1_mst_sem([1 3],:)),'lineprops','r');
        shadedErrorBar(f,nanmean(psd2_mst([1 3],:)), nanmean(psd2_mst_sem([1 3],:)),'lineprops','k');
        xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('PSD (\muV^2/Hz)'); % ylabel('Power spectral density (dB)'); %
        box off; set(gca,'xTick', [10 20 30 40 50],'TickDir', 'out', 'FontSize', 22); %axis square; %ylim([0 0.4]);
        
        % ratio
        figure;hold on;
        plot(f,nanmean(psd1_mst([1 3],:))./nanmean(psd2_mst([1 3],:)), 'LineWidth',2, 'Color','k');
        %         plot(mobile_ratio', 'Color', [0.5 0.5 0.5]); plot(nanmean(mobile_ratio), 'LineWidth',2, 'Color','k');
        axis([1 50 0 3]); hline(1,'k'); xlabel('Frequency (Hz)'); ylabel('PSD');
        box off; set(gca,'TickDir', 'out', 'FontSize', 22); title('motion');
        
        %% PPC
        figure;
        shadedErrorBar(f,mean(psd_ppc),mean(psd_ppc_sem), 'lineprops','k');
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 60])
        set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        title('PPC')
        
        %% PFC
        figure; hold on;
        max_psd = max(max(psd1_pfc([1 3],:)));
        shadedErrorBar(f,nanmean(psd1_pfc([1 3],:))./max_psd, nanmean(psd1_pfc_sem([1 3],:))./max_psd,'lineprops','r');
        shadedErrorBar(f,nanmean(psd2_pfc([1 3],:))./max_psd, nanmean(psd2_pfc_sem([1 3],:))./max_psd,'lineprops','k');
        xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('PSD (\muV^2/Hz)'); % ylabel('Power spectral density (dB)'); %
        box off; set(gca,'xTick', [10 20 30 40 50],'TickDir', 'out', 'FontSize', 22); %axis square; %ylim([0 0.4]);
        
        % ratio
        figure;hold on;
        plot(f,nanmean(psd1_pfc([1 3],:))./nanmean(psd2_pfc([1 3],:)), 'LineWidth',2, 'Color','k');
        %         plot(mobile_ratio', 'Color', [0.5 0.5 0.5]); plot(nanmean(mobile_ratio), 'LineWidth',2, 'Color','k');
        axis([1 50 0 3]); hline(1,'k'); xlabel('Frequency (Hz)'); ylabel('PSD');
        box off; set(gca,'TickDir', 'out', 'FontSize', 22); title('motion');
        
    case 'PSD_move_together_modulation'
        f = monk(1).pw.freq; % frequency
        for nmonk =[1 3] %1:length(monk)
            %             psd1_mst(nmonk,:) = monk(nmonk).pw.area.MST.stationary.mu_sess;
            %             psd1_mst_sem(nmonk,:) = monk(nmonk).pw.area.MST.stationary.std_sess;
            %             psd2_mst(nmonk,:) = monk(nmonk).pw.area.MST.mobile.mu_sess;
            %             psd2_mst_sem(nmonk,:) = monk(nmonk).pw.area.MST.mobile.std_sess;
            
            psd1_ppc(nmonk,:) = monk(nmonk).pw.area.PPC.stationary.mu_sess; psd1_ppc = psd1_ppc([1 3],:);
            psd1_ppc_sem(nmonk,:) = monk(nmonk).pw.area.PPC.stationary.std_sess; psd1_ppc_sem = psd1_ppc_sem([1 3],:);
            psd2_ppc(nmonk,:) = monk(nmonk).pw.area.PPC.mobile.mu_sess; psd2_ppc = psd2_ppc([1 3],:);
            psd2_ppc_sem(nmonk,:) = monk(nmonk).pw.area.PPC.mobile.std_sess;  psd2_ppc_sem = psd2_ppc_sem([1 3],:);
        end
        
        
        %         figure; %  MST
        %         shadedErrorBar(f,mean(psd_mst),mean(psd_mst_sem), 'lineprops','k');
        %         xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('Power spectral density (\muV^2/Hz)'); ylim([0 0.08*10e-4])
        %         set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        %         title('MST')
        
        figure; %  PPC
        plot(f,mean(psd1_ppc./psd2_ppc),'k'); % stationary./mobile   plot(f,mean(psd2_ppc)./mean(psd1_ppc),'k');
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('PSD modulation'); ylim([0.5 2.1])
        set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
        title('PPC'); hline(1,'--k'); axis square;
        
        figure; hold on;
        for ii = 1:size(psd1_ppc,1)
            plot(f,psd1_ppc(ii,:)./psd2_ppc(ii,:),'Color',[0.5 0.5 0.5]);
        end
        xlim([4 50]); xlabel('Frequency (Hz)'); ylabel('PSD modulation'); ylim([0.5 2.1])
        set(gca,'TickDir', 'out', 'FontSize', 22); box off;
        title('PPC'); hline(1,'--k'); axis square;
        
    case 'PSD_eye'
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).pw.area);
            f = monk(1).pw.freq_eye; % frequency
            for narea = 1:length(areas)
                psd1 = monk(nmonk).pw.area.(areas{narea}).eyesfixed.mu_sess; psd1_sem = monk(nmonk).pw.area.(areas{narea}).eyesfixed.std_sess; % stationary
                psd2 = monk(nmonk).pw.area.(areas{narea}).eyesfree.mu_sess; psd2_sem = monk(nmonk).pw.area.(areas{narea}).eyesfree.std_sess; % move
                %plot
                figure; subplot(1,2,1); hold on;
                shadedErrorBar(f,psd1,psd1_sem, 'lineprops','m');
                shadedErrorBar(f,psd2,psd2_sem, 'lineprops','c');
                xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('PSD (dB)'); % ylabel('Power spectral density (\muV^2/Hz)');
                set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
                % ratio
                subplot(1,2,2); hold on; plot(psd1./psd2, 'k'); plot(nanmean(psd1./psd2), 'LineWidth',2, 'Color','k');
                axis([1 50 0 6]); hline(1,'k'); xlabel('Frequency (Hz)'); ylabel('PSD');
                box off; set(gca,'TickDir', 'out', 'FontSize', 22); title('motion');
                title(['Eyes free vs eyes fixed' areas(narea) 'monk ' nmonk])
                
                figure; hold on;
                % theta
                l_lim = 4 ; h_lim = 7.9;
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(1,nanmean(psd1_band),std(psd1_band),'om','MarkerFaceColor', 'm','LineWidth',1,'CapSize',0);
                errorbar(2,nanmean(psd2_band),std(psd2_band),'oc','MarkerFaceColor', 'c','LineWidth',1,'CapSize',0);
                % alpha
                l_lim = 8 ; h_lim = 11.9; % beta band
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(3,nanmean(psd1_band),std(psd1_band),'om','MarkerFaceColor', 'm','LineWidth',1,'CapSize',0);
                errorbar(4,nanmean(psd2_band),std(psd2_band),'oc','MarkerFaceColor', 'c','LineWidth',1,'CapSize',0);
                % Beta
                l_lim = 12 ; h_lim = 19; % beta band
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(5,nanmean(psd1_band),std(psd1_band),'om','MarkerFaceColor', 'm','LineWidth',1,'CapSize',0);
                errorbar(6,nanmean(psd2_band),std(psd2_band),'oc','MarkerFaceColor', 'c','LineWidth',1,'CapSize',0); ylabel('PSD');
                set(gca,'xlim', [0 6],'xTick',[1 3 5 ],'xTickLabel',[{'theta'},{'alpha'}, {'beta'}], 'TickDir', 'out', 'FontSize',22);
                title(['Monk ' areas(narea) num2str(nmonk)]); axis square;
            end
        end
        
    case 'PSD_reward'
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(nmonk).spec.area);
            f = monk(nmonk).pw.freq; % frequency
            for narea = 1:length(areas)
                psd1 = monk(nmonk).pw.area.(areas{narea}).reward(1).mu_sess; psd1_sem = monk(nmonk).pw.area.(areas{narea}).reward(1).std_sess; % stationary
                psd2 = monk(nmonk).pw.area.(areas{narea}).reward(2).mu_sess; psd2_sem = monk(nmonk).pw.area.(areas{narea}).reward(2).std_sess; % move
                if nmonk == 2
                    psd1 = monk(nmonk).pw.area.(areas{narea}).reward(1).mu_sess; psd1_sem = monk(nmonk).pw.area.(areas{narea}).reward(1).std_sess; % stationary
                    psd2 = monk(nmonk).pw.area.(areas{narea}).reward(2).mu_sess; psd2_sem = monk(nmonk).pw.area.(areas{narea}).reward(2).std_sess; % move
                end
                %plot
                figure; subplot(1,2,1); hold on;
                shadedErrorBar(f,psd1,psd1_sem, 'lineprops','k');
                shadedErrorBar(f,psd2,psd2_sem, 'lineprops','g');
                xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('PSD (dB)'); % ylabel('Power spectral density (\muV^2/Hz)');
                set(gca,'TickDir', 'out', 'FontSize', 22); box off; % set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
                % ratio
                subplot(1,2,2); hold on; plot(psd1./psd2, 'k'); plot(nanmean(psd1./psd2), 'LineWidth',2, 'Color','k');
                axis([1 50 0 6]); hline(1,'k'); xlabel('Frequency (Hz)'); ylabel('PSD');
                box off; set(gca,'TickDir', 'out', 'FontSize', 22); title('motion');
                title(['Eyes free vs eyes fixed' areas(narea) 'monk ' nmonk])
                
                figure; hold on;
                % theta
                l_lim = 4 ; h_lim = 7.9;
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(1,nanmean(psd1_band),std(psd1_band),'ok','MarkerFaceColor', 'k','LineWidth',1,'CapSize',0);
                errorbar(2,nanmean(psd2_band),std(psd2_band),'og','MarkerFaceColor', 'g','LineWidth',1,'CapSize',0);
                % alpha
                l_lim = 8 ; h_lim = 11.9; % beta band
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(3,nanmean(psd1_band),std(psd1_band),'ok','MarkerFaceColor', 'k','LineWidth',1,'CapSize',0);
                errorbar(4,nanmean(psd2_band),std(psd2_band),'og','MarkerFaceColor', 'g','LineWidth',1,'CapSize',0);
                % Beta
                l_lim = 12 ; h_lim = 19; % beta band
                psd1_band = psd1(f>l_lim & f<h_lim); psd2_band = psd2(f>l_lim & f<h_lim);
                errorbar(5,nanmean(psd1_band),std(psd1_band),'ok','MarkerFaceColor', 'k','LineWidth',1,'CapSize',0);
                errorbar(6,nanmean(psd2_band),std(psd2_band),'og','MarkerFaceColor', 'g','LineWidth',1,'CapSize',0); ylabel('PSD');
                set(gca,'xlim', [0 6],'xTick',[1 3 5 ],'xTickLabel',[{'theta'},{'alpha'}, {'beta'}], 'TickDir', 'out', 'FontSize',22);
                title(['Monk ' areas(narea) num2str(nmonk)]); axis square;
            end
        end
        
    case 'PSD_reward'
        % cond 2 is correct trials
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).pw.area);
            f = monk(1).pw.freq; % frequency
            for narea = 1:length(areas)
                figure;
                for ncond = 1:length(monk(nmonk).pw.area.(areas{narea}).reward)
                    psd = monk(nmonk).pw.area.(areas{narea}).reward(ncond).mu_sess; psd_sem = monk(nmonk).pw.area.(areas{narea}).reward(ncond).std_sess;
                    %plot
                    shadedErrorBar(f,psd,psd_sem, 'lineprops','k');
                    xlim([2 50]); xlabel('Frequency (Hz)'); ylabel('PSD (dB)'); %ylabel('Power spectral density (\muV^2/Hz)');
                    set(gca,'TickDir', 'out', 'FontSize', 22); box off; %set(gca,'TickDir', 'out', 'FontSize', 22, 'YScale', 'log'); box off;
                    title(['Rewarded' areas(narea) (ncond) 'Monk' nmonk]);
                end
            end
        end
        
    case 'spectrogram_move'
        type = 'reward'; ev = 'move';
        ncond = length(monk(1).erp.PPC.sess(1).lfps(1).trialtype.(type));
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).spec.area);
            for narea = 1:length(areas)
                for cond = 1:ncond
                    %                 trialtype = fieldnames(monk(1).spec.area.MST);
                    
                    freq = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).freq_sess;
                    ts = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).ts_sess;
                    p_spectro = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).mu_sess;
                    % plot
                    figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev 'cond '  num2str(cond)]); colormap(winter);
                    imagesc(ts-1,freq,p_spectro,[0 0.2]); axis xy; colorbar;
                    set(gca,'xlim',[-0.25 0.25], 'ylim',[5 50], 'FontSize', 22)
                    xlabel('time (s)'); ylabel('frequency (Hz)')
                end
            end
        end
        
        %% reward or density
        
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).spec.area);
            for narea = 1:length(areas)
                for cond = 1:ncond
                    %                 trialtype = fieldnames(monk(1).spec.area.MST);
                    freq = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).freq_sess;
                    ts = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).ts_sess;
                    p_spectro = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).mu_sess;
                    % plot
                    figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev]);
                    imagesc(ts-1,freq,p_spectro,[0 0.2]); axis xy; colorbar;
                    set(gca,'xlim',[-0.25 0.25], 'ylim',[4 50], 'FontSize', 22)
                    xlabel('time (s)'); ylabel('frequency (Hz)')
                end
            end
        end
        
        
        
    case 'spectrogram_target'
        ev = 'target';
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).spec.area);
            for narea = 1:length(areas)
                %                 trialtype = fieldnames(monk(1).spec.area.MST);
                freq = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).freq_sess;
                ts = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).ts_sess;
                p_spectro = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).mu_sess;
                % plot
                figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev]);
                imagesc(ts-1,freq,p_spectro,[0 0.9]); axis xy; colorbar;
                set(gca,'xlim',[-0.1 0.25], 'ylim',[4 50], 'FontSize', 22)
                xlabel('time (s)'); ylabel('frequency (Hz)')
            end
        end
        
    case 'spectrogram_stop'
        ev = 'stop';
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(1).spec.area);
            for narea = 1:length(areas)
                %                 trialtype = fieldnames(monk(1).spec.area.MST);
                freq = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).freq_sess;
                ts = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).ts_sess;
                p_spectro = monk(nmonk).spec.area.(areas{narea}).all.events.(ev).mu_sess;
                % plot
                figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev]);
                imagesc(ts-1,freq,p_spectro,[0 0.9]); axis xy; colorbar;
                set(gca,'xlim',[-0.25 0.25], 'ylim',[4 50], 'FontSize', 22)
                xlabel('time (s)'); ylabel('frequency (Hz)')
            end
        end
        
        
    case 'spectrogram_reward_density'
        type = 'reward'
        ev = 'move'
        align_t = 1;  % 1.3 for target (so it's aligned to target offset
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(nmonk).spec.area);
            for narea = 1:length(areas)
                ncond = length(monk(nmonk).spec.area.(areas{narea}).(type));
                for cond = 1:ncond
                    freq = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).freq_sess;
                    ts = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).ts_sess-align_t;
                    p_spectro = monk(nmonk).spec.area.(areas{narea}).(type)(cond).events.(ev).mu_sess;
                    ts_win = ts(ts>-0.45 & ts<0.45);
                    p_spectro_win = p_spectro(freq>5 & freq<51,ts>-0.45 & ts<0.45);
                    max_spectro = max(max(p_spectro_win));
                    % plot
                    figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev ' rew ' num2str(cond)]);
                    colormap(copper);
                    imagesc(ts,freq,p_spectro,[0 0.3]); axis xy; colorbar;
                    set(gca,'xlim',[-0.25 0.25], 'ylim',[5 50], 'FontSize', 22)
                    xlabel('time (s)'); ylabel('frequency (Hz)');
                    
                    % plot win
                    figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev ' rew ' num2str(cond)]);
                    colormap(copper);
                    imagesc(ts_win,freq(freq>3 & freq<51),p_spectro_win/max_spectro); axis xy; colorbar;
                    set(gca,'xlim',[-0.4 0.25], 'ylim',[5 50], 'FontSize', 22)
                    xlabel('time (s)'); ylabel('frequency (Hz)'); %vline(-0.3,'b')
                    
                end
            end
        end
        
        
    case 'spectrogram_reward_density_all_theta'
        area = 'MST'
        type = 'reward' % LOOK FOR TIME OF ALIGN target off is -1.3
        ev = 'move'
        th = [4 12];
        bet = [12 20];
        align_t = 1;  % 1.3 for target (so it's aligned to target offset
        for nmonk = [1 3] % 1:length(monk) %; [1 3]
            freq = monk(nmonk).spec.area.(area).(type)(1).events.(ev).freq_sess;
            ts = monk(nmonk).spec.area.(area).(type)(1).events.(ev).ts_sess-align_t;
            p_spectro_err = monk(nmonk).spec.area.(area).(type)(1).events.(ev).mu_sess;
            p_spectro_corr = monk(nmonk).spec.area.(area).(type)(2).events.(ev).mu_sess;
            ts_win = ts(ts>-0.31 & ts<0.31);
            p_spectro_win = p_spectro_corr(freq>4 & freq<51,ts>-0.31 & ts<0.31);
            % plot average per band across time
            theta_mu_err(nmonk,:) = nanmean(p_spectro_err(th(1):th(2),:)); theta_sem_err(nmonk,:) = nanstd(p_spectro_err(th(1):th(2),:))/sqrt(size(th(1):th(2),2));
            %theta_mu_corr(nmonk,:) = nanmean(p_spectro_corr(th(1):th(2),:)); theta_sem_corr(nmonk,:) = nanstd(p_spectro_corr(th(1):th(2),:))/sqrt(size(th(1):th(2),2));
            theta_mu_corr(nmonk,:) = nanmean(p_spectro_win); theta_sem_corr(nmonk,:) = nanstd(p_spectro_win)/sqrt(size(th(1):th(2),2));
        end
        
        col_mag = [0.3 0 0.3;
            0.45 0 0.139;
            0.128,0,0.128];
        
        figure; hold on;
        for nmonk = 1:length(monk)
            % hold on; s_1 = shadedErrorBar(ts, theta_mu_err(nmonk,:),theta_sem_err(nmonk,:),'lineprops', {'color',[col_mag(nmonk,:)]});
            hold on; s_2 = shadedErrorBar(ts_win, theta_mu_corr(nmonk,:)/max(theta_mu_corr(nmonk,:)),theta_sem_corr(nmonk,:),'lineprops', {'-','color',[col_mag(nmonk,:)]});
            % set(s_1.edge, 'LineStyle', 'none');   set(s_1.patch,'FaceAlpha', 0.1);
            set(s_2.edge, 'LineStyle', 'none'); set(s_2.patch,'FaceAlpha', 0.1);
        end
        set(gca,'xlim',[-0.25 0.25], 'FontSize', 22, 'TickDir', 'out'); axis square;
        ylim([0.6 1]);
        
    case 'spectrogram_reward_density_all_theta_diff'
        area = 'PFC'
        type = 'reward' % LOOK FOR TIME OF ALIGN target off is -1.3
        ev = 'stop'
        th = [4 12];
        bet = [12 20];
        align_t = 1;  % 1.3 for target (so it's aligned to target offset
        for nmonk = 3  % 1:length(monk) %; [1 3]
            freq = monk(nmonk).spec.area.(area).(type)(1).events.(ev).freq_sess;
            ts = monk(nmonk).spec.area.(area).(type)(1).events.(ev).ts_sess;
            p_spectro_err = monk(nmonk).spec.area.(area).(type)(1).events.(ev).mu_sess;
            p_spectro_corr = monk(nmonk).spec.area.(area).(type)(2).events.(ev).mu_sess;
            % plot average per band across time
            theta_mu_err(nmonk,:) = nanmean(p_spectro_err(th(1):th(2),:)); theta_sem_err(nmonk,:) = nanstd(p_spectro_err(th(1):th(2),:))/sqrt(size(th(1):th(2),2));
            theta_mu_corr(nmonk,:) = nanmean(p_spectro_corr(th(1):th(2),:)); theta_sem_corr(nmonk,:) = nanstd(p_spectro_corr(th(1):th(2),:))/sqrt(size(th(1):th(2),2));
        end
        
        col_mag = [0.3 0 0.3;
            0.45 0 0.139;
            0.128,0,0.128];
        
        figure; hold on;
        for nmonk = 1:length(monk)
            plot(ts-align_t, theta_mu_err(nmonk,:)- theta_mu_corr(nmonk,:),'Color',[col_mag(nmonk,:)]);
        end
        set(gca,'xlim',[-0.25 0.25], 'FontSize', 22, 'TickDir', 'out'); axis square;
        hline(0,'--k')
        
    case 'spectrogram_reward_density_all_beta'
        area = 'PFC'
        type = 'reward' % LOOK FOR TIME OF ALIGN target off is -1.3 % reward (separate by rew or unrew) or densities
        ev = 'move'
        th = [4 12];
        bet = [12 20];
        align_t = 1;
        for nmonk = 3  % 1:length(monk) % [1 3] %
            freq = monk(nmonk).spec.area.(area).(type)(1).events.(ev).freq_sess;
            ts = monk(nmonk).spec.area.(area).(type)(1).events.(ev).ts_sess-align_t;
            p_spectro_err = monk(nmonk).spec.area.(area).(type)(1).events.(ev).mu_sess;
            p_spectro_corr = monk(nmonk).spec.area.(area).(type)(2).events.(ev).mu_sess;
            ts_win = ts(ts>-0.31 & ts<0.31);
            p_spectro_win = p_spectro_corr(freq>4 & freq<51,ts>-0.31 & ts<0.31);
            % plot average per band across time
            beta_mu_err(nmonk,:) = nanmean(p_spectro_err(bet(1):bet(2),:)); beta_sem_err(nmonk,:) = nanstd(p_spectro_err(bet(1):bet(2),:))/sqrt(size(bet(1):bet(2),2));
            beta_mu_corr(nmonk,:) = nanmean(p_spectro_win(bet(1):bet(2),:)); beta_sem_corr(nmonk,:) = nanstd(p_spectro_win(bet(1):bet(2),:))/sqrt(size(bet(1):bet(2),2));
        end
        
        col_blue = [0 0.5 0.8;
            0 0.9 1;
            0 0 1];
        
        figure; hold on;
        for nmonk = 1:length(monk)
            % hold on; s_1 = shadedErrorBar(ts, theta_mu_err(nmonk,:),theta_sem_err(nmonk,:),'lineprops', {'color',[col_mag(nmonk,:)]});
            hold on; s_2 = shadedErrorBar(ts_win, beta_mu_corr(nmonk,:)/max(beta_mu_corr(nmonk,:)),beta_sem_corr(nmonk,:),'lineprops', {'-','color',[col_blue(nmonk,:)]});
            % set(s_1.edge, 'LineStyle', 'none');   set(s_1.patch,'FaceAlpha', 0.1);
            set(s_2.edge, 'LineStyle', 'none'); set(s_2.patch,'FaceAlpha', 0.1);
        end
        set(gca,'xlim',[-0.25 0.25], 'FontSize', 22, 'TickDir', 'out'); axis square;
        ylim([0.8 1.01]);
        
    case 'spectrogram_reward_density_all_beta_diff'
        area = 'PFC'
        type = 'reward' % LOOK FOR TIME OF ALIGN target off is -1.3
        ev = 'stop'
        th = [4 12];
        bet = [12 20];
        align_t = 1;
        for nmonk = 3   % 1:length(monk) % [1 3]
            freq = monk(nmonk).spec.area.(area).(type)(1).events.(ev).freq_sess;
            ts = monk(nmonk).spec.area.(area).(type)(1).events.(ev).ts_sess;
            p_spectro_err = monk(nmonk).spec.area.(area).(type)(1).events.(ev).mu_sess;
            p_spectro_corr = monk(nmonk).spec.area.(area).(type)(2).events.(ev).mu_sess;
            % plot average per band across time
            beta_mu_err(nmonk,:) = nanmean(p_spectro_err(bet(1):bet(2),:)); beta_sem_err(nmonk,:) = nanstd(p_spectro_err(bet(1):bet(2),:))/sqrt(size(bet(1):bet(2),2));
            beta_mu_corr(nmonk,:) = nanmean(p_spectro_corr(bet(1):bet(2),:)); beta_sem_corr(nmonk,:) = nanstd(p_spectro_corr(bet(1):bet(2),:))/sqrt(size(bet(1):bet(2),2));
        end
        
        col_blue = [0 0.5 0.8;
            0 0.9 1;
            0 0 1];
        
        figure; hold on;
        for nmonk = 1:length(monk)
            plot(ts-align_t, beta_mu_err(nmonk,:)-beta_mu_corr(nmonk,:),'Color',[col_blue(nmonk,:)]);
        end
        set(gca,'xlim',[-0.25 0.25], 'FontSize', 22, 'TickDir', 'out'); axis square;
        hline(0,'--k')
        
        
    case 'spectrogram_reward_density_diff'
        type = 'reward'
        ev = 'target'
        th = [4 12];
        bet = [12 20];
        win = [-0.25 0.25];
        for nmonk = 1:length(monk)
            areas = fieldnames(monk(nmonk).spec.area);
            for narea = 1:length(areas)
                %                 trialtype = fieldnames(monk(1).spec.area.MST);
                freq = monk(nmonk).spec.area.(areas{narea}).(type)(1).events.(ev).freq_sess;
                ts = monk(nmonk).spec.area.(areas{narea}).(type)(1).events.(ev).ts_sess;
                p_spectro_err = monk(nmonk).spec.area.(areas{narea}).(type)(1).events.(ev).mu_sess;
                p_spectro_corr = monk(nmonk).spec.area.(areas{narea}).(type)(2).events.(ev).mu_sess;
                
                % plot average per band across time
                theta_mu_err = nanmean(p_spectro_err(th(1):th(2),:)); theta_sem_err = nanstd(p_spectro_err(th(1):th(2),:))/sqrt(size(th(1):th(2),2));
                theta_mu_corr = nanmean(p_spectro_corr(th(1):th(2),:)); theta_sem_corr = nanstd(p_spectro_corr(th(1):th(2),:))/sqrt(size(th(1):th(2),2));
                beta_mu_err = nanmean(p_spectro_err(bet(1):bet(2),:)); beta_sem_err = nanstd(p_spectro_err(bet(1):bet(2),:))/sqrt(size(bet(1):bet(2),2));
                beta_mu_corr = nanmean(p_spectro_corr(bet(1):bet(2),:)); beta_sem_corr = nanstd(p_spectro_corr(bet(1):bet(2),:))/sqrt(size(bet(1):bet(2),2));
                %% plot
                %                 figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev]); colormap(winter); hold on;
                %                 %               min_p = min(min(p_spectro_err - p_spectro_corr)); max_p = max(max(p_spectro_err - p_spectro_corr));
                %                 %               imagesc(ts-1,freq,(p_spectro_err - p_spectro_corr),[min_p max_p]); axis xy; colorbar;
                %                 pcolor(ts-1,freq,(p_spectro_err - p_spectro_corr)); shading interp;
                %                 set(gca,'xlim',[-0.25 0.25], 'ylim',[5 50])
                %% plot
                %               plot(ts,p_spectro_err - p_spectro_corr, 'k'); alpha(0.5)
                figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev]); hold on;
                plot(freq,(nanmean(p_spectro_err,2) - nanmean(p_spectro_corr,2)));
                set(gca,'xlim',[4 50],'TickDir','out', 'FontSize', 22)
                xlabel('Freq (Hz)'); ylabel('Diff in power'); title('error - correct')
                %% plot colormap
                figure('Name',['Monk ' num2str(nmonk) ' ' char(areas(narea)) ' ' ev]); colormap(winter); hold on;
                imagesc(ts-1,freq,(p_spectro_err - p_spectro_corr),[-0.03 0.03]); colorbar;   % [-13e-8 6e-8]
                set(gca,'xlim',[-0.25 0.25], 'ylim',[5 50], 'FontSize', 22)
                xlabel('time (s)'); ylabel('frequency (Hz)'); title('error - correct')
                
                % plot diff as a line over time
                
                
                
            end
        end
        
    case 'spectrogram_trl'
        % aligned to target onset
        type = 'reward'
        delta = [0.5 4];
        theta = [4 12];
        low_beta = [12 20];
        high_beta = [20 30];
        for nmonk = 1 %1:length(monk)
            for nsess = 1:length(monk(nmonk).sess)
                areas = fieldnames(monk(nmonk).sess(nsess).trialtype.(type)(1).area);
                for narea = 1:length(areas)
                    for cond = 1:2%length(monk(nmonk).sess(nsess).trialtype.(type))
                        % get trl indx and behavior marker
                        clear trl_end move_on t_sacc t_fix t_rew pad_size_trials ts
                        indx_trls = monk(nmonk).behavior(nsess).stats.trialtype.(type)(cond).trlindx; trl_num = find(indx_trls);
                        for j = 1:length(trl_num)
                            trl_end(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_end;
                            move_on(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_move;
                            t_sacc{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_sac;
                            t_fix{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_fix;
                            t_rew(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_rew;
                            pad_size_trials(j) = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(j).spectro,1);
                        end
                        %% find max trials, pad with Nans and then sort
                        [~,indx_max] = max(pad_size_trials); % find the longest trial
                        pad_size = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(indx_max).spectro,1); % find how big it is
                        freq = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(1).freq;
                        ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(indx_max).ts-1;
                        % Pad with NaNs to compensate for different sized vectors
                        for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl)
                            % find location of t_rew
                            %                             if cond >= 2
                            %                                 this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts;
                            %                                 indx_rew = this_ts(this_ts<=t_rew(ntrl)+0.250); indx_rew = length(indx_rew);
                            %                                 monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro(indx_rew+1:pad_size,:) = NaN;
                            %                             else
                            this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts;
                            indx_end = this_ts(this_ts<=trl_end(ntrl)+0.250); indx_end = length(indx_end);
                            this_spectro = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro;
                            this_spectro(indx_end+1:pad_size,:) = NaN;
                            %                             end
                            padded_spectro(ntrl).spectro = this_spectro;
                        end
                        % sort spectrograms
                        
                        [trl_end_sort, trl_end_sort_indx] = sort(trl_end); % sort the trials
                        spectro_sorted = padded_spectro(trl_end_sort_indx)';
                        % extract different bands, average and normalize by max
                        clear de th be_low be_high max_de max_th max_be_low max_be_high de_norm th_norm be_low_norm be_high_norm
                        for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl)
                            de(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=delta(1) & freq<=delta(2)),2); max_de(ntrl) = max(de(:,ntrl)); de_norm(:,ntrl) = (de(:,ntrl)./max_de(ntrl));
                            th(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=theta(1) & freq<=theta(2)),2); max_th(ntrl) = max(th(:,ntrl)); th_norm(:,ntrl) = (th(:,ntrl)./max_th(ntrl));
                            be_low(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=low_beta(1) & freq<=low_beta(2)),2); max_be_low(ntrl) = max(be_low(:,ntrl)); be_low_norm(:,ntrl) = (be_low(:,ntrl)./max_be_low(ntrl));
                            be_high(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=high_beta(1) & freq<=high_beta(2)),2); max_be_high(ntrl) = max(be_high(:,ntrl)); be_high_norm(:,ntrl) = (be_high(:,ntrl)./max_be_high(ntrl));
                        end
                        max_bands = max([max(max_th) max(max_be_low) max(max_be_high)]);
                        
                        % plot
                        figure('Position',[1838 295 1679 381],'Name',['Session ', num2str(nsess)]); hold on; colormap(winter);
                        %% plot delta
                        % figure; hold on; colormap(winter);
                        subplot(1,4,1); hold on;
                        imagesc(ts,1:size(de_norm,2),de_norm'); colorbar;
                        if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(th,2),(th./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(th,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.25 4],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['delta M' num2str(nmonk) ' area ' (areas{narea}) ' cond ' num2str(cond)])
                        %ylabel('trial number'); xlabel('Time (s)'); axis square
                        %% plot theta
                        %figure; hold on; colormap(winter);
                        subplot(1,4,2);hold on;
                        imagesc(ts,1:size(th_norm,2),th_norm'); colorbar;
                        if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(th,2),(th./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(th,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.25 4],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['Theta M' num2str(nmonk) ' area ' (areas{narea}) ' cond ' num2str(cond)])
                        %ylabel('trial number'); xlabel('Time (s)'); axis square
                        %% plot low beta
                        %figure; hold on; colormap(winter);
                        subplot(1,4,3);hold on;
                        imagesc(ts,1:size(be_low_norm,2),be_low_norm'); colorbar;
                        if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(be_low,2),(be_low./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(be_low,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.25 4],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['Low-B M' num2str(nmonk) ' area ' (areas{narea}) ' cond ' num2str(cond)])
                        %ylabel('trial number'); xlabel('Time (s)');axis square
                        %% plot high beta
                        %figure; hold on; colormap(winter);
                        subplot(1,4,4);hold on;
                        imagesc(ts,1:size(be_high_norm,2),be_high_norm'); colorbar;
                        if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(be_high,2),(be_high./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(be_high,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.25 4],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['High-B M' num2str(nmonk) ' area ' (areas{narea}) ' cond ' num2str(cond)])
                        %ylabel('trial number'); xlabel('Time (s)'); axis square
                        
                        % plot by trial to inspect
                        figure('Name',['Monk ' num2str(nmonk) ' area ' (areas{narea}) ' cond ' num2str(cond) ' sess ' num2str(nsess)]); hold on;
                        %p = numSubplots(50);
                        p = numSubplots(size(th_norm,2));
                        % pick_trls = randperm(size(th_norm,2)); cnt=1;
                        % for ntrl = pick_trls(1:50)
                        for ntrl = 1:size(th_norm,2)
                            freq = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).freq;
                            ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts-1;
                            spectro = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro';
                            subplot(p(1),p(2),ntrl); %colormap(winter);
                            imagesc(ts,freq, spectro); axis xy; %colorbar;
                            set(gca,'xlim',[-0.25 ts(end)],'ylim',[0 50], 'FontSize', 22)
                            %xlabel('time (s)'); ylabel('frequency (Hz)');
                            %title(['Monkey ' num2str(nmonk) ' area ' (areas{narea}) ' trial ' num2str(ntrl) ' cond ' num2str(cond)])
                            hline([4 12 20 30],'-w'); vline(0,'-w'); vline(0.3,'-w'); %vline(move_on(ntrl),'-r');
                            %if cond == 2, vline(t_rew,'-r'), end
                            axis off;
                            % cnt=cnt+1;
                        end
                        
                    end
                end
            end
        end
        
    case 'spectrogram_trl_align_stop'
        % aligned to target onset
        type = 'reward'
        delta = [0.5 4];
        theta = [4 12];
        low_beta = [12 20];
        high_beta = [20 30];
        th_norm_monk_corr=[]; th_norm_monk_err=[]; be_low_norm_monk_corr=[]; be_low_norm_monk_err=[];
        for nmonk = 1:length(monk)  % [1 3]
            areas = {'PPC'} %fieldnames(monk(nmonk).sess(1).trialtype.(type)(1).area);
            for narea = 1:length(areas)
                th_norm_all_err = []; be_low_norm_all_err = []; th_norm_all_corr = []; be_low_norm_all_corr = [];
                for nsess = 1:length(monk(nmonk).sess)
                    for cond = 1:2%length(monk(nmonk).sess(nsess).trialtype.(type))
                        % get trl indx and behavior marker
                        clear trl_end move_on t_sacc t_fix t_rew pad_size_trials ts
                        indx_trls = monk(nmonk).behavior(nsess).stats.trialtype.(type)(cond).trlindx; trl_num = find(indx_trls);
                        for j = 1:length(trl_num)
                            trl_end(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_end;
                            move_on(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_move;
                            t_sacc{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_sac;
                            t_fix{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_fix;
                            t_rew(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_rew;
                            pad_size_trials(j) = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(j).spectro,1);
                        end
                        %% find max trials, pad with Nans and then sort
                        [~,indx_max] = max(pad_size_trials); % find the longest trial
                        pad_size = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(indx_max).spectro,1); % find how big it is
                        freq = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(1).freq;
                        ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(indx_max).ts-1;
                        % Pad with NaNs to compensate for different sized vectors
                        for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl)
                            % find location of t_rew
                            %                             if cond >= 2
                            %                                 this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts;
                            %                                 indx_rew = this_ts(this_ts<=t_rew(ntrl)+0.250); indx_rew = length(indx_rew);
                            %                                 monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro(indx_rew+1:pad_size,:) = NaN;
                            %                             else
                            this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts;
                            indx_end = this_ts(this_ts<=trl_end(ntrl)+0.250); indx_end = length(indx_end);
                            this_spectro = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro;
                            this_spectro(indx_end+1:pad_size,:) = NaN;
                            %                             end
                            padded_spectro(ntrl).spectro = this_spectro;
                        end
                        % sort spectrograms
                        
                        [trl_end_sort, trl_end_sort_indx] = sort(trl_end); % sort the trials
                        spectro_sorted = padded_spectro(trl_end_sort_indx)';
                        % extract different bands, average and normalize by max
                        clear de th be_low be_high max_de max_th max_be_low max_be_high de_norm th_norm be_low_norm be_high_norm
                        for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl)
                            de(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=delta(1) & freq<=delta(2)),2); max_de(ntrl) = max(de(:,ntrl)); de_norm(:,ntrl) = (de(:,ntrl)./max_de(ntrl));
                            th(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=theta(1) & freq<=theta(2)),2); max_th(ntrl) = max(th(:,ntrl)); th_norm(:,ntrl) = (th(:,ntrl)./max_th(ntrl));
                            be_low(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=low_beta(1) & freq<=low_beta(2)),2); max_be_low(ntrl) = max(be_low(:,ntrl)); be_low_norm(:,ntrl) = (be_low(:,ntrl)./max_be_low(ntrl));
                            be_high(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=high_beta(1) & freq<=high_beta(2)),2); max_be_high(ntrl) = max(be_high(:,ntrl)); be_high_norm(:,ntrl) = (be_high(:,ntrl)./max_be_high(ntrl));
                        end
                        max_bands = max([max(max_th) max(max_be_low) max(max_be_high)]);
                        
                        % plot
                        figure('Position',[1838 295 1679 381],'Name',['Session ', num2str(nsess)]); hold on; % colormap(copper);
                        %% plot delta
                        % figure; hold on; colormap(winter);
                        subplot(1,4,1); hold on;
                        %                         pcolor(ts,1:size(de_norm,2),de_norm'); colorbar;
                        %                         shading interp;
                        imagesc(ts,1:size(de_norm,2),de_norm'); colorbar;
                        % if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(th,2),(th./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(th,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['delta M' num2str(nmonk) ' ' (areas{narea}) ' cond ' num2str(cond)])
                        ylabel('trial number'); xlabel('Stop Time (s)'); axis square;vline(0,'-w'); vline(0,'-w');vline([-0.5 0.5], '--w')
                        %                         if cond==2, c = [0 1 0]; else c = [0 0 0]; end
                        %                         plot(ts,nanmean(de_norm'*(size(de_norm',1))), 'LineWidth',2, 'Color', c);
                        %                         set(gca,'xlim',[-0.8 0.8],'ylim',[0.4 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');
                        
                        %% plot theta
                        %figure; hold on; colormap(winter);
                        subplot(1,4,2);hold on;
                        imagesc(ts,1:size(th_norm,2),th_norm'); colorbar;
                        % if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(th,2),(th./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(th,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['theta M' num2str(nmonk) ' ' (areas{narea}) ' cond ' num2str(cond)]); axis square; vline(0,'-w');vline([-0.5 0.5], '--w')
                        
                        %                         plot(ts,nanmean(th_norm'*(size(th_norm',1))), 'LineWidth',1, 'k');
                        %                         set(gca,'xlim',[-0.8 0.8],'ylim',[0.3 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');
                        
                        %ylabel('trial number'); xlabel('Time (s)');
                        %% plot low beta
                        %figure; hold on; colormap(winter);
                        subplot(1,4,3);hold on;
                        imagesc(ts,1:size(be_low_norm,2),be_low_norm'); colorbar;
                        % if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(be_low,2),(be_low./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(be_low,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['Low-B M' num2str(nmonk) ' ' (areas{narea}) ' cond ' num2str(cond)]); axis square;vline(0,'-w');vline([-0.5 0.5], '--w')
                        %ylabel('trial number'); xlabel('Time (s)');axis square
                        
                        %                         plot(ts,nanmean(be_low_norm'*(size(be_low_norm',1))), 'LineWidth',1, 'k');
                        %                         set(gca,'xlim',[-0.8 0.8],'ylim',[0.3 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');
                        
                        %% plot high beta
                        %figure; hold on; colormap(winter);
                        subplot(1,4,4);hold on;
                        imagesc(ts,1:size(be_high_norm,2),be_high_norm'); colorbar;
                        % if cond==2,scatter(sort(t_rew)-0.75,1:size(t_rew,2),3,'k','filled'); end
                        % imagesc(ts,1:size(be_high,2),(be_high./max_bands)'); colorbar;
                        % set(gca,'xlim',[-0.25 trl_end_sort(end)],'ylim',[0 size(be_high,2)], 'FontSize', 22)
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0 size(th,2)], 'FontSize', 22)
                        title(['High-B M' num2str(nmonk) ' ' (areas{narea}) ' cond ' num2str(cond)]); axis square; vline(0,'-w');vline([-0.5 0.5], '--w')
                        %ylabel('trial number'); xlabel('Time (s)'); axis square
                        
                        %                         plot(ts,nanmean(be_low_norm'*(size(be_low_norm',1))), 'LineWidth',1, 'k');
                        %                         set(gca,'xlim',[-0.8 0.8],'ylim',[0.3 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');
                        
                        %% Mean per band
                        
                        figure(6);
                        subplot(1,4,1); hold on; title delta;
                        if cond==2, c = [0 1 0]; else c = [0 0 0]; end
                        % plot(ts,nanmean(de_norm'), 'LineWidth',1, 'Color', c);
                        shadedErrorBar(ts,nanmean(de_norm'),nanstd(de_norm')/sqrt(size(de_norm',1)), 'lineprops',{'Color', c})
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k')
                        
                        subplot(1,4,2); hold on; title theta;
                        if cond==2, c = [0 1 0]; else c = [0 0 0]; end
                        % plot(ts,nanmean(th_norm'), 'LineWidth',1, 'Color', c);
                        shadedErrorBar(ts,nanmean(th_norm'),nanstd(th_norm')/sqrt(size(th_norm',1)), 'lineprops',{'Color', c})
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k')
                        
                        subplot(1,4,3); hold on; title Low-B;
                        if cond==2, c = [0 1 0]; else c = [0 0 0]; end
                        % plot(ts,nanmean(be_low_norm'), 'LineWidth',1, 'Color', c);
                        shadedErrorBar(ts,nanmean(be_low_norm'),nanstd(be_low_norm')/sqrt(size(be_low_norm',1)), 'lineprops',{'Color', c})
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k')
                        
                        subplot(1,4,4); hold on; title High-B;
                        if cond==2, c = [0 1 0]; else c = [0 0 0]; end
                        % plot(ts,nanmean(be_high_norm'), 'LineWidth',1, 'Color', c);
                        shadedErrorBar(ts,nanmean(be_high_norm'),nanstd(be_high_norm')/sqrt(size(be_high_norm',1)), 'lineprops',{'Color', c})
                        set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k')
                        
                        if cond == 1
                            th_norm_all_err = [th_norm' ; th_norm'];
                            be_low_norm_all_err = [be_low_norm' ; be_low_norm'];
                        else
                            th_norm_all_corr = [th_norm' ; th_norm'];
                            be_low_norm_all_corr = [be_low_norm' ; be_low_norm'];
                        end
                        % plot by trial to inspect
                        %                         figure('Name',['Monk ' num2str(nmonk) ' area ' (areas{narea}) ' cond ' num2str(cond) ' sess ' num2str(nsess)]); hold on;
                        %p = numSubplots(50);
                        %                         p = numSubplots(size(th_norm,2));
                        % pick_trls = randperm(size(th_norm,2)); cnt=1;
                        % for ntrl = pick_trls(1:50)
                        %                         for ntrl = 1:size(th_norm,2)
                        %                             freq = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).freq;
                        %                             ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts-1;
                        %                             spectro = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro';
                        %                             subplot(p(1),p(2),ntrl); %colormap(winter);
                        %                             imagesc(ts,freq, spectro); axis xy; %colorbar;
                        %                             set(gca,'xlim',[-0.25 ts(end)],'ylim',[0 50], 'FontSize', 22)
                        %                             %xlabel('time (s)'); ylabel('frequency (Hz)');
                        %                             %title(['Monkey ' num2str(nmonk) ' area ' (areas{narea}) ' trial ' num2str(ntrl) ' cond ' num2str(cond)])
                        %                             hline([4 12 20 30],'-w'); vline(0,'-w'); vline(0.3,'-w'); %vline(move_on(ntrl),'-r');
                        %                             %if cond == 2, vline(t_rew,'-r'), end
                        %                             axis off;
                        %                             % cnt=cnt+1;
                        %                         end
                        
                    end
                end
                % plot all power for all sessions per monkey
                figure; hold on;
                subplot(1,2,1); hold on; title([(areas{narea}) ' theta']);
                shadedErrorBar(ts,nanmean(th_norm_all_err), nanstd(th_norm_all_err)/sqrt(size(th_norm_all_err,1)), 'lineprops', 'k'); 
                shadedErrorBar(ts,nanmean(th_norm_all_corr), nanstd(th_norm_all_corr)/sqrt(size(th_norm_all_corr,1)), 'lineprops', 'g'); 
                set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k'); 
                
                subplot(1,2,2); hold on; title([(areas{narea}) ' beta']);
                shadedErrorBar(ts,nanmean(be_low_norm_all_err), nanstd(be_low_norm_all_err)/sqrt(size(be_low_norm_all_err,1)), 'lineprops', 'k'); 
                shadedErrorBar(ts,nanmean(be_low_norm_all_corr), nanstd(be_low_norm_all_corr)/sqrt(size(be_low_norm_all_corr,1)), 'lineprops', 'g');
                set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k'); 
            end
            close all
            % gather for all monkeys
            th_norm_monk_corr = [th_norm_monk_corr ; th_norm_all_corr]; th_norm_monk_corr_std = [th_norm_monk_corr ; nanstd(th_norm_all_corr)]; 
            th_norm_monk_err = [th_norm_monk_err ; th_norm_all_err];  th_norm_monk_err_std = [th_norm_monk_err ; nanstd(th_norm_all_err)];
            
            be_low_norm_monk_corr = [be_low_norm_monk_corr ; be_low_norm_all_corr];
            be_low_norm_monk_err = [be_low_norm_monk_err ; be_low_norm_all_err]; 
        end
        
        % plot avg for all monkeys
        figure; hold on
        subplot(1,2,1); hold on; title([(areas{narea}) ' theta']);
        shadedErrorBar(ts,nanmean(th_norm_monk_err), nanstd(th_norm_monk_err)/sqrt(size(th_norm_monk_err,1)), 'lineprops', 'k');
        shadedErrorBar(ts,nanmean(th_norm_monk_corr), nanstd(th_norm_monk_corr)/sqrt(size(th_norm_monk_corr,1)), 'lineprops', 'g');
        set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k');

        subplot(1,2,2); hold on; title([(areas{narea}) ' beta']);
        shadedErrorBar(ts,nanmean(be_low_norm_monk_err), nanstd(be_low_norm_monk_err)/sqrt(size(be_low_norm_monk_err,1)), 'lineprops', 'k');
        shadedErrorBar(ts,nanmean(be_low_norm_monk_corr), nanstd(be_low_norm_monk_corr)/sqrt(size(be_low_norm_monk_corr,1)), 'lineprops', 'g');
        set(gca,'xlim',[-0.8 0.8],'ylim',[0.1 0.8],'yTick', [0.1 0.8],'TickDir','out', 'FontSize', 22);  axis square; vline(0,'-k');vline([-0.5 0.5], '--k');
        
        
        
    case 'spectrogram_trl_session_align_target'
        type = 'reward'
        delta = [0.5 4];
        theta = [4 12];
        low_beta = [12 20];
        high_beta = [20 30];
        for nmonk = 2 %1:length(monk)
            for nsess = 1:length(monk(nmonk).sess)
                for cond = 2
                    areas = fieldnames(monk(nmonk).sess(nsess).trialtype.(type)(1).area); % get areas
                    for narea = 1:length(areas)
                        nCh = length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch);
                        if nCh < 17
                            electrode = 'linearprobe16'
                        elseif nCh > 17 & nCh < 31
                            electrode = 'linearprobe24'
                        elseif nCh > 31 & nCh < 36
                            electrode = 'linearprobe32'
                        elseif nCh > 37 & nCh < 90
                            electrode = 'utah2x48'
                        else
                            electrode = 'utah96'
                        end
                        if electrode(1:4) == 'line',[xloc,yloc] = map_linearprobe([],electrode); else [xloc,yloc] = map_utaharray([],electrode);end
                        [channel_id,electrode_id] = MapChannel2Electrode(electrode);
                        [~,indx_ch] = sort(electrode_id); reorderindx = channel_id(indx_ch); if strcmp(electrode,'utah2x48'), reorderindx = reorderindx(1:48);end 
                        figure('Name',['Session ' num2str(nsess) ' area ' areas{narea}]); hold on; colormap(winter);
                        for chNum = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch)
                            % get trl indx and behavior marker
                            clear trl_end move_on t_sacc t_fix t_rew pad_size_trials ts
                            indx_trls = monk(nmonk).behavior(nsess).stats.trialtype.(type)(cond).trlindx; trl_num = find(indx_trls);
                            for j = 1:length(trl_num)
                                trl_end(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_end;
                                move_on(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_move;
                                t_sacc{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_sac;
                                t_fix{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_fix;
                                t_rew(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_rew;
                                
                                pad_size_trials(j) = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(j).spectrogram,1);
                            end
                            
                            %% find max trials, pad with Nans and then sort
                            [~,indx_max] = max(pad_size_trials); % find the longest trial
                            pad_size = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(indx_max).spectrogram,1); % find how big it is
                            freq = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(1).freq_spectrogram;
                            ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(1).ts_spectrogram-1;
                            % Pad with NaNs to compensate for different sized vectors
                            for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl)
                                % find location of t_rew
                                %                             if cond >= 2
                                %                                 this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).ts;
                                %                                 indx_rew = this_ts(this_ts<=t_rew(ntrl)+0.250); indx_rew = length(indx_rew);
                                %                                 monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).pw_trl(ntrl).spectro(indx_rew+1:pad_size,:) = NaN;
                                %                             else
                                this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(ntrl).ts_spectrogram; % CHECK!!
                                indx_end = this_ts(this_ts<=trl_end(ntrl)); indx_end = length(indx_end);
                                this_spectro = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(ntrl).spectrogram;
                                this_spectro(indx_end+1:pad_size,:) = NaN;
                                %                             end
                                padded_spectro(ntrl).spectro = this_spectro;
                            end
                            
                            [trl_end_sort, trl_end_sort_indx] = sort(trl_end); % sort the trials
                            spectro_sorted = padded_spectro(trl_end_sort_indx)';
                            % extract different bands, average and normalize by max
                            clear de th be_low be_high max_de max_th max_be_low max_be_high de_norm th_norm be_low_norm be_high_norm
                            for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl)
                                de(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=delta(1) & freq<=delta(2)),2); max_de(ntrl) = max(de(:,ntrl)); de_norm(:,ntrl) = (de(:,ntrl)./max_de(ntrl));
                                th(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=theta(1) & freq<=theta(2)),2); max_th(ntrl) = max(th(:,ntrl)); th_norm(:,ntrl) = (th(:,ntrl)./max_th(ntrl));
                                be_low(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=low_beta(1) & freq<=low_beta(2)),2); max_be_low(ntrl) = max(be_low(:,ntrl)); be_low_norm(:,ntrl) = (be_low(:,ntrl)./max_be_low(ntrl));
                                be_high(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=high_beta(1) & freq<=high_beta(2)),2); max_be_high(ntrl) = max(be_high(:,ntrl)); be_high_norm(:,ntrl) = (be_high(:,ntrl)./max_be_high(ntrl));
                            end
                            max_bands = max([max(max_th) max(max_be_low) max(max_be_high)]);
                            
                            %% plot
                            if electrode(1:4) == 'line'
                                p = numSubplots(nCh);
                                %                                 subplot(p(1),p(2),chNum)
                                subplot(nCh,1,chNum); hold on;
%                                 imagesc(ts,1:size(de_norm,2),de_norm'); %colorbar;
%                                 imagesc(ts,1:size(th_norm,2),th_norm'); %colorbar;
%                                 imagesc(ts,1:size(be_low_norm,2),be_low_norm'); %colorbar;
                                                                imagesc(ts,1:size(be_high_norm,2),be_high_norm'); %colorbar;
                                set(gca,'xlim',[-0.25 4],'ylim',[0 size(th,2)], 'FontSize', 22)
                                axis off; box off; axis square; vline(0,'w')
                            else
                                subplot(10,10,10*(xloc(chNum)-1) + yloc(chNum)); hold on;
%                                 imagesc(ts,1:size(de_norm,2),de_norm'); %colorbar;
%                                 imagesc(ts,1:size(th_norm,2),th_norm'); %colorbar;
%                                 imagesc(ts,1:size(be_low_norm,2),be_low_norm'); %colorbar;
                                                                imagesc(ts,1:size(be_high_norm,2),be_high_norm'); %colorbar;
                                set(gca,'xlim',[-0.25 4],'ylim',[0 size(th,2)], 'FontSize', 22)
                                axis off; box off; axis square; vline(0,'w')
                            end
                        end
                    end
                end
            end
        end
        
    case 'spectrogram_trl_session_align_stop'
        type = 'reward' % in trialtype, pick rewarded and unrewarded trials (other options are all, density, etc.)
        delta = [0.5 4];
        theta = [4 12];
        low_beta = [12 20];
        high_beta = [20 30];
        for nmonk = 1 %1:length(monk)
            for nsess = 1:length(monk(nmonk).sess)
                for cond = 2
                    areas = fieldnames(monk(nmonk).sess(nsess).trialtype.(type)(1).area); % get areas
                    for narea = 1:length(areas)
                        nCh = length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch);
                        if nCh < 17
                            electrode = 'linearprobe16'
                        elseif nCh > 17 & nCh < 31
                            electrode = 'linearprobe24'
                        elseif nCh > 31 & nCh < 36
                            electrode = 'linearprobe32'
                        elseif nCh > 37 & nCh < 90
                            electrode = 'utah2x48'
                        else
                            electrode = 'utah96'
                        end
                        if electrode(1:4) == 'line',[xloc,yloc] = map_linearprobe([],electrode); else [xloc,yloc] = map_utaharray([],electrode);end
                        [channel_id,electrode_id] = MapChannel2Electrode(electrode);
                        [~,indx_ch] = sort(electrode_id); reorderindx = channel_id(indx_ch); if strcmp(electrode,'utah2x48'), reorderindx = reorderindx(1:48);end 
                        figure('Name',['Session ' num2str(nsess) ' area ' areas{narea}]); hold on; colormap(winter);
                        for chNum = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch)
                            % get trl indx and behavior marker
                            clear trl_end move_on t_sacc t_fix t_rew pad_size_trials ts
                            indx_trls = monk(nmonk).behavior(nsess).stats.trialtype.(type)(cond).trlindx; trl_num = find(indx_trls);
                            for j = 1:length(trl_num)
                                trl_end(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_end;
                                move_on(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_move;
                                t_stop(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_stop;
                                t_sacc{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_sac;
                                t_fix{j} = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_fix;
                                t_rew(j) = monk(nmonk).behavior(nsess).trials(trl_num(j)).events.t_rew;
                                pad_size_trials(j) = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(j).spectrogram_stop,1);
                            end
                            
                            %% find max trials, pad with Nans and then sort
                            [~,indx_max] = max(pad_size_trials); % find the longest trial
                            pad_size = size(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(indx_max).spectrogram_stop,1); % find how big it is
                            freq = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(1).freq_spectrogram_stop;
                            ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(1).ts_spectrogram_stop-1;
                            % Pad with NaNs to compensate for different sized vectors
                            for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl)
                                this_ts = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(ntrl).ts_spectrogram_stop-1; 
                                indx_end = this_ts(this_ts<=trl_end(ntrl)); indx_end = length(indx_end);
                                this_spectro = monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl(ntrl).spectrogram_stop;
                                this_spectro(indx_end+1:pad_size,:) = NaN;
                                padded_spectro(ntrl).spectro = this_spectro;
                            end
                            
                            [trl_end_sort, trl_end_sort_indx] = sort(trl_end); % sort the trials
                            spectro_sorted = padded_spectro(trl_end_sort_indx)';
                            % extract different bands, average and normalize by max
                            clear de th be_low be_high max_de max_th max_be_low max_be_high de_norm th_norm be_low_norm be_high_norm
                            for ntrl = 1:length(monk(nmonk).sess(nsess).trialtype.(type)(cond).area.(areas{narea}).ch(reorderindx(chNum)).trl)
                                de(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=delta(1) & freq<=delta(2)),2); max_de(ntrl) = max(de(:,ntrl)); de_norm(:,ntrl) = (de(:,ntrl)./max_de(ntrl));
                                th(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=theta(1) & freq<=theta(2)),2); max_th(ntrl) = max(th(:,ntrl)); th_norm(:,ntrl) = (th(:,ntrl)./max_th(ntrl));
                                be_low(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=low_beta(1) & freq<=low_beta(2)),2); max_be_low(ntrl) = max(be_low(:,ntrl)); be_low_norm(:,ntrl) = (be_low(:,ntrl)./max_be_low(ntrl));
                                be_high(:,ntrl) = nanmean(spectro_sorted(ntrl).spectro(:,freq>=high_beta(1) & freq<=high_beta(2)),2); max_be_high(ntrl) = max(be_high(:,ntrl)); be_high_norm(:,ntrl) = (be_high(:,ntrl)./max_be_high(ntrl));
                            end
                            max_bands = max([max(max_th) max(max_be_low) max(max_be_high)]);
                            
                            %% plot
                            if electrode(1:4) == 'line'
                                p = numSubplots(nCh);
                                %                                 subplot(p(1),p(2),chNum)
                                subplot(nCh,1,chNum); hold on;
%                                 imagesc(ts,1:size(de_norm,2),de_norm'); %colorbar;
%                                 imagesc(ts,1:size(th_norm,2),th_norm'); %colorbar;
%                                 imagesc(ts,1:size(be_low_norm,2),be_low_norm'); %colorbar;
                                                                imagesc(ts,1:size(be_high_norm,2),be_high_norm'); %colorbar;
                                set(gca,'xlim',[-0.8 0.8],'ylim',[0 size(th,2)], 'FontSize', 22)
                                axis off; box off; axis square; vline(0,'w')
                            else
                                subplot(10,10,10*(xloc(chNum)-1) + yloc(chNum)); hold on;
%                                 imagesc(ts,1:size(de_norm,2),de_norm'); %colorbar;
%                                 imagesc(ts,1:size(th_norm,2),th_norm'); %colorbar;
%                                 imagesc(ts,1:size(be_low_norm,2),be_low_norm'); %colorbar;
                                                                imagesc(ts,1:size(be_high_norm,2),be_high_norm'); %colorbar;
                                set(gca,'xlim',[-0.8 0.8],'ylim',[0 size(th,2)], 'FontSize', 22)
                                axis off; box off; axis square; vline(0,'w')
                            end
                        end
                    end
                end
            end
        end
        
        
    case 'coherogram_move'
        ev = 'move'
        for nmonk = 1:length(monk)
            coher_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_mu';
            freq_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_freq;
            ts_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_ts;
            
            % plot PPC MST
            figure('Name', ['PPC-->MST' nmonk]);colormap(copper);
            imagesc(ts_ppc_mst-1, freq_ppc_mst,coher_ppc_mst, [0 0.6]);
            axis xy; set(gca,'xlim',[-0.25 0.25], 'ylim', [4 50],'FontSize', 22); colorbar;
            
        end
        
        
        
        
    case 'coherogram_target'
        ev = 'target'
        for nmonk = 1:length(monk)
            coher_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_mu;
            freq_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_freq;
            ts_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_ts;
            
            % plot PPC MST
            figure('Name', ['PPC-->MST' nmonk]);colormap(copper);
            imagesc(ts_ppc_mst-1.3, freq_ppc_mst,coher_ppc_mst', [0 0.5]); vline(-0.3,'b'); vline(0,'w');
            axis xy; set(gca,'xlim',[-0.4 0.25], 'ylim', [4 50],'FontSize', 22); colorbar;
            
        end
        
    case 'coherogram_stop'
        ev = 'reward'
        for nmonk = [1 3] %1:length(monk)
            coher_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_mu;
            freq_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_freq;
            ts_ppc_mst = monk(nmonk).coher.trialtype.all.events.(ev).PPCMST.coher_ts;
            
            % plot MST PPC
            figure('Name', ['MST--PPC' nmonk]);colormap(copper);
            imagesc(ts_ppc_mst-1, freq_ppc_mst,coher_ppc_mst', [0 0.5]);
            axis xy; set(gca,'xlim',[-0.25 0.25], 'ylim', [4 50],'FontSize', 22); colorbar;
        end
        
    case 'coherogram_reward_PPC_MST'
        type = 'reward'
        ev = 'stop'
        time_win = [-0.5 0.5];   % [-0.25 0.25]; [-0.45 0.25];
        freq_win = [3 51];
        align_t = 1;  % 1.3 for target (so it's aligned to target offset
        p = numSubplots(12);  % 9 for reward cond, 12 for rest
        cnt=1;
        for nmonk = 1 %[1 3]
            for nsess = 1:length(monk(nmonk).coher.sess)
                 ncond = length(monk(nmonk).coher.sess(nsess).trialtype.(type));
                for cond = 1:ncond
                    freq = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).PPCMST.coher_freq;
                    ts = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).PPCMST.coher_ts-align_t;
                    p_cohero = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).PPCMST.coher_mu';
                    ts_win = ts(ts>time_win(1) & ts<time_win(2));
                    p_cohero_win = p_cohero(freq>freq_win(1) & freq<freq_win(2),ts>time_win(1) & ts<time_win(2));
                    p_cohero_mu_theta = p_cohero(freq>4 & freq<12,ts>time_win(1) & ts<time_win(2));
                    p_cohero_mu_beta = p_cohero(freq>12 & freq<21,ts>time_win(1) & ts<time_win(2));
                    p_cohero_mu_wideband = p_cohero(freq>5 & freq<40,ts>time_win(1) & ts<time_win(2));
                    max_cohero = max(max(p_cohero_win));
                    % plot
                    %figure('Name',['Monk ' num2str(nmonk) ' sess ' num2str(nsess)  ' PPC--MST ' ev ' rew ' num2str(cond)]);
                    %                     subplot(p(1),p(2),cnt)
                    %
                    %                     colormap(copper);
                    %                     imagesc(ts,freq,p_cohero,[0 1]); axis xy; %colorbar;
                    %                     set(gca,'xlim',[-0.25 0.25], 'ylim',[4 50], 'FontSize', 22); vline(0,'w'); %vline(-0.3,'b');
                    %                     axis square
                    
                    % xlabel('time (s)'); ylabel('frequency (Hz)');
                    
                    % plot win
                    %                     figure('Name',['Normalized  Monk ' num2str(nmonk) ' sess ' num2str(nsess)  ' PPC--MST ' ev ' rew ' num2str(cond)]);
                    figure(1); subplot(p(1),p(2),cnt)
                    
                    %colormap(copper);
                    % imagesc(ts_win,freq(freq>freq_win(1) & freq<freq_win(2)),p_cohero_win/max_cohero, [0 1]); axis xy; %colorbar;
                    imagesc(ts_win,freq(freq>freq_win(1) & freq<freq_win(2)),p_cohero_win, [0 0.6]); axis xy; %colorbar;
                    set(gca,'xlim',[-0.5 0.5], 'ylim',[4 50], 'FontSize', 22); vline(0,'w'); axis square; % vline(-0.3,'b');
                    title (['sess ' num2str(nsess) ' cond ' num2str(cond)])
                    % xlabel('time (s)'); ylabel('frequency (Hz)');
                    
                    % theta mu
                    figure(2); set(gcf,'Name','Theta');
                    subplot(p(1),p(2),cnt)
                    plot(ts_win,mean(p_cohero_mu_theta))
                    set(gca,'xlim',[-0.5 0.5],'ylim',[0.28 0.38], 'FontSize', 22); vline(0,'w'); axis square;
                    title (['sess ' num2str(nsess) ' cond ' num2str(cond)])
               
                    % beta mu
                    figure(3); set(gcf,'Name','Beta');
                    subplot(p(1),p(2),cnt)
                    plot(ts_win,mean(p_cohero_mu_beta))
                    set(gca,'xlim',[-0.5 0.5],'ylim',[0.08 0.2], 'FontSize', 22); vline(0,'w'); axis square;
                    title (['sess ' num2str(nsess) ' cond ' num2str(cond)])
              
                    % wideband mu
                    figure(4); set(gcf,'Name','Wideband');
                    subplot(p(1),p(2),cnt)
                    plot(ts_win,mean(p_cohero_mu_wideband))
                    set(gca,'xlim',[-0.5 0.5],'ylim',[0.12 0.18], 'FontSize', 22); vline(0,'w'); axis square;
                    title (['sess ' num2str(nsess) ' cond ' num2str(cond)])
                    cnt = cnt+1;
                    
                end 
            end
        end
        
    case 'coherogram_target_stop_PPC_MST'
        %Contains average per session
        type = 'reward'
        ev = 'stop'
        brain_area = 'MSTPFC'  % 'PFCMST'  % 'PPCMST'
        time_win = [-0.5 0.5];   % [-0.25 0.25]; [-0.45 0.25];
        freq_win = [3 51];
        align_t = 1;  % 1.3 for target (so it's aligned to target offset
        p = numSubplots(12);  % 9 for reward cond, 12 for rest
        cnt=1;
        for nmonk = 2 %[1 3]
            for nsess = 1:length(monk(nmonk).coher.sess)
                 ncond = length(monk(nmonk).coher.sess(nsess).trialtype.(type));
                for cond = 2 % 1:ncond
                    freq = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).(brain_area).coher_freq;
                    ts = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).(brain_area).coher_ts-align_t;
                    p_cohero = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).(brain_area).coher';
                    ts_win = ts(ts>time_win(1) & ts<time_win(2));
                    p_cohero_win = p_cohero(freq>freq_win(1) & freq<freq_win(2),ts>time_win(1) & ts<time_win(2));
                    coher_all(nsess,:,:) = p_cohero_win; 
                    max_cohero = max(max(p_cohero_win));
                    % plot
                    %figure('Name',['Monk ' num2str(nmonk) ' sess ' num2str(nsess)  ' PPC--MST ' ev ' rew ' num2str(cond)]);
                    %                     subplot(p(1),p(2),cnt)
                    %
                    %                     colovrmap(copper);
                    %                     imagesc(ts,freq,p_cohero,[0 1]); axis xy; %colorbar;
                    %                     set(gca,'xlim',[-0.25 0.25], 'ylim',[4 50], 'FontSize', 22); vline(0,'w'); %vline(-0.3,'b');
                    %                     axis square
                    
                    % xlabel('time (s)'); ylabel('frequency (Hz)');
                    
                    % plot win
                    %                     figure('Name',['Normalized  Monk ' num2str(nmonk) ' sess ' num2str(nsess)  ' PPC--MST ' ev ' rew ' num2str(cond)]);
                    figure(1); subplot(p(1),p(2),cnt)
                    
                    %colormap(parula); copper
                    % imagesc(ts_win,freq(freq>freq_win(1) & freq<freq_win(2)),p_cohero_win/max_cohero, [0 1]); axis xy; %colorbar;
                    imagesc(ts_win,freq(freq>freq_win(1) & freq<freq_win(2)),p_cohero_win, [0 max_cohero]); axis xy; colorbar;
                    set(gca,'xlim',[-0.5 0.5], 'ylim',[4 50], 'FontSize', 22); vline(0,'w'); axis square; % vline(-0.3,'b');
                    title (['sess ' num2str(nsess) ' cond ' num2str(cond)]); if align_t == 1.3; vline(-0.3,'--w'); end
                    % xlabel('time (s)'); ylabel('frequency (Hz)');
                    cnt = cnt+1;
                end 
            end
            % Mean for all sessions for one condition 
            coher_mu = mean(coher_all);
            % plot
            figure; 
            imagesc(ts_win,freq(freq>freq_win(1) & freq<freq_win(2)),squeeze(coher_mu),[0 0.6]); axis xy; colorbar;
            set(gca,'xlim',[-0.5 0.5], 'ylim',[4 50], 'FontSize', 22); vline(0,'w'); axis square; if align_t == 1.3; vline(-0.3,'--w'); end
            title (['Cond ' num2str(cond)]); box off
            xlabel('time (s)'); ylabel('frequency (Hz)');
            
            % mu all
            figure;
            plot(ts_win,smooth(mean(squeeze(coher_mu))))
            set(gca, 'ylim',[0.14 0.31], 'xlim',[-0.5 0.5], 'FontSize', 22, 'TickDir', 'out'); vline(0,'k'); axis square;
            title (['cond ' num2str(cond)]); box off; %if align_t == 1.3; vline(-0.3,'--k'); end
            
            % theta mu
            figure; coher_theta = coher_mu(1,freq>3 & freq<12,:);
            plot(ts_win, smooth(squeeze(mean(coher_theta)))); box off
            set(gca, 'ylim',[0.25 0.36], 'xlim',[-0.5 0.5], 'FontSize', 22,'TickDir', 'out'); vline(0,'k'); axis square;
            title(['theta cond ' num2str(cond)])
            
            % beta mu
            figure; coher_beta = coher_mu(1,freq>12 & freq<21,:);
            plot(ts_win, smooth(squeeze(mean(coher_beta)))); box off
            set(gca, 'ylim',[0.18 0.32], 'xlim',[-0.5 0.5], 'FontSize', 22,'TickDir', 'out'); vline(0,'k'); axis square;
            title(['beta cond ' num2str(cond)])
            
            % wideband mu
            figure; coher_wideband = coher_mu(1,freq>5 & freq<40,:);
            plot(ts_win, squeeze(mean(coher_wideband))); box off
            set(gca, 'ylim',[0.13 0.19], 'xlim',[-0.5 0.5], 'FontSize', 22,'TickDir', 'out'); vline(0,'k'); axis square;
            title('wideband') % compute statistics
            
            % save per monk
            save(['monk ' num2str(nmonk) ' cond' num2str(cond)], 'coher_mu', 'coher_theta','coher_beta', 'coher_wideband');
            
        end
        
        case 'coherogram_target_stop_per_band'
        %Contains average per session
        type = 'reward'
        ev = 'stop'
        brain_area = 'MSTPFC'  % 'PFCMST'  % 'PPCMST' $ PFCPPC
        band = 'beta';
        time_win = [-0.5 0.5];   % [-0.25 0.25]; [-0.45 0.25];
        freq_win = [3 51];
        align_t = 1;  % 1.3 for target (so it's aligned to target offset
        p = numSubplots(12);  % 9 for reward cond, 12 for rest
        cnt=1;
        for nmonk = 2 %[1 3]
            for nsess = 1:length(monk(nmonk).coher.sess)
                 ncond = length(monk(nmonk).coher.sess(nsess).trialtype.(type));
                for cond = 1 % 1:ncond
                    % freq = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).(brain_area).coher_freq;
                    ts = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).(brain_area).(band).coher_ts-align_t;
                    p_cohero = monk(nmonk).coher.sess(nsess).trialtype.(type)(cond).events.(ev).(brain_area).(band).coher';
                    ts_win = ts(ts>time_win(1) & ts<time_win(2));
                    p_cohero_win = p_cohero(:,ts>time_win(1) & ts<time_win(2));
                    coher_all(nsess,:,:) = p_cohero_win; 
                    max_cohero = max(max(p_cohero_win));
                   
                    figure(1); subplot(p(1),p(2),cnt)
                    plot(ts_win,mean(p_cohero_win)); box off; 
                    set(gca,'xlim',[-0.5 0.5], 'FontSize', 22); vline(0,'w'); axis square; % vline(-0.3,'b');
                    title (['sess ' num2str(nsess) ' cond ' num2str(cond)]); if align_t == 1.3; vline(-0.3,'--w'); end
                    % xlabel('time (s)'); ylabel('frequency (Hz)');
                    cnt = cnt+1;
                end 
            end
            % Mean for all sessions for one condition 
            coher_mu = mean(coher_all);
            % plot
            figure; 
            plot(ts_win,smooth(squeeze(mean(coher_mu)))); ylim([0.24 0.33])
            set(gca,'xlim',[-0.5 0.5], 'FontSize', 22); vline(0,'k'); axis square; if align_t == 1.3; vline(-0.3,'--w'); end
            title ([ band ' Cond ' num2str(cond)]); box off
            xlabel('time (s)'); ylabel('Coherence');
            
        end
        
        
    case 'coherence_dist'
        
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                if ~isempty(monk(m).coher)
                    %% MST
                    dist_mst = monk(m).pop(sess).trialtype.all.MST.dist; coher_mst = monk(m).pop(sess).trialtype.all.MST.coherByDist;
                    phase_mst = monk(m).pop(sess).trialtype.all.MST.phaseByDist; freq =  monk(m).pop(sess).trialtype.all.crosslfp.freq;
                    %plot coher
                    nch = size(coher_mst,2); cmap = jet(nch);
                    figure; hold on; for k=1:nch, plot(freq,coher_mst(:,k),'Color',cmap(k,:)); end
                    set(gca,'xlim',[2 50],'ylim',[0.6 1], 'TickDir', 'out', 'FontSize', 20); box off
                    title('Coherence MST'); xlabel('frequency');  ylabel('coherence');
                    %plot phase
                    figure; hold on; for k=1:nch, plot(freq,phase_mst(:,k),'Color',cmap(k,:)); end
                    set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off
                    title('Phase MST'); xlabel('frequency');  ylabel('rad');
                    % plot across distance (theta)
                    figure('Name',['Monk ' num2str(m) 'Sess ' num2str(sess)]);hold on;
                    shadedErrorBar(dist_mst,nanmean(coher_mst(4:8,:)),nanstd(coher_mst(4:8,:)), 'lineprops','b');  % plot across distance (theta)
                    shadedErrorBar(dist_mst,nanmean(coher_mst(8:12,:)),nanstd(coher_mst(8:12,:)),'lineprops','c');  % plot across distance (alfa)
                    shadedErrorBar(dist_mst,nanmean(coher_mst(12:30,:)),nanstd(coher_mst(12:30,:)),'lineprops','r'); % plot across distance (beta)
                    set(gca,'xlim',[2 ceil(dist_mst(end))],'ylim',[0.6 1],'xTick',[2 ceil(dist_mst(end))], 'TickDir', 'out', 'FontSize', 20); box off
                    title('Coherence across distance MST'); xlabel('electrode distance'); ylabel('coherence');
                    legend({'theta', 'alpha', 'beta'}, 'box', 'off');
                    
                    % mean coherence in different bands
                    figure; hold on;
                    plot(1,nanmean(coher_mst(4:7,2:end)),'.k'); plot(1,nanmean(nanmean(coher_mst(4:7,2:end))), '.r', 'MarkerSize',20)
                    plot(2,nanmean(coher_mst(8:11,2:end)),'.k'); plot(2,nanmean(nanmean(coher_mst(8:11,2:end))), '.r', 'MarkerSize',20)
                    plot(3,nanmean(coher_mst(12:30,2:end)),'.k'); plot(3,nanmean(nanmean(coher_mst(12:30,2:end))), '.r', 'MarkerSize',20)
                    set(gca,'xlim',[0.5 3.5],'xTick',[1 2 3],'xTickLabel',[{'theta'} {'alfa'} {'beta'}], 'TickDir', 'out', 'FontSize', 20); box off
                    ylim([0.6 1])
                end
                %% PPC
                dist_ppc = monk(m).pop(sess).trialtype.all.PPC.dist; coher_ppc = monk(m).pop(sess).trialtype.all.PPC.coherByDist;
                phase_ppc = monk(m).pop(sess).trialtype.all.PPC.phaseByDist;
                %plot coher
                nch = size(coher_ppc,2); cmap = jet(nch);
                figure; hold on; for k=1:nch, plot(freq,coher_ppc(:,k),'Color',cmap(k,:)); end
                set(gca,'xlim',[2 50],'ylim',[0.6 1], 'TickDir', 'out', 'FontSize', 20); box off
                title('Coherence PPC'); xlabel('frequency');  ylabel('coherence');
                %plot phase
                figure; hold on; for k=1:nch, plot(freq,phase_ppc(:,k),'Color',cmap(k,:)); end
                set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off
                title('Phase PPC'); xlabel('frequency');  ylabel('rad');
                % plot across distance
                figure('Name',['Monk ' num2str(m) 'Sess ' num2str(sess)]);hold on;
                shadedErrorBar(dist_ppc,nanmean(coher_ppc(4:8,:)),nanstd(coher_ppc(4:8,:)), 'lineprops','b');  % plot across distance (theta)
                shadedErrorBar(dist_ppc,nanmean(coher_ppc(8:12,:)),nanstd(coher_ppc(8:12,:)),'lineprops','c');  % plot across distance (alfa)
                shadedErrorBar(dist_ppc,nanmean(coher_ppc(12:30,:)),nanstd(coher_ppc(12:30,:)),'lineprops','r'); % plot across distance (beta)
                set(gca,'xlim',[2 ceil(dist_ppc(end))],'ylim',[0.7 0.9],'xTick',[2 ceil(dist_ppc(end))], 'TickDir', 'out', 'FontSize', 20); box off
                title('Coherence across distance PPC'); xlabel('electrode distance'); ylabel('coherence');
                legend({'theta', 'alpha', 'beta'}, 'box', 'off');
                
                % mean coherence in different bands
                figure; hold on;
                plot(1,nanmean(coher_ppc(4:7,2:end)),'.k'); plot(1,nanmean(nanmean(coher_ppc(4:7,2:end))), '.r', 'MarkerSize',20)
                plot(2,nanmean(coher_ppc(8:11,2:end)),'.k'); plot(2,nanmean(nanmean(coher_ppc(8:11,2:end))), '.r', 'MarkerSize',20)
                plot(3,nanmean(coher_ppc(12:30,2:end)),'.k'); plot(3,nanmean(nanmean(coher_ppc(12:30,2:end))), '.r', 'MarkerSize',20)
                set(gca,'xlim',[0.5 3.5],'xTick',[1 2 3],'xTickLabel',[{'theta'} {'alfa'} {'beta'}], 'TickDir', 'out', 'FontSize', 20); box off
                ylim([0.6 1])
            end
        end
        
    case 'coherence_dist_all'
        for m = 1:length(monk)
            clear coh_mst ph_mst coh_ppc ph_ppc
            % plot per monkey
            for sess = 1:length(monk(m).pop)
                coh_mst(sess,:,:) = monk(m).pop(sess).trialtype.all.MST.coherByDist;
                ph_mst(sess,:,:) = monk(m).pop(sess).trialtype.all.MST.phaseByDist;
                coh_ppc(sess,:,:) = monk(m).pop(sess).trialtype.all.PPC.coherByDist;
                ph_ppc(sess,:,:) = monk(m).pop(sess).trialtype.all.PPC.phaseByDist;
            end
            coher_mst = squeeze(nanmean(coh_mst,1)); phase_mst = squeeze(nanmean(ph_mst,1));
            coher_ppc = squeeze(nanmean(coh_ppc,1)); phase_ppc = squeeze(nanmean(ph_ppc,1));
            dist_mst = monk(m).pop(sess).trialtype.all.MST.dist;
            dist_ppc = monk(m).pop(sess).trialtype.all.PPC.dist;
            freq = monk(m).pop(sess).trialtype.all.crosslfp.freq;
            
            % plot per monkey
            %% MST
            %plot coher
            nch = size(coher_mst,2);
            cmap = jet(nch);
            figure; hold on; for k=1:nch, plot(freq,coher_mst(:,k),'Color',cmap(k,:)); end
            set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off; ylim([0.6 1]);
            title('Coherence MST'); xlabel('frequency');  ylabel('coherence');
            vline([4 8],'k'); vline([8 12],'k'); vline([12 30],'k');
            %plot phase
            figure; hold on; for k=1:nch, plot(freq,phase_mst(:,k),'Color',cmap(k,:)); end
            set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off
            title('Phase MST'); xlabel('frequency');  ylabel('rad');
            % plot across distance
            figure('Name',['Monk ' num2str(m)]);hold on;
            shadedErrorBar(dist_mst,nanmean(coher_mst(4:8,:)),nanstd(coher_mst(4:8,:)),'lineprops','b');  % plot across distance (theta)
            shadedErrorBar(dist_mst,nanmean(coher_mst(8:12,:)),nanstd(coher_mst(8:12,:)),'lineprops','c');  % plot across distance (alfa)
            shadedErrorBar(dist_mst,nanmean(coher_mst(12:30,:)),nanstd(coher_mst(12:30,:)),'lineprops','r'); % plot across distance (beta)
            set(gca,'xlim',[2 24],'ylim',[0.7 1],'xTick',[2 23], 'TickDir', 'out', 'FontSize', 20); box off
            xlabel('electrode distance'); ylabel('coherence');
            legend({'theta', 'alpha', 'beta'}, 'box', 'off'); title(['Coh across distance MST Monk ' num2str(m)]);
            
            % mean coherence in different bands for all sessions in one monkey
            c_mst_theta = []; c_mst_alfa = []; c_mst_beta = [];
            for s = 1:length(coh_mst(:,1,1))
                c_mst_theta = [c_mst_theta ; nanmean(squeeze(coh_mst(s,4:8,2:end)))'];
                c_mst_alfa = [c_mst_alfa ; nanmean(squeeze(coh_mst(s,8:12,2:end)))'];
                c_mst_beta = [c_mst_beta ; nanmean(squeeze(coh_mst(s,12:30,2:end)))'];
            end
            
            figure; hold on;
            plot(1,c_mst_theta,'.k'); plot(1,nanmean(c_mst_theta), '.r', 'MarkerSize',20);
            plot(2,c_mst_alfa,'.k'); plot(2,nanmean(c_mst_alfa), '.r', 'MarkerSize',20);
            plot(3,c_mst_beta,'.k'); plot(3,nanmean(c_mst_beta), '.r', 'MarkerSize',20);
            set(gca,'xlim',[0.5 3.5],'xTick',[1 2 3],'xTickLabel',[{'theta'} {'alfa'} {'beta'}], 'TickDir', 'out', 'FontSize', 20); box off
            ylim([0.6 1]); ylabel('Coherence'); title(['MST Monk ' num2str(m)]);
            
            %% PPC
            nch = size(coher_ppc,2);
            %plot coher
            cmap = jet(nch);
            figure; hold on; for k=1:nch, plot(freq,coher_ppc(:,k),'Color',cmap(k,:)); end
            set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off; ylim([0.6 1]);
            title(['Coherence PPC ' num2str(m)]); xlabel('frequency');  ylabel('coherence');
            vline([4 8],'k'); vline([8 12],'k'); vline([12 30],'k');
            %plot phase
            figure; hold on; for k=1:nch, plot(freq,phase_ppc(:,k),'Color',cmap(k,:)); end
            set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off
            title(['Phase PPC ' num2str(m)]); xlabel('frequency');  ylabel('rad');
            % plot across distance
            figure('Name',['Monk ' num2str(m)]);hold on;
            shadedErrorBar(dist_ppc,nanmean(coher_ppc(4:8,:)),nanstd(coher_ppc(4:8,:)),'lineprops','b');  % plot across distance (theta)
            shadedErrorBar(dist_ppc,nanmean(coher_ppc(8:12,:)),nanstd(coher_ppc(8:12,:)),'lineprops','c');  % plot across distance (alfa)
            shadedErrorBar(dist_ppc,nanmean(coher_ppc(12:30,:)),nanstd(coher_ppc(12:30,:)),'lineprops','r'); % plot across distance (beta)
            set(gca,'xlim',[2 12],'ylim',[0.7 1],'xTick',[2 12], 'TickDir', 'out', 'FontSize', 20); box off
            xlabel('electrode distance'); ylabel('coherence');
            legend({'theta', 'alpha', 'beta'}, 'box', 'off'); title(['Coh across distance PPC Monk ' num2str(m)]);
            
            % mean coherence in different bands for all sessions in one monkey
            c_ppc_theta = []; c_ppc_alfa = []; c_ppc_beta = [];
            for s = 1:length(coh_mst(:,1,1))
                c_ppc_theta = [c_ppc_theta ; nanmean(squeeze(coh_mst(1,4:7,2:end)))'];
                c_ppc_alfa = [c_ppc_alfa ; nanmean(squeeze(coh_mst(1,8:11,2:end)))'];
                c_ppc_beta = [c_ppc_beta ; nanmean(squeeze(coh_mst(1,12:30,2:end)))'];
            end
            
            figure; hold on;
            plot(1,c_ppc_theta,'.k'); plot(1,nanmean(c_ppc_theta), '.r', 'MarkerSize',20);
            plot(2,c_ppc_alfa,'.k'); plot(2,nanmean(c_ppc_alfa), '.r', 'MarkerSize',20);
            plot(3,c_ppc_beta,'.k'); plot(3,nanmean(c_ppc_beta), '.r', 'MarkerSize',20);
            set(gca,'xlim',[0.5 3.5],'xTick',[1 2 3],'xTickLabel',[{'theta'} {'alfa'} {'beta'}], 'TickDir', 'out', 'FontSize', 20); box off
            ylim([0.6 1]); ylabel('Coherence'); title(['PPC Monk ' num2str(m)]);
        end
        
    case 'coherence_dist_reward_density'
        type = 'reward'
        area = 'PPC'
        for m = 1:length(monk)
            clear coh_err coh_corr ph_err ph_corr
            % plot per monkey
            for sess = 1:length(monk(m).pop)
                coh_err(sess,:,:) = monk(m).pop(sess).trialtype.(type)(1).PPC.coherByDist;
                coh_corr(sess,:,:) = monk(m).pop(sess).trialtype.(type)(2).PPC.coherByDist;
                ph_err(sess,:,:) = monk(m).pop(sess).trialtype.(type)(1).PPC.phaseByDist;
                ph_corr(sess,:,:) = monk(m).pop(sess).trialtype.(type)(2).PPC.phaseByDist;
            end
            coher_err = squeeze(nanmean(coh_err,1)); phase_err = squeeze(nanmean(ph_err,1));
            coher_err_mu = nanmean(coher_err,2);
            coher_corr = squeeze(nanmean(coh_corr,1)); phase_corr = squeeze(nanmean(ph_corr,1));
            coher_corr_mu = nanmean(coher_corr,2);
            
            dist_ppc = monk(m).pop(sess).trialtype.all.PPC.dist;
            freq = monk(m).pop(sess).trialtype.all.crosslfp.freq;
            
            % plot per monkey
            %% PPC
            nch = size(coher_corr,2);
            %plot coher
            figure; hold on; plot(freq,coher_err,'k');
            plot(freq,coher_corr,'g'); legend({'unrewarded', 'rewarded'},'box', 'off');
            set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off; ylim([0.6 1]);
            title(['Coherence monk ' num2str(m) ' area ' area]); xlabel('frequency');  ylabel('coherence');
            vline([4 8],'k'); vline([8 12],'k'); vline([12 30],'k');
            
            % ratio
            figure; hold on;
            plot(freq, coher_err - coher_corr, 'Color', [0.5 0.5 0.5]);
            plot(freq, nanmean(coher_err,2)- nanmean(coher_corr,2), 'k', 'LineWidth',2);
            set(gca,'xlim',[4 50], 'TickDir', 'out', 'FontSize', 20); box off; ylim([-0.005 0.04]);
            ylabel('unrew - rew'); xlabel('frequency');
            
            %             %plot phase
            %             figure; hold on; for k=1:nch, plot(freq,phase_ppc(:,k),'Color',cmap(k,:)); end
            %             set(gca,'xlim',[2 50], 'TickDir', 'out', 'FontSize', 20); box off
            %             title(['Phase PPC ' num2str(m)]); xlabel('frequency');  ylabel('rad');
            % plot across distance
            figure('Name',['Monk ' num2str(m)]); hold on;
            plot(dist_ppc,nanmean(coher_err(4:8,:)),'b');  % plot across distance (theta)
            plot(dist_ppc,nanmean(coher_err(8:12,:)),'c');  % plot across distance (alfa)
            plot(dist_ppc,nanmean(coher_err(12:30,:)),'r'); % plot across distance (beta)
            
            plot(dist_ppc,nanmean(coher_corr(4:8,:)),'--b');  % plot across distance (theta)
            plot(dist_ppc,nanmean(coher_corr(8:12,:)),'--c');  % plot across distance (alfa)
            plot(dist_ppc,nanmean(coher_corr(12:30,:)),'--r'); % plot across distance (beta)
            set(gca,'xlim',[2 max(round(dist_ppc))],'xTick',[2 max(round(dist_ppc))], 'TickDir', 'out', 'FontSize', 20); box off; ylim([0.75 0.95]);
            xlabel('electrode distance'); ylabel('coherence');
            legend({'theta', 'alpha', 'beta'}, 'box', 'off'); title(['Coh across distance PPC Monk ' num2str(m)]);
            
            figure('Name',['Monk ' num2str(m)]);hold on;
            shadedErrorBar(dist_ppc,nanmean(coher_err(4:8,:)),nanstd(coher_err(4:8,:)),'lineprops','b');  % plot across distance (theta)
            shadedErrorBar(dist_ppc,nanmean(coher_err(8:12,:)),nanstd(coher_err(8:12,:)),'lineprops','c');  % plot across distance (alfa)
            shadedErrorBar(dist_ppc,nanmean(coher_err(12:30,:)),nanstd(coher_err(12:30,:)),'lineprops','r'); % plot across distance (beta)
            
            shadedErrorBar(dist_ppc,nanmean(coher_corr(4:8,:)),nanstd(coher_corr(4:8,:)),'lineprops','--b');  % plot across distance (theta)
            shadedErrorBar(dist_ppc,nanmean(coher_corr(8:12,:)),nanstd(coher_corr(8:12,:)),'lineprops','--c');  % plot across distance (alfa)
            shadedErrorBar(dist_ppc,nanmean(coher_corr(12:30,:)),nanstd(coher_corr(12:30,:)),'lineprops','--r'); % plot across distance (beta)
            
            set(gca,'xlim',[2 12],'ylim',[0.7 1],'xTick',[2 12], 'TickDir', 'out', 'FontSize', 20); box off
            xlabel('electrode distance'); ylabel('coherence');
            legend({'theta', 'alpha', 'beta'}, 'box', 'off'); title(['Coh across distance PPC Monk ' num2str(m)]);
            
            % mean coherence in different bands for all sessions in one monkey
            c_theta = []; c_alfa = []; c_beta = [];
            for s = 1:length(coh_err(:,1,1))
                c_theta = [c_theta ; nanmean(squeeze(coh_err(1,4:7,2:end)))'];
                c_alfa = [c_alfa ; nanmean(squeeze(coh_err(1,8:11,2:end)))'];
                c_beta = [c_beta ; nanmean(squeeze(coh_err(1,12:30,2:end)))'];
            end
            
            figure; hold on;
            plot(1,c_theta,'.k'); plot(1,nanmean(c_theta), '.r', 'MarkerSize',20);
            plot(2,c_alfa,'.k'); plot(2,nanmean(c_alfa), '.r', 'MarkerSize',20);
            plot(3,c_beta,'.k'); plot(3,nanmean(c_beta), '.r', 'MarkerSize',20);
            set(gca,'xlim',[0.5 3.5],'xTick',[1 2 3],'xTickLabel',[{'theta'} {'alfa'} {'beta'}], 'TickDir', 'out', 'FontSize', 20); box off
            ylim([0.6 1]); ylabel('Coherence'); title(['PPC Monk ' num2str(m)]);
        end
        
        
        
    case 'sim_coherence_between_areas'
        %% MST to PPC
        for m = 1:length(monk)
            % plot per monkey
            clear mst_ppc_sess ppc_mst_sess
            for sess = 1:length(monk(m).pop)
                mst_ppc_sess(sess,:,:) = monk(m).pop(sess).trialtype.all.MSTPPC.coherByElectrode;
                ppc_mst_sess(sess,:,:) = monk(m).pop(sess).trialtype.all.PPCMST.coherByElectrode;
            end
            mst_ppc = squeeze(nanmean(mst_ppc_sess)); mst_ppc_std = squeeze(nanstd(mst_ppc_sess));
            ppc_mst = squeeze(nanmean(ppc_mst_sess)); ppc_mst_std = squeeze(nanstd(ppc_mst_sess));
            
            freq = monk(1).pop(1).trialtype.all.crosslfp.freq;
            % plot each session to check
            
            
            %% MST to PPC
            figure('Name',['Monk ' num2str(m)]); hold on; cmap = jet(24);
            for k=1:24, plot(freq,mst_ppc(:,k),'Color',cmap(k,:)); end
            set(gca,'xlim',[2 50],'ylim',[0.5 0.7], 'TickDir', 'out', 'FontSize', 20); box off
            title('MST --> PPC'); xlabel('frequency');  ylabel('coherence');
            
            %% PPC to MST
            if size(ppc_mst,2) == 48
                figure('Name',['Monk ' num2str(m)]); hold on; cmap = jet(48);
                for k=1:48, plot(freq,ppc_mst(:,k),'Color',cmap(k,:)); end
                set(gca,'xlim',[2 50],'ylim',[0.5 0.7], 'TickDir', 'out', 'FontSize', 22); box off
                title('PPC --> MST'); xlabel('frequency');  ylabel('coherence');
            else
                figure('Name',['Monk ' num2str(m)]); hold on; cmap = jet(96);
                for k=1:96, plot(freq,ppc_mst(:,k),'Color',cmap(k,:)); end
                set(gca,'xlim',[2 50],'ylim',[0.5 0.7], 'TickDir', 'out', 'FontSize', 22); box off
                title('PPC --> MST'); xlabel('frequency');  ylabel('coherence');
            end
            % plot mean for both
            figure; hold on;
            plot(freq,mean(mst_ppc,2),'b', 'LineWidth',2);
            plot(freq,mean(ppc_mst,2),'m','LineWidth',2);
            set(gca,'xlim',[2 50],'ylim',[0.5 0.7], 'TickDir', 'out', 'FontSize', 20); box off
            legend({'MST->PPC', 'PPC->MST'}, 'box', 'off'); title(['Coherence across areas Monk ' num2str(m)]);
            
        end
        
    case 'speed_dependent_LFP_activity_MST'
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for ch = 1:length(monk(1).cont.MST.sess(sess).lfps)
                    monk(m).sess(sess).theta_v(ch,:) = monk(m).cont.MST.sess(sess).lfps(ch).trialtype.all.continuous.v.thetafreq.tuning.rate.mu;
                    monk(m).sess(sess).theta_w(ch,:) = monk(m).cont.MST.sess(sess).lfps(ch).trialtype.all.continuous.w.thetafreq.tuning.rate.mu;
                    
                    monk(m).sess(sess).alpha_v(ch,:) = monk(m).cont.MST.sess(sess).lfps(ch).trialtype.all.continuous.v.alphafreq.tuning.rate.mu;
                    monk(m).sess(sess).alpha_w(ch,:) = monk(m).cont.MST.sess(sess).lfps(ch).trialtype.all.continuous.w.alphafreq.tuning.rate.mu;
                    
                    monk(m).sess(sess).beta_v(ch,:) = monk(m).cont.MST.sess(sess).lfps(ch).trialtype.all.continuous.v.betafreq.tuning.rate.mu;
                    monk(m).sess(sess).beta_w(ch,:) = monk(m).cont.MST.sess(sess).lfps(ch).trialtype.all.continuous.w.betafreq.tuning.rate.mu;
                    
                end
            end
        end
        % average across channels
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                monk(m).sess(sess).theta_v_mu = nanmean(monk(m).sess(sess).theta_v);
                monk(m).sess(sess).theta_w_mu = nanmean(monk(m).sess(sess).theta_w);
                monk(m).sess(sess).alpha_v_mu = nanmean(monk(m).sess(sess).alpha_v);
                monk(m).sess(sess).alpha_w_mu = nanmean(monk(m).sess(sess).alpha_w);
                monk(m).sess(sess).beta_v_mu = nanmean(monk(m).sess(sess).beta_v);
                monk(m).sess(sess).beta_w_mu = nanmean(monk(m).sess(sess).beta_w);
            end
        end
        
        % average across sessions
        for m = 1:length(monk)
            clear th_v th_w bet_v bet_w
            for sess = 1:length(monk(m).sess)
                th_v(sess,:) = monk(m).sess(sess).theta_v_mu;
                th_w(sess,:) = monk(m).sess(sess).theta_w_mu;
                al_v(sess,:) = monk(m).sess(sess).alpha_v_mu;
                al_w(sess,:) = monk(m).sess(sess).alpha_w_mu;
                bet_v(sess,:) = monk(m).sess(sess).beta_v_mu;
                bet_w(sess,:) = monk(m).sess(sess).beta_w_mu;
            end
            %
            monk(m).speed.theta_v = nanmean(th_v);
            monk(m).speed.theta_w = nanmean(th_w);
            monk(m).speed.alpha_v = nanmean(al_v);
            monk(m).speed.alpha_w = nanmean(al_w);
            monk(m).speed.beta_v = nanmean(bet_v);
            monk(m).speed.beta_w = nanmean(bet_w);
        end
        
        v = monk(1).cont.MST.sess(1).lfps(1).trialtype.all.continuous.v.thetafreq.tuning.stim.mu;
        w = monk(1).cont.MST.sess(1).lfps(1).trialtype.all.continuous.w.thetafreq.tuning.stim.mu;
        
        % plot monk separate -- theta
        for m = 1:length(monk)
            figure('Position',[1797 490 824 339]); subplot(1,2,1); hold on; plot(w,monk(m).speed.theta_w,'.','color','k','MarkerSize', 16);
            xlabel('Angular velocity (deg/s)'); ylabel('\theta - frequency (Hz)');
            set(gca,'xlim', [-80 80],'xTick',[-75 75],'TickDir', 'out', 'FontSize', 22);
            
            subplot(1,2,2); hold on; plot(v,monk(m).speed.theta_v,'.','color','k','MarkerSize', 16);
            xlabel('Linear velocit(cm/s)'); ylabel('\theta - frequency (Hz)');
            set(gca, 'xlim', [0 200], 'TickDir', 'out', 'FontSize', 22);
        end
        
        for m = 1:length(monk) % -- alpha
            figure('Position',[1797 490 824 339]); subplot(1,2,1); hold on; plot(w,monk(m).speed.alpha_w,'.','color','k','MarkerSize', 16);
            xlabel('Angular velocity (deg/s)'); ylabel('\alpha - frequency (Hz)');
            set(gca,'xlim', [-80 80],'xTick',[-75 75],'TickDir', 'out', 'FontSize', 22);
            
            subplot(1,2,2); hold on; plot(v,monk(m).speed.alpha_v,'.','color','k','MarkerSize', 16);
            xlabel('Linear velocit(cm/s)'); ylabel('\alpha - frequency (Hz)');
            set(gca, 'xlim', [0 200], 'TickDir', 'out', 'FontSize', 22);
        end
        
        for m = 1:length(monk) % -- beta
            figure('Position',[1797 490 824 339]); subplot(1,2,1); hold on; plot(w,monk(m).speed.beta_w,'.','color','k','MarkerSize', 16);
            xlabel('Angular velocity (deg/s)'); ylabel('\beta - frequency (Hz)');
            set(gca,'xlim', [-80 80],'xTick',[-75 75],'TickDir', 'out', 'FontSize', 22);
            
            subplot(1,2,2); hold on; plot(v,monk(m).speed.beta_v,'.','color','k','MarkerSize', 16);
            xlabel('Linear velocit(cm/s)'); ylabel('\beta - frequency (Hz)');
            set(gca, 'xlim', [0 200], 'TickDir', 'out', 'FontSize', 22);
        end
        
    case 'speed_dependent_LFP_activity_PPC'
        for m = 1:length(monk)
            for sess = 1:length(monk(m).pop)
                for ch = 1:length(monk(1).cont.PPC.sess(sess).lfps)
                    monk(m).sess(sess).theta_v(ch,:) = monk(m).cont.PPC.sess(sess).lfps(ch).trialtype.all.continuous.v.thetafreq.tuning.rate.mu;
                    monk(m).sess(sess).theta_w(ch,:) = monk(m).cont.PPC.sess(sess).lfps(ch).trialtype.all.continuous.w.thetafreq.tuning.rate.mu;
                    
                    monk(m).sess(sess).alpha_v(ch,:) = monk(m).cont.PPC.sess(sess).lfps(ch).trialtype.all.continuous.v.alphafreq.tuning.rate.mu;
                    monk(m).sess(sess).alpha_w(ch,:) = monk(m).cont.PPC.sess(sess).lfps(ch).trialtype.all.continuous.w.alphafreq.tuning.rate.mu;
                    
                    monk(m).sess(sess).beta_v(ch,:) = monk(m).cont.PPC.sess(sess).lfps(ch).trialtype.all.continuous.v.betafreq.tuning.rate.mu;
                    monk(m).sess(sess).beta_w(ch,:) = monk(m).cont.PPC.sess(sess).lfps(ch).trialtype.all.continuous.w.betafreq.tuning.rate.mu;
                    
                end
            end
        end
        % average across channels
        for m = 1:length(monk)
            for sess = 1:length(monk(m).sess)
                monk(m).sess(sess).theta_v_mu = nanmean(monk(m).sess(sess).theta_v);
                monk(m).sess(sess).theta_w_mu = nanmean(monk(m).sess(sess).theta_w);
                monk(m).sess(sess).alpha_v_mu = nanmean(monk(m).sess(sess).alpha_v);
                monk(m).sess(sess).alpha_w_mu = nanmean(monk(m).sess(sess).alpha_w);
                monk(m).sess(sess).beta_v_mu = nanmean(monk(m).sess(sess).beta_v);
                monk(m).sess(sess).beta_w_mu = nanmean(monk(m).sess(sess).beta_w);
            end
        end
        
        % average across sessions
        for m = 1:length(monk)
            clear th_v th_w bet_v bet_w
            for sess = 1:length(monk(m).sess)
                th_v(sess,:) = monk(m).sess(sess).theta_v_mu;
                th_w(sess,:) = monk(m).sess(sess).theta_w_mu;
                al_v(sess,:) = monk(m).sess(sess).alpha_v_mu;
                al_w(sess,:) = monk(m).sess(sess).alpha_w_mu;
                bet_v(sess,:) = monk(m).sess(sess).beta_v_mu;
                bet_w(sess,:) = monk(m).sess(sess).beta_w_mu;
            end
            %
            monk(m).speed.theta_v = nanmean(th_v);
            monk(m).speed.theta_w = nanmean(th_w);
            monk(m).speed.alpha_v = nanmean(al_v);
            monk(m).speed.alpha_w = nanmean(al_w);
            monk(m).speed.beta_v = nanmean(bet_v);
            monk(m).speed.beta_w = nanmean(bet_w);
        end
        
        v = monk(1).cont.PPC.sess(1).lfps(1).trialtype.all.continuous.v.thetafreq.tuning.stim.mu;
        w = monk(1).cont.PPC.sess(1).lfps(1).trialtype.all.continuous.w.thetafreq.tuning.stim.mu;
        
        % plot monk separate -- theta
        for m = 1:length(monk)
            figure('Position',[1797 490 824 339]); subplot(1,2,1); hold on; plot(w,monk(m).speed.theta_w,'.','color','k','MarkerSize', 16);
            xlabel('Angular velocity (deg/s)'); ylabel('\theta - frequency (Hz)');
            set(gca,'xlim', [-80 80],'xTick',[-75 75],'TickDir', 'out', 'FontSize', 22);
            
            subplot(1,2,2); hold on; plot(v,monk(m).speed.theta_v,'.','color','k','MarkerSize', 16);
            xlabel('Linear velocit(cm/s)'); ylabel('\theta - frequency (Hz)');
            set(gca, 'xlim', [0 200], 'TickDir', 'out', 'FontSize', 22);
        end
        
        for m = 1:length(monk) % -- alpha
            figure('Position',[1797 490 824 339]); subplot(1,2,1); hold on; plot(w,monk(m).speed.alpha_w,'.','color','k','MarkerSize', 16);
            xlabel('Angular velocity (deg/s)'); ylabel('\alpha - frequency (Hz)');
            set(gca,'xlim', [-80 80],'xTick',[-75 75],'TickDir', 'out', 'FontSize', 22);
            
            subplot(1,2,2); hold on; plot(v,monk(m).speed.alpha_v,'.','color','k','MarkerSize', 16);
            xlabel('Linear velocit(cm/s)'); ylabel('\alpha - frequency (Hz)');
            set(gca, 'xlim', [0 200], 'TickDir', 'out', 'FontSize', 22);
        end
        
        for m = 1:length(monk) % -- beta
            figure('Position',[1797 490 824 339]); subplot(1,2,1); hold on; plot(w,monk(m).speed.beta_w,'.','color','k','MarkerSize', 16);
            xlabel('Angular velocity (deg/s)'); ylabel('\beta - frequency (Hz)');
            set(gca,'xlim', [-80 80],'xTick',[-75 75],'TickDir', 'out', 'FontSize', 22);
            
            subplot(1,2,2); hold on; plot(v,monk(m).speed.beta_v,'.','color','k','MarkerSize', 16);
            xlabel('Linear velocit(cm/s)'); ylabel('\beta - frequency (Hz)');
            set(gca, 'xlim', [0 200], 'TickDir', 'out', 'FontSize', 22);
        end
        
        
    case 'trial_band_passed'
        area = 'PPC'     % MST PPC PFC
        
        for m = 3  % 1:length(monk)
            nsess = 1:length(monk(m).session);
            for sess = 1:length(nsess)
                corr = monk(m).behavior(sess).stats.trialtype.reward(strcmp({monk(m).behavior(sess).stats.trialtype.reward.val},'rewarded')).trlindx;
                incorr = monk(m).behavior(sess).stats.trialtype.reward(strcmp({monk(m).behavior(sess).stats.trialtype.reward.val},'unrewarded')).trlindx;
                behv_corr = monk(m).behavior(sess).trials(corr);
                behv_incorr = monk(m).behavior(sess).trials(incorr);
                pick_ch = randperm(length(monk(m).session(sess).area.(area).lfp)); 
                %p = numSubplots(6); 
                cnt=1;
                for nlfp = pick_ch(1:5) %1:length(monk(m).session(sess).area.(area).lfp)
                    trl_corr = monk(m).session(sess).area.(area).lfp(nlfp).trials(corr);
                    trl_incorr = monk(m).session(sess).area.(area).lfp(nlfp).trials(incorr);
                    % extract t_stop corr
                    clear trl_end_corr trl_end_incorr
                    for j = 1:length(behv_corr), trl_end_corr(j) = behv_corr(j).events.t_end; end
                    for j = 1:length(behv_incorr), trl_end_incorr(j) = behv_incorr(j).events.t_end; end
                    %% corr
                    % align to t_stop and extract around stop time
                    for trl = 1:length(behv_corr)
                        this_ts = behv_corr(trl).continuous.ts-trl_end_corr(trl);
                        if this_ts(end)>0.5
                            this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                            this_lfp_corr = trl_corr(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_corr(trl,:) = this_lfp_corr(1:333,:);  %% Choose band
                            this_ninety = prctile(real(beta_corr(trl,:)),95); beta_corr_indx_90(trl,:) = real(beta_corr(trl,:))>this_ninety;
                            trial_corr(trl).tspk = ts(trl,beta_corr_indx_90(trl,:))';
                        else
                            beta_corr(trl,:) = NaN(1,333); 
                        end
                    end
                     %% take a psth of the 95th prcnt
                    [rate_90_corr, ts_90_corr] = Spiketimes2Rate(trial_corr, [-1:0.02:1], 0.02);
                    figure; hold on; colormap(jet);title('Correct')
                    %subplot(p(1),p(2),cnt); hold on; colormap(jet)
                    %subplot(1,5,cnt); hold on; colormap(jet)
                    for trl = 16:36 % 1:length(behv_corr)
                        % plot(ts(trl,:),real(beta_corr(trl,:)))
                        %plot(ts(trl,:),real(beta_corr(trl,:))/max(abs(real(beta_corr(trl,:))))) % normalized
                        %                         waterfall(ts(trl,:),trl,real(beta_corr(trl,:))/max(abs(real(beta_corr(trl,:)))));
                        %                         waterfall(ts(trl,:),trl,abs(beta_corr(trl,:))/max(abs(real(beta_corr(trl,:)))));
                       
                        waterfall(ts(trl,:),trl,real(beta_corr(trl,:)));
                        waterfall(ts(trl,:),trl,abs(beta_corr(trl,:)));
                        set(gca,'xTick',[-1 0 1], 'yTick', [], 'zTick',[])
                    end
                    view(gca,[0 64.4]); grid off; axis square; vline(0,'k')
                    cnt=cnt+1;
                    % plot raster of time points above 95
                    figure; hold on; 
                    for trl = 1:size(beta_corr_indx_90,1)
                    plot(ts(trl,beta_corr_indx_90(trl,:)),trl,'.k', 'MarkerSize',4)
                    end
                    vline(0,'-r')
                    
                    %% incorrect 
                     for trl = 1:length(behv_incorr)
                        this_ts = behv_incorr(trl).continuous.ts-trl_end_incorr(trl);
                        if this_ts(end)>0.5
                            this_ts_win = this_ts(this_ts>-1 & this_ts<1); ts(trl,:) = this_ts_win(1:333,:);
                            this_lfp_incorr = trl_incorr(trl).lfp_beta(this_ts>-1 & this_ts<1); beta_incorr(trl,:) = this_lfp_incorr(1:333,:);  %% Choose band
                            this_ninety = prctile(real(beta_incorr(trl,:)),95); beta_incorr_indx_90(trl,:) = real(beta_incorr(trl,:))>this_ninety;
                            trial_incorr(trl).tspk = ts(trl,beta_incorr_indx_90(trl,:))';
                        else
                            beta_incorr(trl,:) = NaN(1,333); 
                        end
                     end
                    %% take a psth of the 95th prcnt
                    [rate_90_incorr, ts_90_incorr] = Spiketimes2Rate(trial_incorr, [-1:0.02:1], 0.02);
                    %%
                    figure; hold on; colormap(jet); title('Incorrect')
                    %subplot(p(1),p(2),cnt); hold on; colormap(jet)
                    %subplot(1,5,cnt); hold on; colormap(jet)
                    for trl = 16:36 % 1:length(behv_corr)
                        % plot(ts(trl,:),real(beta_corr(trl,:)))
                        %plot(ts(trl,:),real(beta_corr(trl,:))/max(abs(real(beta_corr(trl,:))))) % normalized
                        %                         waterfall(ts(trl,:),trl,real(beta_corr(trl,:))/max(abs(real(beta_corr(trl,:)))));
                        %                         waterfall(ts(trl,:),trl,abs(beta_corr(trl,:))/max(abs(real(beta_corr(trl,:)))));
                       
                        waterfall(ts(trl,:),trl,real(beta_incorr(trl,:)));
                        waterfall(ts(trl,:),trl,abs(beta_incorr(trl,:)));
                        set(gca,'xTick',[-1 0 1], 'yTick', [], 'zTick',[])
                    end
                    view(gca,[0 64.4]); grid off; axis square; vline(0,'k')
                    % plot raster of time points above 95
                    figure; hold on; 
                    for trl = 1:size(beta_corr_indx_90,2)
                    plot(ts(trl,beta_incorr_indx_90(trl,:)),trl,'.k', 'MarkerSize',4)
                    end
                    vline(0,'-r')
                    %plot psths
                    figure; hold on
                    plot(ts_90_corr,smooth(rate_90_corr),'b')
                    plot(ts_90_incorr,smooth(rate_90_incorr),'r')
                    set(gca,'TickDir', 'out', 'FontSize', 22); axis square; box off
                    title([area ' Channel ' num2str(nlfp)])
                    cnt=cnt+1;
                end
            end
        end
        
        
        
        
        
        
    case 'pop_psth_band_passed'
        % load behavioral file if saved separately
        areaToLoad = 'PPC'
        monkey = '1'; 
        band = 'beta'
        win = [-0.5 0.5];
        %% load and extract each file per area
        fprintf(['Time:  ' num2str(clock) '\n']);
        fnames = dir(['monk ' (monkey) ' * ' (areaToLoad) ' *.mat']);
        r_corr = []; r_incorr = []; peak_t_corr = []; peak_t_incorr = [];
        for sess = 1:length(fnames)
            load(fnames(sess).name)
            monk = dataToSave;
            clear dataToSave
            ts = monk.area.(areaToLoad).band_passed.(band).lfp(1).psth.ts_95_corr; ts_win_indx = ts>=win(1) & ts<=win(2); ts_win = ts(ts>=win(1) & ts<=win(2)); 
            for ch = 1:length(monk.area.(areaToLoad).band_passed.(band).lfp)
                r_corr(ch,:,sess) = monk.area.(areaToLoad).band_passed.(band).lfp(ch).psth.rate_95_corr; [~,indx_max_corr] = max(r_corr(ch,ts_win_indx)); peak_t_corr(ch,sess) = ts_win(indx_max_corr);
                r_incorr(ch,:,sess) = monk.area.(areaToLoad).band_passed.(band).lfp(ch).psth.rate_95_incorr; [~,indx_max_incorr] = max(r_incorr(ch,ts_win_indx)); peak_t_incorr(ch,sess) = ts_win(indx_max_incorr);
            end
            
            % mean and sem
            r_corr_mu(sess,:) = mean(r_corr(:,:,sess));
            r_corr_sem(sess,:) = std(r_corr(:,:,sess))/sqrt(size(r_corr(:,:,sess),1));
            r_incorr_mu(sess,:) = mean(r_incorr(:,:,sess));
            r_incorr_sem(sess,:) = std(r_incorr(:,:,sess))/sqrt(size(r_incorr(:,:,sess),1));
            %% plot all channels
            figure(1); hold on
            subplot(1,length(fnames),sess); hold on
            for ii = 1:size(r_corr(:,:,sess),1), plot(ts,smooth(r_corr(ii,:,sess)), 'g'); end
            for ii = 1:size(r_incorr(:,:,sess),1)plot(ts,smooth(r_incorr(ii,:,sess)), 'k'); end
            set(gca, 'xlim', [-0.95 0.95], 'TickDir', 'out', 'FontSize', 22); axis square; box off
            vline(0,'-r'); title([(areaToLoad) ' ' num2str(sess)])
            
            % plot mean +/- sem psth per area for corr and incorr
            figure(2); hold on
            subplot(1,length(fnames),sess); hold on
            shadedErrorBar(ts,smooth(r_corr_mu(sess,:)),r_corr_sem(sess,:), 'lineprops', 'g')
            shadedErrorBar(ts,smooth(r_incorr_mu(sess,:)),r_incorr_sem(sess,:), 'lineprops', 'k')
            set(gca, 'xlim', [-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); axis square; box off
            title([(areaToLoad) ' ' num2str(sess)])
            xlabel('Time(s)'); vline(0,'-r');
            
%             figure(3); hold on
%             shadedErrorBar(ts,smooth(r_corr_mu(sess,:)),r_corr_sem(sess,:), 'lineprops', 'g')
%             shadedErrorBar(ts,smooth(r_incorr_mu(sess,:)),r_incorr_sem(sess,:), 'lineprops', 'k')
%             set(gca, 'xlim', [-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); axis square; box off
%             title(['MST ' num2str(sess)])
%             xlabel('Time(s)'); vline(0,'-r');
            %% bar mean plots
            figure(4); hold on
            errorbar(sess, mean(peak_t_corr(:,sess)),std(peak_t_corr(:,sess)), 'og', 'MarkerSize', 20)
            errorbar(sess+0.25, mean(peak_t_incorr(:,sess)),std(peak_t_incorr(:,sess)), 'ok', 'MarkerSize', 20)
            set(gca, 'xlim', [0.5 sess+0.5], 'TickDir', 'out', 'FontSize', 22); axis square; box off
           
            %% Ratio plot
            figure(5); hold on
            subplot(1,length(fnames),sess); hold on
            plot(ts,smooth(r_corr_mu(sess,:))./smooth(r_incorr_mu(sess,:)), 'k', 'LineWidth',1)
            set(gca, 'xlim', [-0.5 0.5],'ylim',[0 4],'yTick',[0 1 4], 'TickDir', 'out', 'FontSize', 22); axis square; box off
            ylabel('Correct / Incorrect'); xlabel('Time(s)')
            hline(1,'--r')
            
            % clear variable and load next one
        end
        % plot summary for all sessions for this monkey 
        % just bar plot and ratio
        %% ratios
        figure(6); hold on
        for sess = 1:length(fnames)
            plot(ts,smooth(r_corr_mu(sess,:))./smooth(r_incorr_mu(sess,:)),'Color', [0.5 0.5 0.5], 'LineWidth',1);
        end
        corr_mu = mean(r_corr_mu); incorr_mu = mean(r_incorr_mu); 
        plot(ts,smooth(corr_mu./incorr_mu), 'k', 'LineWidth',2);
        set(gca, 'xlim', [-0.5 0.5],'ylim',[0.5 4], 'TickDir', 'out', 'FontSize', 22); axis square; box off
        ylabel('Correct / Incorrect'); xlabel('Time(s)'); %hline(1,'--r')
        % save pop data for this monkey
        save(['monkey_' num2str(monkey) '_' areaToLoad '_band_' band '_r_pop'],'r_corr_mu','r_incorr_mu','peak_t_corr', 'peak_t_incorr')

        
    case 'psth_peak_pop'
        areaToLoad = 'PPC'
        band = 'beta'
        %% load and extract each file per area
        fprintf(['Time:  ' num2str(clock) '\n']);
        fnames = dir(['monkey*_' (areaToLoad) '_band_' (band) '_*.mat']);
        for m = 1:length(fnames)
            load(fnames(m).name) % load
            % extract
            corr_mu(m,:) = mean(r_corr_mu); corr_sem(m,:) = std(r_corr_mu)/sqrt(size(r_corr_mu,1)); 
            peak_corr_mu(m,:) = mean(mean(peak_t_corr)); peak_corr_sem(m,:) = mean(std(peak_t_corr)/sqrt(size(peak_t_corr,1)));
            
            incorr_mu(m,:) = mean(r_incorr_mu); incorr_sem(m,:) = std(r_incorr_mu)/sqrt(size(r_incorr_mu,1));
            peak_incorr_mu(m,:) = mean(mean(peak_t_incorr)); peak_incorr_sem(m,:) = mean(std(peak_t_incorr)/sqrt(size(peak_t_incorr,1)));
        end
        
        % plot ratios for all monkeys in one area
        load('ts_psth')
        figure(1); hold on;
        for m = 1:length(fnames)
            plot(ts,smooth(corr_mu(m,:)./incorr_mu(m,:)), 'Color', [0.5 0.5 0.5], 'LineWidth',0.5);
        end
        plot(ts,smooth(mean(corr_mu)./mean(incorr_mu)), 'k', 'LineWidth',2); hline(1,'r')
        set(gca,'ylim',[0 3], 'TickDir', 'out', 'FontSize', 22); axis square; box off
        ylabel('Correct / Incorrect'); xlabel('Time(s)'); title(band)
        
        % plot max peaks for all monkeys in one area
        figure(2); hold on; cnt=1;
        col = {'ks' 'ko' 'kv' }
        for m = 1:length(fnames)
           errorbar(cnt,peak_corr_mu(m),peak_corr_sem(m), col{m}, 'MarkerSize', 14)
           cnt=cnt+0.15;
        end
        set(gca,'xlim',[0.5 2],'ylim',[-0.3 0.3], 'TickDir', 'out', 'FontSize', 22); axis square; box off
        ylabel('Time (s)'); xlabel(''); title(band)
        
    case 'band_passed_vs_accuracy'
        % needs single session files "...band_passed area..."
        areaToLoad = 'PPC'
        monkey = '3';
        band = 'beta'
        win = [-0.5 0.5];
        %% load and extract each file per area
        fprintf(['Time:  ' num2str(clock) '\n']);
        fnames = dir(['monk ' (monkey) ' * ' (areaToLoad) ' *.mat']);
        r_corr_low = []; r_corr_middle = []; r_corr_upper = [];  peak_t_corr_low = []; peak_t_corr_middle = []; peak_t_corr_upper = [];
        for sess = 1:length(fnames)
            load(fnames(sess).name)
            monk = dataToSave;
            clear dataToSave
            ts = monk.area.(areaToLoad).band_passed.(band).lfp(1).psth.middle.ts_95_corr; ts_win_indx = ts>=win(1) & ts<=win(2); ts_win = ts(ts>=win(1) & ts<=win(2));
            for ch = 1:length(monk.area.(areaToLoad).band_passed.(band).lfp)
                r_corr_low(ch,:,sess) = monk.area.(areaToLoad).band_passed.(band).lfp(ch).psth.low.rate_95_corr_low; [~,indx_max_corr_low] = max(r_corr_low(ch,ts_win_indx)); peak_t_corr_low(ch,sess) = ts_win(indx_max_corr_low);
                r_corr_middle(ch,:,sess) = monk.area.(areaToLoad).band_passed.(band).lfp(ch).psth.middle.rate_95_corr; [~,indx_max_corr_middle] = max(r_corr_middle(ch,ts_win_indx)); peak_t_corr_middle(ch,sess) = ts_win(indx_max_corr_middle);
                r_corr_upper(ch,:,sess) = monk.area.(areaToLoad).band_passed.(band).lfp(ch).psth.upper.rate_95_corr; [~,indx_max_corr_upper] = max(r_corr_upper(ch,ts_win_indx)); peak_t_corr_upper(ch,sess) = ts_win(indx_max_corr_upper);
            end
            % mean and sem
            r_corr_low_mu(sess,:) = mean(r_corr_low(:,:,sess)); r_corr_low_sem(sess,:) = std(r_corr_low(:,:,sess))/sqrt(size(r_corr_low(:,:,sess),1));
            r_corr_middle_mu(sess,:) = mean(r_corr_middle(:,:,sess)); r_corr_middle_sem(sess,:) = std(r_corr_middle(:,:,sess))/sqrt(size(r_corr_middle(:,:,sess),1));
            r_corr_upper_mu(sess,:) = mean(r_corr_upper(:,:,sess)); r_corr_upper_sem(sess,:) = std(r_corr_upper(:,:,sess))/sqrt(size(r_corr_upper(:,:,sess),1));
            %% plot all channels
            
            figure(1); hold on
            subplot(1,length(fnames),sess); hold on
            for ii = 1:size(r_corr_low(:,:,sess),1), plot(ts,smooth(r_corr_low(ii,:,sess)), 'Color', [0 1 0]); end
            for ii = 1:size(r_corr_middle(:,:,sess),1), plot(ts,smooth(r_corr_middle(ii,:,sess)), 'Color', [0 0.65 0]); end
            %for ii = 1:size(r_corr_upper(:,:,sess),1), plot(ts,smooth(r_corr_upper(ii,:,sess)), 'Color', [0 0.25 0]); end
            set(gca, 'xlim', [-0.95 0.95], 'TickDir', 'out', 'FontSize', 22); axis square; box off
            vline(0,'-r'); title([(areaToLoad) ' ' num2str(sess)])
            
            % plot mean +/- sem psth per area for corr and incorr
            figure(2); hold on
            subplot(1,length(fnames),sess); hold on
            shadedErrorBar(ts,smooth(r_corr_low_mu(sess,:),5),r_corr_low_sem(sess,:), 'lineprops', {'Color',[0 1 0]})
            shadedErrorBar(ts,smooth(r_corr_middle_mu(sess,:)),r_corr_middle_sem(sess,:), 'lineprops', {'Color',[0 0.65 0]})
            %shadedErrorBar(ts,smooth(r_corr_upper_mu(sess,:)),r_corr_upper_sem(sess,:), 'lineprops', {'Color',[0 0.25 0]})
            set(gca, 'xlim', [-0.5 0.5], 'TickDir', 'out', 'FontSize', 22); axis square; box off
            title([(areaToLoad) ' ' num2str(sess)])
            xlabel('Time(s)'); vline(0,'-r');
            
            %% bar mean plots (calculate latency to peak)
            figure(4); hold on
            errorbar(sess, mean(peak_t_corr_low(:,sess)),std(peak_t_corr_low(:,sess)),'o', 'Color',[0 1 0], 'MarkerSize', 20)
            errorbar(sess+0.25, mean(peak_t_corr_middle(:,sess)),std(peak_t_corr_middle(:,sess)),'o', 'Color',[0 0.65 0], 'MarkerSize', 20)
            %errorbar(sess+0.55, mean(peak_t_corr_upper(:,sess)),std(peak_t_corr_upper(:,sess)), 'ok', 'MarkerSize', 20)
            set(gca, 'xlim', [0.5 sess+0.5], 'TickDir', 'out', 'FontSize', 22); axis square; box off
        end
         z=1; 
   
     case 'distance_vs_reward'
        figure; hold on;
        r_targ = behv.stats.pos_final.r_targ(incorrect);
        r_monk = behv.stats.pos_final.r_monk(incorrect);
        if ntrls > maxtrls
            trl_indx = randperm(ntrls);
            trl_indx = trl_indx(1:maxtrls);
            plot(r_targ(trl_indx), r_monk(trl_indx), '.r','markersize',4);
        else
            plot(r_targ, r_monk, '.r','markersize',4);
        end
        
        r_targ = behv.stats.pos_final.r_targ(correct);
        r_monk = behv.stats.pos_final.r_monk(correct);
        if ntrls > maxtrls
            trl_indx = randperm(ntrls);
            trl_indx = trl_indx(1:maxtrls);
            plot(r_targ(trl_indx), r_monk(trl_indx), '.b','markersize',4);
        else
            plot(r_targ, r_monk, '.b','markersize',4);
        end
        
        axis([0 400 0 400]);
        plot(0:400,0:400,'--k','Linewidth',1);
        set(gca, 'XTick', [0 200 400], 'XTickLabel', [0 2 4], 'Fontsize',14);
        xlabel('Target, r(m)');
        set(gca, 'YTick', [0 200 400], 'YTickLabel', [0 2 4]);
        ylabel('Response, r(m)');

end