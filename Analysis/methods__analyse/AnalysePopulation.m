function stats = AnalysePopulation(units,trials_behv,behv_stats,lfps,prs)

nunits = length(units);
dt = prs.dt; % sampling resolution (s)
fs = 500;
%% which analayses to do
fitGAM_coupled = prs.fitGAM_coupled;
compute_canoncorr = prs.compute_canoncorr;
regress_popreadout = prs.regress_popreadout;
simulate_population = prs.simulate_population;
compute_coherencyLFP = prs.compute_coherencyLFP;

%% load cases
trialtypes = fields(behv_stats.trialtype);
events = cell2mat({trials_behv.events});
continuous = cell2mat({trials_behv.continuous});

%% fit GAM with cross-neuronal coupling
if fitGAM_coupled
    stats = [];
    GAM_prs.varname = prs.GAM_varname; varname = GAM_prs.varname;
    GAM_prs.vartype = prs.GAM_vartype; vartype = GAM_prs.vartype;
    GAM_prs.nbins = prs.GAM_nbins;
    GAM_prs.binrange = [];
    GAM_prs.nfolds = prs.nfolds;
    GAM_prs.dt = dt;
    GAM_prs.filtwidth = prs.neuralfiltwidth;
    GAM_prs.linkfunc = prs.GAM_linkfunc;
    GAM_prs.lambda = prs.GAM_lambda;
    GAM_prs.alpha = prs.GAM_alpha;
    GAM_prs.varchoose = prs.GAM_varchoose;
    GAM_prs.method = prs.GAM_method;
    for i=1% if i=1, fit model using data from all trials rather than separately to data from each condition
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        fprintf(['.........fitting GAM model :: trialtype: ' (trialtypes{i}) '\n']);
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).models.(GAM_prs.linkfunc) = stats.trialtype.all.models.(GAM_prs.linkfunc);
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                if ~isempty(lfps)
                    lfps_temp = lfps([units.channel_id]);
                    trials_lfps_temp = cell(1,nunits);
                    for l = 1:nunits, trials_lfps_temp{l} = lfps_temp(l).trials(trlindx); end
                end
                %% select variables of interest and load their details
                vars = cell(length(varname),1);
                GAM_prs.binrange = cell(1,length(varname));
                for k=1:length(varname)
                    if isfield(continuous_temp,varname(k)), vars{k} = {continuous_temp.(varname{k})};
                    elseif isfield(behv_stats.pos_rel,varname(k)), vars{k} = behv_stats.pos_rel.(varname{k})(trlindx);
                    elseif strcmp(varname(k),'d')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'phi')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'eye_ver')
                        isnan_le = all(isnan(cell2mat({continuous_temp.zle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.zre}')));
                        if isnan_le, vars{k} = {continuous_temp.zre};
                        elseif isnan_re, vars{k} = {continuous_temp.zle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'eye_hor')
                        isnan_le = all(isnan(cell2mat({continuous_temp.yle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.yre}')));
                        if isnan_le, vars{k} = {continuous_temp.yre};
                        elseif isnan_re, vars{k} = {continuous_temp.yle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'phase')
                        vars{k} = [];
                        var_phase = cell(1,nunits);
                        for l = 1:nunits, var_phase{l} = cellfun(@(x) angle(hilbert(x)), {trials_lfps_temp{l}.lfp},'UniformOutput',false); end
                    elseif strcmp(vartype(k),'event')
                        vars{k} = [events_temp.(prs.varlookup(varname{k}))];
                        if strcmp(varname(k),'target_OFF'), vars{k} = vars{k} + prs.fly_ONduration; end % target_OFF = t_targ + fly_ONduration
                    end
                    GAM_prs.binrange{k} = prs.binrange.(varname{k});
                end
                %% define time windows for computing tuning
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                timewindow_full = [min([events_temp.t_move],[events_temp.t_targ]) - prs.pretrial ;... % from "min(move,targ) - pretrial_buffer"
                    [events_temp.t_end] + prs.posttrial]'; % till "end + posttrial_buffer"
                %% concatenate stimulus data from all trials
                xt = []; yt = [];
                trials_spks_temp = units(1).trials(trlindx); % dummy
                for k=1:length(vars)
                    if strcmp(varname(k),'phase')
                        xt(:,k) = nan; % works as long as 'phase' is not the first entry in varname
                    elseif ~strcmp(vartype(k),'event')
                        [xt(:,k),~,yt] = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    elseif ~strcmp(varname(k),'spikehist')
                        [~,xt(:,k),yt] = ConcatenateTrials([],mat2cell(vars{k}',ones(length(events_temp),1)),{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    end
                end
                if any(strcmp(varname,'spikehist')), xt(:,strcmp(varname,'spikehist')) = yt; end % pass spike train back as an input to fit spike-history kernel
                %% concatenate units
                Yt = zeros(size(xt,1),nunits);
                for k=1:nunits
                    trials_spks_temp = units(k).trials(trlindx);
                    [~,~,Yt(:,k)] = ConcatenateTrials(vars{1},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                end
                %% fit fully coupled GAM model to each unit
                for k=1:nunits
                    % replace xt(:,'phase') with the unit-specific LFP channel
                    xt(:,strcmp(varname,'phase')) = ConcatenateTrials(var_phase{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    xt_k = mat2cell(xt,size(xt,1),ones(1,size(xt,2))); % convert to cell
                    models = BuildGAMCoupled(xt_k,Yt(:,k),Yt(:,[1:k-1 k+1:nunits]),GAM_prs);
                    stats.trialtype.(trialtypes{i})(j).models.(GAM_prs.linkfunc).units(k) = models;
                end
            end
        end
    end
end

%% cannonical correlation analysis
if compute_canoncorr
    varname = prs.canoncorr_varname;
    filtwidth = prs.neuralfiltwidth;
    for i=1% if i=1, fit model using data from all trials rather than separately to data from each condition
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        fprintf(['.........computing canonical correlations :: trialtype: ' (trialtypes{i}) '\n']);
        for j=1:nconds
            stats.trialtype.(trialtypes{i})(j).canoncorr.vars = varname;
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).canoncorr = stats.trialtype.all.canoncorr;
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                %% select variables of interest and load their details
                vars = cell(length(varname),1);
                for k=1:length(varname)
                    if isfield(continuous_temp,varname{k}), vars{k} = {continuous_temp.(varname{k})};
                    elseif isfield(behv_stats.pos_rel,varname{k}), vars{k} = behv_stats.pos_rel.(varname{k})(trlindx);
                    elseif strcmp(varname{k},'d')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname{k},'dv')
                        vars{k} = cellfun(@(x) [0 ; diff(x)/dt],{continuous_temp.v},'UniformOutput',false);
                    elseif strcmp(varname{k},'dw')
                        vars{k} = cellfun(@(x) [0 ; diff(x)/dt],{continuous_temp.w},'UniformOutput',false);
                    elseif strcmp(varname{k},'phi')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'eye_ver')
                        isnan_le = all(isnan(cell2mat({continuous_temp.zle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.zre}')));
                        if isnan_le, vars{k} = {continuous_temp.zre};
                        elseif isnan_re, vars{k} = {continuous_temp.zle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'eye_hor')
                        isnan_le = all(isnan(cell2mat({continuous_temp.yle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.yre}')));
                        if isnan_le, vars{k} = {continuous_temp.yre};
                        elseif isnan_re, vars{k} = {continuous_temp.yle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false);
                        end
                    end
                end
                %% define time windows for computing tuning
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                timewindow_full = [min([events_temp.t_move],[events_temp.t_targ]) - prs.pretrial ;... % from "min(move,targ) - pretrial_buffer"
                    [events_temp.t_end] + prs.posttrial]'; % till "end + posttrial_buffer"
                %% concatenate stimulus data from all trials
                trials_spks_temp = units(1).trials(trlindx);
                xt = [];
                for k=1:length(vars)
                    xt(:,k) = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    xt(isnan(xt(:,k)),k) = 0;
                end
                %% concatenate units
                Yt = zeros(size(xt,1),nunits);
                for k=1:nunits
                    trials_spks_temp = units(k).trials(trlindx);
                    [~,~,Yt(:,k)] = ConcatenateTrials(vars{1},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                end
                Yt_smooth = SmoothSpikes(Yt, 3*filtwidth);
                %% compute canonical correlation
                [X,Y,R,~,~,pstats] = canoncorr(xt,Yt_smooth/dt);
                stats.trialtype.(trialtypes{i})(j).canoncorr.stim = X;
                stats.trialtype.(trialtypes{i})(j).canoncorr.resp = Y;
                stats.trialtype.(trialtypes{i})(j).canoncorr.coeff = R;
                stats.trialtype.(trialtypes{i})(j).canoncorr.pval = pstats;
                %% canonical task dimensionality
                Dmax = numel(R); taskcov = zeros(1,Dmax);
                for k=1:Dmax
                    wx = X(:,k)/sqrt(X(:,k)'*X(:,k));
                    wy = Y(:,k)/sqrt(Y(:,k)'*Y(:,k));
                    xt_proj = xt*wx; yt_proj = (Yt_smooth/dt)*wy;
                    cov_temp = cov([xt_proj yt_proj]);
                    taskcov(k) = cov_temp(1,2); % off-diagonal entry
                end
                stats.trialtype.(trialtypes{i})(j).canoncorr.dimensionality = sum(taskcov)^2/sum(taskcov.^2); % defined analogously to participation ratio
                %% compute pairwise correlations
                stats.trialtype.(trialtypes{i})(j).responsecorr_raw = corr(Yt);
                stats.trialtype.(trialtypes{i})(j).responsecorr_smooth = corr(Yt_smooth);
                %% compute pairwise noise correlations
                for k=1:numel(varname)
                    binrange{k} = prs.binrange.(varname{k});
                    nbins{k} = prs.GAM_nbins{strcmp(prs.GAM_varname,varname{k})};
                end
                [stats.trialtype.(trialtypes{i})(j).noisecorr_raw,...
                    stats.trialtype.(trialtypes{i})(j).noisecorr_smooth] = ComputeNoisecorr(Yt,xt,binrange,nbins);
            end
        end
    end
end

%% compute population readout weights via ordinary-least-squares
if regress_popreadout
    varname = prs.readout_varname;
    decodertype = prs.decodertype;
    filtwidth = prs.neuralfiltwidth;
    for i=1 % i=1 means compute only for all trials together
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        for j=1:nconds
            trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
            events_temp = events(trlindx);
            continuous_temp = continuous(trlindx);
            %% define time windows for decoding
            timewindow_move = [[events_temp.t_move]' [events_temp.t_stop]']; % only consider data when the subject is integrating path
            timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % only consider data when the subject is integrating path
            timewindow_full = [min([events_temp.t_move],[events_temp.t_targ]) - prs.pretrial ;... % from "min(move,targ) - pretrial_buffer"
                [events_temp.t_end] + prs.posttrial]'; % till "end + posttrial_buffer"
            nunits = length(units);
            %% gather spikes from all units
            Yt = [];
            for k=1:nunits
                trials_spks_temp = units(k).trials(trlindx);
                [~,~,Yt(:,k)] = ConcatenateTrials({continuous_temp.v},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
            end
            %% build decoder for each variable
            vars = cell(length(varname),1);
            trials_spks_temp = units(k).trials(trlindx);
            for k=1:length(varname)
                fprintf(['Building decoder for ' varname{k} '...\n']);
                if isfield(continuous_temp,varname{k}), vars{k} = {continuous_temp.(varname{k})};
                elseif isfield(behv_stats.pos_rel,varname{k}), vars{k} = behv_stats.pos_rel.(varname{k})(trlindx);
                elseif strcmp(varname{k},'d')
                    vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                elseif strcmp(varname{k},'dv')
                    vars{k} = cellfun(@(x) [0 ; diff(x)/dt],{continuous_temp.v},'UniformOutput',false);
                elseif strcmp(varname{k},'dw')
                    vars{k} = cellfun(@(x) [0 ; diff(x)/dt],{continuous_temp.w},'UniformOutput',false);
                elseif strcmp(varname{k},'phi')
                    vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                elseif strcmp(varname(k),'eye_ver')
                    isnan_le = all(isnan(cell2mat({continuous_temp.zle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.zre}')));
                    if isnan_le, vars{k} = {continuous_temp.zre};
                    elseif isnan_re, vars{k} = {continuous_temp.zle};
                    else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false);
                    end
                elseif strcmp(varname(k),'eye_hor')
                    isnan_le = all(isnan(cell2mat({continuous_temp.yle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.yre}')));
                    if isnan_le, vars{k} = {continuous_temp.yre};
                    elseif isnan_re, vars{k} = {continuous_temp.yle};
                    else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false);
                    end
                end
                xt = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                % filter dv and dw
                t = linspace(-2*filtwidth,2*filtwidth,4*filtwidth + 1); h = exp(-t.^2/(2*filtwidth^2)); h = h/sum(h);
                if any(strcmp(varname{k},{'dv','dw'})), xt = conv(xt,h,'same'); end
                xt(isnan(xt)) = 0;
                %% fit smoothing window
                if prs.lineardecoder_fitkernelwidth
                    filtwidths=1:5:100; decodingerror = nan(1,length(filtwidths));
                    fprintf('...optimising hyperparameter\n');
                    for l=1:length(filtwidths)
                        Yt_temp = SmoothSpikes(Yt, filtwidths(l)); % smooth spiketrains before fitting model
                        wts = (Yt_temp'*Yt_temp)\(Yt_temp'*xt); % analytical
                        decodingerror(l) = sqrt(sum((Yt_temp*wts - xt).^2));
                    end
                    [~,bestindx] = min(decodingerror); bestfiltwidth = filtwidths(bestindx);
                    stats.(decodertype).(varname{k}).bestfiltwidth = bestfiltwidth;
                    fprintf('**********decoding**********\n');
                    Yt_temp = SmoothSpikes(Yt, bestfiltwidth); % smooth spiketrains before fitting model
                    stats.(decodertype).(varname{k}).wts = (Yt_temp'*Yt_temp)\(Yt_temp'*xt); % analytical
                    stats.(decodertype).(varname{k}).true = xt;
                    stats.(decodertype).(varname{k}).pred = (Yt_temp*stats.(decodertype).(varname{k}).wts);
                    stats.(decodertype).(varname{k}).corr = corr(stats.(decodertype).(varname{k}).true,stats.(decodertype).(varname{k}).pred);
                end
                %% subsample neurons
                if prs.lineardecoder_subsample
                    N_neurons = prs.N_neurons; N_samples = prs.N_neuralsamples; Nt = size(Yt,1);
                    for l=1:numel(N_neurons)
                        fprintf(['.........decoding ' num2str(N_neurons(l)) ' neuron(s) \n']);
                        if N_neurons(l)< nunits, sampleindx = cell2mat(arrayfun(@(x) randperm(nunits,N_neurons(l))',1:N_samples,'UniformOutput',false));
                        else, sampleindx = [repmat((1:nunits)',1,N_samples) ; randi(nunits,[N_neurons(l)-nunits N_samples])]; end
                        Yt_temp = reshape(Yt(:,sampleindx),[Nt N_neurons(l) N_samples]);
                        Yt_temp = SmoothSpikes(Yt_temp, stats.(decodertype).(varname{k}).bestfiltwidth);
                        for m=1:N_samples
                            stats.(decodertype).(varname{k}).corr_subsample(l,m) = ...
                                corr(xt,squeeze(Yt_temp(:,:,m))*((squeeze(Yt_temp(:,:,m))'*squeeze(Yt_temp(:,:,m)))\(squeeze(Yt_temp(:,:,m))'*xt)));
                            stats.(decodertype).(varname{k}).popsize_subsample(l,m) = N_neurons(l);
                        end
                    end
                end
                %% fixed filtwidth
                Yt_temp = SmoothSpikes(Yt, 3*filtwidth); % smooth spiketrains before fitting model
                stats.(decodertype).(varname{k}).wts2 = (Yt_temp'*Yt_temp)\(Yt_temp'*xt); % analytical
                stats.(decodertype).(varname{k}).true2 = xt;
                stats.(decodertype).(varname{k}).pred2 = (Yt_temp*stats.(decodertype).(varname{k}).wts2);
                stats.(decodertype).(varname{k}).corr2 = corr(stats.(decodertype).(varname{k}).true2,stats.(decodertype).(varname{k}).pred2);
            end
        end
    end
end

%% evaluate model responses for coupled and uncoupled models
if simulate_population && exist('stats','var') && isfield(stats.trialtype.all,'models')
    varname = prs.simulate_varname;
    vartype = prs.simulate_vartype;
    filtwidth = prs.neuralfiltwidth;
    for i=1% if i=1, fit model using data from all trials rather than separately to data from each condition
        nconds = length(behv_stats.trialtype.(trialtypes{i}));
        if ~strcmp((trialtypes{i}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        fprintf(['.........computing canonical correlations for actual and model-simulated responses :: trialtype: ' (trialtypes{i}) '\n']);
        for j=1:nconds
            if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
                stats.trialtype.(trialtypes{i})(j).canoncorr = stats.trialtype.all.canoncorr;
            else
                trlindx = behv_stats.trialtype.(trialtypes{i})(j).trlindx;
                events_temp = events(trlindx);
                continuous_temp = continuous(trlindx);
                %% select variables of interest and load their details
                vars = cell(length(varname),1); binrange = cell(length(varname),1); nbins = cell(length(varname),1);
                for k=1:length(varname)
                    if isfield(continuous_temp,varname{k}), vars{k} = {continuous_temp.(varname{k})};
                    elseif isfield(behv_stats.pos_rel,varname{k}), vars{k} = behv_stats.pos_rel.(varname{k})(trlindx);
                    elseif strcmp(varname{k},'d')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.v},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname{k},'phi')
                        vars{k} = cellfun(@(x,y) [zeros(sum(y<=0),1) ; cumsum(x(y>0)*dt)],{continuous_temp.w},{continuous_temp.ts},'UniformOutput',false);
                    elseif strcmp(varname(k),'eye_ver')
                        isnan_le = all(isnan(cell2mat({continuous_temp.zle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.zre}')));
                        if isnan_le, vars{k} = {continuous_temp.zre};
                        elseif isnan_re, vars{k} = {continuous_temp.zle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.zle},{continuous_temp.zre},'UniformOutput',false);
                        end
                    elseif strcmp(varname(k),'eye_hor')
                        isnan_le = all(isnan(cell2mat({continuous_temp.yle}'))); isnan_re = all(isnan(cell2mat({continuous_temp.yre}')));
                        if isnan_le, vars{k} = {continuous_temp.yre};
                        elseif isnan_re, vars{k} = {continuous_temp.yle};
                        else, vars{k} = cellfun(@(x,y) 0.5*(x + y),{continuous_temp.yle},{continuous_temp.yre},'UniformOutput',false);
                        end
                    end
                    binrange{k} = prs.binrange.(varname{k});
                    nbins{k} = prs.GAM_nbins{strcmp(prs.GAM_varname,varname{k})};
                end
                %% define time windows for computing tuning
                timewindow_path = [[events_temp.t_targ]' [events_temp.t_stop]']; % when the subject is integrating path
                timewindow_full = [min([events_temp.t_move],[events_temp.t_targ]) - prs.pretrial ;... % from "min(move,targ) - pretrial_buffer"
                    [events_temp.t_end] + prs.posttrial]'; % till "end + posttrial_buffer"
                %% concatenate stimulus data from all trials
                trials_spks_temp = units(1).trials(trlindx);
                xt = [];
                for k=1:length(vars)
                    xt(:,k) = ConcatenateTrials(vars{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
                    xt(isnan(xt(:,k)),k) = 0;
                end
                %% encode stimulus as one-hot variables
                x_1hot = [];
                for k=1:length(vars), x_1hot(:,:,k) = Encode1hot(xt,vartype{k},binrange{k},nbins{k}); end
                %% simulate uncoupled model
                Yt_uncoupled = zeros(size(xt,1),nunits);
                for k=1:nunits
                    y_temp = zeros(length(vars),size(xt,1));
                    uncoupledmodel = stats.trialtype.all.models.log.units(k).Uncoupledmodel;
                    if ~isnan(uncoupledmodel.bestmodel), wts = uncoupledmodel.wts{uncoupledmodel.bestmodel};
                    else, wts = uncoupledmodel.wts{1}; end
                    for l = 1:length(vars)
                        % simulated response to each variable before exp
                        y_temp(l,:) = sum(repmat(wts{strcmp(prs.GAM_varname,varname{l})},[size(x_1hot,1),1]).*squeeze(x_1hot(:,:,l)),2);
                    end
                    y_temp = exp(sum(y_temp));
                    Yt_uncoupled(:,k) = y_temp;
                end
                Yt_uncoupled_smooth = SmoothSpikes(Yt_uncoupled, 3*filtwidth);
                %% compute canonical correlation
                [X_uncoupled,Y_uncoupled,R_uncoupled,~,~,pstats_uncoupled] = canoncorr(xt,Yt_uncoupled_smooth/dt);
                stats.trialtype.(trialtypes{i})(j).canoncorr.uncoupled_stim = X_uncoupled;
                stats.trialtype.(trialtypes{i})(j).canoncorr.uncoupled_resp = Y_uncoupled;
                stats.trialtype.(trialtypes{i})(j).canoncorr.uncoupled_coeff = R_uncoupled;
                stats.trialtype.(trialtypes{i})(j).canoncorr.uncoupled_pval = pstats_uncoupled;
                %% canonical task dimensionality
                Dmax = numel(R_uncoupled); taskcov = zeros(1,Dmax);
                for k=1:Dmax
                    wx = X_uncoupled(:,k)/sqrt(X_uncoupled(:,k)'*X_uncoupled(:,k));
                    wy = Y_uncoupled(:,k)/sqrt(Y_uncoupled(:,k)'*Y_uncoupled(:,k));
                    xt_proj = xt*wx; yt_proj = (Yt_uncoupled_smooth/dt)*wy;
                    cov_temp = cov([xt_proj yt_proj]);
                    taskcov(k) = cov_temp(1,2); % off-diagonal entry
                end
                stats.trialtype.(trialtypes{i})(j).canoncorr.uncoupled_dimensionality = sum(taskcov)^2/sum(taskcov.^2); % defined analogously to participation ratio
                %% compute pairwise correlations
                stats.trialtype.(trialtypes{i})(j).responsecorr_uncoupled_raw = corr(Yt_uncoupled);
                stats.trialtype.(trialtypes{i})(j).responsecorr_uncoupled_smooth = corr(SmoothSpikes(Yt_uncoupled, filtwidth));
                %% compute pairwise noise correlations
                for k=1:numel(varname)
                    binrange{k} = prs.binrange.(varname{k});
                    nbins{k} = prs.GAM_nbins{strcmp(prs.GAM_varname,varname{k})};
                end
                [stats.trialtype.(trialtypes{i})(j).noisecorr_uncoupled_raw,...
                    stats.trialtype.(trialtypes{i})(j).noisecorr_uncoupled_smooth] = ComputeNoisecorr(Yt_uncoupled,xt,binrange,nbins);
                %% simulate coupled model
                Yt_coupled = zeros(size(xt,1),nunits);
                for k=1:nunits
                    y_temp = zeros(length(vars),size(xt,1));
                    coupledmodel = stats.trialtype.all.models.log.units(k).Coupledmodel;
                    wts = coupledmodel.wts;
                    for l = 1:length(vars)
                        % simulated response to each variable before exp
                        y_temp(l,:) = sum(repmat(wts{strcmp(prs.GAM_varname,varname{l})},[size(x_1hot,1),1]).*squeeze(x_1hot(:,:,l)),2);
                    end
                    Yt_coupled(:,k) = exp(sum(y_temp)' + ...
                        sum(Yt(:,1:nunits~=k).*repmat(stats.trialtype.all.models.log.units(k).Coupledmodel.wts{end},[size(Yt,1), 1]),2));
                end
                Yt_coupled_smooth = SmoothSpikes(Yt_coupled, 3*filtwidth);
                %% compute canonical correlation
                [X_coupled,Y_coupled,R_coupled,~,~,pstats_coupled] = canoncorr(xt,Yt_coupled_smooth/dt);
                stats.trialtype.(trialtypes{i})(j).canoncorr.coupled_stim = X_coupled;
                stats.trialtype.(trialtypes{i})(j).canoncorr.coupled_resp = Y_coupled;
                stats.trialtype.(trialtypes{i})(j).canoncorr.coupled_coeff = R_coupled;
                stats.trialtype.(trialtypes{i})(j).canoncorr.coupled_pval = pstats_coupled;
                %% canonical task dimensionality
                Dmax = numel(R_coupled); taskcov = zeros(1,Dmax);
                for k=1:Dmax
                    wx = X_coupled(:,k)/sqrt(X_coupled(:,k)'*X_coupled(:,k));
                    wy = Y_coupled(:,k)/sqrt(Y_coupled(:,k)'*Y_coupled(:,k));
                    xt_proj = xt*wx; yt_proj = (Yt_coupled_smooth/dt)*wy;
                    cov_temp = cov([xt_proj yt_proj]);
                    taskcov(k) = cov_temp(1,2); % off-diagonal entry
                end
                stats.trialtype.(trialtypes{i})(j).canoncorr.coupled_dimensionality = sum(taskcov)^2/sum(taskcov.^2); % defined analogously to participation ratio
                %% compute pairwise correlations
                stats.trialtype.(trialtypes{i})(j).responsecorr_coupled_raw = corr(Yt_coupled);
                stats.trialtype.(trialtypes{i})(j).responsecorr_coupled_smooth = corr(SmoothSpikes(Yt_coupled, filtwidth));
                %% compute pairwise noise correlations
                for k=1:numel(varname)
                    binrange{k} = prs.binrange.(varname{k});
                    nbins{k} = prs.GAM_nbins{strcmp(prs.GAM_varname,varname{k})};
                end
                [stats.trialtype.(trialtypes{i})(j).noisecorr_coupled_raw,...
                    stats.trialtype.(trialtypes{i})(j).noisecorr_coupled_smooth] = ComputeNoisecorr(Yt_coupled,xt,binrange,nbins);
            end
        end
    end
end


%% coherence between LFPs
if compute_coherencyLFP
    for type = 1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{type}));
        %if ~strcmp((trialtypes{type}),'all') && nconds==1, copystats = true; else, copystats = false; end % only one condition means variable was not manipulated
        for cond = 1:nconds
            %             if copystats % if only one condition present, no need to recompute stats --- simply copy them from 'all' trials
            %                 stats.trialtype.(trialtypes{i})(j) = stats.trialtype.all;
            %             else
            clear trlindx lfp_concat triallen sMarkers
            trlindx = behv_stats.trialtype.(trialtypes{type})(cond).trlindx;
            lfp_concat = nan(length(cell2mat({units(1).trials(trlindx).lfp}')),nunits);
            % params
            spectralparams.tapers = prs.spectrum_tapers;
            spectralparams.Fs = 1/dt;
            spectralparams.trialave = prs.spectrum_trialave;
            % data
            for i=1:nunits, lfp_concat(:,i) = cell2mat({units(i).trials(trlindx).lfp}'); end
            triallen = cellfun(@(x) length(x), {units(1).trials(trlindx).lfp});
            sMarkers(:,1) = cumsum([1 triallen(1:end-1)]); sMarkers(:,2) = cumsum(triallen); % demarcate trial onset and end
            fprintf('**********Computing pairwise coherence between LFPs********** \n');
            fprintf(['Time:  ' num2str(clock) '\n']);
            [stats.trialtype.(trialtypes{type})(cond).crosslfp.coher , stats.trialtype.(trialtypes{type})(cond).crosslfp.phase, ~, ~, stats.trialtype.(trialtypes{type})(cond).crosslfp.freq] = ...
                coherencyc_unequal_length_trials(lfp_concat, prs.spectrum_movingwin , spectralparams, sMarkers); % needs http://chronux.org/
            % store output of coherencyc_unequal_length_trials to 3-D matrix (freq x electrode_site x electrode_site)
            ind2row = @(i,j) min(i,j) + (max(i,j)-1)*(max(i,j)-2)/2; % to read the output of "coherencyc_unequal_length_trials" function from Chronux
            spatial_coher = []; spatial_phase = [];
            for i=2:nunits
                for j=1:i-1
                    stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,i,j) = stats.trialtype.(trialtypes{type})(cond).crosslfp.coher(:,ind2row(i,j));
                    stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,i,j) = stats.trialtype.(trialtypes{type})(cond).crosslfp.phase(:,ind2row(i,j));
                end
            end
            stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,end,end+1) = 0; % pad column of zeros to squarify the matrix
            stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,end,end+1) = 0;
            
            %% split electrode pairs within same area based on brain area
            unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
            for p = 1:num_brain_areas
                unitindx = strcmp({units.brain_area}, unique_brain_areas{p});
                coherMat = stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,unitindx,unitindx);
                phaseMat = stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,unitindx,unitindx);
                stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).electrode_type = cell2mat(unique({units(unitindx).electrode_type}));
                stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).electrode_ids = [units(unitindx).electrode_id];
                coherMatFull = coherMat + permute(coherMat,[1 3 2]); coherMatFull(coherMatFull == 0) = NaN; % symmetrify
                stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).coherByElectrode = nanmean(coherMatFull,3);  % save freq x elec
                phaseMatFull = coherMat + permute(phaseMat,[1 3 2]); phaseMatFull(phaseMatFull == 0) = NaN; % symmetrify
                stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).phaseByElectrode = nanmean(phaseMatFull,3);  % save freq x elec
                stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).coherMatFull = coherMatFull; % save freq x elec x elec
                stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).phaseMatFull = phaseMatFull; % save freq x elec x elec
                [stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).coherByDist,stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).phaseByDist,stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).dist] = ...
                    ComputeCoherenceByDistance(coherMat,phaseMat,stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).electrode_ids,stats.trialtype.(trialtypes{type})(cond).(unique_brain_areas{p}).electrode_type);  % save freq x dist
            end
            %% analyse cross-area
            if num_brain_areas==2
                unitindx1 = strcmp({units.brain_area}, unique_brain_areas(1));
                unitindx2 = strcmp({units.brain_area}, unique_brain_areas(2));
                
                if max(find(unitindx1)) > min(find(unitindx2)), stats.trialtype.(trialtypes{type})(cond).crossarea.coher12 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,unitindx1,unitindx2),3));
                else, stats.trialtype.(trialtypes{type})(cond).crossarea.coher12 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,unitindx2,unitindx1),2)); end
                
                if max(find(unitindx2)) > min(find(unitindx1)), stats.trialtype.(trialtypes{type})(cond).crossarea.coher21 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,unitindx2,unitindx1),3));
                else, stats.trialtype.(trialtypes{type})(cond).crossarea.coher21 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher(:,unitindx1,unitindx2),2)); end
                
                if max(find(unitindx1)) > min(find(unitindx2)), stats.trialtype.(trialtypes{type})(cond).crossarea.phase12 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,unitindx1,unitindx2),3));
                else, stats.trialtype.(trialtypes{type})(cond).crossarea.phase12 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,unitindx2,unitindx1),2)); end
                
                if max(find(unitindx2)) > min(find(unitindx1)), stats.trialtype.(trialtypes{type})(cond).crossarea.phase21 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,unitindx2,unitindx1),3));
                else, stats.trialtype.(trialtypes{type})(cond).crossarea.phase21 = squeeze(mean(stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_phase(:,unitindx1,unitindx2),2)); end
                
            end
            %% split electrode pairs across two areas based on brain areas
            coherMat = stats.trialtype.(trialtypes{type})(cond).crosslfp.spatial_coher;
            coherMatFull = coherMat + permute(coherMat,[1 3 2]); coherMatFull(coherMatFull == 0) = NaN; % symmetrify
            for area1 = 1:num_brain_areas
                for area2 = find(1:num_brain_areas ~= area1)
                    unitindx1 = strcmp({units.brain_area}, unique_brain_areas{area1});
                    unitindx2 = strcmp({units.brain_area}, unique_brain_areas{area2});
                    
                    stats.trialtype.(trialtypes{type})(cond).([unique_brain_areas{area2} unique_brain_areas{area1}]).electrode_ids = [units(unitindx2).electrode_id];
                    stats.trialtype.(trialtypes{type})(cond).([unique_brain_areas{area2} unique_brain_areas{area1}]).electrode_type = cell2mat(unique({units(unitindx2).electrode_type}));
                    stats.trialtype.(trialtypes{type})(cond).([unique_brain_areas{area2} unique_brain_areas{area1}]).coherByElectrode = ...
                        squeeze(mean(coherMatFull(:,unitindx1,unitindx2),2));
                end
            end
        end
        
        
        
    end
    
end

%% compute spectrogram for each trial for all channels (average all channels per trial) and per channel aligned to target/move/stop/reward onset
if prs.compute_spectrum_whole_trial
    fprintf('**********Computing spectrogram for the whole trial between LFPs********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    gettuning = prs.tuning_events; unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    for type = 2 %1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{type}));
        for cond = 1:nconds
            for ev = 1:length(gettuning)
                for area = 1:num_brain_areas
                    unitindx = strcmp({units.brain_area}, unique_brain_areas{area});
                    %% extract LFP trace for each trial for all channels in each brain area
                    if ~isempty(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align)
                        for j = 1:length(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align(1,:)) % trials
                            ar = find(unitindx); clear lfp_trl_area
                            for n = 1:length(ar)  % first 24 ch for MST if applicable
                                lfp_trl_area(n,:) =  lfps(ar(n)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align(:,j);  % extract lfp for all ch per trial  % get # of trials (samp x ch) 334xch
                                % compute spectrogram for each trial for each channel
                                %                                 [stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).ch(n).trl(j).spectrogram, stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).ch(n).trl(j).ts_spectrogram,...
                                %                                     stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).ch(n).trl(j).freq_spectrogram] = ...
                                %                                     mtspecgramc(lfp_trl_area(n,:),prs.spectrogram_movingwin,spectralparams);
                            end
                            % compute spectrogram for all channels for each trial
                            [stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).spectrogram,...
                                stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).ts_spectrogram,...
                                stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).freq_spectrogram] = ...
                                mtspecgramc(lfp_trl_area',prs.spectrogram_movingwin,spectralparams);
                            
                            stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).ts_spectrogram = ...
                                stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).ts_spectrogram-...
                                abs(lfps(ar(n)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.ts_lfp_align(end));
                            %                                                         figure; imagesc((stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).ts_spectrogram),stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).freq_spectrogram,...
                            %                                                             stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).events.(gettuning{ev}).pop_trl(j).spectrogram'); axis xy;
                            %                                                         set(gca, 'ylim',[0 50])
                        end
                    end
                end
            end
        end
    end
end


%% compute spectrogram for each trial for all channels (average all channels per trial) and per channel aligned to movement stop
if prs.compute_spectrum_whole_trial_align_stop
    fprintf('**********Computing spectrogram for the whole trial between LFPs********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    gettuning = prs.tuning_events; unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    for type = 2 %1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{type}));
        for cond = 1:nconds %1:2 %nconds
            %             for ev = 3 % ('stop') %1:length(gettuning)
            for area = 1:num_brain_areas
                unitindx = strcmp({units.brain_area}, unique_brain_areas{area});
                %% extract LFP trace for each trial for all channels in each brain area
                if ~isempty(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.stop.all_freq.lfp_align)
                    for j = 1:length(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.stop.all_freq.lfp_align(1,:)) % trials
                        ar = find(unitindx); clear lfp_trl_area
                        for n = 1:length(ar)  % first 24 ch for MST if applicable
                            lfp_trl_area(n,:) =  lfps(ar(n)).stats.trialtype.(trialtypes{type})(cond).events.stop.all_freq.lfp_align(:,j); % extract lfp for all ch per trial  % get # of trials (samp x ch) 334xch
                            % compute spectrogram for each trial for each channel
                            %                             [stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).ch(n).trl(j).spectrogram_stop,...
                            %                                 stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).ch(n).trl(j).ts_spectrogram_stop,...
                            %                                 stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).ch(n).trl(j).freq_spectrogram_stop] = ...
                            %                                 mtspecgramc(lfp_trl_area(n,:),prs.spectrogram_movingwin,spectralparams);
                            %                             stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).ch(n).trl(j).ts_spectrogram_stop = stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).ch(n).trl(j).ts_spectrogram_stop-...
                            %                                 (lfps(ar(n)).stats.trialtype.(trialtypes{type})(cond).events.stop.all_freq.ts_lfp_align(1));
                        end
                        % compute spectrogram
                        [stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).spectrogram_stop,...
                            stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).ts_spectrogram_stop,...
                            stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).freq_spectrogram_stop] = ...
                            mtspecgramc(lfp_trl_area',prs.spectrogram_movingwin,spectralparams);
                        
                        stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).ts_spectrogram_stop = stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).ts_spectrogram_stop-...
                            abs(lfps(ar(n)).stats.trialtype.(trialtypes{type})(cond).events.stop.all_freq.ts_lfp_align(1));
                        
                        %                          figure; imagesc(stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).ts_spectrogram_stop,...
                        %                                              stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).freq_spectrogram_stop,...
                        %                                              stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).pop_trl(j).spectrogram_stop'); axis xy;
                        %                                          set(gca,'xlim',[-1.5 1.5])
                    end
                end
            end
            %             end
        end
    end
end

%% compute spectrogram for each trial for all channels (average all channels per trial) and per channel aligned to movement stop for each band separately
if prs.compute_spectrum_whole_trial_align_stop_per_band
    fprintf('**********Computing spectrogram for the whole trial between LFPs********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    gettuning = prs.tuning_events; unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    for type = 2 %1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{type}));
        for cond = 1:nconds %2 %nconds
            %             for ev = 3 % ('stop') %1:length(gettuning)
            for area = 1:num_brain_areas
                unitindx = strcmp({units.brain_area}, unique_brain_areas{area});
                %% extract LFP trace for each trial for all channels in each brain area
                bands = {'wideband'} %{'theta' 'beta' 'wideband'}
                for b = 1:length(bands)
                    if ~isempty(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.stop.(bands{b}).lfp_align)
                        for j = 1:length(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.stop.(bands{b}).lfp_align(1,:)) % trials
                            ar = find(unitindx); clear lfp_trl_area
                            for n = 1:length(ar)  % first 24 ch for MST if applicable
                                lfp_trl_area(n,:) =  lfps(ar(n)).stats.trialtype.(trialtypes{type})(cond).events.stop.(bands{b}).lfp_align(:,j); % extract lfp for all ch per trial  % get # of trials (samp x ch) 334xch
                                % compute spectrogram for each trial for each channel
                                [stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).(bands{b}).ch(n).trl(j).spectrogram_stop, stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).(bands{b}).ch(n).trl(j).ts_spectrogram_stop,...
                                    stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).(bands{b}).ch(n).trl(j).freq_spectrogram_stop] = ...
                                    mtspecgramc(lfp_trl_area(n,:),prs.spectrogram_movingwin,spectralparams);
                            end
                            % compute spectrogram
                            [stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).(bands{b}).pop_trl(j).spectrogram_stop, stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).(bands{b}).pop_trl(j).ts_spectrogram_stop,...
                                stats.trialtype.(trialtypes{type})(cond).area.(unique_brain_areas{area}).(bands{b}).pop_trl(j).freq_spectrogram_stop] = ...
                                mtspecgramc(lfp_trl_area',prs.spectrogram_movingwin,spectralparams);
                        end
                    end
                end
            end
            %             end
        end
    end
end


%% compute coherence over time
if prs.compute_coherogram
    fprintf('**********Computing coherogram between LFPs********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    gettuning = prs.tuning_events; unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    for type = 2%1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{type}));
        for cond = 1:nconds
            for ev = [2 3] %1:length(gettuning)
                for area1 = 1:num_brain_areas
                    for area2 = find(1:num_brain_areas ~= area1)
                        unitindx1 = strcmp({units.brain_area}, unique_brain_areas{area1});
                        unitindx2 = strcmp({units.brain_area}, unique_brain_areas{area2});
                        clear lfp_trl_a1 lfp_trl_a2 a1_lfp a2_lfp
                        %% extract LFP trace for each trial for all channels in each brain area
                        if ~isempty(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align)
                            for j = 1:length(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align(1,:)) % trials
                                
                                a1 = find(unitindx1);
                                for n = 1:length(a1)  % first 24 ch for MST
                                    lfp_trl_a1(n,j,:) = lfps(a1(n)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align(:,j)';  % extract lfp for all ch per trial  % get # of trials (samp x trl) 334xch
                                end
                                
                                a2 = find(unitindx2);
                                
                                for n2 = 1:length(a2) % ch
                                    lfp_trl_a2(n2,j,:) = lfps(a2(n2)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_align(:,j)';  % extract lfp for all ch per trial  % get # of trials (samp x trl) 334xch
                                end
                            end
                            % average over all channels
                            for j = 1:length(lfp_trl_a1(1,:,1))
                                a1_lfp(j,:) = squeeze(nanmean(lfp_trl_a1(:,j,:)))';
                                a2_lfp(j,:) = squeeze(nanmean(lfp_trl_a2(:,j,:)))';
                            end
                            % compute coherogram MST-PPC
                            [stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).coher,...
                                stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).coherPhi,...
                                stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).crossSpec, ~,~, ...
                                stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_ts, ...
                                stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_freq] = cohgramc(a1_lfp',a2_lfp',prs.spectrogram_movingwin,spectralparams);
                            
                            
                            stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_ts = ...
                                stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_ts-...
                                (lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.ts_lfp_align(1));
                            
                            %                                                     figure; imagesc((stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_ts)-1,...
                            %                                                         stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_freq,...
                            %                                                         stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).([unique_brain_areas{area2} unique_brain_areas{area1}]).coher', [0 1]);
                            %                                                     axis xy; set(gca,'xlim',[-0.7 0.7], 'ylim', [4 50]); colorbar;
                        end
                    end
                end
            end
        end
    end
end

if prs.compute_coherogram_band_passed
    fprintf('**********Computing coherogram between LFPs********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    spectralparams.tapers = prs.spectrum_tapers;
    spectralparams.Fs = 1/dt;
    spectralparams.trialave = prs.spectrum_trialave;
    gettuning = prs.tuning_events; unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    for type = 2%1:length(trialtypes)
        nconds = length(behv_stats.trialtype.(trialtypes{type}));
        for cond = 1:nconds
            for ev = 1:length(gettuning)
                for area1 = 1:num_brain_areas
                    for area2 = find(1:num_brain_areas ~= area1)
                        unitindx1 = strcmp({units.brain_area}, unique_brain_areas{area1});
                        unitindx2 = strcmp({units.brain_area}, unique_brain_areas{area2});
                        %% extract LFP trace for each trial for all channels in each brain area per band
                        bands = {'all_freq','theta', 'beta', 'wideband'};
                        %%%%%%%%%%%%%% Add loop to calculate coherence per band
                        for b = 1:length(bands)
                            clear lfp_trl_a1 lfp_trl_a2 a1_lfp a2_lfp
                            if ~isempty(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).(bands{b}).lfp_align)
                                for j = 1:length(lfps(1).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).(bands{b}).lfp_align(1,:)) % trials
                                    
                                    a1 = find(unitindx1);
                                    for n = 1:length(a1)  % first 24 ch for MST
                                        lfp_trl_a1(n,j,:) = lfps(a1(n)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).(bands{b}).lfp_align(:,j)';  % extract lfp for all ch per trial  % get # of trials (samp x trl) 334xch
                                    end
                                    
                                    a2 = find(unitindx2);
                                    
                                    for n2 = 1:length(a2) % ch
                                        lfp_trl_a2(n2,j,:) = lfps(a2(n2)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).(bands{b}).lfp_align(:,j)';  % extract lfp for all ch per trial  % get # of trials (samp x trl) 334xch
                                    end
                                end
                                % average over all channels
                                for j = 1:length(lfp_trl_a1(1,:,1))
                                    a1_lfp(j,:) = squeeze(nanmean(lfp_trl_a1(:,j,:)))';
                                    a2_lfp(j,:) = squeeze(nanmean(lfp_trl_a2(:,j,:)))';
                                end
                                % compute coherogram MST-PPC
                                [stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).(bands{b}).coher,...
                                    stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).(bands{b}).coherPhi,...
                                    stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).(bands{b}).crossSpec, ~,~, ...
                                    stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).(bands{b}).coher_ts, ...
                                    stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).coherogram.([unique_brain_areas{area2} unique_brain_areas{area1}]).(bands{b}).coher_freq] = cohgramc(a1_lfp',a2_lfp',prs.spectrogram_movingwin,spectralparams);
                                
                                
                                %                                                     figure; imagesc((stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_ts)-1,...
                                %                                                         stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).([unique_brain_areas{area2} unique_brain_areas{area1}]).coher_freq,...
                                %                                                         stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).([unique_brain_areas{area2} unique_brain_areas{area1}]).coher', [0 1]);
                                %                                                     axis xy; set(gca,'xlim',[-0.7 0.7], 'ylim', [4 50]); colorbar;
                            end
                        end
                    end
                end
            end
        end
    end
end



if prs.analyse_band_passed
    fprintf('**********Analyses on band passed LFPs********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    %% all of these are aligned to stop
    unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    % extract correct vs incorrect trials
    corr_indx = behv_stats.trialtype.reward(strcmp({behv_stats.trialtype.reward.val},'rewarded')).trlindx;
    incorr_indx = behv_stats.trialtype.reward(strcmp({behv_stats.trialtype.reward.val},'unrewarded')).trlindx;
    
    % for band passed vs accuracy
    % get monkey and target position only for correct trials
    monk_abs_x = behv_stats.pos_abs.x_monk(corr_indx); ntrls_correct = length(monk_abs_x);
    monk_abs_y = behv_stats.pos_abs.y_monk(corr_indx);
    r_targ = behv_stats.pos_final.r_targ(corr_indx);
    theta_targ = behv_stats.pos_final.theta_targ(corr_indx)*pi/180;
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
    stats.indx_accuracy.low_half_indx = monk2targ < prs.band_pass_acc_thresh_low_third;
    stats.indx_accuracy.upper_half_indx = monk2targ > prs.band_pass_acc_thresh_upper_third;
    
    %% sanity check plots
    % plot histogram
    %             histogram(monk2targ,100);
    %             title('End point to Target')
    % plot monk2targ vs target linear and angular distance.
    % figure; hold on; plot(theta_targ, monk2targ, '.k', 'MarkerSize', 14); xlabel('Target position'); ylabel('Distance to target'); ylim([0 65]); title('Angular distance')
    % figure; hold on; plot(r_targ, monk2targ, '.k', 'MarkerSize', 14); xlabel('Target position'); ylabel('Distance to target'); ylim([0 65]); title('Linear distance')
    gettuning = prs.tuning_events;
    %% add variable with all 95th pct timings per channel per condition
    for area = 1:num_brain_areas
        for ev = 2:length(gettuning)
            unitindx = strcmp({units.brain_area}, unique_brain_areas{area});
            % extract all 95th pct timings per channel
            ar = find(unitindx); corr_trl = find(corr_indx); incorr_trl = find(incorr_indx); low_half_indx = find(stats.indx_accuracy.low_half_indx); upper_half_indx = find(stats.indx_accuracy.upper_half_indx);
            theta_corr.low_half=[]; beta_corr.low_half = []; theta_corr.upper_half=[]; beta_corr.upper_half=[]; corr_theta = []; incorr_theta = []; corr_beta = []; incorr_beta = [];
            corr_theta_bef_stop = []; corr_theta_af_stop = []; incorr_theta_bef_stop = []; incorr_theta_af_stop = []; corr_beta_bef_stop = []; corr_beta_af_stop = []; incorr_beta_bef_stop = []; incorr_beta_af_stop = [];
            %% corr
            for ch = 1:length(ar)  % first 24 ch for MST if applicable
                %% move
                for ntrl = 1:sum(corr_indx)
                    theta.chan(ch).corr.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.trl(ntrl).theta_95_indx)'; % extract timing of value (this is the important val)
                    beta.chan(ch).corr.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.trl(ntrl).beta_95_indx)'; % extract timing of value(this is the important val)
                end
                % add psth per channel
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(theta.chan(ch).corr.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(beta.chan(ch).corr.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                
                %% compute for band passed vs accuracy low half
                for ntrl = 1:sum(stats.indx_accuracy.low_half_indx)
                    theta_corr.low_half.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.trl(low_half_indx(ntrl)).theta_95_indx)'; % extract timing of value (this is the important val)
                    beta_corr.low_half.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.trl(low_half_indx(ntrl)).beta_95_indx)'; % extract timing of value(this is the important val)
                end
                % add psth per channel
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.low_half.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.low_half.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(theta_corr.low_half.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.low_half.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.low_half.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(beta_corr.low_half.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                
                % compute for band passed vs accuracy upper half
                for ntrl = 1:sum(stats.indx_accuracy.upper_half_indx)
                    theta_corr.upper_half.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.trl(upper_half_indx(ntrl)).theta_95_indx)'; % extract timing of value (this is the important val)
                    beta_corr.upper_half.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).corr.trl(upper_half_indx(ntrl)).beta_95_indx)'; % extract timing of value(this is the important val)
                end
                % add psth per channel
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.upper_half.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.upper_half.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(theta_corr.upper_half.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.upper_half.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.upper_half.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(beta_corr.upper_half.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                
                %% incorr
                for ntrl = 1:sum(incorr_indx)
                    theta.chan(ch).incorr.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).err.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).err.trl(ntrl).theta_95_indx)';
                    beta.chan(ch).incorr.trl(ntrl).ch_95th_ts = lfps(ar(ch)).stats.band_passed.(gettuning{ev}).err.ts(lfps(ar(ch)).stats.band_passed.(gettuning{ev}).err.trl(ntrl).beta_95_indx)'; % extract timing of value
                end
                % add psth calc per channel
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(theta.chan(ch).incorr.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.incorr.rate_95(ch,:),stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.incorr.ts_rate_95(ch,:)] = ...
                    BandPassedAmp2Rate(beta.chan(ch).incorr.trl, [-2:prs.temporal_binwidth:2], prs.temporal_binwidth);
                %% gather
                if ev == 2
                    ts = mean(stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.ts_rate_95);
                    % theta
                    corr_theta = [corr_theta ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.rate_95(ch,ts>=0 & ts<=0.3)'];
                    incorr_theta = [incorr_theta ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.rate_95(ch,ts>=0 & ts<=0.3)'];
                    % beta
                    corr_beta = [corr_beta ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.rate_95(ch,ts>=0 & ts<=0.3)'];
                    incorr_beta = [incorr_beta ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.incorr.rate_95(ch,ts>=0 & ts<=0.3)'];
                elseif ev == 3
                    ts = mean(stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.ts_rate_95);
                    % theta
                    corr_theta_bef_stop = [corr_theta_bef_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.rate_95(ch,ts>=-1.5 & ts<=0)'];
                    corr_theta_af_stop = [corr_theta_af_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.corr.rate_95(ch,ts>=0 & ts<=1.5)'];
                    incorr_theta_bef_stop = [incorr_theta_bef_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.rate_95(ch,ts>=-1.5 & ts<=0)'];
                    incorr_theta_af_stop = [incorr_theta_af_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.incorr.rate_95(ch,ts>=0 & ts<=1.5)'];
                    % beta
                    corr_beta_bef_stop = [corr_beta_bef_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.rate_95(ch,ts>=-1.5 & ts<=0)'];
                    corr_beta_af_stop = [corr_beta_af_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.corr.rate_95(ch,ts>=0 & ts<=1.5)'];
                    incorr_beta_bef_stop = [incorr_beta_bef_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.incorr.rate_95(ch,ts>=-1.5 & ts<=0)'];
                    incorr_beta_af_stop = [incorr_beta_af_stop ; stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.incorr.rate_95(ch,ts>=0 & ts<=1.5)'];
                end
            end
            % compute stats between correct vs incorrect responses in each
            % area and using the window of interest
            if ev == 2
                % theta
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.p_val, stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.p_val_flag] = ranksum(corr_theta,incorr_theta);
                % beta
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.p_val, stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.p_val_flag] = ranksum(corr_beta,incorr_beta);
                
            elseif ev == 3
                %theta
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.p_val_bef_stop,stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.p_val_flag_bef_stop] = ranksum(corr_theta_bef_stop, incorr_theta_bef_stop);
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.p_val_af_stop,stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).theta.p_val_flag_af_stop] = ranksum(corr_theta_af_stop, incorr_theta_af_stop);
                %beta
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.p_val_bef_stop,stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.p_val_flag_bef_stop] = ranksum(corr_beta_bef_stop, incorr_beta_bef_stop);
                [stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.p_val_af_stop,stats.area.(unique_brain_areas{area}).band_pass.(gettuning{ev}).beta.p_val_flag_af_stop] = ranksum(corr_beta_af_stop, incorr_beta_af_stop);
            end
            
        end
    end
end

if prs.analyse_phase
    fprintf('********** Analyses on LFP phases ********** \n');
    fprintf(['Time:  ' num2str(clock) '\n']);
    unique_brain_areas = unique({units.brain_area}); num_brain_areas = numel(unique_brain_areas);
    gettuning = prs.tuning_events;
    for area = 1:num_brain_areas
        for type = 2 %1:length(trialtypes)
            for ev = 1:length(gettuning)
                unitindx = strcmp({units.brain_area}, unique_brain_areas{area});
                ar = find(unitindx);
                % extract phase for each time point and compute phase clustering (phase locking)
                nconds = length(behv_stats.trialtype.reward);
                for cond = 1:nconds
                    clear theta theta_angle beta beta_angle plv_all
                    for ch = 1:length(ar)
                        f = lfps(1).stats.trialtype.reward(2).events.move.all_freq.freq_spectrogram;
                        % gather phase for each freq for all channels
                        if prs.analyse_phase_within_area, plv_all(ch,:,:) = lfps(ar(ch)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.plv;
                            
                            stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_angle(ch,:,:,:)...
                                = lfps(ar(ch)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.lfp_angle; end
                        
                        %% gather angles at each time point for each electrode for band passed signal
                        % theta
                        if prs.analyse_theta
                            t_temp_theta = lfps(ar(ch)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ts_lfp_align;
                            if prs.analyse_move_before_after
                                for all_cond = 1:nconds, num_trls(all_cond) = size(lfps(ar(ch)).stats.trialtype.(trialtypes{type})(all_cond).events.(gettuning{ev}).theta.lfp_align,2);end
                                ntrls = min(num_trls);
                            else
                                ntrls = size(lfps(ar(ch)).stats.trialtype.(trialtypes{type})(1).events.(gettuning{ev}).theta.lfp_align,2); % % match trial number. PLV is sensitive to trial num
                            end
                        end
                        theta_angle = angle(lfps(ar(ch)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.lfp_align);
                        % store to later compute plv across areas
                        stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.angle_mu(ch,:,:) = theta_angle;
                        %% Compute phase clustering for all trials in one timepoint for each electrode: abs(mean(exp(1i*angles_at_one_time_point_across_trials)))
                        if ~isempty(theta_angle)
                            for tmp_indx = 1:length(t_temp_theta)
                                theta(ch).itpc(tmp_indx,:) = abs ( nanmean ( exp ( 1i * theta_angle(tmp_indx,1:ntrls) ) ) ) ;
                                % compute p-value of the observed ITPC p = exp(-trl*ITPC^2)
                                theta(ch).itpc_pval = exp(-ntrls * (theta(ch).itpc(tmp_indx,:))^2);
                            end
                        end
                        
                        %% compute Rayleigh test for non-uniformiuty of circular
                        % data for each time point
                        if ~isempty(theta_angle)
                            for tmp_indx = 1:length(t_temp_theta)
                                [stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).chan(ch).theta.pval_circ(tmp_indx), stats.area.(unique_brain_areas{area}).trialtype.reward(cond).events.(gettuning{ev}).chan(ch).theta.z_val_circ(tmp_indx)]...
                                    = circ_rtest(theta_angle(tmp_indx,1:ntrls));
                            end
                        end
                        
                        % beta
                        if prs.analyse_beta
                            t_temp_beta = lfps(ar(ch)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ts_lfp_align;
                            if prs.analyse_move_before_after
                                for all_cond = 1:nconds, num_trls(all_cond) = size(lfps(ar(ch)).stats.trialtype.(trialtypes{type})(all_cond).events.(gettuning{ev}).theta.lfp_align,2);end
                                ntrls = min(num_trls);
                            else
                                ntrls = size(lfps(ar(ch)).stats.trialtype.(trialtypes{type})(1).events.(gettuning{ev}).theta.lfp_align,2); % % match trial number. PLV is sensitive to trial num
                            end
                        end
                        beta_angle = angle(lfps(ar(ch)).stats.trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.lfp_align);
                        % store to later compute plv across areas
                        stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.angle_mu(ch,:,:) = beta_angle;
                        %%%%%%%%%   f = instantaneous_frequency( amplitude (abs(x)), Fs ); % f contains "instantaneous frequency"
                        
                        %% Compute phase clustering for all trials in one timepoint for each electrode: abs(mean(exp(1i*angles_at_one_time_point_across_trials)))
                        if ~isempty(beta_angle)
                            for tmp_indx = 1:length(t_temp_beta)
                                beta(ch).itpc(tmp_indx,:) = abs ( nanmean ( exp ( 1i * beta_angle(tmp_indx,1:ntrls) ) ) ) ;
                                % compute p-value of the observed ITPC p = exp(-trl*ITPC^2)
                                beta(ch).itpc_pval = exp(-ntrls * (beta(ch).itpc(tmp_indx,:))^2);
                            end
                        end
                        % plot sanity check
                        %                     figure; hold on
                        %                     plot(t_temp_theta, theta(ch).itpc);
                        %                     plot(t_temp_beta,theta(ch).itpc); xlim([-1.5 1.5]); ylim([0 5])
                        
                        %% compute Rayleigh test for non-uniformiuty of circular
                        if ~isempty(beta_angle)
                            for tmp_indx = 1:length(t_temp_beta)
                                [stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).chan(ch).beta.pval_circ(tmp_indx), stats.area.(unique_brain_areas{area}).trialtype.reward(cond).events.(gettuning{ev}).chan(ch).beta.z_val_circ(tmp_indx)]...
                                    = circ_rtest(beta_angle(tmp_indx,1:ntrls));
                            end
                        end
                        
                    end
                    %                     % plot polar plot and histogram
                    %                     figure(1); hold on; set(gcf, 'Position',[1 704 1919 401]); title(['channel ' num2str(ch)])
                    %                     for i = 1:length(timepoints)
                    %                         subplot(4,10,i)
                    %                         polarplot([zeros(1,ntrl);squeeze(stats.area.(unique_brain_areas{area}).band.beta.reward(cond).angle(1,:,i))],[zeros(1,ntrl); ones(1,ntrl)],'k', 'LineWidth', 0.01);
                    %                         thetaticks(0:90:315)
                    %
                    %                         hold on; subplot(4,10,i+10)
                    %                         histogram(rad2deg(squeeze(stats.area.(unique_brain_areas{area}).band.beta.reward(cond).angle(1,:,i))),20)
                    %                         set(gca, 'xlim', [-180 180], 'xTick', []); box off
                    %                         xlabel([ num2str(timepoints(i)) 's'])
                    %                     end
                    
                    %% average itpc for all channels -1.5 to 1.5s store
                    if prs.analyse_phase_within_area, stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).all_freq.plv_all = squeeze(nanmean(plv_all)); end
                    
                    % for theta and beta only
                    if prs.analyse_theta | prs.analyse_beta
                        for ch = 1:length(ar)
                            if ~isempty(theta_angle)
                                stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ang_itpc(ch,:) = theta(ch).itpc;
                                stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ts = t_temp_theta;
                            end
                            if ~isempty(beta_angle)
                                stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ang_itpc(ch,:) = beta(ch).itpc;
                                stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ts = t_temp_beta;
                            end
                        end
                        % theta
                        if ~isempty(theta_angle)
                            stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ang_itpc_mu = ...
                                nanmean(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ang_itpc);
                            stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ang_itpc_sem = ...
                                nanstd(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ang_itpc)/sqrt(size(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).theta.ang_itpc,2));
                        end
                        % beta
                        if ~isempty(beta_angle)
                            stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ang_itpc_mu = ...
                                nanmean(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ang_itpc);
                            stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ang_itpc_sem = ...
                                nanstd(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ang_itpc)/sqrt(size(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).beta.ang_itpc,2));
                        end
                    end
                    
                    % plot for all channels for one area
                    %                     try
                    %                         figure; hold on; title([(unique_brain_areas{area}) ' ' (gettuning{ev}) ' cond ' num2str(cond)])
                    %                         for ch = 1:length(ar)
                    %                             plot(t_temp_beta,stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).chan(ch).beta.itpc)
                    %                             plot(t_temp_beta,(stats.area.(unique_brain_areas{area}).trialtype.(trialtypes{type})(cond).events.(gettuning{ev}).chan(ch).beta.pval_circ)<0.05)
                    %                             xlabel('time to event'); ylabel('ITPC'); xlim([-1.5 1.5])
                    %                         end
                    %                     end
                    
                    
                end
            end
        end
    end
    
    
    % compute phase locking value across areas (plv, Lachaux 1999)
    for area1 = 1:num_brain_areas
        for area2 = find(1:num_brain_areas ~= area1)
            for ev = 1:length(gettuning)
                theta_all_corr = []; theta_all_incorr = []; beta_all_corr = []; beta_all_incorr = []; theta_all_pli_corr = []; theta_all_pli_incorr = []; beta_all_pli_corr = []; beta_all_pli_incorr = [];
                for cond = 1:2 %nconds
                    if prs.analyse_move_before_after
                        clear num_trls
                        for all_cond = 1:nconds, num_trls(all_cond) = sum(behv_stats.trialtype.reward(all_cond).trlindx);end
                        ntrls = min(num_trls);
                    else
                        ntrls = sum(behv_stats.trialtype.reward(1).trlindx); % match trial number. PLV is sensitive to trial num
                    end
                    
                    if (ev == 4 && cond == 1) %||(ev == 4 && cond == 3)||(ev == 4 && cond == 4)  % no reward in cond==1
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).(['ch' num2str(ch_area1) '_to_' num2str(ch_area2)]) = NaN;
                    else
                        clear e_theta theta_plv e_beta beta_plv phase_lag_indx_theta phase_lag_indx_beta all_freq_plv phase_lag_indx_theta phase_lag_indx_theta_dir phase_lag_indx_beta phase_lag_indx_beta_dir
                        cnt = 1;
                        
                        for ch_area1 = 1:length(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).chan)   % ch in one area
                            for ch_area2 = 1:length(stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).chan)
                                %% uncomment when running plv all_freq
                                %                                 for ch_area1 = 1:length(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).all_freq.lfp_angle(:,1,1,1))  % ch in one area
                                %                                     for ch_area2 = 1:length(stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).all_freq.lfp_angle(:,1,1,1))
                                %                                         %% compute PLV for all freq
                                %                                         freq = [prs.lfp_freqmin:prs.lfp_freqmax];
                                %                                         for f_indx = 1:length(freq)
                                %                                             if f_indx ~= length(freq)
                                %                                                 e_all_freq = squeeze(exp( 1i*(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).all_freq.lfp_angle(ch_area1,f_indx,:,1:ntrls)...
                                %                                                     - stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).all_freq.lfp_angle(ch_area2,f_indx,:,1:ntrls))));
                                %                                                 all_freq_plv(cnt,f_indx,:) = abs(sum(e_all_freq,2)) / ntrls;
                                %                                             end
                                %                                         end
                                %% theta
                                if prs.analyse_theta
                                    %% phase locking value
                                    e_theta = squeeze(exp( 1i*(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).theta.angle_mu(ch_area1,:,1:ntrls)...
                                        - stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).theta.angle_mu(ch_area2,:,1:ntrls))));
                                    theta_plv(cnt,:) = abs(sum(e_theta,2)) / ntrls;
                                    
                                    %% compute phase lag-index as described in Nolte et al 2008 and Cohen, MX (2014) phase_lag_indx = abs(mean(sign(imag(phase_area_1 - phase_area_2))));
                                    diff_theta = squeeze(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).theta.angle_mu(ch_area1,:,1:ntrls)...
                                        - stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).theta.angle_mu(ch_area2,:,1:ntrls));
                                    cdd_theta = exp(1i*diff_theta);
                                    phase_lag_indx_theta(cnt,:) = abs(nanmean( abs(imag(cdd_theta)).*sign(imag(cdd_theta)),2) )./nanmean(abs(imag(cdd_theta)),2); % Weighted PLI as described in Vinck et al 2018
                                    phase_lag_indx_theta_dir(cnt,:) = nanmean( abs(imag(cdd_theta)).*sign(imag(cdd_theta)),2)./nanmean(abs(imag(cdd_theta)),2);  % same but remove abs to get directionality
                                end
                                
                                %% beta
                                if prs.analyse_beta
                                    %% phase locking value
                                    e_beta = squeeze(exp( 1i*(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).beta.angle_mu(ch_area1,:,1:ntrls)...
                                        - stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).beta.angle_mu(ch_area2,:,1:ntrls))));
                                    beta_plv(cnt,:) = abs(sum(e_beta,2)) / ntrls;
                                    
                                    %% compute phase lag-index as described in Stam et al 2007 and Cohen, MX (2014)
                                    diff_beta = squeeze(stats.area.(unique_brain_areas{area1}).trialtype.reward(cond).events.(gettuning{ev}).beta.angle_mu(ch_area1,:,1:ntrls)...
                                        - stats.area.(unique_brain_areas{area2}).trialtype.reward(cond).events.(gettuning{ev}).beta.angle_mu(ch_area2,:,1:ntrls));
                                    cdd_beta = exp(1i*diff_beta);
                                    phase_lag_indx_beta(cnt,:) = abs(nanmean( abs(imag(cdd_beta)).*sign(imag(cdd_beta)),2) )./nanmean(abs(imag(cdd_beta)),2); % Weighted PLI as described in Vinck et al 2018
                                    phase_lag_indx_beta_dir(cnt,:) = nanmean( abs(imag(cdd_beta)).*sign(imag(cdd_beta)),2)./nanmean(abs(imag(cdd_beta)),2);  % same but remove abs to get directionality
                                end
                                
                                
                                cnt=cnt+1;
                            end
                        end
                        %% average for all rows
                        %                         % all freq
                        %                         for f_indx = 1:length(all_freq_plv(1,:,1))
                        %                             stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).all_freq.PLV_mu(f_indx,:) = squeeze(nanmean(all_freq_plv(:,f_indx,:)));
                        %                             stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).all_freq.PLV_mu(f_indx,:) = squeeze(nanstd(all_freq_plv(:,f_indx,:)))./ ...
                        %                                 sqrt(length(all_freq_plv(:,1,1)));
                        %                         end
                        % theta
                        if prs.analyse_theta
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).theta.PLV_mu = nanmean(theta_plv);
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).theta.PLV_sem = nanstd(theta_plv)./sqrt(size(theta_plv,1));
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(cond).events.(gettuning{ev}).theta.PLI_mu = nanmean(phase_lag_indx_theta);
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(cond).events.(gettuning{ev}).theta.PLI_sem = nanstd(phase_lag_indx_theta)./sqrt(size(theta_plv,1));
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(cond).events.(gettuning{ev}).theta.PLI_dir_mu = nanmean(phase_lag_indx_theta_dir);
                            % save correct
                            if cond == 2, theta_all_corr = [theta_all_corr ; theta_plv ]; theta_all_pli_corr = [theta_all_pli_corr ; phase_lag_indx_theta ]; ...
                            else theta_all_incorr = [theta_all_incorr; theta_plv]; theta_all_pli_incorr = [theta_all_pli_incorr; phase_lag_indx_theta]; end
                        
                        end
                        
                        if prs.analyse_beta
                            % beta
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).beta.PLV_mu = nanmean(beta_plv);
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(cond).events.(gettuning{ev}).beta.PLV_sem = nanstd(beta_plv)./sqrt(size(beta_plv,1));
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(cond).events.(gettuning{ev}).beta.PLI_mu = nanmean(phase_lag_indx_beta);
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(cond).events.(gettuning{ev}).beta.PLI_sem = nanstd(phase_lag_indx_beta)./sqrt(size(beta_plv,1));
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(cond).events.(gettuning{ev}).beta.PLI_dir_mu = nanmean(phase_lag_indx_beta_dir);
                            % save correct
                            if cond == 2, beta_all_corr = [beta_all_corr ; beta_plv ]; beta_all_pli_corr = [beta_all_pli_corr ; phase_lag_indx_beta ]; ...
                            else beta_all_incorr = [beta_all_incorr; beta_plv]; beta_all_pli_incorr = [beta_all_pli_incorr; phase_lag_indx_beta]; end
                        end
                    end
                end
                
                clear n_boot_indx shuffled_theta shuffled_beta Z_theta Z_beta shuffled_theta_pli shuffled_beta_pli
                if prs.analyse_theta | prs.analyse_beta
                    %% create a shuffled distr to get signif
                    % find min between the two
                    if (ev == 4 && cond == 2), min_size_theta = min([size(theta_all_corr,2)]);  min_size_beta = min([size(beta_all_corr,2)]); else ...
                            min_size_theta = min([size(theta_all_corr,2) size(theta_all_incorr,2)]); min_size_beta = min([size(beta_all_corr,2) size(beta_all_incorr,2)]); end
                    % concatenate correct and incorrect
                    if (ev == 4 && cond == 2),  theta_all = theta_all_corr(:,1:min_size_theta); theta_all_pli = theta_all_pli_corr(:,1:min_size_theta);
                        beta_all = beta_all_corr(:,1:min_size_beta); beta_all_pli = beta_all_pli_corr(:,1:min_size_beta); else ...
                            theta_all = [theta_all_corr(:,1:min_size_theta) ; theta_all_incorr(:,1:min_size_theta)]; theta_all_pli = [theta_all_pli_corr(:,1:min_size_theta) ; theta_all_pli_incorr(:,1:min_size_theta)];
                        beta_all = [beta_all_corr(:,1:min_size_beta) ; beta_all_incorr(:,1:min_size_beta)]; beta_all_pli = [beta_all_pli_corr(:,1:min_size_beta) ; beta_all_pli_incorr(:,1:min_size_beta)]; end
                    
                    n_boot = 1000;
                    n_boot_indx_theta = randi(size(theta_all,1),n_boot,1);
                    n_boot_indx_beta = randi(size(beta_all,1),n_boot,1);
                    
                    for i = 1:length(n_boot_indx_theta),shuffled_theta(i,:) = randsample(theta_all(n_boot_indx_theta(i),:),length(theta_all(n_boot_indx_theta(i),:)),'true');end
                    for i = 1:length(n_boot_indx_theta),shuffled_theta_pli(i,:) = randsample(theta_all_pli(n_boot_indx_theta(i),:),length(theta_all_pli(n_boot_indx_theta(i),:)),'true');end
                    
                    for i = 1:length(n_boot_indx_beta),shuffled_beta(i,:) = randsample(beta_all(n_boot_indx_beta(i),:),length(beta_all(n_boot_indx_beta(i),:)),'true');end
                    for i = 1:length(n_boot_indx_beta),shuffled_beta_pli(i,:) = randsample(beta_all_pli(n_boot_indx_beta(i),:),length(beta_all_pli(n_boot_indx_beta(i),:)),'true');end
                    
                    % convert plv values to std dev units of the shuffled (null) distribution Z = ( PLV - mean(PLVShuff) )/ std(PLVShuff)
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.shuffled_theta = nanmean(shuffled_theta);
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.shuffled_theta_pli = nanmean(shuffled_theta_pli);
                    
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.shuffled_theta = nanmean(shuffled_beta);
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.shuffled_theta_pli = nanmean(shuffled_beta_pli);
                    
                    % theta
                    if strcmp((gettuning{ev}),'reward'), stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).theta.PLV_mu_Z = NaN; else ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).theta.PLV_mu_Z = ...
                            (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).theta.PLV_mu-nanmean(nanmean(shuffled_theta)))./...
                            nanstd(nanstd(shuffled_theta)); end
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.PLV_mu_Z = ...
                        (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.PLV_mu-nanmean(nanmean(shuffled_theta)))./...
                        nanstd(nanstd(shuffled_theta));
                    
                    if strcmp((gettuning{ev}),'reward'), stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).theta.PLI_mu_Z = NaN; else ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).theta.PLI_mu_Z = ...
                            (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).theta.PLI_mu-nanmean(nanmean(shuffled_theta_pli)))./...
                            nanstd(nanstd(shuffled_theta_pli)); end
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.PLI_mu_Z = ...
                        (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.PLI_mu-nanmean(nanmean(shuffled_theta_pli)))./...
                        nanstd(nanstd(shuffled_theta_pli));
                    
                    % beta
                    if strcmp((gettuning{ev}),'reward'), stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).beta.PLV_mu_Z = NaN; else ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).beta.PLV_mu_Z = ...
                            (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).beta.PLV_mu-nanmean(nanmean(shuffled_beta)))./...
                            nanstd(nanstd(shuffled_beta)); end
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.PLV_mu_Z = ...
                        (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.PLV_mu-nanmean(nanmean(shuffled_beta)))./...
                        nanstd(nanstd(shuffled_beta));
                    
                    if strcmp((gettuning{ev}),'reward'), stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).beta.PLI_mu_Z = NaN; else ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).beta.PLI_mu_Z = ...
                            (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).beta.PLI_mu-nanmean(nanmean(shuffled_beta_pli)))./...
                            nanstd(nanstd(shuffled_beta_pli)); end
                    stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.PLI_mu_Z = ...
                        (stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.PLI_mu-nanmean(nanmean(shuffled_beta_pli)))./...
                        nanstd(nanstd(shuffled_beta_pli));
                    
                    
                    p_val_theta = prctile(nanmean(shuffled_theta),[2.5 97.5]);
                    p_val_beta = prctile(nanmean(shuffled_beta),[2.5 97.5]);
                    
                    p_val_theta_pli = prctile(nanmean(shuffled_theta_pli),[2.5 97.5]);
                    p_val_beta_pli = prctile(nanmean(shuffled_beta_pli),[2.5 97.5]);
                    
                    % are correct and incorrect significant compared to permuted vals?
                    if strcmp((gettuning{ev}),'reward')
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.indx_signif = NaN;
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.indx_signif = NaN;
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.indx_signif = NaN;
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.indx_signif = NaN;
                    else
                        % theta plv
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.indx_signif = ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.PLV_mu(1:min_size_theta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).theta.PLV_mu(1:min_size_theta) > p_val_theta(2) | ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).theta.PLV_mu(1:min_size_theta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).theta.PLV_mu(1:min_size_theta) < p_val_theta(1);
                        % theta pli
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.indx_signif = ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.PLI_mu(1:min_size_theta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).theta.PLI_mu(1:min_size_theta) > p_val_theta_pli(2) | ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).theta.PLI_mu(1:min_size_theta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).theta.PLI_mu(1:min_size_theta) < p_val_theta_pli(1);
                        % beta plv
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.indx_signif = ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.PLV_mu(1:min_size_beta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).beta.PLV_mu(1:min_size_beta) > p_val_beta(2) | ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(2).events.(gettuning{ev}).beta.PLV_mu(1:min_size_beta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLV']).trialtype.reward(1).events.(gettuning{ev}).beta.PLV_mu(1:min_size_beta) < p_val_beta(1);
                        
                        stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.indx_signif = ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.PLI_mu(1:min_size_beta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).beta.PLI_mu(1:min_size_beta) > p_val_beta_pli(2) | ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(2).events.(gettuning{ev}).beta.PLI_mu(1:min_size_beta) - ...
                            stats.area.([unique_brain_areas{area1} '_' unique_brain_areas{area2} '_PLI']).trialtype.reward(1).events.(gettuning{ev}).beta.PLI_mu(1:min_size_beta) < p_val_beta_pli(1);
                    end
                end
            end
        end
    end
    
    %% compute phase slope index
    %     if prs.analyse_theta | prs.analyse_beta
    %         for ev = 1:length(gettuning)
    %             for ncond = 1:length(lfps(1).stats.trialtype.reward)
    %                 clear lfp_data lfp_data_ep1 lfp_data_ep2 lfp_data_ep3 lfp_data_ep4
    %                 % times to check psi: -1 to -0.5s, -0.5 to 0,0 to 0.5, 0.5 to 1
    %                 if ncond == 1 && strcmp((gettuning{ev}), 'reward')
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep1 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep1_Z = NaN;
    %                     %
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep2 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep2_Z = NaN;
    %                     %
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep3 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep3_Z = NaN;
    %                     %
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep4 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep4_Z = NaN;
    %                     %% beta
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep1 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep1_Z = NaN;
    %                     %
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep2 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep2_Z = NaN;
    %                     %
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep3 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep3_Z = NaN;
    %                     %
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep4 = NaN;
    %                     stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep4_Z = NaN;
    %                 else
    %                     for ch = 1:length(lfps)
    %                         ts = lfps(ch).stats.trialtype.reward(ncond).events.(gettuning{ev}).all_freq.ts_lfp_align;
    %                         lfp_data_ep1(ch,:,:) = lfps(ch).stats.trialtype.reward(ncond).events.(gettuning{ev}).all_freq.lfp_align(ts>-1 & ts<-0.5,:);
    %                         lfp_data_ep2(ch,:,:) = lfps(ch).stats.trialtype.reward(ncond).events.(gettuning{ev}).all_freq.lfp_align(ts>-0.5 & ts<0,:);
    %                         lfp_data_ep3(ch,:,:) = lfps(ch).stats.trialtype.reward(ncond).events.(gettuning{ev}).all_freq.lfp_align(ts>0 & ts<0.5,:);
    %                         lfp_data_ep4(ch,:,:) = lfps(ch).stats.trialtype.reward(ncond).events.(gettuning{ev}).all_freq.lfp_align(ts>0.5 & ts<1,:);
    %                     end
    %                     %% theta
    %                     if ~isempty(lfp_data_ep1)
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep1 = data2psiX(lfp_data_ep1,150,[4 12],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep1_Z = data2psiX(lfp_data_ep1,150,[4 12],1);
    %
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep1 = data2psiX(lfp_data_ep1,150,[12 20],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep1_Z = data2psiX(lfp_data_ep1,150,[12 20],1);
    %                     end
    %                     %
    %                     if ~isempty(lfp_data_ep2)
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep2 = data2psiX(lfp_data_ep2,150,[4 12],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep2_Z = data2psiX(lfp_data_ep2,150,[4 12],1);
    %
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep2 = data2psiX(lfp_data_ep2,150,[12 20],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep2_Z = data2psiX(lfp_data_ep2,150,[12 20],1);
    %                     end
    %                     %
    %                     if ~isempty(lfp_data_ep3)
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep3 = data2psiX(lfp_data_ep3,150,[4 12],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep3_Z = data2psiX(lfp_data_ep3,150,[4 12],1);
    %
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep3 = data2psiX(lfp_data_ep3,150,[12 20],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep3_Z = data2psiX(lfp_data_ep3,150,[12 20],1);
    %                     end
    %                     %
    %                     if ~isempty(lfp_data_ep4)
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep4 = data2psiX(lfp_data_ep4,150,[4 12],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).theta.psi_theta_ep4_Z = data2psiX(lfp_data_ep4,150,[4 12],1);
    %
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep4 = data2psiX(lfp_data_ep4,150,[12 20],0);
    %                         stats.trialtype.reward(ncond).events.(gettuning{ev}).beta.psi_beta_ep4_Z = data2psiX(lfp_data_ep4,150,[12 20],1);
    %                     end
    %                 end
    %             end
    %         end
    %     end
end



fprintf('**********End of LFP Pop Analyses********** \n');
fprintf(['Time:  ' num2str(clock) '\n']);

