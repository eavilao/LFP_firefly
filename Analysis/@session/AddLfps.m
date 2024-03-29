%% add lfps
function AddLfps(this,prs)
cd(prs.filepath_neur);
% determine type of electrode
linearprobe_type = []; utaharray_type = [];
for k=1:length(prs.electrode_type)
    linearprobe_type = [linearprobe_type find(cellfun(@(electrode_type) strcmp(prs.electrode_type{k},electrode_type), prs.linearprobe.types),1)];
    utaharray_type = [utaharray_type find(cellfun(@(electrode_type) strcmp(prs.electrode_type{k},electrode_type), prs.utaharray.types),1)];
end

if ~isempty(linearprobe_type) % assume linearprobe is recorded using Plexon
    brain_area = prs.area{strcmp(prs.electrode_type,prs.linearprobe.types{linearprobe_type})};
    file_ead=dir('*_ead.plx'); file_lfp=dir('*_lfp.plx'); prs.neur_filetype = 'plx';
    % read events
    fprintf(['... reading events from ' file_ead.name '\n']);
    [events_plx, prs.fs_spk] = GetEvents_plx(file_ead.name);
    % read lfp
    if length(this.behaviours.trials)==length(events_plx.t_end)
        fprintf(['... reading ' file_lfp.name '\n']);
        [ch_id,electrode_id] = MapChannel2Electrode(prs.linearprobe.types{linearprobe_type});
        for j=1:prs.linearprobe.channelcount(linearprobe_type)
            fprintf(['...... channel ' num2str(j) '/' num2str(prs.linearprobe.channelcount(linearprobe_type)) '\n']);
            [adfreq, n, ~, fn, ad] = plx_ad_v(file_lfp.name, j-1);
            if n == fn
                if adfreq > prs.fs_lfp, N = round(adfreq/prs.fs_lfp); ad = downsample(ad,N); end
                channel_id = j;
                fprintf(['Segmenting LFP :: channel ' num2str(channel_id) '\n']);
                fprintf(['Time:  ' num2str(clock) '\n']);
                this.lfps(end+1) = lfp(channel_id,electrode_id(ch_id == channel_id),prs.linearprobe.types{linearprobe_type});
                this.lfps(end).brain_area = brain_area;
                this.lfps(end).AddTrials(ad,adfreq/N,events_plx,this.behaviours,prs);
            else
                fprintf('...... LFP is fragmented. Use a machine with more RAM or contact KL\n');
            end
        end
    else
        fprintf('Cannot segment LFP: Trial counts in smr and plx files do not match \n');
        fprintf(['Trial end events: PLX file - ' num2str(length(events_plx.t_end)) ...
            ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
        fprintf('Debug and try again! \n');
    end
else
    fprintf('No Plexon neural data files in the specified path \n');
end

if prs.isRipple
    if ~isempty(utaharray_type)
        file_nev=dir('*.nev'); file_ns2=dir('*.ns2'); prs.neur_filetype = 'nev';
        fprintf(['... reading events from ' file_nev.name '\n']);
        fprintf(['Time:  ' num2str(clock)])
        [events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK
        if length(this.behaviours.trials)~=length(events_nev.t_end)
            events_nev = FixEvents_nev(events_nev,this.behaviours.trials);
        end
        if length(this.behaviours.trials)==length(events_nev.t_end)
            NS2 = openNSx(['/' file_ns2.name],'report','read', 'uV');
            if NS2.MetaTags.ChannelCount ~= prs.utaharray.channelcount(utaharray_type), warning('Unexpected channel count in the file \n'); end
            [ch_id,electrode_id] = MapChannel2Electrode(prs.utaharray.types{utaharray_type});
            brain_area = prs.area{strcmp(prs.electrode_type,prs.utaharray.types{utaharray_type})};
            for j=[1:64 129:192] %1:prs.utaharray.channelcount(utaharray_type) 
                channel_id = NS2.MetaTags.ChannelID(j);
                fprintf(['Segmenting LFP :: channel ' num2str(channel_id) '\n']);
                fprintf(['Time:  ' num2str(clock) '\n']);
                this.lfps(end+1) = lfp(channel_id,electrode_id(ch_id == channel_id),prs.utaharray.types{utaharray_type});
                this.lfps(end).brain_area = prs.MapDualArray2BrainArea_128(brain_area, this.lfps(end).electrode_id);
                this.lfps(end).AddTrials(NS2.Data(j,:),NS2.MetaTags.SamplingFreq,events_nev,this.behaviours,prs);
            end
        else
            fprintf('Cannot segment LFP: Trial counts in smr and nev files do not match \n');
            fprintf(['Trial end events: NEV file - ' num2str(length(events_nev.t_end)) ...
                ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
            fprintf('Debug and try again! \n');
        end
    else
        fprintf('No Ripple neural data files in the specified path \n');
    end
    
elseif utaharray_type~=3
    if ~isempty(utaharray_type) % assume utaharray is recorded using Cereplex
        file_nev=dir('*.nev'); file_ns1=dir('*.ns1'); prs.neur_filetype = 'nev';
        fprintf(['... reading events from ' file_nev.name '\n']);
        fprintf(['Time:  ' num2str(clock)])
        [events_nev,prs] = GetEvents_nev(file_nev.name,prs); % requires package from Blackrock Microsystems: https://github.com/BlackrockMicrosystems/NPMK
        if length(this.behaviours.trials)~=length(events_nev.t_end)
            events_nev = FixEvents_nev(events_nev,this.behaviours.trials);
        end
        if length(this.behaviours.trials)==length(events_nev.t_end)
            NS1 = openNSx(['/' file_ns1.name],'report','read', 'uV');
            if NS1.MetaTags.ChannelCount ~= prs.utaharray.channelcount(utaharray_type), warning('Unexpected channel count in the file \n'); end
            [ch_id,electrode_id] = MapChannel2Electrode(prs.utaharray.types{utaharray_type});
            brain_area = prs.area{strcmp(prs.electrode_type,prs.utaharray.types{utaharray_type})};
            for j=1:prs.utaharray.channelcount(utaharray_type) % [2 41 54 96]
                channel_id = NS1.MetaTags.ChannelID(j);
                fprintf(['Segmenting LFP :: channel ' num2str(channel_id) '\n']);
                fprintf(['Time:  ' num2str(clock) '\n']);
                this.lfps(end+1) = lfp(channel_id,electrode_id(ch_id == channel_id),prs.utaharray.types{utaharray_type});
                if strcmp(prs.utaharray.types{utaharray_type},'utah96'), this.lfps(end).brain_area = brain_area;
                else, this.lfps(end).brain_area = prs.MapDualArray2BrainArea(brain_area, this.lfps(end).electrode_id); end
                this.lfps(end).AddTrials(NS1.Data(j,:),NS1.MetaTags.SamplingFreq,events_nev,this.behaviours,prs);
            end
        else
            fprintf('Cannot segment LFP: Trial counts in smr and nev files do not match \n');
            fprintf(['Trial end events: NEV file - ' num2str(length(events_nev.t_end)) ...
                ' , SMR file - ' num2str(length(this.behaviours.trials)) '\n']);
            fprintf('Debug and try again! \n');
        end
    else
        fprintf('No Cereplex neural data files in the specified path \n');
    end
    
    
end

end


