function trials = AddTrials2Behaviour(prs)

trials = []; % initialise

%% list all files to read
flist_log=dir('*.log');
%% remove files starting with ._ in macosx
% remove all folders
isBadFile = cat(1,flist_log.isdir); %# all directories are bad
% loop to identify hidden files
for iFile = find(~isBadFile)' %'# loop only non-dirs
    % on OSX, hidden files start with a dot
    isBadFile(iFile) = strcmp(flist_log(iFile).name(1),'.');
    if ~isBadFile(iFile) && ispc
        % check for hidden Windows files - only works on Windows
        [~,stats] = fileattrib(fullfile(folder,flist_log(iFile).name));
        if stats.hidden
            isBadFile(iFile) = true;
        end
    end
end
% remove bad files
flist_log(isBadFile) = [];
%%
for i=1:length(flist_log), fnum_log(i) = str2num(flist_log(i).name(end-6:end-4)); end
flist_smr=dir('*.smr');
%%
%% remove files starting with ._ in macosx
% remove all folders
isBadFile = cat(1,flist_smr.isdir); %# all directories are bad
% loop to identify hidden files
for iFile = find(~isBadFile)' %'# loop only non-dirs
    % on OSX, hidden files start with a dot
    isBadFile(iFile) = strcmp(flist_smr(iFile).name(1),'.');
    if ~isBadFile(iFile) && ispc
        % check for hidden Windows files - only works on Windows
        [~,stats] = fileattrib(fullfile(folder,flist_smr(iFile).name));
        if stats.hidden
            isBadFile(iFile) = true;
        end
    end
end
% remove bad files
flist_smr(isBadFile) = [];
%%
for i=1:length(flist_smr), fnum_smr(i) = str2num(flist_smr(i).name(end-6:end-4)); end
% flist_mat=dir('*.mat');
% for i=1:length(flist_mat), fnum_mat(i) = str2num(flist_mat(i).name(end-6:end-4)); end
nfiles = length(flist_log);
%% read files
for i=1:nfiles
    fprintf(['... reading ' flist_log(i).name '\n']);
    % read .log file
    if prs.monk_Ody
         trials_log = AddLOGData_Ody(flist_log(i).name,prs);
    elseif prs.isRipple
        trials_log = AddLOGData_NYU(flist_log(i).name,prs);
    else    
        trials_log = AddLOGData(flist_log(i).name);
    end
    % read all .smr files associated with this log file
    if i<nfiles, indx_smr = find(fnum_smr >= fnum_log(i) & fnum_smr < fnum_log(i+1));
    else indx_smr = find(fnum_smr >= fnum_log(i)); end
    trials_smr = [];
    for j = indx_smr
        data_smr = ImportSMR(flist_smr(j).name);
        if prs.isRipple
            trials_smr = [trials_smr AddSMRData_NYU(data_smr,prs)];
        else
            trials_smr = [trials_smr AddSMRData(data_smr,prs)];
        end
        
    end
    % merge contents of .log and .smr files
    ntrls_log = length(trials_log); ntrls_smr = length(trials_smr);
    if ntrls_smr <= ntrls_log
        %       for j=1:length(trials_smr), trials_temp(j) = catstruct(trials_smr(j),trials_log(j)) ; end
        for j=1:length(trials_smr), trials_temp(j) = concat_smr_log(trials_smr(j),trials_log(j)) ; end
    else  % apply a very dirty fix if spike2 was not "stopped" on time (can happen when replaying stimulus movie)
        %         for j=1:ntrls_log, trials_temp(j) = catstruct(trials_smr(j),trials_log(j)) ; end
        for j=1:ntrls_log, trials_temp(j) = concat_smr_log(trials_smr(j),trials_log(j)) ; end
        dummy_trials_log = trials_log(1:ntrls_smr-ntrls_log);
        %         for j=1:(ntrls_smr-ntrls_log); trials_temp(ntrls_log+j) = catstruct(trials_smr(ntrls_log+j),dummy_trials_log(j)); end
        for j=1:(ntrls_smr-ntrls_log); trials_temp(ntrls_log+j) = concat_smr_log(trials_smr(ntrls_log+j),dummy_trials_log(j)); end
    end
    % add contents of .mat file
    %     trials_temp = AddMATData(flist_mat(i).name,trials_temp);
    trials = [trials trials_temp];
    clear trials_temp;
    fprintf(['... total trials = ' num2str(length(trials)) '\n']);
end