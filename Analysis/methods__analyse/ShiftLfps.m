function [shiftedtrials_lfp, shiftedtrials_ts] = ShiftLfps(trials,continuous,eventtimes,fieldname)
% shifts spike trains by eventtimes - used to align spike trains to events

%% check if there as as many event times as there are trials
if length(trials)~=length(eventtimes)
    fprintf('error: event times should be a vector of same length as trials \n');
    return;
end

%% shift spike train on each trial by the event time on that trial
ntrls = length(trials);
t_min = min(cell2mat({continuous.ts}'));
t_max = max(cell2mat({continuous.ts}'));
dt = median(diff(continuous(1).ts));
nt = ceil((t_max - t_min)/dt);
ts = linspace(4*t_min, 4*t_max, 4*nt);
shiftedtrials(ntrls) = struct();
for i=1:length(trials)
    indx = find(ts > continuous(i).ts(1),1)-1;
    shiftlen = round(eventtimes(i)/dt);
    shiftedtrials(i).(fieldname) = nan(4*nt,1);
    try 
        if ~isnan(shiftlen), shiftedtrials(i).(fieldname)(indx-shiftlen : indx+length(trials(i).(fieldname))-1 - shiftlen) = trials(i).(fieldname); end 
    catch
    end
end
shiftedtrials_lfp = cell2mat({shiftedtrials.(fieldname)});
shiftedtrials_ts = ts;