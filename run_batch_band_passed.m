
%%
% rerun_analyseLFP = false; 
% rerun_AddPop = true; 
% load('experiments_lfp_Schro_113_spectro_all_freq_cohero')
% prs = default_prs(53,113);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Schro_113_cohero', 'experiments', '-v7.3')
clear
%%
rerun_analyseLFP = false; 
rerun_AddPop = true; 
load('experiments_lfp_Schro_107_spectro_all_freq_cohero')
prs = default_prs(53,107);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Schro_107_cohero', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = false; 
rerun_AddPop = true;  
load('experiments_lfp_Schro_86_spectro_all_freq_cohero')
prs = default_prs(53,86);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Schro_86_cohero', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = false; 
rerun_AddPop = true;  
load('experiments_lfp_Quigley_207_spectro_all_freq_cohero')
prs = default_prs(44,207);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Quigley_207_cohero', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = false; 
rerun_AddPop = true; 
load('experiments_lfp_Quigley_188_spectro_all_freq_cohero')
prs = default_prs(44,188);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Quigley_188_cohero', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = false; 
rerun_AddPop = true; 
load('experiments_lfp_Quigley_185_spectro_all_freq_cohero')
prs = default_prs(44,185);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Quigley_185_cohero', 'experiments', '-v7.3')
clear experiments
% %%
% rerun_analyseLFP = true; 
% rerun_AddPop = true;  
% load('experiments_lfp_Bruno_43_phase')
% prs = default_prs(51,43);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_43_phase', 'experiments', '-v7.3')
% clear experiments
% %%
% rerun_analyseLFP = true; 
% rerun_AddPop = true;  
% load('experiments_lfp_Bruno_42_phase')
% prs = default_prs(51,42);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_42_phase', 'experiments', '-v7.3')
% clear experiments
% %%
% rerun_analyseLFP = true; 
% rerun_AddPop = true;  
% load('experiments_lfp_Bruno_41_phase')
% prs = default_prs(51,41);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_41_phase', 'experiments', '-v7.3')
% clear experiments
% %%
% rerun_analyseLFP = true; 
% rerun_AddPop = true;  
% load('experiments_lfp_Bruno_38_phase')
% prs = default_prs(51,38);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_38_phase', 'experiments', '-v7.3')
% clear experiments
% 
