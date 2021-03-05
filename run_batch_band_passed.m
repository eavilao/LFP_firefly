rerun_analyseLFP = true; 
rerun_AddPop = true; 

%%    
load('experiments_lfp_Schro_113_band_passed')
prs = default_prs(53,113);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Schro_113_band_passed', 'experiments', '-v7.3')
clear
%%
load('experiments_lfp_Schro_107_band_passed')
prs = default_prs(53,107);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
save('experiments_lfp_Schro_107_band_passed', 'experiments', '-v7.3')
clear experiments
%%
% load('experiments_lfp_Schro_86_band_passed')
% prs = default_prs(53,86);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Schro_86_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Quigley_207_band_passed')
% prs = default_prs(44,207);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Quigley_207_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Quigley_188_band_passed')
% prs = default_prs(44,188);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Quigley_188_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Quigley_185_band_passed')
% prs = default_prs(44,185);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Quigley_185_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Bruno_43_band_passed')
% prs = default_prs(51,43);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_43_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Bruno_42_band_passed')
% prs = default_prs(51,42);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_42_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Bruno_41_band_passed')
% prs = default_prs(51,41);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_41_band_passed', 'experiments', '-v7.3')
% clear experiments
% %%
% load('experiments_lfp_Bruno_38_band_passed')
% prs = default_prs(51,38);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% save('experiments_lfp_Bruno_38_band_passed', 'experiments', '-v7.3')
% clear experiments

