
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 1... . . ')
load('experiments_lfp_Schro_113_phase')
prs = default_prs(53,113);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Schro_113_phase', 'experiments', '-v7.3')
clear
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;
disp('Loading file 2... . . ')
load('experiments_lfp_Schro_107_phase')
prs = default_prs(53,107);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Schro_107_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 3... . . ')
load('experiments_lfp_Schro_86_phase')
prs = default_prs(53,86);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Schro_86_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 4... . . ')
load('experiments_lfp_Quigley_207_phase')
prs = default_prs(44,207);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Quigley_207_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 5... . . ')
load('experiments_lfp_Quigley_188_phase')
prs = default_prs(44,188);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Quigley_188_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 6... . . ')
load('experiments_lfp_Quigley_185_phase')
prs = default_prs(44,185);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Quigley_185_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 7... . . ')
load('experiments_lfp_Bruno_43_phase')
prs = default_prs(51,43);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_43_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;
disp('Loading file 8... . . ')
load('experiments_lfp_Bruno_42_phase')
prs = default_prs(51,42);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_42_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 9... . . ')
load('experiments_lfp_Bruno_41_phase')
prs = default_prs(51,41);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_41_phase', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 10... . . ')
load('experiments_lfp_Bruno_38_phase')
prs = default_prs(51,38);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_38_phase', 'experiments', '-v7.3')
clear experiments

