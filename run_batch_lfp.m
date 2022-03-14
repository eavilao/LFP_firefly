name_load_exp = 'spectro_move_ba';
name_save_exp = 'phase_move_ba';
% %%
% rerun_analyseLFP = true; 
% rerun_AddPop = true; 
% disp('Loading file 1... . . ')
% disp(['experiments_lfp_Schro_113_' name_load_exp])
% load(['experiments_lfp_Schro_113_' name_load_exp])
% prs = default_prs(53,113);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Schro_113_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;
disp('Loading file 2... . . ')
load(['experiments_lfp_Schro_107_' name_load_exp])
prs = default_prs(53,107);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Schro_107_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 3... . . ')
load(['experiments_lfp_Schro_86_' name_load_exp])
prs = default_prs(53,86);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Schro_86_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 4... . . ')
load(['experiments_lfp_Quigley_207_' name_load_exp])
prs = default_prs(44,207);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Quigley_207_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 5... . . ')
load(['experiments_lfp_Quigley_188_' name_load_exp])
prs = default_prs(44,188);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Quigley_188_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 6... . . ')
load(['experiments_lfp_Quigley_185_' name_load_exp])
prs = default_prs(44,185);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Quigley_185_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 7... . . ')
load(['experiments_lfp_Bruno_43_' name_load_exp])
prs = default_prs(51,43);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Bruno_43_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;
disp('Loading file 8... . . ')
load(['experiments_lfp_Bruno_42_' name_load_exp])
prs = default_prs(51,42);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Bruno_42_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 9... . . ')
load(['experiments_lfp_Bruno_41_' name_load_exp])
prs = default_prs(51,41);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Bruno_41_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 10... . . ')
load(['experiments_lfp_Bruno_38_' name_load_exp])
prs = default_prs(51,38);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Bruno_38_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 11... . . ')
load(['experiments_lfp_Vik_1_' name_load_exp])
prs = default_prs(71,1);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Viktor_1_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 11... . . ')
load(['experiments_lfp_Vik_2_' name_load_exp])
prs = default_prs(71,1);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Viktor_2_' name_save_exp], 'experiments', '-v7.3')
clear experiments
%%
% rerun_analyseLFP = true; 
% rerun_AddPop = true;  
% disp('Loading file 12... . . ')
% load(['experiments_lfp_Vik_3_'name_load_exp])
% prs = default_prs(71,3);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% disp('Saving... . .' )
% save(['experiments_lfp_Viktor_3_' name_save_exp], 'experiments', '-v7.3')
% clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 13... . . ')
load(['experiments_lfp_Vik_4_' name_load_exp])
prs = default_prs(71,4);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save(['experiments_lfp_Viktor_4_' name_save_exp], 'experiments', '-v7.3')
clear experiments
