
% %%
% rerun_analyseLFP = true; 
% rerun_AddPop = true; 
% disp('Loading file 1... . . ')
% load('experiments_lfp_Schro_113_spectro_move_ba')
% prs = default_prs(53,113);
% if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
% if rerun_AddPop
% experiments.sessions(1).populations(1)=[];
% experiments.sessions(1).AddPopulation('lfps',prs)  
% end
% disp('Saving... . .' )
% save('experiments_lfp_Schro_113_phase_mba', 'experiments', '-v7.3')
% clear
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;
disp('Loading file 2... . . ')
load('experiments_lfp_Schro_107_spectro_move_ba')
prs = default_prs(53,107);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Schro_107_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 3... . . ')
load('experiments_lfp_Schro_86_spectro_move_ba')
prs = default_prs(53,86);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Schro_86_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 4... . . ')
load('experiments_lfp_Quigley_207_spectro_move_ba')
prs = default_prs(44,207);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Quigley_207_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 5... . . ')
load('experiments_lfp_Quigley_188_spectro_move_ba')
prs = default_prs(44,188);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Quigley_188_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true; 
disp('Loading file 6... . . ')
load('experiments_lfp_Quigley_185_spectro_move_ba')
prs = default_prs(44,185);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Quigley_185_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 7... . . ')
load('experiments_lfp_Bruno_43_spectro_move_ba')
prs = default_prs(51,43);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_43_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;
disp('Loading file 8... . . ')
load('experiments_lfp_Bruno_42_spectro_move_ba')
prs = default_prs(51,42);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_42_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 9... . . ')
load('experiments_lfp_Bruno_41_spectro_move_ba')
prs = default_prs(51,41);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_41_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 10... . . ')
load('experiments_lfp_Bruno_38_spectro_move_ba')
prs = default_prs(51,38);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Bruno_38_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 11... . . ')
load('experiments_lfp_Vik_1_spectro_move_ba')
prs = default_prs(71,1);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Viktor_1_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 11... . . ')
load('experiments_lfp_Vik_2_band_passed')
prs = default_prs(71,1);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Viktor_2_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 12... . . ')
load('experiments_lfp_Vik_3_spectro_move_ba')
prs = default_prs(71,3);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Viktor_3_phase_mba', 'experiments', '-v7.3')
clear experiments
%%
rerun_analyseLFP = true; 
rerun_AddPop = true;  
disp('Loading file 13... . . ')
load('experiments_lfp_Vik_4_spectro_move_ba')
prs = default_prs(71,4);
if rerun_analyseLFP, experiments.sessions.AnalyseLfps(prs); end
if rerun_AddPop
experiments.sessions(1).populations(1)=[];
experiments.sessions(1).AddPopulation('lfps',prs)  
end
disp('Saving... . .' )
save('experiments_lfp_Viktor_4_phase_mba', 'experiments', '-v7.3')
clear experiments
