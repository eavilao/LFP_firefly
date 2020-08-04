% remove daw data lfp

for i = 1:length(experiments.sessions)
    for j = 1:length(experiments.sessions(i).lfps)
    experiments.sessions(i).lfps(1).trials = []; 
    experiments.sessions(i).lfps(1).stationary = [];
    experiments.sessions(i).lfps(1).mobile = [];
    experiments.sessions(i).lfps(1).eyesfixed = [];
    experiments.sessions(i).lfps(1).eyesfree = []; 
    end 
end 
    
save('experiments_lfp_Quig_182_Schro_113_114_noTrials','experiments', '-v7.3')