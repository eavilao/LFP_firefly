%Run batch of files to analyze LFP

name_exp = 'spectro_move_ba';
extract_behv_only = false;


if extract_behv_only
%     experiments = experiment('firefly-monkey');
%     experiments.AddSessions(44,185,{'behv'})
%     experiments.AddSessions(44,188,{'behv'})
%     experiments.AddSessions(44,207,{'behv'})
%     cd('/Volumes/Extreme SSD/temp_data/Behavior')
%     disp('   Saving ....')
%     save('experiments_lfp_Quigley_185_188_207_only_behv', 'experiments', '-v7.3')
%     clear experiments
    
    experiments = experiment('firefly-monkey');
%     experiments.AddSessions(51,41,{'behv'})
    experiments.AddSessions(51,42,{'behv'})
    experiments.AddSessions(51,43,{'behv'})
    cd('/Volumes/Extreme SSD/temp_data/Behavior')
    disp('   Saving ....')
    save('experiments_lfp_Bruno_38_41_42_43', 'experiments', '-v7.3')
    clear experiments
    
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(53,86,{'behv'})
    experiments.AddSessions(53,107,{'behv'})
    experiments.AddSessions(53,113,{'behv'})
    cd('/Volumes/Extreme SSD/temp_data/Behavior')
    disp('   Saving ....')
    save('experiments_lfp_Schro_86_107_113', 'experiments', '-v7.3')
    clear experiments
    
else
    
    %% Quigley
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,182,{'behv','lfps','pop'})
    % save(['experiments_lfp_Quigley_182_' name_exp], 'experiments', '-v7.3')
    % clear experiments
    %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,183,{'behv','lfps','pop'})
    % save(['experiments_lfp_Quigley_183_' name_exp], 'experiments', '-v7.3')
    % clear experiments
    %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,184,{'behv','lfps','pop'})
    % save(['experiments_lfp_Quigley_184_' name_exp], 'experiments', '-v7.3')
    % clear experiments
    % %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(44,185,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Quigley_185_' name_exp], 'experiments', '-v7.3')
    clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,186,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_186_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,187,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_187_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(44,188,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Quigley_188_' name_exp], 'experiments', '-v7.3')
    clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,189,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_189_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,190,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_190_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,191,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_191_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,192,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_192_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,193,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_193_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,194,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_194_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,195,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_196_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,196,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_196_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,197,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_197_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,198,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_198_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,199,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_199_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,200,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_201_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,201,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_201_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,202,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_202_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,203,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_203_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,204,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_204_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,205,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_205_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,206,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_206_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % % % %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(44,207,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Quigley_207_' name_exp], 'experiments', '-v7.3')
    clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,208,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_208_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,209,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_209_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,210,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_210_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,211,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_211_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,212,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_212_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,213,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_213_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,214,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_214_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,215,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_215_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,216,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_216_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,217,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_217_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,218,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_218_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,219,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_219_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,220,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_220_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,221,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_221_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(44,222,{'behv','lfps','pop'})
    % save('experiments_lfp_Quigley_222_' name_exp, 'experiments', '-v7.3')
    % clear experiments
    % %
    
    %% Bruno
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(51,38,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Bruno_38_' name_exp, 'experiments'], '-v7.3')
    clear experiments
    %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(51,41,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Bruno_41_' name_exp], 'experiments', '-v7.3')
    clear experiments
    
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(51,42,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Bruno_42_' name_exp], 'experiments', '-v7.3')
    clear experiments
    %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(51,43,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Bruno_43_' name_exp], 'experiments', '-v7.3')
    clear experiments
    
    
    %% Schro
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(53,86,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Schro_86_' name_exp], 'experiments', '-v7.3')
    clear experiments
    %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(53,93,{'behv','lfps','pop'})
    % save(['experiments_lfp_Schro_93_' name_exp], 'experiments', '-v7.3')
    % clear experiments
    % %
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(53,95,{'behv','lfps','pop'})
    % save(['experiments_lfp_Schro_95_' name_exp], 'experiments', '-v7.3')
    % clear experiments
    %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(53,107,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Schro_107_' name_exp], 'experiments', '-v7.3')
    clear experiments
    
    % experiments = experiment('firefly-monkey')
    % experiments.AddSessions(53,108,{'behv','lfps','pop'})
    % save(['experiments_lfp_Schro_108_' name_exp], 'experiments', '-v7.3')
    % clear experiments
    %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(53,113,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Schro_113_' name_exp], 'experiments', '-v7.3')
    clear experiments
    %
    
    %% Viktor
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(71,1,{'behv','lfps','pop'})
    cd('E:\Output\band_passed')
    save(['experiments_lfp_Vik_1_' name_exp], 'experiments', '-v7.3')
    clear experiments
    %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(71,2,{'behv','lfps','pop'})
    cd('E:\Output\spectrograms\move_before_after')
    save(['experiments_lfp_Vik_2_' name_exp], 'experiments', '-v7.3')
    clear experiments
    %
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(71,4,{'behv','lfps','pop'})
    cd('E:\Output\band_passed')
    save(['experiments_lfp_Vik_4_' name_exp], 'experiments', '-v7.3')
    clear experiments
end