function behv2 = CombineBehv(behv)

nsessions = length(behv);

%% combine trials
trials = [];
for i=1:nsessions
    trials = [trials behv(i).trials];
end
behv2.trials = trials;

%% combine stats
trlindx.correct = []; trlindx.incorrect = []; trlindx.crazy = [];
pos_final.r_monk = []; pos_final.theta_monk = []; pos_final.r_targ = []; pos_final.theta_targ = [];
pos_abs.x_monk = []; pos_abs.y_monk = []; 
pos_abs.z_leye = []; pos_abs.y_leye = []; pos_abs.z_reye = []; pos_abs.y_reye = []; 
pos_rel.x_targ = []; pos_rel.y_targ = []; pos_rel.r_targ = []; pos_rel.theta_targ = []; 
pos_rel.x_leye = []; pos_rel.y_leye = []; pos_rel.x_reye = []; pos_rel.y_reye = []; 
rewardwin = []; pCorrect = []; pcorrect_shuffled_mu = [];

for i=1:nsessions
    trlindx.correct = [trlindx.correct behv(i).stats.trialtype.reward(2).trlindx];
    trlindx.incorrect = [trlindx.incorrect behv(i).stats.trialtype.reward(1).trlindx];
    trlindx.crazy = ~(trlindx.correct | trlindx.incorrect); ntrls = sum(~trlindx.crazy);
    
    pos_final.r_monk = [pos_final.r_monk behv(i).stats.pos_final.r_monk];
    pos_final.theta_monk = [pos_final.theta_monk behv(i).stats.pos_final.theta_monk];
    pos_final.r_targ = [pos_final.r_targ behv(i).stats.pos_final.r_targ];
    pos_final.theta_targ = [pos_final.theta_targ behv(i).stats.pos_final.theta_targ];
    
    pos_abs.x_monk = [pos_abs.x_monk behv(i).stats.pos_abs.x_monk];
    pos_abs.y_monk = [pos_abs.y_monk behv(i).stats.pos_abs.y_monk];
%     pos_abs.z_leye = [pos_abs.z_leye behv(i).stats.pos_abs.z_leye];
%     pos_abs.y_leye = [pos_abs.y_leye behv(i).stats.pos_abs.y_leye];
%     pos_abs.z_reye = [pos_abs.z_reye behv(i).stats.pos_abs.z_reye];
%     pos_abs.y_reye = [pos_abs.y_reye behv(i).stats.pos_abs.y_reye];
    
    pos_rel.x_targ = [pos_rel.x_targ behv(i).stats.pos_rel.x_targ];
    pos_rel.y_targ = [pos_rel.y_targ behv(i).stats.pos_rel.y_targ];
    pos_rel.r_targ = [pos_rel.r_targ behv(i).stats.pos_rel.r_targ];
    pos_rel.theta_targ = [pos_rel.theta_targ behv(i).stats.pos_rel.theta_targ];
    pos_rel.x_leye = [pos_rel.x_leye behv(i).stats.pos_rel.x_leye];
    pos_rel.y_leye = [pos_rel.y_leye behv(i).stats.pos_rel.y_leye];
    pos_rel.x_reye = [pos_rel.x_reye behv(i).stats.pos_rel.x_reye];
    pos_rel.y_reye = [pos_rel.y_reye behv(i).stats.pos_rel.y_reye];
    
    rewardwin = [rewardwin ; behv(i).stats.trialtype.all.accuracy.rewardwin];
    pCorrect = [pCorrect ; behv(i).stats.trialtype.all.accuracy.pCorrect];
    pcorrect_shuffled_mu = [pcorrect_shuffled_mu ; behv(i).stats.trialtype.all.accuracy.pcorrect_shuffled_mu];
end
behv2.stats.trlindx = trlindx;
behv2.stats.pos_final = pos_final;
behv2.stats.pos_abs = pos_abs;
behv2.stats.pos_rel = pos_rel;
behv2.stats.accuracy.rewardwin = mean(rewardwin);
behv2.stats.accuracy.pCorrect = mean(pCorrect);
behv2.stats.accuracy.pCorrect_shuffled_mu = mean(pcorrect_shuffled_mu);