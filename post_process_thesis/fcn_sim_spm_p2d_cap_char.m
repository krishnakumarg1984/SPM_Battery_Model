% function fcn_sim_spm_p2d_cap_char()

%% Simulate cnst current SPM followed by p2d
t_start = cputime;
run('user_inputs_for_sim_cap_char');
run('disc_spm_core_sim_cnst_curr');   % spm simulation loop is in this script
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script
t_total = cputime-t_start;
return;