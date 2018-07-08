%% Run polynomial model
run('poly_ce_model_loop.m');
run('postprocess_pp2d_results');
extradata_for_solver = rmfield(extradata_for_solver,{'pp2d_model_params'});

%% Setup the ML model, run, post-process and plot
run('setup_ML_model_parameters.m');
run('ML_ce_model_loop.m');
run('postprocess_ml_results_essential.m');
run('ce_movie_creation.m');
run('investigating_Qe');
% clear initial_soc_pct;

