% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
addpath('../spm_results');
load('cts_sim_Northrop_cnst_dischg_initial_soc_100pct_May_26_2018_13_39_12');

run('plots_spm.m');
clear;

load('disc_sim_Northrop_cnst_dischg_initial_soc_100pct_May_26_2018_13_39_24');
figure(1); subplot(211); hold on; subplot(212); hold on;
figure(2); subplot(211); hold on; subplot(212); hold on;

run('plots_spm.m');

figure(1);
shg;