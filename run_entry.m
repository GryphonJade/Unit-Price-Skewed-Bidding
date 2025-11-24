%% Counterfactual Simlation
clear; clc;


addpath(genpath('config'));
addpath(genpath('utils'));
addpath(genpath('strategy'));
addpath(genpath('sim'));
addpath(genpath('stats'));
addpath(genpath('counterfactual'));


fprintf('Counterfactual Simulation...\n');
cf_run_fixed_entry(10);  % n_target = 5

fprintf('Completed. Results saved in output/ directory\n');
