function [trial_target_numbers, trial_type, prescribed_PT, trial_ret_period, trial_ITI, trial_stim_wait] = generate_trial_table_E1retention_fMRI_v1(subject_code)
% subject_code = 'S025';
%% setup block structure
trial_count = 24 + 27*6;
trial_target_numbers = nan(trial_count, 1);
trial_type = nan(trial_count, 1);
prescribed_PT = nan(trial_count, 1);
trial_ret_period = nan(trial_count, 1);
trial_ITI = nan(trial_count, 1);
trial_stim_wait = nan(trial_count, 1);

%% Block0 - test and practice
type = nan(24, 1);
targ_dir = nan(24, 1);
pt = nan(24,1);
retention = nan(24, 1);
iti = nan(24, 1);
ctch = zeros(24, 1);

type(1:8) = 0;
type(9:16) = 1;
type(17:24) = 2;

targ_dir_seed = [1:4, 1:4];
targ_dir(1:8) = targ_dir_seed(randperm(8))';
targ_dir(9:16) = targ_dir_seed(randperm(8))';
targ_dir(17:24) = targ_dir_seed(randperm(8))';

pt_seed = linspace(.100, .600, 8);
pt(1:8) = pt_seed(randperm(8))';
pt(9:16) = pt_seed(randperm(8))';
pt(17:24) = pt_seed(randperm(8))';

ret_seed = linspace(6,8,8);
retention(1:8) = ret_seed(randperm(8))';
retention(9:16) = ret_seed(randperm(8))';
retention(17:24) = ret_seed(randperm(8))';

iti_seed = linspace(2,4,8);
iti(1:8) = iti_seed(randperm(8))';
iti(9:16) = iti_seed(randperm(8))';
iti(17:24) = iti_seed(randperm(8))';

stim_wait_seed = linspace(2,4,8);
stim_wait(1:8) = stim_wait_seed(randperm(8))';
stim_wait(9:16) = stim_wait_seed(randperm(8))';
stim_wait(17:24) = stim_wait_seed(randperm(8))';

block0_inds = 1:24;
trial_target_numbers(block0_inds) = targ_dir;
trial_type(block0_inds) = type;
prescribed_PT(block0_inds) = pt;
trial_ret_period(block0_inds) = retention;
trial_ITI(block0_inds) = iti;
trial_stim_wait(block0_inds) = stim_wait;

%% Blocks 1 - 4
block_inds_ = 24 + (1:(4*27));
block_inds = reshape(block_inds_, 27, 4);
for i_block = 1:4
    block_mat = randomize_full_block;

    trial_type(block_inds(:, i_block)) = block_mat(:,1);
    trial_target_numbers(block_inds(:, i_block)) = block_mat(:,2);
    prescribed_PT(block_inds(:, i_block)) = block_mat(:,3);
    trial_ret_period(block_inds(:, i_block)) = block_mat(:,4);
    trial_ITI(block_inds(:, i_block)) = block_mat(:,5);
    trial_stim_wait(block_inds(:, i_block)) = block_mat(:, 6);
end

%% save subject data
save(['trial_parameters_', subject_code], 'trial_target_numbers', 'trial_type', 'prescribed_PT', 'trial_ret_period', 'trial_ITI', 'trial_stim_wait');

function block_mat = randomize_full_block

type = nan(27, 1);
targ_dir = nan(27, 1);
pt = nan(27,1);
retention = nan(27, 1);
iti = nan(27, 1);
stim_wait = nan(27, 1);

% designate catch trials, and randomize o1:8rder of trials
block_inds_ = 9:27; % exclude first 8 trials from catch trials
block_inds_rand = block_inds_(randperm(length(block_inds_)));
ctch_inds = block_inds_rand(1:3);
reg_inds_ = setdiff([1:8, block_inds_rand], ctch_inds);
reg_inds = reg_inds_(randperm(length(reg_inds_)));

% assign trial type and target direction together; equal # of trial_ret_perioddirections
% per type
type_seed = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2];
type(reg_inds) = type_seed;

targ_dir_seed = [1:4, 1:4, 1:4, 1:4, 1:4, 1:4];
targ_dir(reg_inds) = targ_dir_seed;

% assign PT, retention, and iti randomly without regard to type or target
pt_seed = linspace(.100, .600, 8);
pt(reg_inds(1:8)) = pt_seed(randperm(8))';
pt(reg_inds(9:16)) = pt_seed(randperm(8))';
pt(reg_inds(17:24)) = pt_seed(randperm(8))';

ret_seed = linspace(6,8,8);
retention(reg_inds(1:8)) = ret_seed(randperm(8))';
retention(reg_inds(9:16)) = ret_seed(randperm(8))';
retention(reg_inds(17:24)) = ret_seed(randperm(8))';

iti_seed = linspace(2,4,8);
iti(reg_inds(1:8)) = iti_seed(randperm(8))';
iti(reg_inds(9:16)) = iti_seed(randperm(8))';
iti(reg_inds(17:24)) = iti_seed(randperm(8))';

stim_wait_seed = linspace(2,4,8);
stim_wait(reg_inds(1:8)) = stim_wait_seed(randperm(8))';
stim_wait(reg_inds(9:16)) = stim_wait_seed(randperm(8))';
stim_wait(reg_inds(17:24)) = stim_wait_seed(randperm(8))';

ctch_pt_seed = .2*rand(1,3) + .5; % ensure catch pt is 500 < pt < 700
ctch_ret_seed = 2*rand(1,3) + 6; % same distribution as non-ctch
ctch_iti_seed = 2*rand(1,3) + 2; % same distribution as non-ctch
ctch_stim_wait_seed = 2*rand(1,3) + 2; % same distribution as non-ctch
ctch_type_seed_1 = [3 3 4 4];
ctch_type_seed = ctch_type_seed_1(randperm(length(ctch_type_seed_1)));
ctch_dir_seed_1 = [1 2 3 4];
ctch_dir_seed = ctch_dir_seed_1(randperm(length(ctch_dir_seed_1)));
type(ctch_inds) = ctch_type_seed(1:3);
targ_dir(ctch_inds) = ctch_dir_seed(1:3);
pt(ctch_inds) = ctch_pt_seed;
retention(ctch_inds) = ctch_ret_seed;
iti(ctch_inds) = ctch_iti_seed;
stim_wait(ctch_inds) = ctch_stim_wait_seed;

block_mat = [type, targ_dir, pt, retention, iti, stim_wait];
