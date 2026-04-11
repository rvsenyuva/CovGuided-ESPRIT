% =========================================================
% Script: setup_path.m
% Purpose: Add all repository subfolders to MATLAB path.
%          Run this once per MATLAB session before any script.
%          (Entry-point scripts in experiments/ are self-contained
%           and call addpath internally; this file is provided as
%           a convenience for interactive use from the repo root.)
% Paper  : Sensors-100615-2025
% Version: v1.0 (updated for modernized repo structure)
% =========================================================
%
% USAGE
%   >> cd /path/to/CovGuided-ESPRIT
%   >> setup_path
%   >> cd experiments
%   >> runMainExperiment   % or any other entry-point script

repo_root = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(repo_root, 'src')));   % +covguided package
addpath(fullfile(repo_root, 'config'));          % getSimParams.m
addpath(fullfile(repo_root, 'validation'));      % validateReproducibility.m
addpath(fullfile(repo_root, 'experiments'));     % entry-point scripts

fprintf('[setup_path] Repository root: %s\n', repo_root);
fprintf('[setup_path] Added to path: src/+covguided, config, validation, experiments\n');
fprintf('[setup_path] Call functions as: covguided.tlsEsprit(...) etc.\n');
fprintf('[setup_path] Run ''which covguided.tlsEsprit'' to confirm.\n');
