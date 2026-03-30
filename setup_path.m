% =========================================================
% Script: setup_path.m
% Purpose: Add all repository subfolders to MATLAB path.
%          Run this once per MATLAB session before any script.
% Paper  : Sensors-100615-2025
% =========================================================
%
% USAGE
%   >> cd /path/to/CovGuided-ESPRIT
%   >> setup_path
%   >> paperCodeMain        % or any other entry-point script

repo_root = fileparts(mfilename('fullpath'));

addpath(fullfile(repo_root, 'functions'));
addpath(fullfile(repo_root, 'params'));
addpath(fullfile(repo_root, 'scripts'));
addpath(fullfile(repo_root, 'validation'));

fprintf('[setup_path] Repository root: %s\n', repo_root);
fprintf('[setup_path] Added to path: functions, params, scripts, validation\n');
fprintf('[setup_path] Run ''which getSimParams'' to confirm.\n');
