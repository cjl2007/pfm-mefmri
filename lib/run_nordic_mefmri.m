function run_nordic_mefmri(fn_magn_in, fn_phase_in, fn_out_prefix, nordic_code_dir, log_path)
%RUN_NORDIC_MEFMRI Thin pipeline wrapper around the official NIFTI_NORDIC.

if nargin < 4 || isempty(nordic_code_dir)
    nordic_code_dir = '/home/charleslynch/Desktop/NORDIC-Test/code';
end
if nargin < 5
    log_path = '';
end

if ~isempty(log_path)
    log_dir = fileparts(log_path);
    if ~isempty(log_dir) && ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end
    diary(log_path);
    cleanup_obj = onCleanup(@() diary('off')); %#ok<NASGU>
end

fprintf('NORDIC wrapper start: %s\n', datestr(now, 31));
fprintf('Magnitude input: %s\n', fn_magn_in);
if isempty(fn_phase_in)
    fprintf('Phase input: <none; magnitude-only NORDIC>\n');
else
    fprintf('Phase input: %s\n', fn_phase_in);
end
fprintf('Output prefix: %s\n', fn_out_prefix);
fprintf('NORDIC code dir: %s\n', nordic_code_dir);

if ~exist(fn_magn_in, 'file')
    error('Missing magnitude input: %s', fn_magn_in);
end
if ~isempty(fn_phase_in) && ~exist(fn_phase_in, 'file')
    error('Missing phase input: %s', fn_phase_in);
end
if ~exist(fullfile(nordic_code_dir, 'NIFTI_NORDIC.m'), 'file')
    error('Missing NIFTI_NORDIC.m under: %s', nordic_code_dir);
end

addpath(nordic_code_dir);

[out_dir, out_base, out_ext] = fileparts(fn_out_prefix);
out_base = [out_base out_ext];
if isempty(out_dir)
    out_dir = pwd;
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

ARG = [];
ARG.DIROUT = [out_dir filesep];
ARG.noise_volume_last = 0;
ARG.NORDIC = 1;
ARG.MP = 0;
ARG.temporal_phase = 1;
ARG.phase_filter_width = 10;
ARG.write_gzipped_niftis = 1;
if isempty(fn_phase_in)
    ARG.magnitude_only = 1;
end

fprintf('ARG.noise_volume_last = %d\n', ARG.noise_volume_last);
fprintf('ARG.NORDIC = %d\n', ARG.NORDIC);
fprintf('ARG.MP = %d\n', ARG.MP);
fprintf('ARG.temporal_phase = %d\n', ARG.temporal_phase);
fprintf('ARG.phase_filter_width = %d\n', ARG.phase_filter_width);
fprintf('ARG.write_gzipped_niftis = %d\n', ARG.write_gzipped_niftis);
if isfield(ARG, 'magnitude_only')
    fprintf('ARG.magnitude_only = %d\n', ARG.magnitude_only);
end

NIFTI_NORDIC(fn_magn_in, fn_phase_in, out_base, ARG);

expected = fullfile(out_dir, [out_base '.nii.gz']);
if ~exist(expected, 'file')
    error('Expected NORDIC output was not created: %s', expected);
end

fprintf('NORDIC wrapper complete: %s\n', expected);
fprintf('NORDIC wrapper end: %s\n', datestr(now, 31));
end
