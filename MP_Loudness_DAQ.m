%% Panoramic ECAP Data Collection Script for Advanced Bionics
%
% Garcia et al 2021 described a method for estimating neural activation
% patterns and separating current spread and neural health variation along
% the length of an electrode array for individual Cochlear Implant (CI) 
% users. This method requires measuring Electrically-Evoked Compound
% Action-Potentials (ECAPs) using the Forward-Masking Artefact-Reduction
% Technique for every combination of masker and probe electrodes at the
% Maximum Comfortable Level (MCL).
% This script allows the user to collect the standard Panoramic ECAP
% (PECAP) matrix. It requires the running of BEDCS 1.8 in the background
% and runs a Graphical User Interface (GUI) that allows the user to
% loudness-scale multiple electrodes along the implant array for an
% Advanced Bionics CI user. These loudness levels are then used to first
% record ECAPs for every electrode at MCL using the forward-masking
% artefact-reduction technique where the probe and the masker are
% co-located on the same electrode. The user is able to then adjust the
% volume if it is uncomfortably loud for the CI user before collecting the
% entire PECAP Matrix. Time is saved throughout by skipping recordings of
% redundant ECAP recording frames within the forward-masking paradigm. 
%
% This script MP_Loudness_DAQ.m is 1 of 2 and includes the loudness-scaling
% portion of the experiment. It must be run before MP_PECAP_DAQ.m. 
% 
% Required Software: 
%       - Bionic Ear Data Collection System (BEDCS) version 1.8
%       - Note: this code has only been debugged with versions 1.18.321 and
%               1.18.337
%       - MATLAB 2018a
%         Note: this script may work with other versions of MATLAB but the
%               only version that has been debugged / checked is 2018a
% Required Hardware:
%       - Advanced Bionics Research Hardware required for measuring ECAPs
%         (i.e. a Clarion Programming Cable, Programming Interface, 
%         associated Power Supply, Advanced Bionics
%       - Optional: Advanced Bionics Load board (for testing prior to use
%         with a patient)
% 
% Note: This script does not include checks to confirm that stimulation
% levels are presented to the CI patient within compliance. Please check
% this using separate software prior to use.
%
% Usage Disclaimer: Users of this script accept responsibility for safety
% checks undertaken whilst testing with CI patients. The developers of this
% software and their institutions (The MRC Cognition & Brain Sciences Unit
% and the University of Cambridge) accept no responsibility for the safety
% of the application of this script outside their institution. 
% 
% Reference:
% (If you use this software, please cite the following publication that
% contains details of how the updated PECAP2 method works)
% Garcia, C., Goehring, T., Cosentino, S. et al. The Panoramic ECAP Method: 
% Estimating Patient-Specific Patterns of Current Spread and Neural Health 
% in Cochlear Implant Users. JARO (2021). 
% https://doi.org/10.1007/s10162-021-00795-2
%
% We hope you will find the software useful, but neither the authors nor
% their employers accept any responsibility for the consequencies of its
% use. The USER is responsible for ensuring the saftety of any subjects
% during testing.
%
% legal disclaimer from the University of Cambridge:
% The Software is the result work conducted within the MRC Cognition & 
% Brain Sciences unit at the University of Cambridge (the “University”)
% "The Software shall be used for non-commercial research purposes only.  
% The USER will not reverse engineer, manufacture, sell or sublicense for 
% manufacture and sale upon a commercial basis the Software, incorporate it 
% into other software or products or use the Software other than herein 
% expressly specified. The USER agrees that it will not translate, alter, 
% decompile, disassemble, reverse engineer, reverse compile, attempt to 
% derive, or reproduce source code from the Software. The USER also agrees 
% that it will not remove any copyright or other proprietary or product 
% identification notices from the Software. The Software is provided 
% without warranty of merchantability or fitness for a particular purpose 
% or any other warranty, express or implied, and without any representation 
% or warranty that the use or supply of the Software will not infringe any 
% patent, copyright, trademark or other right. In no event shall the 
% University or their staff and students be liable for any use by the USER 
% of the Software. The supply of the Software does not grant the USER any 
% title or interest in and to the Software, the related source code, and 
% any and all associated patents, copyrights, and other intellectual 
% property rights. The University and their staff and students reserve the 
% right to revise the Software and to make changes therein from time to 
% time without obligation to notify any person or organisation of such 
% revision or changes. While the University will make every effort to 
% ensure the accuracy of Software and the data contained therein, however 
% neither the University nor its employees or agents may be held 
% responsible for errors, omissions or other inaccuracies or any 
% consequences thereof."

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Developed by Charlotte Garcia, 2021                             %
% Loudness-Scaling GUI provided by Francois Guerit                       %
% MRC Cognition & Brain Sciences Unit, Cambridge, UK                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add source files to path
[current_path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(current_path, 'source'))
clear current_path


%% Enter Parameters

% Participant parameters
fprintf(['\nInitiating Standard Monopolar PECAP Loudness Scaling ' ...
    'Sequence\n\n']);
param.condition = 'MP';         % 'MP' for monopolar mode
if ~isfield(param,'ID')
    fprintf('Please enter your research participant ID\n');
    fprintf('\tNote: do this in ''string'' form\n');
    fprintf('\ti.e. ''AB01'' instead of AB01\n\t');
    param.ID = 0; 
    param.ID = input('Research Participant ID: ');
    while ~ischar(param.ID)
        fprintf('\n\tID must be a string.');
        fprintf('\n\tPlease put single quotes around text.\n')
        param.ID = input('\nResearch Participant ID: ');
    end
end

% if pTP Loudness has not already been run and the basic recording
% parameters have not been selected from there, then enter basic recording
% parameters now
if isfield(param,'pTP_electrodes')
    fprintf(['Note that recording parameters are being imported to ' ...
        'be \n\tconsistent with pTP stimulation']);
else
    % basic parameters for Monopolar PECAP
    param.debug = 0;
    param.extraplot = 1;        % switch to 1 to see ABCD frames plotted
    param.BEDCSVisible = 1;     % switch to 1 to make BEDCS visible 
    param.available_sweeps = [1 5 10 25 50 75 100 125 150 175 200];

    % Recording parameters

    % enter vector of electrodes to scale loudness (rest will be interpolated)
    %       if electrodes are turned off at the end of the array, simply
    %       exclude them from this vector (i.e. for a patient for whom
    %       electrodes 1, 2, and 16 are switched off, consider scaling a 
    %       loudness vector of [3 6 9 12 15] or [3 5 7 9 11 13 15]
    %       NOTE: This script currently does not support data collection 
    %       for patients who have an electrode switched off in the middle 
    %       of their array. To handle this ad-hoc, exclude the mid-array
    %       electrode from the loudness_electrodes vector below, and
    %       manually set the current level to 0 before collecting ECAPs for 
    %       the diagonal.
    %       This vector must contain integers between 1 and 16.
    param.loudness_electrodes = [1 4 7 10 13 16];
    param.electrodes = ...
        param.loudness_electrodes(1):param.loudness_electrodes(end);

    % enter the loudness scaling options
    %       The loudness GUI will not allow you to click the 'done' button
    %       until you have recorded responses for each of the following 
    %       levels. It will also enable you to collect PECAP at any of the 
    %       specified loudness levels and will interpolate the current 
    %       level for the un-scaled electrodes.
    %       This vector must contain integers between 1 and 10.
    %       Suggested: [5 6 7]
    param.loudness_options = [5 6 7 8];

    % enter other recording parameters
    param.gain = 1000;          % gain of the amplifier [1,3,100,300,1000]
    param.fs_kHz = 56;          % amplifier sampling frequency [56,28,9kHz]
    param.rel_rec_EL = -3;      % -3 if recording apically, +3 if basally
    param.sweeps = 10;          % number of sweeps (restricted as in ln125)

    % enter stimulus properties
    param.phase_duration_us = 10.776*4;     % has to be multiples of 10.776
    param.masker_probe_delay_us = 600;      % MPI (400us in Cochlear)
    param.level_start_masker = 100;         % these will be adjusted in GUI
    param.level_start_probe = 100;          % these will be adjusted in GUI
end

% assign BEDCS input filename based on set parameters
% While 50 should be used to record, 10 sweeps are sufficient for loudness
param.exp_file = ...
    sprintf('ForwardMaskingMP_gain_%d_%dkHz_%d_repeats.bExp', ...
    param.gain, param.fs_kHz, param.sweeps); 

% print out recording parameters
fprintf(['\n\nPECAP Recording Parameters for Participant ' param.ID ...
    ':\n']);
fprintf(['\tAmplifier Gain: \t\t\t\t' num2str(param.gain) '\n']);
fprintf(['\tSampling Frequency: \t\t\t' num2str(param.fs_kHz) ' kHz\n']);
if param.rel_rec_EL ~= 3 && param.rel_rec_EL ~= -3
    error(['Recording electrode must be 3 electrodes away from the ' ...
        'probe. Please enter 3 or -3 for param.rel_rec_EL and restart.']);
elseif param.rel_rec_EL == 3
    fprintf('\tRecording Electrode:\t\t\t3 basal to the probe\n');
elseif param.rel_rec_EL == -3
    fprintf('\tRecording Electrode:\t\t\t3 apical to the probe\n');
end
if mod(param.phase_duration_us,10.776)~=0
    error(['Phase Duration must be a multiple of 10.776 microseconds. ' ...
        'Please enter a multiple of 10.776 for param.phase_duration_us' ...
        ' and restart.']);
else
    fprintf(['\tPhase Duration: \t\t\t\t' ...
        num2str(param.phase_duration_us) ' us\n']);
end
fprintf(['\tMasker-Probe Interval (MPI): \t' ...
    num2str(param.masker_probe_delay_us) ' us\n']);
fprintf(['\tLoudness Scaling sweeps: \t\t' num2str(param.sweeps) '\n']);
if nnz(param.available_sweeps - param.sweeps) == ...
        length(param.available_sweeps)
    fprintf(['Error. Please enter a number of sweeps from the' ...
        '\n\tfollowing and restart: ']);
    for ii = 1:length(param.available_sweeps)
        fprintf([num2str(param.available_sweeps(ii)) ' ']);
    end
    fprintf('\n\t');
    error(' ');
end 
fprintf(['\tLoudness Scaling will be done for the following' ...
    '\n\t\telectrodes: ']);
for ii = 1:length(param.loudness_electrodes)
    if ii == length(param.loudness_electrodes)
        fprintf(['and ' num2str(param.loudness_electrodes(ii))]);
    else
        fprintf([num2str(param.loudness_electrodes(ii)) ', ']);
    end
end
fprintf(['\n\tLoudness Scaling will required for the following' ...
    '\n\t\tlevels: ']);
for ii = 1:length(param.loudness_options)
    if ii == length(param.loudness_options)
        fprintf(['and ' num2str(param.loudness_options(ii))]);
    else
        fprintf([num2str(param.loudness_options(ii)) ', ']);
    end
end

% pause to allow user to cancel and adjust parameters if desired
fprintf('\n\nAre these the parameters you would like to continue with?\n')
fprintf('\tIf no, press ctrl+c, and adjust parameters in the script.\n\t')
nth = input('If yes, press enter to initiate the Loudness Scaling GUI');


%% Call the GUI function & save loudness results & ECAPs in ECAP_Struct

% initiate loudness ECAP storage structure
MP_Loudness = struct();

% set parameters for loudness scaling for each electrode
for ii = 1:length(param.loudness_electrodes)
    param.electrode_masker = param.loudness_electrodes(ii);
    param.electrode_probe = param.loudness_electrodes(ii);
    if param.loudness_electrodes(ii) + param.rel_rec_EL >= 1 && ...
            param.loudness_electrodes(ii) + param.rel_rec_EL <= 16
        param.rec_EL = param.loudness_electrodes(ii) + param.rel_rec_EL;
    else
        param.rec_EL = param.loudness_electrodes(ii) - param.rel_rec_EL;
    end
    % call ECAP recording & loudness GUI
    MP_Loudness(ii).ECAP_Struct = run_loudness_and_ecap(param);
end


%% Interpolate loudness levels for all electrodes

% initiate loudness vectors
electrodes = param.loudness_electrodes(1):param.loudness_electrodes(end);
MCL_levels = zeros(1,length(electrodes));

% request loudness level to collect ECAPs at
loudnesslevel = -1;
fprintf('\nAt what Loudness Level would you like to record PECAP?');
fprintf('\nRecomended: 6 (MCL)');
while nnz(param.loudness_options - loudnesslevel) == ...
        length(param.loudness_options)
    fprintf('\nPlease enter a value from ');
    for ii = 1:length(param.loudness_options)
        if ii == length(param.loudness_options)
            fprintf(['and ' num2str(param.loudness_options(ii))]);
        else
            fprintf([num2str(param.loudness_options(ii)) ', ']);
        end
    end
    fprintf('\n');
    loudnesslevel = input('Loudness Recording Level: ');
end

% enter scaled loudness levels
for ii = 1:length(MP_Loudness)
    for jj = 1:length(MP_Loudness(ii).ECAP_Struct)
        if MP_Loudness(ii).ECAP_Struct(jj).LoudnessLevel == loudnesslevel
            % converts device units (linear) to dB re 1 uA to interpolate
            MCL_levels(MP_Loudness(ii).ECAP_Struct(jj).Probe - ...
                min(electrodes) + 1) = ...
                20*log10(MP_Loudness(ii).ECAP_Struct(jj).Probe_lvl);
        end
    end
end

% interpolate missing values by interpolating linearly on a dB scale
MCL_idx = find(MCL_levels);
MCL_idx_val = MCL_levels(MCL_idx);
MCL_levels_ext = MCL_levels;
for ii = 1:(length(param.loudness_electrodes)-1)
    % if there are non-scaled electrodes between the scaled ones
    if MCL_idx(ii+1) - MCL_idx(ii) ~= 1
        % determine the number of electrodes between them
        interps = MCL_idx(ii+1) - MCL_idx(ii) - 1;
        % determine the difference between the two scaled values
        span = MCL_idx_val(ii+1) - MCL_idx_val(ii);
        % determine the step size based on the difference & spacing
        step = span/(interps+1);
        % create vector sized to fill in the interpolated electrodes
        interp_vals = zeros(1,interps);
        % enter the
        interp_vals(1) = MCL_idx_val(ii) + step;
        % fill in the rest of the vector
        if interps > 1
            for jj = 2:interps
                interp_vals(jj) = interp_vals(jj-1) + step;
            end
        end
        % put the interpolated values in the MCL vector
        MCL_levels_ext((MCL_idx(ii)+1):(MCL_idx(ii+1)-1)) = interp_vals;
    end
end
MCL_levels = round(10.^(MCL_levels_ext/20));

% store MCL levels in params structure
param.MCL_levels = MCL_levels;
param.electrodes = electrodes;
param.MCL_loudnesslevel = loudnesslevel;

% save PECAP data to file
datetimestr = datestr(datetime);
save(['data/Subj' param.ID '_cond' param.condition '_LoudnessData_' date ...
    '_' datetimestr(13:14) '-' datetimestr(16:17) '-', ...
    datetimestr(19:20)], 'param', 'MP_Loudness');

% clear workspace
clear start_val MCL_idx secnd_val idx jj interps interp_vals step span ii
clear electrodes MCL_levels datetimestr MCL_idx_val loudnesslevel MCL_idx


%% Display the interpolated MCL levels

% display loudness parameters
fprintf(['\nLogarithmically Extrapolated MCL (' ...
    num2str(param.MCL_loudnesslevel) ') Levels\n']);
fprintf('El.\t\tMCL (Device Units)\tdB re 1 uA\n');
for ii = 1:length(param.electrodes)
    % prints an asterisk for the electrodes that were loudness scaled
    if nnz(param.loudness_electrodes - param.electrodes(ii)) ~= ...
            length(param.loudness_electrodes)
        fprintf([num2str(param.electrodes(ii)) '*\t\t' ...
            num2str(param.MCL_levels(ii)) '\t\t\t\t\t' ...
            num2str(20*log10(param.MCL_levels(ii)),3) '\n']);
    % leaves out the asterisk for the interpolated electrodes
    else
        fprintf([num2str(param.electrodes(ii)) '\t\t' ...
            num2str(param.MCL_levels(ii)) '\t\t\t\t\t' ...
            num2str(20*log10(param.MCL_levels(ii)),3) '\n']);
    end
end
fprintf('* indicates loudness-scaled electrodes\n');

clear ii nth MCL_levels_ext