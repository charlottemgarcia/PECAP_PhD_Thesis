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
% This script MP_PECAP_DAQ.m is 2 of 2 and includes the data acquisition
% portion of the PECAP experiment that collects the Diagonal of the PECAP
% matrix and subsequently the entire PECAP matrix. MP_Loudness_DAQ.m must
% be run before this script.
%
% Required Software: 
%       - MP_Loudness_DAQ.m (must be run before this script)
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


%% Collect ECAPs for the Diagonal of M

% Best to record ECAPs at at least 50 sweeps
param.sweeps = 50;                  % number of sweeps 

% hide or show BEDCS and other things
param.extraplot = 1;            % switch to 1 to see ABCD frames plotted
param.BEDCSVisible = 1;         % switch to 1 to make BEDCS visible 

% check to make sure there is a param structure to pull from
if ~isfield(param,'MCL_levels')
    error('You must run MP_Loudness_DAQ.m prior to collecting PECAP')
end

% set recording condition
param.condition = 'MP';

% set base filename with the correct sweeps
param.exp_file = ...
    sprintf('ForwardMaskingMP_gain_%d_%dkHz_%d_repeats.bExp',...
    param.gain, param.fs_kHz, param.sweeps); 

% print out recording parameters
fprintf(['\nInitiating Standard Monopolar PECAP DAQ Sequence\n\n']);
fprintf(['PECAP Recording Parameters for Participant ' param.ID ':\n']);
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
fprintf(['\tLoudness Scaling was done for the following' ...
    '\n\t\telectrodes: ']);
for ii = 1:length(param.loudness_electrodes)
    if ii == length(param.loudness_electrodes)
        fprintf(['and ' num2str(param.loudness_electrodes(ii))]);
    else
        fprintf([num2str(param.loudness_electrodes(ii)) ', ']);
    end
end
fprintf(['\n\tLoudness Scaling was done for the following' ...
    '\n\t\tlevels: ']);
for ii = 1:length(param.loudness_options)
    if ii == length(param.loudness_options)
        fprintf(['and ' num2str(param.loudness_options(ii))]);
    else
        fprintf([num2str(param.loudness_options(ii)) ', ']);
    end
end

% ask the user if they want to record at a different level
loudnesslevel = -1;
fprintf('\n\nAt what Loudness Level would you like to record PECAP?');
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

% initiate loudness vectors
electrodes = param.loudness_electrodes(1):param.loudness_electrodes(end);
MCL_levels = zeros(1,length(electrodes));

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

% check parameters before commencing diagonal recording
cont = 1;
fprintf('\nWould you like to manually enter new MCL levels?\n')
while cont ~= 0
    cont = input('Yes (1) or No (0): ');
    % enters this loop of the user wants to manually adjust MCL levels
    if cont == 1
        fprintf('\nEnter the desired MCL levels below\n');
        % updates the MCL level vector (in Device Units - linear scale)
        param.MCL_levels = input('MCL Vector (Device Units): ');
        % displays the new MCL levels 
        if length(param.electrodes) == length(param.MCL_levels)
            fprintf('\nManually Entered MCL Levels:\n')
            fprintf('El.\t\tMCL (Device Units)\tdB re 1 uA\n');
            for ii = 1:length(param.electrodes)
                fprintf([num2str(param.electrodes(ii)) '\t\t' ...
                    num2str(param.MCL_levels(ii)) '\t\t\t\t\t' ...
                    num2str(20*log10(param.MCL_levels(ii)),3) '\n']);
            end
            % asks if the new manually entered levels are correct
            fprintf('\nRe-enter MCL Levels?\n');
        % or throws an error if the wrong number of levels is entered
        else
            fprintf('\nMismatch # of electrodes & MCL levels.\n');
            param.MCL_levels = round(10.^(MCL_levels_ext/20));
            fprintf('\nRe-enter MCL Levels?\n')
            fprintf('(Answering No will revert to automatic levels)\n')
        end
    end
end 
clear MCL_levels_ext

% allow user to set the number of sweeps to collect
fprintf('\nECAP Recording is currently set to %d sweeps\n',param.sweeps);
fprintf('\tWould you like to change this?\n');
fprintf('\tIf you select No, the PECAP Diagonal will commence.\n\t');
cont = input('Yes (1) or No (0): ');
while cont ~= 0
    % enters this loop if the user wants to adjust the number of sweeps
    if cont == 1
        fprintf('Please enter a number of sweeps from the following: ')
        for ii = 1:length(param.available_sweeps)
            fprintf([num2str(param.available_sweeps(ii)) ' ']);
        end
        fprintf('\n\t');
        param.sweeps = input('Number of Sweeps: ');
        % throws an error and requires a new sweeps parameter if the
        % enetred value is not included in the available number of sweeps
        % parameters
        while nnz(param.available_sweeps - param.sweeps) == ...
                length(param.available_sweeps)
            fprintf(['Error. Please enter a number of sweeps from the' ...
                'following: ']);
            for ii = 1:length(param.available_sweeps)
                fprintf([num2str(param.available_sweeps(ii)) ' ']);
            end
            fprintf('\n\t');
            param.sweeps = input('Number of Sweeps: ');
        end
        % displays updated sweep parameter before commencing ECAP recording
        fprintf('\nECAP Recording is now set to %d sweeps\n',param.sweeps);
        fprintf('\tWould you like to change this?\n\t');
        fprintf('If you select No, the PECAP Diagonal will commence.\n\t');
        cont = input('Yes (1) or No (0): ');
    end
end

% update the base file based on the specified number of sweeps
param.exp_file = ...
    sprintf('ForwardMasking_AP_MP_gain_%d_%dkHz_%d_repeats.bExp',...
    param.gain, param.fs_kHz, param.sweeps); 

% collect diagonal ECAPs & record time
Diagonal_ECAP_Data = run_diagonal_ecap(param);

% plot diagonal ECAPs
plot_Diagonal(Diagonal_ECAP_Data,param);

% save diagonal data to file
datetimestr = datestr(datetime);
save(['data/Subj' param.ID '_cond' param.condition '_DiagData_' date ...
    '_' datetimestr(13:14) '-' datetimestr(16:17) '-', ...
    datetimestr(19:20)], 'param', 'Diagonal_ECAP_Data');
clear datetimestr


%% Confirm that levels are as desired for PECAP (maybe re-record diagonal)

% check that the levels are comfortable for the participant. If not, then
% change them so that they are.
cont = 0;
while cont ~= 1
    fprintf('\nDo you want to record PECAP at this volume?\n');
    fprintf('\tIf you select yes, PECAP will commence.\n');
    fprintf('\tIf you select no, you can change the volume.\n\t');
    cont = input('Yes (1) or No (0): ');
    % enters this loop if the user needs to change the volume
    if cont == 0
        % displays current loudness level parameters
        fprintf(['\nMCL (' num2str(param.MCL_loudnesslevel) ') Levels\n']);
        fprintf('El.\t\tMCL (Device Units)\tdB re 1 uA\n');
        for ii = 1:length(param.electrodes)
            fprintf([num2str(param.electrodes(ii)) '\t\t' ...
                num2str(param.MCL_levels(ii)) '\t\t\t\t\t' ...
                num2str(20*log10(param.MCL_levels(ii)),3) '\n']);
        end
        % prompts the user to enter the reduction in volume
        fprintf('What Loudness Level would you like to switch to?\n')
        loudnesslevel = -1;
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
        
        % initiate loudness vectors
        electrodes = ...
            param.loudness_electrodes(1):param.loudness_electrodes(end);
        MCL_levels = zeros(1,length(electrodes));
        
        % enter scaled loudness levels
        for ii = 1:length(MP_Loudness)
            for jj = 1:length(MP_Loudness(ii).ECAP_Struct)
                if MP_Loudness(ii).ECAP_Struct(jj).LoudnessLevel == ...
                        loudnesslevel
                    % converts device units (linear) to dB re 1 uA to 
                    % interpolate
                    MCL_levels(MP_Loudness(ii).ECAP_Struct(jj).Probe - ...
                        min(electrodes) + 1) = ...
                        20*log10(MP_Loudness(ii).ECAP_Struct(jj).Probe_lvl);
                end
            end
        end

        % interpolate missing values by interpolating linearly on dB scale
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
                % create vector sized to fill in  interpolated electrodes
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
                MCL_levels_ext((MCL_idx(ii)+1):(MCL_idx(ii+1)-1)) = ...
                    interp_vals;
            end
        end
        MCL_levels = round(10.^(MCL_levels_ext/20));

        % store MCL levels in params structure
        param.MCL_levels = MCL_levels;
        param.electrodes = electrodes;
        param.MCL_loudnesslevel = loudnesslevel;
        
        % prints out new loudness levels
        fprintf(['\nMCL (' num2str(param.MCL_loudnesslevel) ') Levels\n']);
        fprintf('El.\t\tMCL (Device Units)\tdB re 1 uA\n');
        for ii = 1:length(param.electrodes)
            fprintf([num2str(param.electrodes(ii)) '\t\t' ...
                num2str(param.MCL_levels(ii)) '\t\t\t\t\t' ...
                num2str(20*log10(param.MCL_levels(ii)),3) '\n']);
        end
        
        % pauses to allow the user to view the new levels before commencing
        % the recording sequence
        nth = input('\nClick enter to collect Diagonal again');
        
        % collect diagonal ECAPs
        Diagonal_ECAP_Data = run_diagonal_ecap(param);
        
        % save diagonal data to file
        datetimestr = datestr(datetime);
        save(['data/Subj' param.ID '_cond' param.condition '_DiagData_' ...
            date '_' datetimestr(13:14) '-' datetimestr(16:17) '-', ...
            datetimestr(19:20)], 'param', 'Diagonal_ECAP_Data');
        clear datetimestr
        
        % plot diagonal ECAPs (this updates previous )
        plot_Diagonal(Diagonal_ECAP_Data,param);
    end
end
% clear workspace
clear nth ii cont change loudnesslevel MCL_levels electrodes interps 
clear interp_vals step span ii MCL_idx_val MCL_idx

%% Collect Monopolar PECAP

% hide BEDCS visible while collecting PECAP or the BC frames will cause
% errors in smoothing within BEDCS
param.BEDCSVisible = 0;

% if they have indicated that the levels are fine, collect normal PECAP
PECAP_Data = run_PECAP(param);

% plot M matrix for collected PECAP Data
PECAP_Data = plot_PECAP(PECAP_Data, param);

% plot M matrix in SOE curves
plot_SOEs(PECAP_Data, param);

% set marker to indicate that MP PECAP has been collected
param.MP_PECAP_Collected = 1;

% save PECAP data to file
datetimestr = datestr(datetime);
save(['data/Subj' param.ID '_cond' param.condition '_PECAPData_' date ...
    '_' datetimestr(13:14) '-' datetimestr(16:17) '-', ...
    datetimestr(19:20)], 'param', 'PECAP_Data');
clear datetimestr ans jj