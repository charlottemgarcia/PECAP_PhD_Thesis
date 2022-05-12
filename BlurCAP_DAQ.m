%% Blurred Panoramic ECAP Data Collection Script for Advanced Bionics
% 
% Garcia et al 2021 described a method for estimating neural activation
% patterns and separating current spread and neural health variation along
% the length of an electrode array for individual Cochlear Implant (CI) 
% users. This method requires measuring Electrically-Evoked Compound
% Action-Potentials (ECAPs) using the Forward-Masking Artefact-Reduction
% Technique for every combination of masker and probe electrodes at the
% Maximum Comfortable Level (MCL).
% Goehring et al 2020 described the effect of introducing spectral
% blurring on speech perception in cochlear implant users. SRTs began to
% rise above 3 blurred electrodes. This script intends to use the PECAP
% method to determine if increased current spread assumed to exist with
% blurred compared to monopolar stimulation is observable in the current
% spread estimate of the PECAP model. Prior to running this script,
% MP_PECAP_DAQ.m must be run in order to collect the standard, monopolar
% PECAP matrix. This script first calls a loudness GUI that allows the user
% to scale the loudness of a selected electrode to blur, at two blurring
% factors: 3 blurred electrodes, and 5 blurred electrodes. Once the
% loudness levels have been determined, one row and column of the PECAP
% matrix are recorded for each blurring level for the specified blurred
% electrode. These ECAPs are then combined with the Monopolar PECAP matrix
% in order to create a full 'BlurCAP' matrix for each of the two blurring
% factors.
% 
% Required Software: 
%       - MP_PECAP_DAQ.m (must be run before this script)
%       - Bionic Ear Data Collection System (BEDCS) version 1.8
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
% References:
% (If you use this software, please cite the following publications that
% contains details of how the updated PECAP2 method works)
% Garcia, C., Goehring, T., Cosentino, S. et al. The Panoramic ECAP Method: 
% Estimating Patient-Specific Patterns of Current Spread and Neural Health 
% in Cochlear Implant Users. JARO (2021). 
% https://doi.org/10.1007/s10162-021-00795-2
% Goehring, T., Arenberg, J. G., & Carlyon, R. P. (2020). Using Spectral 
% Blurring to Assess Effects of Channel Interaction on Speech-in-Noise 
% Perception with Cochlear Implants. Journal of the Association for 
% Research in Otolaryngology : JARO, 21(4), 353–371. 
% https://doi.org/10.1007/s10162-020-00758-z
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

%% Blurring Parameters

% set blurring condition
param.condition = 'BL';

% check to see if monopolar PECAP has been recorded already
if ~isfield(param,'MP_PECAP_Collected')
    error('Please record Monopolar PECAP before initiating BlurCAP');
end

% start parameter checking
fprintf('\nInitiating Blurred PECAP (BlurCAP) Sequence\n\n');

% specify which electrode to blur
fprintf('What electrode would you like to blur?\n');
fprintf('\t(Must be between e%d and e%d)\n\t', min(param.electrodes) + 4, ...
    max(param.electrodes) - 4); % must be at least 4 away from the edge 
param.blurred_electrode = input('Blurred Electrode: ');
if ~isnumeric(param.blurred_electrode)
    while ~isnumeric(param.blurred_electrode)
        fprintf('\nPlease enter a number for a blurring electrode\n\t');
        param.blurred_electrode = input('Blurred Electrode: ');
    end
end
if isnumeric(param.blurred_electrode)
    while param.blurred_electrode < min(param.electrodes) + 4 || ...
            param.blurred_electrode > max(param.electrodes) - 4
        fprintf(['\nPlease select an electrode to blur between e%d and' ...
            ' e%d\n\t'],min(param.electrodes) + 4, ...
            max(param.electrodes) - 4);
        param.blurred_electrode = input('Blurred Electrode: ');
    end
end

% set blurring parameters
param.sweeps = 10;                   % this is just for loudness testing
param.blurring_conditions = [3,5];   % 3, 5, and 7 (with holes) are goal

% print out recording parameters
fprintf(['\n\nBlurCAP Recording Parameters for Participant ' param.ID ...
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
fprintf(['\tBlurred Electrode: \t\t\t\t' ...
    num2str(param.blurred_electrode)]);
fprintf('\n\tBlurring Factors: \t\t\t\t');
for ii = 1:length(param.blurring_conditions)
    if ii == length(param.blurring_conditions)
        fprintf(['and ' num2str(param.blurring_conditions(ii))]);
    else
        fprintf([num2str(param.blurring_conditions(ii)) ' ']);
    end
end

% pause to allow user to cancel and adjust parameters if desired
fprintf('\n\nAre these the parameters you would like to continue with?\n')
fprintf('\tIf no, press ctrl+c, and adjust parameters in the script.\n\t')
nth = input('If yes, press enter to initiate the Loudness Scaling GUI');

%% Blurring Loudness

% initiate loudness ECAP storage structure
Blurring_Loudness = struct();

% cycle through loudness GUI for each of the 3 blurring conditions
for ii = 1:length(param.blurring_conditions)
    param.electrode_masker = param.blurred_electrode;
    param.electrode_probe = param.blurred_electrode;
    if param.blurred_electrode + param.rel_rec_EL >= 1 && ...
            param.blurred_electrode + param.rel_rec_EL <= 16
        param.rec_EL = param.blurred_electrode + param.rel_rec_EL;
    else
        param.rec_EL = param.blurred_electrode - param.rel_rec_EL;
    end
    while param.blurred_electrode < min(param.electrodes) + 3 || ...
            param.blurred_electrode > max(param.electrodes) - 3
        fprintf('\nCannot Blur chosen electrodes.\n')
        fprintf('Please enter a blurring electrode between e4-13.\n');
        fprintf('\t(At least 3 elecs away from the edge of the evaluated PECAP elecs)\n)')
    end
    % rename file for blurring level
    param.blur_factor = param.blurring_conditions(ii);
    param.exp_file = sprintf('ForwardMaskingBlur%d_gain_%d_%dkHz_%d_repeats.bExp',...
        param.blur_factor, param.gain, param.fs_kHz, ...
        param.sweeps); 
    % call ECAP recording & loudness GUI
    Blurring_Loudness(ii).ECAP_Struct = run_loudness_and_ecap(param);
    for jj = 1:length(Blurring_Loudness(ii).ECAP_Struct)
        Blurring_Loudness(ii).ECAP_Struct(jj).BlurFactor = param.blurring_conditions(ii);
    end
end

% add MCL levels for blurring conditions to param structure
param.Blurred_levels = zeros(1,length(Blurring_Loudness));
for ii = 1:length(Blurring_Loudness)
    for jj = 1:length(Blurring_Loudness(ii).ECAP_Struct)
        if Blurring_Loudness(ii).ECAP_Struct(jj).LoudnessLevel == ...
                param.MCL_loudnesslevel
            param.Blurred_levels(ii) = ...
                Blurring_Loudness(ii).ECAP_Struct(jj).Probe_lvl;
        end
    end
end

% save Blurring Loudness data to file
datetimestr = datestr(datetime);
save(['data/Subj' param.ID '_cond' param.condition '_LoudnessData_' ...
    date '_' datetimestr(13:14) '-' datetimestr(16:17) '-', ...
    datetimestr(19:20)], 'param', 'Blurring_Loudness');
clear datetimestr

%% BlurCAP Data Collection

% initiate structure & reset sweeps to data collection level
BlurCAP_Data = struct();
param.sweeps = 50;

% allow user to set the number of sweeps to collect
fprintf('\n\nBlurCAP Recording is currently set to %d sweeps\n',param.sweeps);
fprintf('It is recommended to use the same # of sweeps as Monopolar PECAP.\n');
fprintf('\tWould you like to change this?\n');
fprintf('\tIf you select No, the BlurCAP will commence.\n\t');
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
        fprintf('If you select No, BlurCAP will commence.\n\t');
        cont = input('Yes (1) or No (0): ');
    end
end

% Cycle through data collection
for ii = 1:length(param.blurring_conditions)
    % update blurring condition
    param.blur_factor = param.blurring_conditions(ii);
    % print out the blurring factor being collected
    fprintf('\nCommencing Data Collection for Blurring\n');
    fprintf(['\tBlurred Electrode: \t' num2str(param.blurred_electrode)]);
    fprintf(['\n\tBlurring Factor: \t' num2str(param.blur_factor)]);
    fprintf(['\n\tBlur' num2str(param.blurred_electrode) ' MCL (' ...
        num2str(param.MCL_loudnesslevel) ') Level: \t' ...
        num2str(param.Blurred_levels(find(~(param.blurring_conditions - ...
        param.blur_factor)))) ' Device Units (' ...
        num2str(round(20*log10(param.Blurred_levels(find(~(param.blurring_conditions ...
        - param.blur_factor)))),1)) ' dB re 1 uA)']);
    fprintf('\n\n');
    BlurCAP_Data(ii).Blurred_Electrode = param.blurred_electrode;
    BlurCAP_Data(ii).Blur_Factor = param.blur_factor;
    BlurCAP_Data(ii).Blurred_ECAPs = run_BlurCAP(param);
end

% plot M matrix for collected PECAP Data
BlurCAP_Data = plot_BlurCAP(PECAP_Data, BlurCAP_Data, param);

% save PECAP data to file
datetimestr = datestr(datetime);
save(['data/Subj' param.ID '_cond' param.condition '_PECAPData_' date ...
    '_' datetimestr(13:14) '-' datetimestr(16:17) '-', ...
    datetimestr(19:20)], 'param', 'BlurCAP_Data');
clear datetimestr nth ii jj cont
