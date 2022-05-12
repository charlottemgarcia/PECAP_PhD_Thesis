function pecap_daq(strOutputFilename, Rec_Elec, Diag, Gain_dB, Delay_us, PRate_Hz, sweeps, elecs, e_levels)

%--------------------------------------------------------------------------
% DISCLAIMER
%
% Software for collection of panoramic-ECAP (PECAP) data from Cochlear CI
% devices
%
% Cosentino et al (2015) described a method that allows reconstruction of
% the underlying neural activation pattern of individual channels from ECAP
% amplitudes. The method estimates the joint neural activation produced by
% two electrodes at a time (one acting as the probe and the other as the
% masker), and in order to estimate neural activation patterns along the
% length of the cochlea, it utilizes every combination of masker and probe
% available in the subject's MAP. The software generates .csv files that
% can be used to gather the data required for PECAP as inputs to Custom
% Sound EP Software. It is written in MATLAB and should run on any platform
% supporting that language, although we have only tested it on Windows 10
% using MATLAB R2018a.
%
% If you use the software, please cite the following publication, which
% contains details of how the PECAP method words:
% Cosentino, S., Gaudrain, E., Deeks, J.M. & Carlyon, R.P. (2015) 
% Multistage nonlinear optimization to recover neural activation patterns 
% from evoked compound action potentials of cochlear implant users. 
% Biomedical Engineering, IEEE Transactions on, 63, 833-840. 
%
% We hope you will find the software useful, but neither the authors nor
% their employers accept any responsibility for the consequencies of it's
% use. The USER is responsible for ensuring the saftety of any subjects
% during testing.
%
% disclaimer:
% The Software is the result work conducted within the MRC Cognition & Brain Sciences unit at the University of Cambridge (the “University”)
% "The Software shall be used for non-commercial research purposes only.  The USER will not reverse engineer, manufacture, sell or sublicense for manufacture and sale upon a commercial basis the Software, incorporate it into other software or products or use the Software other than herein expressly specified.
% The USER agrees that it will not translate, alter, decompile, disassemble, reverse engineer, reverse compile, attempt to derive, or reproduce source code from the Software. The USER also agrees that it will not remove any copyright or other proprietary or product identification notices from the Software.
% The Software is provided without warranty of merchantability or fitness for a particular purpose or any other warranty, express or implied, and without any representation or warranty that the use or supply of the Software will not infringe any patent, copyright, trademark or other right.
% In no event shall the University or their staff and students be liable for any use by the USER of the Software.
% The supply of the Software does not grant the USER any title or interest in and to the Software, the related source code, and any and all associated patents, copyrights, and other intellectual property rights.
% The University and their staff and students reserve the right to revise the Software and to make changes therein from time to time without obligation to notify any person or organisation of such revision or changes.
% While the University will make every effort to ensure the accuracy of Software and the data contained therein, however neither the University nor its employees or agents may be held responsible for errors, omissions or other inaccuracies or any consequences thereof."
%--------------------------------------------------------------------------

% pecap_daq(strOutputFilename, Rec_Elec, Diag, Gain_dB, Delay_us,
%       PRate_Hz, sweeps, elecs, e_levels)
%
%   Produces a csv file in the format that is accepted by Custom Sound EP.
%   It also produces a verification csv file that allows the researcher to
%   confirm all electrodes are within comfortable loudness ranges for the
%   subject prior to starting the entire PECAP data collection. The
%   recording elctrode can be -2 or +2 from the probe electrode, but is
%   adjusted by +/- 1 when necessary so that it's never the same as either 
%   the probe or the masker electrode.
%
%   Inputs are as follows:
%       - strOutputFilename: Name of the csv file the script will create
%                            i.e. 'Subject1' will create two files:
%                            'Subject1_compliance.csv' and 'Subject1.csv'
%       - Rec_Elec: number that that records the default offset of the
%                   recording electrode from the probe, +2 (basal) or -2 
%                   (apical)
%       - Diag: turns on/off multiple recording electrodes along diagonal
%               'on' indicates that ECAPs will be recorded both apically
%               and basally along the diagonal (i.e. when probe electrode
%               = masker electrode)
%               'off' indicates that ECAPs will only be recorded from the
%               Rec_Elec direction along the diagonal
%       - Gain_dB: sets the gain of the amplifier in dB
%                  Possible values: 40, 50, 60, 70
%                  Recommended: set this value after running the
%                  optimization routine for maximum ECAP Amplitude
%       - Delay_us: delay after the probe before recording in us
%                   Possible values: 72, 98, 122
%                   Recommended: set this value after running the
%                   optimization routine for maximum ECAP Amplitude
%       - PRate_Hz: probe rate in Hz
%                   Recommended value: 80
%       - sweeps: the number of sweeps to be done for each ECAP
%                 Recommended value: 50
%       - elecs: vector that defines the electrodes from which to stimulate 
%                i.e. [1:22] is all electrodes, and [4:7 10:20] is
%                electrodes 4-7 and 10-20.
%       - e_levels: vector that indicates the electrode current levels in
%                   CU between 0 and 255. The standard is to set these to 
%                   MCL
%
%   Note: This routine will only use electrodes specified in the elecs
%   vector with the corresponding e_levels to stimulate for both the probe
%   and masker conditions, but will use all available electrodes [1:22] in
%   range to record from.

% Example:
% pecap_daq('Subject1', +2, 'on', 40, 98, 80, 50, [4:7 10:20], [204 205 206 207 210 211 212 213 214 215 216 217 218 219 220]);

%--------------------------------------------------------
% Developed: John Deeks
% Ammended: Charlotte Garcia, April 2019
% MRC Cognition and Brain Sciences Unit, Cambridge, UK
%--------------------------------------------------------

%% Set defaults
strProgramName = 'pecap_daq';

ELECTRODE_MIN = 1;          ELECTRODE_MAX = 22;
LEVEL_CU_MIN = 0;           LEVEL_CU_MAX = 255;
PHASE_DUR_US_MIN = 25;      PHASE_DUR_US_MAX = 150;
PHASEGAP_US_MIN = 7;        PHASEGAP_US_MAX = 58;
MIN_DELAY_US = 44;          MAX_DELAY_US = 1607;
GAINS = ([40 50 60 70]);    REC_EL =([-2 2]);  

%% User might want to change these parameters
stim_mode = 'MP1';          rec_mode = 'MP2';
phase_duration_us = 25;     interphase_gap_us = 8;
n_samples = 32;             mpi_us = 400; % must be in range 13 - 200000 us
probe_rate_Hz = PRate_Hz;   n_sweeps = sweeps;

no_electrodes = length(elecs);
electrodes = elecs;
min_spec_electrode = min(electrodes);
max_spec_electrode = max(electrodes);

%% Parameter checks and warnings 
if length(e_levels) ~= no_electrodes
    fprintf(1, '%s: mismatch of no_electrodes and levels - exiting\n', strProgramName);
    return;
end
if (min(electrodes) < ELECTRODE_MIN || max(electrodes) > ELECTRODE_MAX)
    fprintf(1, '%s: Electrode must be integer in range %d to %d - exiting.\n', strProgramName, ELECTRODE_MIN, ELECTRODE_MAX);
    return
end
if (phase_duration_us < PHASE_DUR_US_MIN || phase_duration_us > PHASE_DUR_US_MAX)
    fprintf(1, '%s: Phase duration must be in range %.2f µs to %.2f µs - exiting.\n', strProgramName, PHASE_DUR_US_MIN, PHASE_DUR_US_MAX);
    return
end
if (min(e_levels) < LEVEL_CU_MIN || max(e_levels) > LEVEL_CU_MAX)
    fprintf(1, '%s: Level must be integers in range %d to %d - exiting.\n', strProgramName, LEVEL_CU_MIN, LEVEL_CU_MAX);
    return
end
if length(Diag) == 2
    if min(Diag == 'on') == 1
        diag_rec = 1;
    else
        fprintf(1, '%s: Diag input must either be specified as on or off - exiting.\n', strProgramName);
        return
    end
elseif length(Diag) == 3
    if min(Diag == 'off') == 1
        diag_rec = 0;
    else
        fprintf(1, '%s: Diag input must either be specified as on or off - exiting.\n', strProgramName);
        return
    end
else
    fprintf(1, '%s: Diag input must either be specified as on or off - exiting.\n', strProgramName);
    return
end

% Cochlear charge safety check
max_current_uA = 17.5*(100^(max(e_levels)/255));
max_charge_nC = phase_duration_us*max_current_uA*10^-3;
if max_charge_nC > 212
    fprintf(1, '%s: Charge must not exceed safety limit of 212 nC per phase - exiting\n', strProgramName);
    return;
end
if (interphase_gap_us < PHASEGAP_US_MIN || interphase_gap_us > PHASEGAP_US_MAX)
    fprintf(1, '%s: Interphase gap duration must be in range %.2f µs to %.2f µs - exiting.\n', strProgramName, PHASEGAP_US_MIN, PHASEGAP_US_MAX);
    return
end
if ~( (Rec_Elec == REC_EL(1)) || (Rec_Elec == REC_EL(2)) )
    fprintf(1, '%s: Rec_Elec must be %d or %d re: Probe Electrode\n', strProgramName, REC_EL(1), REC_EL(2));
    return;
end
if ( ((max(electrodes)+Rec_Elec) > ELECTRODE_MAX) || ((max(electrodes)+Rec_Elec) > max_spec_electrode) )
    fprintf(1, '%s: Recording electrode will be reset to -2 where PROBE electrode > %d\n', strProgramName, min(ELECTRODE_MAX, max_spec_electrode)-Rec_Elec);    
end
if ( ((min(electrodes)+Rec_Elec) < ELECTRODE_MIN) || ((min(electrodes)+Rec_Elec) < min_spec_electrode) )
    fprintf(1, '%s: Recording electrode will be reset to +2 where PROBE electrode < %d\n', strProgramName, max(ELECTRODE_MIN, min_spec_electrode)-Rec_Elec);
end
if ~( (Gain_dB == GAINS(1)) || (Gain_dB == GAINS(2)) || (Gain_dB == GAINS(3)) || (Gain_dB == GAINS(4)) )
    fprintf(1, '%s: Gain_dB must be 40, 50, 60 or 70 dB\n', strProgramName);
    return;
end
if (Delay_us < MIN_DELAY_US || Delay_us > MAX_DELAY_US)
    fprintf(1, '%s: Recording Delay must be in range %.2f µs to %.2f µs - exiting.\n', strProgramName, MIN_DELAY_US, MAX_DELAY_US);
    return
end

%% Summarise
fprintf('\n-----\nNOTE: PLEASE READ THE DISCLAIMER IN THE HEADING OF THIS PROGRAM PRIOR TO USE.\n-----\n\n')
fprintf(1, 'Generating sequence for %d electrodes between %d and %d:\nEl.\tLevel CUs\n', no_electrodes, elecs(1), elecs(end));
for i=1:no_electrodes
    fprintf(1, '%d\t%d\n', electrodes(i), e_levels(i));
end
fprintf(1, 'Stimulating Probe and Masker in %s mode\n', stim_mode); 
fprintf(1, 'Phase duration = %2.0f us\t\tInterphase gap = %2.0f us\n', phase_duration_us, interphase_gap_us);
fprintf(1, 'Probe rate = %2.0f Hz\t\t\tMasker-Probe Interval = %2.0f us\n', probe_rate_Hz, mpi_us);
fprintf(1, 'Gain is %d dB\t\t\t\tDelay is %d us\n', Gain_dB, Delay_us);
fprintf(1, 'Number of sweeps = %d\t\tNumber of samples = %d\n', n_sweeps, n_samples);
fprintf(1, 'Recording Active Electrode is %d relative to PROBE electrode, %s mode\n\n', Rec_Elec, rec_mode);

% check to confirm that the filename doesn't already exist, open file
fid = fopen([strOutputFilename '.csv'],'r');
if fid == 3
    fprintf(1, '%s: runfile already exists - aborting\n', strProgramName);
    fclose(fid);
    return;
end
fid = fopen([strOutputFilename '.csv'],'w');

%% Write the PECAP .csv file

% Note: only NRT4.0 works, e.g. NRT4.3 'not supported'
fprintf(fid, '[NRT4.0]\n');

% Generates two vector arrays, one "t" that contains the headers for the 
% .csv file, and one "a" that sets the defaults for the content of each of
% the columns. Some value(s) in "a" will be changed each time a row is
% written to the .csv file
t{1}='NRT Number'; a{1}='1';
t{2}='Probe Active Electrode'; a{2}='11';
t{3}='Probe Indifferent Electrode'; a{3}=stim_mode;
t{4}='Probe Current Level'; a{4}='20';
t{5}='Probe Pulse Width'; a{5}=num2str(phase_duration_us);
t{6}='Probe Inter Phase Gap'; a{6}=num2str(interphase_gap_us);
t{7}='Probe Rate'; a{7}=num2str(probe_rate_Hz);
t{8}='Recording Active Electrode'; a{8}='13';
t{9}='Recording Indifferent Electrode'; a{9}=rec_mode;
t{10}='Gain'; a{10}=([num2str(Gain_dB) ' dB']);
t{11}='Delay'; a{11}=num2str(Delay_us);
t{12}='Artefact Cancellation Technique'; a{12}='Forward Masking';
t{13}='Artefact Reduction'; a{13}='Off';
t{14}='Nr of Sweeps'; a{14}=num2str(n_sweeps);
t{15}='Nr of Samples'; a{15}=num2str(n_samples);
t{16}='Resolution'; a{16}='x1';
t{17}='Artefact Rejection Level'; a{17}='Off';
t{18}='Discard first N Recordings'; a{18}='0';
t{19}='Masker Active Electrode'; a{19}='11';
t{20}='Masker Indifferent Electrode'; a{20}=stim_mode;
t{21}='Masker Current Level'; a{21}='20';
t{22}='Masker Pulse Width'; a{22}=num2str(phase_duration_us);
t{23}='Masker Inter Phase Gap'; a{23}=num2str(interphase_gap_us);
t{24}='Nr of Maskers'; a{24}='1';
t{25}='Masker Rate'; a{25}='100';
t{26}='Masker Probe Interval'; a{26}=num2str(mpi_us);
t{27}='Reference MPI'; a{27}='400';
t{28}='Lock Masker Polarity'; a{28}='Off';
t{29}='Artefact Template Current Level'; a{29}='1';
t{30}='Scaling Factor'; a{30}='1';
t{31}='Nr of Sweeps for Artef Templ'; a{31}='500';
t{32}='Artef Red Current Level'; a{32}='1';
t{33}='Artef Red Pulse Width'; a{33}='10';
t{34}='Nr of Averages for MFD'; a{34}='10';
t{35}='Step Ratio'; a{35}='1000';
t{36}='Targ Offset'; a{36}='0';
t{37}='Tolerance'; a{37}='Normal';
t{38}='Time Out'; a{38}='100';
t{39}='OK Measurements'; a{39}='4';
t{40}='Error Determination'; a{40}='Last minus First';

% Writes t vector array to the second row of the csv file, 40 fields
formatSpec = [repmat('%s,',1,39) '%s\n'];
fprintf(fid,formatSpec,t{:});

% Initiate NRT and probe electrode counter
NRT=1; p=1;

% Start while loop for each NRT probe electrode segment
while p <= no_electrodes
    a{2}=num2str(electrodes(p));                % Probe Active Electrode
    a{4}=num2str(e_levels(p));                  % Probe Current Level    
    
    % Checks to make sure that the set recording electrode is not outside
    % available electrodes - if the default is, then it sets the recording
    % electrode to the other side
    if (electrodes(p) + Rec_Elec) > ELECTRODE_MAX
        rec_el = -2;
    elseif (electrodes(p) + Rec_Elec) < ELECTRODE_MIN
        rec_el = 2;
    else
        rec_el = Rec_Elec;
    end
    
    % Initiate masker electrode counter 
    n=1;
    
    % Start while loop for each NRT segment within each probe segment
    while n <= no_electrodes
        a{19}=num2str(electrodes(n));           % Masker Active Electrode
        a{21}=num2str(e_levels(n));             % Masker Current Level
        
        % If the masker and probe are on the same electrode (diagonal)
        if a{2}==a{19}
            % If we are recording from both sides along the diagonal
            if diag_rec == 1
                % record in the direction specified (within e1-22)
                if (electrodes(p) + rec_el) >= ELECTRODE_MIN && (electrodes(p) + rec_el) <= ELECTRODE_MAX
                    a{8}=num2str(electrodes(p) + rec_el);
                    a{1}=num2str(NRT);
                    fprintf(fid,formatSpec,a{:});
                    NRT = NRT + 1;
                end
                % and also in the other direction (within e1-22)
                if (electrodes(p) - rec_el) >= ELECTRODE_MIN && (electrodes(p) - rec_el) <= ELECTRODE_MAX
                    a{8}=num2str(electrodes(p) - rec_el);
                    a{1}=num2str(NRT);
                    fprintf(fid,formatSpec,a{:});
                    NRT = NRT + 1;
                end
            % If we are only recording along the diagonal once
            else
                % record in the direction specified (within e1-22)
                if (electrodes(p) + rec_el) >= ELECTRODE_MIN && (electrodes(p) + rec_el) <= ELECTRODE_MAX
                    a{8}=num2str(electrodes(p) + rec_el);
                    a{1}=num2str(NRT);
                    fprintf(fid,formatSpec,a{:});
                    NRT = NRT + 1;
                end
            end
        % If we are not on the diagonal    
        else
            a{1}=num2str(NRT);
            % check to make sure that recording and masker electrode
            % are not the same - if so, move the recording electrode
            % between the probe and the masker
            if(electrodes(n) == (electrodes(p)+rec_el))
                a{8}=num2str((electrodes(p)+electrodes(n))/2);
            else
                a{8}=num2str(electrodes(p)+rec_el);            
            end              
            fprintf(fid,formatSpec,a{:});
            NRT = NRT + 1;
        end
        n = n + 1;
    end
    p = p + 1;
end

%  Close first output file
fclose(fid);

%% Write the Loudness Verifciation Check .csv file

% check to confirm that the filename doesn't already exist, open file
fid = fopen([strOutputFilename '_verification.csv'],'r');
if fid == 3
    fprintf(1, '%s: runfile already exists - aborting\n', strProgramName);
    fclose(fid);
    return;
end
fid = fopen([strOutputFilename '_verification.csv'],'w');

% Note: only NRT4.0 works, e.g. NRT4.3 'not supported'
fprintf(fid, '[NRT4.0]\n');

% Writes t vector array to the second row of the csv file, 40 fields
formatSpec = [repmat('%s,',1,39) '%s\n'];
fprintf(fid,formatSpec,t{:});

% Set Nr of Sweeps to 50 so subject can hear MCL levels for each electrode
% Reset NRT & probe counters
a{14} = '50'; NRT = 1; p = 1;

% Write diagonal to the csv file
while p <= no_electrodes
    % set probe and level
    a{2}=num2str(electrodes(p));                % Probe Active Electrode
    a{4}=num2str(e_levels(p));                  % Probe Current Level    
    
    % set masker and level equal to probe
    a{19}=num2str(electrodes(p));               % Masker Active Electrode
    a{21}=num2str(e_levels(p));                 % Masker Current Level
        
    % Checks to make sure that the set recording electrode is not outside
    % available electrodes - if the default is, then it sets the recording
    % electrode to the other side
    if (electrodes(p) + Rec_Elec) > ELECTRODE_MAX
        rec_el = -2;
    elseif (electrodes(p) + Rec_Elec) < ELECTRODE_MIN
        rec_el = 2;
    else
        rec_el = Rec_Elec;
    end
    
    % set recording electrode and NRT, print row to csv file, iterate NRT
    a{8}=num2str(electrodes(p) + rec_el);
    a{1}=num2str(NRT);
    fprintf(fid,formatSpec,a{:});
    NRT = NRT + 1;

    % iterate probe counter
    p = p + 1;
end

%  Close loudess verification check output file
fclose(fid);