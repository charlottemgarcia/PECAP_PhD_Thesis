#!/usr/bin/env python

# This script is for collecting PECAP using NIC2 for Cochlear subjects
# It can be used as an alternative to Custom Sound EP, and allows the flexibility to record multiple repetitions of the
#   same parameters (up to 4 at a Probe Rate of 80Hz, or 5 at 350Hz) with minimal additional time. Recording 12 sweeps
#   4 times in Custom Sound EP requires approximately 4 seconds per NRT (i.e. 4 *4 = 16), whereas in NIC2 this takes
#   approximately 5 seconds.
# While this script can be run in the Python Shell, it is suggested that it is run in PyCharm 2016.1.5. Note that NIC2
#   is only compatible with python 2.4.4. It is also suggested that the script be copied and pasted in to the python
#   console and run from there, as the program functions to store all ECAP information internally in a dictionary called
#   NRTs which it only writes to a csv file after the entirety of PECAP has been run. Should the PECAP sequence be
#   interrupted or a hardware failure occur for any reason, running on the console will allow the researcher to save the
#   data acquired *so far* to a csv file after the fact.
# Method for usage:
#       1. Enter parameters for current subject in the 'PECAP Parameters' section below
#       2. Define the three functions required for the script by copying and pasting them into the console
#       3. Define the parameters by copying and pasting them into the console
#       4. Generate the NRT Order by copying and pasting the 'Generate NRT Order' section of the script into the
#          console and confirm that all parameters are as desired before moving onto step 5. Check to ensure that you
#          have commented out the correct line here. It is recommended to do the diagonal first to ensure the subject
#          is comfortable with the levels for all electrodes for the longer listening time that Full PECAP takes, as
#          well as check to confirm that the ECAPs being recorded are visible with the recording parameters selected.
#       5. Run 'PECAP Acquisition' section of the code by copying and pasting into the console
#          Note: This will tell you the sweep number, the Probe/Masker/Record settings, and the elapsed time for each
#                NRT condition as it's executing. You can calculate total time to go from this information
#          Note: Running this script in the console allows for the ability to extract partially completed PECAP data to
#                a csv file after the fact, should testing be interrupted. Simply execute the gen_csv function in the
#                last line of this script in the console, but replace the 'nrtno' input with len(NRTs).
#          Note: Should testing be interrupted, you can also simply restart the NIC2 script loop from
#                i in range(len(NRTs)+1,nrtno+1) to continue adding to the NRTs dictionary where right after the last
#                recorded value.

########################################
# Charlotte Garcia, March 2019
# MRC - Cogniton & Brain Sciences Unit
########################################


## Define PECAP Functions
# function for getting one ECAP (with reps) from NIC2 using forward masking
def NIC2_ECAP_FM(probe_elec, probe_level, masker_elec, masker_level, rec_elec, params):

    # probe_elec takes an integer input 1-22
    # probe_level takes an integer input 1-255 (Current Units)
    # masker_elec takes an integer input 1-22
    # masker_level takes an integer input 1-255 (Current Units)
    # rec_elec takes an integer input 1-22
    # params takes a dictionary with a number of fields that can be changed:
    #       'gain' = 40.0, 50.0, 60.0, or 70.0 dB
    #       'delay' = 44.0, 73.0, 98.0, 122.0 us
    #       'ipg' = i.e. 7.0 or 8.0 us
    #       'reps' = 1 to 5, indicates the number of recordings to extract from data
    #       'tsweeps' = total number of sweeps to approximate, i.e. 50

    import datetime, time
    import nucleus.nre.nic as nic
    import nucleus.nre.stream.nic as nics

    time_begin = datetime.datetime.fromtimestamp(time.time())

    sequence = nic.Sequence(1)

    # Create the power sub-sequence 1000 null stimuli at 5 kHz
    powerupSequence = nic.Sequence(1000)
    nullStimulus = nic.NullStimulus(200.0)
    nullCommand = nic.StimulusCommand(nullStimulus)
    powerupSequence.append(nullCommand)
    sequence.append(powerupSequence)

    # channel defines the recording electrode (ICE4 = e4), and the second input is ground (ECE2 = MP2)
    rec_channel = nic.Channel(rec_elec, nic.Channel.ECE2)
    nrt_gain = params['gain']		    # 50.0 dB
    num_samples = 32	                # number of samples
    s_freq = 20000.0	                # Sampling always at 20 kHz for CIC4
    delay = params['delay']		        # 122.0 us
    pw = params['PD']                   # Pulse Width aka Phase Duration (generally set to 25us)
    ipg = params['ipg']                 # Interphase gap set to 8us
    reps = params['reps']               # CMG: this is the number of repetitions (recorded ABCD traces)
    sweeps = params['tsweeps']/reps     # CMG: this is the number of sweeps that it averages before recording things

    print 'Probe:', probe_elec, 'Level:', probe_level, 'Masker:', masker_elec, 'Level:', masker_level, 'Record:', rec_elec

    # masker fp = 400 us
    # probe fp = 2000 us
    # together with pufs, Probe Period (@ 80 Hz) should take up 12500 us:  12500 - 400 - 2000 = 10100 us
    # so for instance 20 pufs with fp = 505 us (10100 = 505*20)
    # note: limit of the SP12 processor buffer means that at 80 Hz (with 20 pufs at 505us each), the maximum number of
    #       repeats we can record (defined as reps above) is 4. However, if we increase the Probe Rate, the Probe Period
    #       (@ 350 Hz) should take up 2857 us: 2857 - 400 - 2000 = 457 us
    if params['reps'] == 5:
        nrt_puf_seq = nic.Sequence(1)
        nrt_puf_seq.append(nic.StimulusCommand(nic.NullStimulus(457)))
    else:
        nrt_puf_seq = nic.Sequence(20)
        nrt_puf_seq.append(nic.StimulusCommand(nic.NullStimulus(505)))

    # the below defines the active electrode as e2, the ground reference as MP1, etc)
    # nic.BiphasicStimulus(stimulating_channel,ground_channel,CU_level,pulse_width,interphase_gap,duration?)
    stim_masker0  = nic.BiphasicStimulus(masker_elec, nic.Channel.ECE1, 0,            pw, ipg,  120.0)		         # 120 = 400 - 280 (280 taken by 4 * 70 us)
    stim_probe    = nic.BiphasicStimulus(probe_elec,  nic.Channel.ECE1, probe_level,  pw, ipg, 2000.0)
    stim_masker   = nic.BiphasicStimulus(masker_elec, nic.Channel.ECE1, masker_level, pw, ipg,  120.0)               # changing masker to equal probe level (instead of probe+10)
    stim_probe0   = nic.BiphasicStimulus(probe_elec,  nic.Channel.ECE1, 0,            pw, ipg, 2000.0)

    measurementSeq = nic.Sequence(sweeps)

    # will repeat the ABCD sequence in measurementSeq for the number of sweeps desired
    for i in range(0,reps):
        # define measurement on A
        measurementSeq.append(nrt_puf_seq)		# lasts 10100 us
        measurementSeq.append(nic.StimulusCommand(stim_masker0))			# lasts 120.0 us + 4 * 70.0 = 400 us
        # nic.NRTCommand(biphasicstimulus, recording channel, gain, num_samples to record (generally 32), sample frequency (always 20kz), delay)
        measurementSeq.append(nic.NRTCommand(stim_probe, rec_channel, nrt_gain, num_samples, s_freq, delay))	# lasts 2000 us
        # total: 12500us = 12.5 ms = 1/80 Hz

        # define measurement on B
        measurementSeq.append(nrt_puf_seq)		# lasts 10100 us
        measurementSeq.append(nic.StimulusCommand(stim_masker))		        # lasts 120.0 us + 4 * 70.0 = 400 us
        measurementSeq.append(nic.NRTCommand(stim_probe, rec_channel, nrt_gain, num_samples, s_freq, delay))	# lasts 2000 us
        # total: 12500us = 12.5 ms = 1/80 Hz

        # define measurement on C
        measurementSeq.append(nrt_puf_seq)		# lasts 10100 us
        measurementSeq.append(nic.StimulusCommand(stim_masker))		        # lasts 120.0 us + 4 * 70.0 = 400 us
        measurementSeq.append(nic.NRTCommand(stim_probe0, rec_channel, nrt_gain, num_samples, s_freq, delay))	# lasts 2000 us
        # total: 12500us = 12.5 ms = 1/80 Hz

        # define measurement on D
        measurementSeq.append(nrt_puf_seq)		# lasts 10100 us
        measurementSeq.append(nic.StimulusCommand(stim_masker0))			# lasts 120.0 us + 4 * 70.0 = 400 us
        measurementSeq.append(nic.NRTCommand(stim_probe0, rec_channel, nrt_gain, num_samples, s_freq, delay))	# lasts 2000 us
        # total: 12500us = 12.5 ms = 1/80 Hz

    sequence.append(measurementSeq)

    # Trailing frames are required otherwise the last measurement it not performed.
    trailingSequence = nic.Sequence(100)
    trailingSequence.append(nic.StimulusCommand(nic.NullStimulus(200.0)))
    sequence.append(trailingSequence)

    # Create the client object.
    client = nics.NICStreamClient("sp12-cic4-1")

    # Specify the sequence to the client.
    client.sendData(sequence)

    # Start streaming.
    client.startType = nic.trigger.IMMEDIATE
    client.startStream()

    # Wait until the streaming has finished, by checking the stream status.
    #   Status progression is: Stopped -> Streaming -> Idle
    while (nics.status.IDLE != client.streamStatus()):
        time.sleep(0.2)

    # The system needs to be stopped, even though it's finished streaming.
    client.stopStream()

    # Retrieve the telemetry data and print it to screen.
    data = client.receiveData()

    # create dictionary in which to story ECAP reps
    ECAP = {}
    for i in range(0,reps):
        A = []
        B = []
        C = []
        D = []
        ABCD = []
        for j in range (0,num_samples):
            A.append(data.at(i*4+0).samples[j])
            B.append(data.at(i*4+1).samples[j])
            C.append(data.at(i*4+2).samples[j])
            D.append(data.at(i*4+3).samples[j])
            ABCD.append(data.at(i*4+0).samples[j]-(data.at(i*4+1).samples[j]-data.at(i*4+2).samples[j])-data.at(i*4+3).samples[j])
        ECAP[i+1] = {'Probe':probe_elec,
                   'Masker':masker_elec,
                   'Recording': rec_elec,
                   'A':A, 'B':B, 'C':C, 'D':D, 'ABCD':ABCD}

    # Cease communications with the hardware and clean up the memory.
    del client

    time_end = datetime.datetime.fromtimestamp(time.time())
    print 'Elapsed Time:', str(time_end - time_begin)

    return ECAP
    # end function
# function for generating NRT order for the subject as per PECAP parameters
def generate_nrt_order_fm(elecs, e_levels, record_from, super_elecs, diagrecelec, params):

    # set defaults
    min_elec = 1
    max_elec = 22
    level_cu_min = 0
    level_cu_max = 255
    phase_dur_us_min = 25.0
    phase_dur_us_max = 180.0
    phasegap_us_min = 7.0
    phasegap_us_max = 58.0
    min_delay_us = 44.0
    max_delay_us = 1607.0
    gains = [40.0, 50.0, 60.0, 70.0]
    rec_elec_options = [-3, -2, 2, 3]

    # parameter checks and warnings
    if (elecs[1] - elecs[0] + 1) != len(e_levels):
        print 'Mismatch of number of electrodes and electrode levels - exiting.'
        return
    if min(elecs) < min_elec or max(elecs) > max_elec:
        print 'Electrodes must be within range %d to %d - exiting.' % (min_elec,max_elec)
        return
    if params['PD'] < phase_dur_us_min or params['PD'] > phase_dur_us_max:
        print 'Phase duration must be within range %d to %d - existing.' % (phase_dur_us_min,phase_dur_us_max)
        return
    if min(e_levels) < level_cu_min or max(e_levels) > level_cu_max:
        print 'Level must be integers in range %d to %d - exiting.' % (level_cu_min, level_cu_max)
        return
    if params['reps'] > 5:
        print 'The SP12 Processor does not have space in the buffer to run more than 5 reps - existing.'
        return
    elif params['reps'] < 0:
        print 'Reps must be a positive number between 1 and 5 - exiting'
        return
    if len(super_elecs) > 0:
        if min(super_elecs) < min(elecs) or max(super_elecs) > max(elecs):
            print 'Super Elecs must be within the range of used electrodes, %d to %d - exiting.' % (min(elecs), max(elecs))
            return

    # cochlear charge safety check
    max_current_uA = 17.5*pow(100.0,(max(e_levels)/255.0))
    max_charge_nC = params['PD']*max_current_uA*pow(10.0,-3.0)
    if max_charge_nC > 212:
        print 'Charge must not exceed safety limit of 212 nC per phase - exiting.'
        return
    if params['ipg'] < phasegap_us_min or params['ipg'] > phasegap_us_max:
        print 'Interphase gap duration must be in range %d to %d - exiting.' % (phasegap_us_min, phasegap_us_max)
        return
    if record_from not in rec_elec_options:
        print 'Exiting. Recording electrode must be one of the following:', rec_elec_options
        return
    if params['gain'] not in gains:
        print 'Gain must be 40.0, 50.0, 60.0 or 70.0 dB - exiting.'
        return
    if params['delay'] < min_delay_us or params['delay'] > max_delay_us:
        print 'Recording Delay must be in range of %d to %d us - exiting.' % (min_delay_us, max_delay_us)
        return

    # initiate nrtno counter and NRT dictionary with all masker/probe/rec combinations
    nrtno = 0
    NRT_run_order = {1: {'Probe Electrode':     0,
                         'Probe Level':         0,
                         'Masker Electrode':    0,
                         'Masker Level':        0,
                         'Recording Electrode': 0}}

    for p in range(elecs[0], elecs[1] + 1):
        for m in range(elecs[0], elecs[1] + 1):
            # if statements and loops for along the diagonal
            if p == m:
                # for loop for every recording electrode along the diagonal, super electrodes
                if p in super_elecs:
                    for i in range(min_elec, max_elec + 1):
                        if p != i:
                            nrtno = nrtno + 1
                            r = i
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p-min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m-min(elecs)],
                                                    'Recording Electrode': r}
                # for loop for multiple recording electrodes along the diagonal, non-super electrodes
                else:
                    for i in range(0, len(diagrecelec)):
                        if (p + diagrecelec[i]) >= min_elec and (p + diagrecelec[i]) < max_elec + 1:
                            nrtno = nrtno + 1
                            r = p + diagrecelec[i]
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p-min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m-min(elecs)],
                                                    'Recording Electrode': r}
            # if statements and loops for the rest of the matrix
            else:
                # if probe and masker are on super electrodes
                if (p in super_elecs) and (m in super_elecs):
                    # loop through all possible recording electrodes
                    for ridx in range(min_elec, max_elec + 1):
                        # if the recording index isn't the probe or the masker
                        if (ridx != p) and (ridx != m):
                            nrtno = nrtno + 1
                            r = ridx
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                # if we are recording apically (+2 or +3)
                elif record_from > 0:
                    # if probe is on electrodes 1:20 (r=+2) or 1:19 (r=+3)
                    if p <= (max_elec - record_from):
                        if (p + record_from) != m:
                            nrtno = nrtno + 1
                            r = p + record_from
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                        else:
                            nrtno = nrtno + 1
                            r = p + record_from - 1
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                    # if probe is on electrodes 21:22 (r=+2) or 20:22 (r=+3)
                    else:
                        if (p - record_from) != m:
                            nrtno = nrtno + 1
                            r = p - record_from
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                        else:
                            nrtno = nrtno + 1
                            r = p - record_from + 1
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                # if we are recording basally (-2 or -3)
                elif record_from < 0:
                    # if probe is on electrodes 3:22 (r=-2) or 4:22 (r=-3)
                    if p >= (min_elec - record_from):
                        if (p + record_from) != m:
                            nrtno = nrtno + 1
                            r = p + record_from
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                        else:
                            nrtno = nrtno + 1
                            r = p + record_from + 1
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                    # if probe is on electrodes 1:2 (r=-2) or 1:3 (r=-3)
                    else:
                        if (p - record_from) != m:
                            nrtno = nrtno + 1
                            r = p - record_from
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}
                        else:
                            nrtno = nrtno + 1
                            r = p - record_from - 1
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}


    # Summarize PECAP run order
    print 'Generating PECAP sequence for %d electrodes %d to %d:' % (elecs[1] - elecs[0] + 1, elecs[0], elecs[1])
    print 'El. \tLevel CUs'
    for i in range(0, elecs[1] - elecs[0] + 1):
        print str(i + min(elecs)), '\t\t', e_levels[i]
    print 'Stimulating Probe and Masker in MP1 Mode.'
    print 'Recording Electrode is %d relative to Probe, MP2 Mode.' % record_from
    print 'Phase Duration = %d us\tInterphase gap = %d us' % (params['PD'], params['ipg'])
    if params['reps'] == 5:
        print 'Probe Rate = 350 Hz\t\tMasker-Probe Interval = 400 us'
    else:
        print 'Probe Rate = 80 Hz\t\tMasker-Probe Interval = 400 us'
    print 'Gain is %d dB\t\t\tDelay is %d us' % (params['gain'],params['delay'])
    print 'Number of Reps = %d\t\tTotal number of Sweeps = %d x %d = %d' % (params['reps'], params['reps'],
                                                                            params['tsweeps'] / params['reps'],
                                                                            params['reps'] * (
                                                                            params['tsweeps'] / params['reps']))
    print 'Number of Samples = 32\tNumber of NRTs = %d\n' %nrtno

    return NRT_run_order
# function for generating NRT order for the DIAGNOAL as per PECAP parameters
def generate_nrt_order_fm_diag(elecs, e_levels, record_from, super_elecs, diagrecelec, params):

    # set defaults
    min_elec = 1
    max_elec = 22
    level_cu_min = 0
    level_cu_max = 255
    phase_dur_us_min = 25.0
    phase_dur_us_max = 180.0
    phasegap_us_min = 7.0
    phasegap_us_max = 58.0
    min_delay_us = 44.0
    max_delay_us = 1607.0
    gains = [40.0, 50.0, 60.0, 70.0]
    rec_elec_options = [-3, -2, 2, 3]

    # parameter checks and warnings
    if (elecs[1] - elecs[0] + 1) != len(e_levels):
        print 'Mismatch of number of electrodes and electrode levels - exiting.'
        return
    if min(elecs) < min_elec or max(elecs) > max_elec:
        print 'Electrodes must be within range %d to %d - exiting.' % (min_elec,max_elec)
        return
    if params['PD'] < phase_dur_us_min or params['PD'] > phase_dur_us_max:
        print 'Phase duration must be within range %d to %d - existing.' % (phase_dur_us_min,phase_dur_us_max)
        return
    if min(e_levels) < level_cu_min or max(e_levels) > level_cu_max:
        print 'Level must be integers in range %d to %d - exiting.' % (level_cu_min, level_cu_max)
        return
    if params['reps'] > 5:
        print 'The SP12 Processor does not have space in the buffer to run more than 5 reps - existing.'
        return
    elif params['reps'] < 0:
        print 'Reps must be a positive number between 1 and 5 - exiting'
        return
    if len(super_elecs) > 0:
        if min(super_elecs) < min(elecs) or max(super_elecs) > max(elecs):
            print 'Super Elecs must be within the range of used electrodes, %d to %d - exiting.' % (min(elecs), max(elecs))
            return

    # cochlear charge safety check
    max_current_uA = 17.5*pow(100.0,(max(e_levels)/255.0))
    max_charge_nC = params['PD']*max_current_uA*pow(10.0,-3.0)
    if max_charge_nC > 212:
        print 'Charge must not exceed safety limit of 212 nC per phase - exiting.'
        return
    if params['ipg'] < phasegap_us_min or params['ipg'] > phasegap_us_max:
        print 'Interphase gap duration byst be in range %d to %d - exiting.' % (phasegap_us_min, phasegap_us_max)
        return
    if record_from not in rec_elec_options:
        print 'Exiting. Recording electrode must be one of the following:', rec_elec_options
        return
    if params['gain'] not in gains:
        print 'Gain must be 40.0, 50.0, 60.0 or 70.0 dB - exiting.'
        return
    if params['delay'] < min_delay_us or params['delay'] > max_delay_us:
        print 'Recording Delay must be in range of %d to %d us - exiting.' % (min_delay_us, max_delay_us)
        return

    # initiate nrtno counter and NRT dictionary with all masker/probe/rec combinations
    nrtno = 0
    NRT_run_order = {1: {'Probe Electrode':     0,
                         'Probe Level':         0,
                         'Masker Electrode':    0,
                         'Masker Level':        0,
                         'Recording Electrode': 0}}

    for p in range(elecs[0], elecs[1] + 1):
        for m in range(elecs[0], elecs[1] + 1):
            # if statements and loops for along the diagonal
            if p == m:
                # for loop for every recording electrode along the diagonal, super electrodes
                if p in super_elecs:
                    for i in range(min_elec, max_elec + 1):
                        if p != i:
                            nrtno = nrtno + 1
                            r = i
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p-min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m-min(elecs)],
                                                    'Recording Electrode': r}
                # for loop for multiple recording electrodes along the diagonal, non-super electrodes
                else:
                    for i in range(0, len(diagrecelec)):
                        if (p + diagrecelec[i]) >= min_elec and (p + diagrecelec[i]) < max_elec + 1:
                            nrtno = nrtno + 1
                            r = p + diagrecelec[i]
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p-min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m-min(elecs)],
                                                    'Recording Electrode': r}
                        elif len(diagrecelec)==1 and (diagrecelec[i] + p > 22 or diagrecelec[i] + p < 1):
                            # if we are only using 1 recording electrode, this switches to the other side at ends
                            nrtno = nrtno + 1
                            r = p - diagrecelec[i]
                            NRT_run_order[nrtno] = {'Probe Electrode': p, 'Probe Level': e_levels[p - min(elecs)],
                                                    'Masker Electrode': m, 'Masker Level': e_levels[m - min(elecs)],
                                                    'Recording Electrode': r}


    # Summarize PECAP run order
    print 'Generating Diagonal PECAP sequence for %d electrodes %d to %d:' % (elecs[1] - elecs[0] + 1, elecs[0], elecs[1])
    print 'El. \tLevel CUs'
    for i in range(0, elecs[1] - elecs[0] + 1):
        print str(i + min(elecs)), '\t\t', e_levels[i]
    print 'Stimulating Probe and Masker in MP1 Mode.'
    print 'Recording Electrode is %d relative to Probe, MP2 Mode.' % record_from
    print 'Phase Duration = %d us\tInterphase gap = %d us' % (params['PD'], params['ipg'])
    if params['reps'] == 5:
        print 'Probe Rate = 350 Hz\t\tMasker-Probe Interval = 400 us'
    else:
        print 'Probe Rate = 80 Hz\t\tMasker-Probe Interval = 400 us'
    print 'Gain is %d dB\t\t\tDelay is %d us' % (params['gain'],params['delay'])
    print 'Number of Reps = %d\t\tTotal number of Sweeps = %d x %d = %d' % (params['reps'], params['reps'],
                                                                            params['tsweeps'] / params['reps'],
                                                                            params['reps'] * (
                                                                            params['tsweeps'] / params['reps']))
    print 'Number of Samples = 32\tNumber of NRTs = %d\n' %nrtno

    return NRT_run_order
# function for generating the .csv file in which all the NRTs will be stored - FM
def gen_csv_fm(filename,NRTs,nrtno,params,elecs,e_levels):

    # open csv file
    fp = open(filename, 'w')		# to '.csv' file, for import in Excel

    # write headers for information that applies to all recorded values
    fp.write('Parameters\nGain,Delay,Interphase Gap,Phase Duration,Repetitions,Sweeps,Probe Indifferent Electrode,')
    fp.write('Masker Indifferent Electrode,Recording Indifferent Electrode,Probe Rate,Artefact Cancellation Technique\n')

    # write content for information that applies to all recorded values
    fp.write(str(params['gain']));                      fp.write(' dB,')
    fp.write(str(params['delay']));                     fp.write(' us,')
    fp.write(str(params['ipg']));                       fp.write(' us,')
    fp.write(str(params['PD']));                        fp.write(' us,')
    fp.write(str(params['reps']));                      fp.write(',')
    fp.write(str(params['tsweeps']/params['reps']));    fp.write(',')
    fp.write('MP1,')
    fp.write('MP1,')
    fp.write('MP2,')
    if params['reps']==5:
        fp.write('350 Hz,')
    else:
        fp.write('80 Hz,')
    fp.write('Forward Masking,')
    fp.write('\n\n')

    #write electrodes and MCL levels
    fp.write('Electrode No:,')
    for i in range(elecs[0],elecs[1]+1):
        fp.write(str(i))
        fp.write(',')
    fp.write('\nLevel (CU):,')
    for i in range(0,elecs[1]-elecs[0]+1):
        fp.write(str(e_levels[i]))
        fp.write(',')
    fp.write('\n\n')

    # write NRT Information
    fp.write('NRTs\n')
    for nrt in range(0,nrtno):
        for rep in range(0,params['reps']):
            fp.write('NRT,Repetition,Probe Electrode,Masker Electrode,Recording Electrode\n')
            fp.write(str(nrt+1));                               fp.write(',')
            fp.write(str(rep+1));                               fp.write(',')
            fp.write(str(NRTs[nrt+1][rep+1]['Probe']));         fp.write(',')
            fp.write(str(NRTs[nrt+1][rep+1]['Masker']));        fp.write(',')
            fp.write(str(NRTs[nrt+1][rep+1]['Recording']));     fp.write('\n')
            # write A
            fp.write('A,')
            for n in range(0,32):                               # assumes 32 samples default hasn't been changed
                fp.write(str(NRTs[nrt + 1][rep + 1]['A'][n]))
                fp.write(',')
            #write B
            fp.write('\nB,')
            for n in range(0, 32):                              # assumes 32 samples default hasn't been changed
                fp.write(str(NRTs[nrt + 1][rep + 1]['B'][n]))
                fp.write(',')
            # write C
            fp.write('\nC,')
            for n in range(0, 32):                              # assumes 32 samples default hasn't been changed
                fp.write(str(NRTs[nrt + 1][rep + 1]['C'][n]))
                fp.write(',')
            # write D
            fp.write('\nD,')
            for n in range(0, 32):  # assumes 32 samples default hasn't been changed
                fp.write(str(NRTs[nrt + 1][rep + 1]['D'][n]))
                fp.write(',')
            # write ECAP
            fp.write('\nECAP,')
            for n in range(0, 32):  # assumes 32 samples default hasn't been changed
                fp.write(str(NRTs[nrt + 1][rep + 1]['ABCD'][n]))
                fp.write(',')
            fp.write('\n')

    # now write the csv file will all the data in it

    # close csv file
    fp.close()


## PECAP Parameters
# enter parameters according to this subject
elecs = [1,22]
e_levels = [188, 180, 172, 176, 174, 174, 172, 172, 172, 172, 164, 164, 164, 154, 153, 152, 150, 150, 148, 147, 147, 146]
record_from = 2                    # i.e. '2' = record 2 apical from probe electrode, and -2 is record 2 basal
super_elecs = []                    # vector containing the stimulating electrodes you want to record from all rec elecs
diagrecelec = [-2, 2]                # electrodes to record from along the diagonal (i.e. p = m), i.e. [-6, -4, -2, 2, 4, 6]
params = {'gain':       50,     # in dB i.e. 40, 50, 60, or 70
          'delay':      98,    # in us i.e. 73, 98, or 122
          'ipg':        8.0,    # in us i.e. 7.0, 8.0, 42.0 - generally 8.0 is used
          'reps':       4,      # note that at present, if reps = 5, Probe Rate will be set at 350 Hz. Else, 80 Hz (FM)
          'tsweeps':    50,     # note that this is the aim of the total sweeps, including reps (5*10=50, or 4*12=48)
          'PD':         37.0}   # phase duration in us
filename_fm = 'JD_PECAP_01Nov2021.csv'    # name of csv filename that you'd like to save this data to for this subject

## Generate NRT Order
# check this has desired parameters before continuing to the next section
nrt_order = generate_nrt_order_fm(elecs, e_levels, record_from, super_elecs, diagrecelec, params)
#nrt_order = generate_nrt_order_fm_diag(elecs, e_levels, record_from, super_elecs, diagrecelec, params)


## PECAP Acquisition Script
# import time and date modules to allow use of elapsed time counters
import datetime, time
pecap_timestart = datetime.datetime.fromtimestamp(time.time())

# initiate NRT dictionary and counter
NRTs = {}                   # this dictionary will hold all ECAP responses and details
nrtno = len(nrt_order)      # this will count how many NRTs we have

# call NIC2 script for delivering and recording an ECAP, looped over the nrt_order
for i in range(1,nrtno+1):
    print 'NRT Number:', i, 'of', nrtno
    ECAPs = NIC2_ECAP_FM(nrt_order[i]['Probe Electrode'],nrt_order[i]['Probe Level'],nrt_order[i]['Masker Electrode'],
                      nrt_order[i]['Masker Level'],nrt_order[i]['Recording Electrode'],params)
    # check if the device was disconnected and recorded nothing
    if sum(ECAPs[1]['ABCD']) == 0.0:
        while sum(ECAPs[1]['ABCD']) == 0.0:
            print '\nConnection Bad - trying recording again in 3 seconds.\n'
            time.sleep(3)
            ECAPs = NIC2_ECAP_FM(nrt_order[i]['Probe Electrode'], nrt_order[i]['Probe Level'],
                                 nrt_order[i]['Masker Electrode'],
                                 nrt_order[i]['Masker Level'], nrt_order[i]['Recording Electrode'], params)
        NRTs[i] = ECAPs
    else:
        NRTs[i] = ECAPs

# timer end
pecap_timeend = datetime.datetime.fromtimestamp(time.time())
print "PECAP Total Elapsed Time: ", str(pecap_timeend - pecap_timestart)

# save NRTs stored internally in python dictionary to .csv file
gen_csv_fm(filename_fm,NRTs,nrtno,params,elecs,e_levels)

