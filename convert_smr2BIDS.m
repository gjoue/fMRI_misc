%% This script converts CED Spike2 .smr files for physiological data 
%% collected for noise correction in fMRI data to BIDS format 
%%
%% REQUIRES:
%%     Matlab code exchange tool Spike-smr-reader modified 
%%     with import data in seconds (not ticks)
%% INPUT: config structure with the following fields
%%     cfg.file.smr        the Spike2 smr file to convert
%%     cfg.keys.smr2Ignore any channel labels in .smr log to ignore
%%     cfg.debugPlots      plot channels for debugging? 0 or 1
%%     cfg.file.debugPlots where to save debug plots
%%     cfg.subjID          subject ID (for output file naming)
%%     cfg.task            experiment task label (for output file naming)
%%     cfg.dir.out         directory to output physio files -- by BIDS 
%%                         convention this should be the same folder as
%%                         the fMRI data you collected the physio data for
%%                         This script will create the run directories
%%                         and assume that the folder naming convention for
%%                         run directories within your task is "runN"
%%     cfg.run0            the 1st run number recorded by the smr logfile 
%%                           (for when several .smr files are saved for the
%%                            diff runs of the same task)
%%     cfg.MRI.TR          TR for your fMRI acquisition
%%     cfg.MRI.nVols       number of total volumes for your run 
%%     cfg.MRI.ndummies    number of dummy MRI scans
%% It assumes that if Spike2 was left running across multiple runs 
%% that no other non-task trigger-producing MRI sequences were run.
%% It does not handle cases when Spike2 was started too late (after the
%% experiment was started AND stopped too late so that any non-task-relevant 
%% MRI triggers were also recorded in the Spike2 logs.
%% For continuous Spike2 recordings (many experimental runs for one task in
%% one Spike2 log, no Spike2 data is thrown out but spliced so that the
%% first run has all the data from when Spike 2 started, even when well 
%% before the experiment starts, and subsequent runs have the data from
%% the 1st MR trigger (and has all the extra data when that run finishes until
%% right before the 1st MR trigger of subsequent run.
%% This also ignores "Keyboard" channel.
%%
%%
%% _____ NOTES _____
%% In addition the .smr data files, Spike2 also has the following files:
%% 1) saves .s2r "resource" files which have screen appearance 
%% settings for when the .smr files were last closed (screen layout, 
%% cursor positions)
%% 2) allows you to load s2c "configuration" files (what channels, sample rate)
%%
%% There are several ways to import CED's proprietary smr format -- 
%%  1) save Matlab file directly from Spike2, 
%%  2) use pSPM's file import or 
%%  3) Spike-smr-reader::ImportSMR
%%        which you can download from Matlab's file exchange website
%%  2+3 use the SON library which you can see in the pdf manuals that is
%%  bundled with Spike-smr-reader
%%
%% Some notes on the 3 imports are in ~/estropharm3/NOTES/smr_dataImport_tests.txt
%% Here, opting for option 3: Spike-smr-reader::ImportSMR which I have
%% modified to import .smr in seconds rather than ticks
%%
%% 2020 Gina Joue

function convert_smr2BIDS( cfg )

% %% make smr names match BIDS suggestions (note: BIDS suggests "trigger")
cfg0.keys.smr      = {'PULS','Resp','scanner'};
cfg0.keys.BIDS     = {'cardiac','respiratory','triggerMR'};
cfg0.tol.TR        = 0.1;   % MRI TR tolerance
cfg0.tol.rounding  = 0.001; % rounding error tolerance

cfg0.ext_j = ".json";
cfg0.ext_t = ".tsv";

fout_base = sprintf('%s_%s_physio_run', cfg.subjID, cfg.task);


tStart = tic;
fprintf('________ Importing .smr for %s, task %s [%s] ________\n', cfg.subjID, cfg.task, datestr(now));
data =ImportSMR( cfg.file.smr );

hdrs = [data.hdr];
cols = {hdrs.title};

colNames2use = regexprep(cols,cfg0.keys.smr,cfg0.keys.BIDS);
colNames2use(ismember(colNames2use, cfg.keys.smr2Ignore)) = [];

for cc=1:length(colNames2use)
    key = colNames2use{cc};
    i.(key) = find(strcmp(colNames2use,key)); 
end

%% _____  TIMING  ______
%% 
%% By default, Spike-smr-reader::ImportSMR is set to save timing in terms
%%   of ticks (easily changed to sec/msec).
%% It also reads in .smr headers and saves
%%      hdr.tim.scale = F.usPerTime          MICROSEC/TICK  e.g.  e3: 50
%%      hdr.tim.Units = F.dTimeBase                         e.g.  e3:  1.0000e-06
%%
%% Time stamp in smr files is the count of the number of underlying clock 
%% ticks from the start of the file to the event. 
%% The time of the event in seconds is:
%%     time in sec = timestamp * underlying clock tick period
%%   (Spike-smr-reader::SONTicksToSeconds.m)
%%     time in sec = timestamp * usePerTime * dTimeBase 

%%
%% According to MATLAB SON Library v.2.0 (2005),
%% When the source file is derived from frame based acquisition software 
%% (such as CED's Signal), the epochs in each row of data should be aligned 
%% with the trigger. However, this will not generally be the case with Spike2 
%% derived files. The triggered sampling mode introduced in Spike 2 Version 5 
%% produces sweeps of variable length and the trigger will occur at variable 
%% distances into each row of data (the sweep time and pre-trigger time set 
%% in Spike2 are treated as minimum values for each sweep (not constant 
%% values from sweep to sweep). 
%%
%% in imported Spike2 smr:
%% data(1).hdr.adc
%% ADC (Analogue to Digital Converter)
%%   - SampleInterval: hdr.adc.SampleInterval(1) 
%%     Prob corresponds to Spike2's "Microseconds per time unit" field 
%%     sets the time units for a new data file in the range 1 to 10000 microseconds.
%%     which is the sampling rate?
%%   - input in volts              data(:).hdr.adc.Units = 'volts'
%%   - which can be amplified      data(:).hdr.adc.Multiplex
%%
%% sample interval  = time between two samples in sec
%% sample rate/freq = and the reciprocal of the sample interval
%%                  = # samples/sec
%% Frequency resolution f = sample rate / n 
%%                        = 1 / (n* sample interval) 
%%                        = 1 / (input duration)

sampleInterv_sec = data(i.cardiac).hdr.adc.SampleInterval(1);
%sampleRate_sec   = 1/sampleInterv_sec;

%% input duration might be incorrect
%inputDur_sec     = data(i.cardiac).hdr.adc.Npoints / sampleInterv_sec;

%% According to the Spike2 manual, Spike2 often uses a sampling rate of 20 kHz/channel.
%% This matches the sampleRate_MR = 20kHz
Dt_MRtrig_sec    = diff(data(i.triggerMR).imp.tim);  
TR_MR_sec        = median(Dt_MRtrig_sec);    % = 30800 ticks = 1.54 sec
MRtrig_last      = length(data(i.triggerMR).imp.tim);

if TR_MR_sec - cfg.MRI.TR > cfg0.tol.rounding %% don't use ~= because of rounding diffs
    %% check whether we have times in ticks rather than sec
    if ( TR_MR_sec * data(i.triggerMR).hdr.tim.Scale * data(i.triggerMR).hdr.tim.Units - cfg.MRI.TR < cfg0.tol.rounding )
        error('ImportSMR.m appears to be importing in ticks, not sec as needed by this script. \nFix in ImportSMR l.28 with \nSONGetChannel(fid, chan,''progress'',''seconds'')');
    else
        warning('!!! Specified TR = %.2f sec but calculated %.2f from Spike2 logs. Going by Spike2 log calculations', cfg.MRI.TR, TR_MR_sec);
    end
end
        
%sampleRate_MR    = T_MR_ticks / cfg.MRI.TR;      % e3: 20000 = clock tick interval 

%%... find gap between runs = end of runs
[itrig_run_A, itrig_run_Z, nMRgot] = get_itrigXrun(data(i.triggerMR).imp.tim, TR_MR_sec, MRtrig_last);


%%___________________ fix ghost triggers __________________________________

itrig_ghost  = find( Dt_MRtrig_sec < TR_MR_sec - cfg0.tol.TR );
itrig_missed = find( Dt_MRtrig_sec >= 2 * TR_MR_sec - cfg0.tol.TR & Dt_MRtrig_sec <= 2 * TR_MR_sec + cfg0.tol.TR );

if ~isempty( itrig_ghost )
    
    D = diff([0;diff(itrig_ghost)==1; 0]);
    i_ghostA = itrig_ghost(D>0);
    i_ghostZ = itrig_ghost(D<0);
    
    i_ghosts2del = [];
    
    for ii=1:numel(i_ghostA)
        ii_ghostA = i_ghostA(ii);
        ii_ghostZ = i_ghostZ(ii);
        
        %% fix 1 MR trigger was split into several triggers
        %%   - consec missing triggers add up to what TR should be
          %|| ...
                %( sum( Dt_MRtrig_sec(ii_ghostA:ii_ghostZ) ) < TR_MR_sec + cfg0.tol.TR && ...
                %  ~isempty(find(Dt_MRtrig_sec(ii_ghostA:ii_ghostZ)==0, 1)) ...
                %) ...
        if (...
                ( sum( Dt_MRtrig_sec(ii_ghostA:ii_ghostZ) ) < TR_MR_sec + cfg0.tol.TR && ...
                  sum( Dt_MRtrig_sec(ii_ghostA:ii_ghostZ) ) > TR_MR_sec - cfg0.tol.TR ...
                )... 
            )
            i_ghosts2del = [i_ghosts2del; ii_ghostA+1:ii_ghostZ];
        else
            warning('Received ghost triggers but didn''t do anything about it');
        end
    end
    
    data(i.triggerMR).imp.tim( i_ghosts2del ) = []; % delete ghost(s)
    MRtrig_last = MRtrig_last - numel( i_ghosts2del );

                
    [itrig_run_A, itrig_run_Z, nMRgot] = get_itrigXrun(data(i.triggerMR).imp.tim, TR_MR_sec, MRtrig_last);
            
end

if ~isempty( itrig_missed )
    for ii=1:numel(itrig_missed)
        i_missed = itrig_missed(ii);
        
        if i_missed == numel(data(i.triggerMR).imp.tim)
             data(i.triggerMR).imp.tim = [ data(i.triggerMR).imp.tim;...
                data(i.triggerMR).imp.tim(missed) + TR_MR_sec/sampleInterv_sec]; % insert a trigger
        else
            data(i.triggerMR).imp.tim = [ data(i.triggerMR).imp.tim(1:i_missed);...
                data(i.triggerMR).imp.tim(missed) + TR_MR_sec/sampleInterv_sec ...
                data(i.triggerMR).imp.tim(missed+1:end) ]; % insert a trigger
        end
        MRtrig_last = MRtrig_last + 1;
        
        [itrig_run_A, itrig_run_Z, nMRgot] = get_itrigXrun(data(i.triggerMR).imp.tim, TR_MR_sec, MRtrig_last);
    end
end

%%______________________  __________________________________

nMRdiff  = nMRgot - repmat(cfg.MRI.nVols, length(itrig_run_Z),1);
nRuns    = length(nMRgot);    

%%____________ find diff in MR triggers received across runs ______________
idiff    = find( diff( nMRgot ) ~= 0);

%% diff takes pairs so results in n-1 vec so have to shift indices to match 
%% orig vec (except if the 1st index of the orig vector is what is diff) -- 
%% this might return indices out of bounds but is fixed in
%% the if block afterwards
for dd=1:numel(idiff)
    if idiff(dd) > 1
        idiff(dd) = idiff(dd) + 1;
    end
end

%% idiff returns pairs between which there are differences -- remove the 
%% second of the pair of indices because we had shifted by 1 before
if ~isempty(idiff)
    idiff = idiff([true;diff(idiff)~=true]);
end

%%____________ ID false starts (runs that are way too short)
blocks2skip  = find( nMRdiff ~= 0 & abs(nMRdiff) > cfg.MRI.nVols/2 );
nFalseStarts = numel(blocks2skip);

nRuns_real   = nRuns - nFalseStarts;

fprintf('\tDetected %d run(s), %d of which = false start(s)\n', nRuns, nFalseStarts);

%%____________ remove false start blocks __________________________________
itrig_run_A( blocks2skip) = [];
itrig_run_Z( blocks2skip ) = [];

           

                
%%____________ calculate offset from initial/theoretical starts ___________
%%... init when Spike2 started with experiment start time wrt
%%...    timestamps as logged in our files
startEXP_sec    = zeros(nRuns_real,1);
relstartSPK_sec = zeros(nRuns_real,1);  % start of Spike2 wrt experiment start (in sec)

startEXP_sec(1:nRuns_real)      = data(i.triggerMR).imp.tim( itrig_run_A + cfg.MRI.ndummies );  % experiment start time wrt current timestamps in sec

%% Spike2 start time wrt startEXP_sec (keep the extra physio signal before actual start of exp in case)
%% Note TR_MR_sec * cfg.MRI.ndummies for e3 is 1.54 sec/vol * 5 dummy vols = 7.7 sec
%% For the 1st run, the 1st MR trigger is not when Spike2 started (prob started well before that any triggers),
%% but for any subseq runs concatenated in the same Spike2 log, set start of physio at 1st MR trigger in that run.

relstartSPK_sec(1)   =  0-startEXP_sec(1);

if nRuns_real > 1
    relstartSPK_sec(2:nRuns_real)   = data(i.triggerMR).imp.tim( itrig_run_A(2:nRuns_real) ) - startEXP_sec(2:nRuns_real);  
end                
                
%% idiff at start and end of .smr logs are usually when Spike2 was started
%% or stopped too late. If in middle, then false starts (scanner started and
%% stopped)
if isempty(idiff) % no need to offset other than by dummy scans
    if ( nMRgot(1) - cfg.MRI.nVols > cfg0.tol.rounding )
        warning("\n\t*** Specified INCORRECT number of MR triggers (%d) -- prob should be %d. \n\tWill continue with %d and offset by %d MRI dummy scans **", cfg.MRI.nVols, nMRgot(1), nMRgot(1), cfg.MRI.ndummies);
    end  
    
else
    for dd=1:length(idiff) %% don't really need the loop for script now but leave in for more complicated scenarios
        %% If rec'd less triggers (nMRdiff<0) then started Spike2 too late.
        %% If rec'd more triggers than expected, assume that forgot to stop 
        %% Spike2 before another MR sequence started (for last run).
        irun = idiff(dd);
        
        if ~ismember(irun,blocks2skip) % not a false start
            if irun == 1 %% if 1st run then started Spike2 too late. Only have to adjust 1st run logged.
                startEXP_sec(irun) = startEXP_sec(irun) + nMRdiff(irun) * TR_MR_sec; % nMRdiff returns negative so add to shift earlier (so negative times here mean exp started before Spike2 logging)

                relstartSPK_sec(irun)   = data(i.triggerMR).imp.tim( itrig_run_A(irun) ) - startEXP_sec(irun);   % leaving irun rather than hardcoding 1 to facilitate future mods
            elseif irun == nRuns || irun == nRuns_real %% if last run then stopped Spike2 too late -- don't have to worry too much about this
                %% don't do anything for now    
            else
                figure; 
                sgtitle(sprintf('ERROR with run %d in .smr logs for %s', irun, cfg.subjID),'Color','red')
                plot( data(i.triggerMR).imp.tim);
                error('unhandled situation: missing %d MR triggers not in 1st/last run but in run %d', nMRdiff(irun), irun);            
            end  
        end
    end    
end


%%________________ shift timestamps/MR triggers within runs _______________
series_MRtrig               = zeros(data(i.cardiac).hdr.adc.Npoints, 1); 
ismpl_trig                  = int64( data(i.triggerMR).imp.tim * sampleInterv_sec );
series_MRtrig(ismpl_trig')  = 1;  % populate with trigger 1 corresp to a trigger received

mats2write                  = [];%zeros(1, length(colNames2use), nRuns_real); 

%% build matrix and split at 1st MR trigger of a run to just before 1st MR trigger for subseq run

for r=1:nRuns_real
    r_start = 0;
    r_end   = 0;
    
    if nRuns_real > 1
       %% grab times for 1st MR trigger of current run 
       r_start    = int64( data(i.triggerMR).imp.tim( itrig_run_A(r) ) * sampleInterv_sec );
       
       if r ~= nRuns_real
           r_end    = int64( data(i.triggerMR).imp.tim( itrig_run_Z(r) ) * sampleInterv_sec ); % lose 1 TR of physio signal
       end
    end
    
    %% also covers when have only one run (then one matrix)
    if r==1
        r_start    = int64(1);
    end
    
    if r == nRuns_real
        r_end    = data(i.respiratory).hdr.adc.Npoints; % sometimes resp + cardiac have diff Npoints! take smaller of the two..resp slower? 
    end
    
    for cc=1:length(colNames2use)
        src = colNames2use{cc};

        if ~isempty( find(strcmp(src, cfg.keys.smr2Ignore), 1) )
            continue;
        end

        if strcmp(src,'triggerMR')
            mats2write{r}(:,cc) = series_MRtrig(r_start:r_end);
        else
            mats2write{r}(:,cc) = data(i.(src)).imp.adc(r_start:r_end);
        end
    end    
end


ntrigs = itrig_run_Z - itrig_run_A + 1;	

%% PLOT to CHECK mats2write
if cfg.debugPlots   
    if ~isfield(cfg.file, 'debugPlots')
        cfg.file.debugPlots = pwd;
    end
    
    figID = figure('Visible','off'); 
    sgtitle(sprintf('%s: detected %d run(s) with %d false-start run(s), %d missing/%d ghost MR triggers', cfg.subjID, nRuns, nFalseStarts, numel(itrig_missed), numel(itrig_ghost) ));
    subplot(2,1,1);
    plot( data(i.triggerMR).imp.tim);
    title('MR triggers received');
    xlabel('MR trigger number');
    ylabel('TR #');

    %% annotate with number of triggers/itrig start and end
    ntrigs = itrig_run_Z - itrig_run_A + 1;
w    
    for ii=1:length(itrig_run_A)
        text(itrig_run_A(ii), 10, num2str(itrig_run_A(ii)));
        text(itrig_run_Z(ii), 500, num2str(itrig_run_Z(ii)));
        
        text( itrig_run_A(ii) + ntrigs(ii)/2, 1000, sprintf('#=%d',ntrigs(ii)) );
    end
    
    subplot(2,1,2);
    plot( data(i.cardiac).imp.adc,'Color',[0.6350, 0.0780, 0.1840]); hold on; plot( data(i.respiratory).imp.adc,'Color',[0.9290, 0.6940, 0.1250]); plot(series_MRtrig*1e4,'.','Color',[0, 0.4470, 0.7410]'); hold off
    title(sprintf('skipped run %d. ', blocks2skip));
    xlabel('time (sec)');
    ylabel('volt');
    legend({'cardiac','respiratory','MR trigger'});

    print(figID,'-dpsc2','-append',cfg.file.fplots);  %% NOTE APPEND BREAKS IN PARALLEL MODE -- better save as separate pdfs if want to run in parfor loop

    % figure; plot(mats2write{1}(:,1)); hold on; plot(mats2write{1}(:,2)); plot(mats2write{1}(:,3)*1e4, 'r.'); hold off
    % figure; plot(mats2write{2}(:,1)); hold on; plot(mats2write{2}(:,2)); plot(mats2write{2}(:,3)*1e4, 'r.'); hold off
    % figure; plot(mats2write{3}(:,1)); hold on; plot(mats2write{3}(:,2)); plot(mats2write{3}(:,3)*1e4, 'r.'); hold off
    % figure; plot(mats2write{4}(:,1)); hold on; plot(mats2write{4}(:,2)); plot(mats2write{4}(:,3)*1e4, 'r.'); hold off
    % figure; plot(mats2write{5}(:,1)); hold on; plot(mats2write{5}(:,2)); plot(mats2write{5}(:,3)*1e4, 'r.'); hold off
    % figure; plot(mats2write{6}(:,1)); hold on; plot(mats2write{6}(:,2)); plot(mats2write{6}(:,3)*1e4, 'r.'); hold off
end

%%________________ create files _______________
if ~isfield(cfg,'run0') 
    cfg.run0 = 1; % assume default starting run, if not set is 1
end

if isstring( cfg.run0 )
    cfg.run0 = str2num(cfg.run0);
end

%%__________ print JSON _______
for r=1:nRuns_real
    dout_base = fullfile(cfg.dir.out, ['run' num2str(r)]);
   
    runNo = r + cfg.run0 - 1;
    
    jID = fopen( fullfile(dout_base, sprintf('%s%d%s',fout_base, runNo, cfg0.ext_j) ), 'w');
    
cols2print = sprintf('\"%s\",',colNames2use{:});
    cols2print = strip(cols2print,'right',',');
    
    fprintf(jID, '{\n');
    fprintf(jID, '   \"SamplingFrequency\": %.1f,\n',    sampleInterv_sec);
    fprintf(jID, '   \"StartTime\": %.1f,\n',            relstartSPK_sec(r) );
    fprintf(jID, '   \"Columns\": [%s],\n',              cols2print);
    fprintf(jID, '   \"SourceFile\": \"%s\",\n',         cfg.file.smr); %% this is not part of BIDS
    fprintf(jID, '   \"ntriggerMR\": %d\n',              ntrigs(r));    %% this is not part of BIDS
    fprintf(jID, '}');
    
    fclose(jID);
end

%%__________ print tsv's _______
for r=1:nRuns_real

    runNo = r + cfg.run0 - 1;

    
    dout_base = fullfile(cfg.dir.out, ['run' num2str(runNo)]);
    fout      = fullfile(dout_base, sprintf('%s%d%s',fout_base, runNo, cfg0.ext_t));
    %tID       = fopen(fout, 'w');
    writematrix(mats2write{r},fout, 'FileType','text', 'Delimiter','\t');
    %fclose(tID);
    
    fprintf('\tcreated file %s and json file [%s].......\n', fout, datestr(now));
end

tEnd = toc(tStart);

fprintf('........ FINISHED %d pairs of json/tsv files %s [Took %.1f sec].......\n', nRuns_real, fout_base, tEnd);

end

%%%%%%%%%%%%
%%... find gap between runs = end of runs
function [itrigA, itrigZ, nMRtrigs] = get_itrigXrun(trig_times, TR_sec, lastTR)
        itrigZ    = find(diff(trig_times) > (TR_sec + TR_sec/2) );  % trigger # marking end of runs up to nruns-1
        itrigA    = [ 1; itrigZ + 1 ];  % initialize trigger # marking beginning of all runs
        itrigZ    = [itrigZ; lastTR];   % Do we have the same number of triggers for all runs?
        nMRtrigs  = diff( [0; itrigZ] );
end



%%%%%%%%%%%%%%%  NOTES ON RETROICOR and TAPAS' PhysIO  %%%%%%%%%%%%%%%%%%%%
%% The Fourier expansion of cardiac and respiratory phases was introduced as
%% RETROICOR (RETROspective Image CORrection, Glover2000, see also
%% Josephs1997). These Fourier Terms can enter a General Linear Model (GLM)
%% as nuisance regressors, analogous to movement parameters. As the
%% physiological noise regressors augment the GLM and explain variance in
%% the time series, they increase sensitivity in all contrasts of interest.
%% 
%% The core reference for PhysIO is: The PhysIO Toolbox for Modeling
%% Physiological Noise in fMRI Data
%% (http://dx.doi.org/10.1016/j.jneumeth.2016.10.019) Please cite this paper
%% if you use PhysIO in your work. Moreover, this paper is also a good
%% source for more information on PhysIO (see next question). A standard
%% snippet to include in your method section could look like the following,
%% assuming you use our specific implementation of RETROICOR, which uses
%% Fourier expansions of different order for the estimated phases of cardiac
%% pulsation (3rd order), respiration (4th order) and cardio-­?respiratory
%% interactions (1st order) following (Harvey et al., 2008)

%%%%%%%%%%%%%%%  NOTES ON NOISE HANDLING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  low-frequency (< 0.15 Hz) fluctuations of the BOLD fMRI signal in the 
%%  resting brain have consistently revealed significant correlations 
%%  between distinct brain regions giving rise to a number of functional 
%%  networks, termed resting-state networks (RSNs) (Biswal et al., 1995; 
%%  Smith et al., 2009). 
%% 
%% ___ FLUCTUATIONS IN BOLD__
%% 
%% * cardiac pulsation 
%%    - pushes brainstem into the surrounding brain tissue, causing 
%%      deformation and cerebrospinal fluid movement (Dagli et al., 1999) 
%% 
%% * respiration-induced fluctuations result
%%    - bulk movement of the head (Hu et al., 1995). 
%% 
%% varying rate or depth of respiration => change in arterial tension of 
%% CO2 => vasodilation => changes in CBF (Birn et al., 2006; 
%% Wise et al., 2004). 
%% 
%% 
%% global CBF changes => low-frequency (~0.1 Hz) fluctuations in the BOLD 
%%  signal (Birn et al., 2006, 2008)
%% 
%% spontaneous fluctuations in heart rate => fluctuations in the BOLD 
%%  (Napadow et al., 2008; Shmueli et al., 2007)
%% 
%% 
%% ___MODEL SOLUTIONS FOR FMRI PHYSIO NOISE CORREX__
%% 
%% RETROICOR (Glover et al. 2000)
%% * one of the most widely used methods 
%% * method: physiological regressors are estimated as a linear combination 
%%   of sinusoidal signals coinciding with the cardiac and respiratory cycles 
%%   using concurrent cardiac and respiratory measurements and subsequently 
%%   regressed out.
%% * results: can effectively remove the high-frequency cardiac (~1 Hz) and 
%%   respiratory (~0.3 Hz) artifacts despite aliasing due to low sampling 
%%   rate/TR in typical fMRI acquisition (Jones et al., 2008).
%% 
%% Chang et al. (2009)
%% * physiological response function (PRF) model
%% * cardiac convolved with cardiac response function (CRF).
%% * in physiological noise modelling (PNM) toolbox of FSL and (Jenkinson 
%%    et al., 2012) and the PhysIO SPM toolbox (Kasper et al., 2017)
%% * does not account for between-subject PRF variability
%% 
%% Birn et al. (2006, 2008)
%% * physiological response function (PRF) model
%% * respiration volume per time (RVT) (proportional to breathing rate and 
%%    depth estimated based on respiratory measurements) convolved with the 
%%    respiration response function (RRF) 
%% * in physiological noise modelling (PNM) toolbox of FSL and (Jenkinson 
%%    et al., 2012) and the PhysIO SPM toolbox (Kasper et al., 2017)
%% * does not account for between-subject PRF variability
%% 
%% Falahpour et al. (2013) 
%% * subject-specific PRF curves based on global signal (mean BOLD across 
%%     all voxels in the brain) from each scan
%% * accounts for a considerably more variance in fMRI timeseries than 
%%     standard PRF <- but may be overfitting (Falahpour et al., 2013).
%% 
%% 
%% 3 ways suggested for modeling respiration changes
%% 
%% 1) breathing changes convolved with BOLD-HRF (Abbott et al., 2005; 
%%      Thomason et al., 2005; Thomason et al., 2007); 
%% 
%% 2) time-shifting a boxcar waveform representing the cue for breathing 
%%      changes (e.g. cues for breath-holding) (Kastrup et al., 1999; 
%%      Kastrup et al., 1999; Li et al., 1999);
%% 
%% 3) time-shifting an estimate of the respiration volume per unit time 
%%      (Birn et al., 2006). 
%% 
%% 
%% ___DATA-DRIVEN SOLUTIONS FOR FMRI PHYSIO NOISE CORREX__
%% 
%% global signal regression (GSR)
%% * implicitly assumes that processes that globally affect the fMRI BOLD 
%%      signal are mostly uncorrelated to neural activity (Fox et al., 2005; 
%%      Greicius et al., 2003; Qing et al., 2015)
%% * disadv: also regresses out neuronal signal (Liu et al., 2017; 
%%      Murphy and Fox, 2017)
%% 
%% ICA/PCA
%% * identify physiological/?noisy? components based on temporal, spatial 
%%      and spectral features (Churchill and Strother, 2013; Kay et al., 
%%      2013; Pruim et al., 2015; Salimi-Khorshidi et al., 2014)
%% * e.g. FIX (?FMRIB?s ICA-based X-noisefier?)
%%     - good perf so used in HCP rs-fMRI connectomes (Salimi-Khorshidi 
%%      et al., 2014). 
%% * disadv: may not be enough correction (Chang and Glover, 2009a; 
%%      Falahpour et al., 2013). 

%%%%%%%%%%%%%%% NOTES ON BIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Physiological recordings specified using two files: 
%% 1. a gzip compressed TSV file with data (without header line) 
%%        time 
%% 2. JSON for storing start time (sec, wrt time of 1st vol -- neg vals okay)
%%                     sampling frequency (Hz), and 
%%                     name of the columns from the TSV.
%%                            naming conventions:
%%                                - "cardiac": cts pulse measurement
%%                                - "respiratory": cts breathing measurement
%%                                - "trigger": cts scanner trigger signal
%%  e.g.
%%       sub­control01/
%%                     func/
%%                          sub­control01_task­nback_physio.tsv.gz
%%                        ...................................
%%                               34    110     0
%%                               44    112     0
%%                               23    100     1
%%                        ...................................
%%       sub­control01/
%%                     func/
%%                          sub­control01_task­nback_physio.json
%%                        ...................................
%%                               {
%%                                  "SamplingFrequency": 100.0,                          in Hz
%%                                  "StartTime": ­22.3,                                  in sec
%%                                  "Columns": ["cardiac", "respiratory", "trigger"]     headers not in tsv
%%                               }
%%                        ...................................
%% FILE NAME CONVENTIONS:
%%  *  Physiological recordings 
%%        (including eye tracking + scnner motion correx) should use the 
%%        suffix "_physio"
%%  *  stimulus: "_stim"
%%
