function apply_audio_p_coherence_batch(SUBJECT,CLEANMETHOD)
% apply_xcorr_audio_p
%
% This protocol performs a cross correlation analysis of all channels
% against the produced audio

bml_defaults
CRITERIA = 'C';

%defining paths
%CLEANMETHOD='TesPCDwithNotchConnectors'


PATH_PROTOCOL = fullfile('Z:\Commits\Vibration_artifacts\denoising_method\coherence-analysis\AllSubjects',SUBJECT);
PATH_DATA = 'Z:\DBS';
%PATH_PROTOCOL = '/Volumes/Nexus/Users/busha/Analysis/2020-08-05-audio_p-coherence';
%PATH_DATA = '/Volumes/Nexus/DBS';
PATH_SUBJECT = [PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT '/Preprocessed Data/Sync'];
PATH_CLEAN = ['Z:\Commits\Vibration_artifacts\denoising_method\Data\' CLEANMETHOD '\Data\' SUBJECT];
%PATH_CLEAN_v2 = ['Z:\Commits\Vibration_artifacts\denoising_method\Data\' CLEANMETHOD '\Outputs\Data\' SUBJECT];

PATH_ANALYSIS = PATH_PROTOCOL;  
PATH_FIG = [PATH_ANALYSIS '/figures'];
DATE=datestr(now,'yyyymmdd');

fprintf('audio_p coherence analysis for subject %s \n',SUBJECT)
  
cd(PATH_SYNC)
%loading annotation tables
%session = bml_annot_read(['annot/' SUBJECT '_session.txt']);
%tf_epoch = bml_annot_read(['annot/' SUBJECT '_trial_epoch.txt']);
%electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
artifact = bml_annot_read(['annot/' SUBJECT '_artifact_criteria_' CRITERIA '.txt']);
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);

%loading behavioral epochs
produced_syllable = bml_annot_read(['annot/' SUBJECT '_produced_syllable.txt']);
produced_triplet = bml_annot_read(['annot/' SUBJECT '_produced_triplet.txt']);

if isfile(['annot/' SUBJECT '_stimulus_syllable.txt'])
    stimulus_syllable = bml_annot_read(['annot/' SUBJECT '_stimulus_syllable.txt']);
else
	stimulus_syllable = bml_annot_read(['annot/' SUBJECT '_cue_presentation.txt']);
end


if isfile(['annot/' SUBJECT '_inter_trial_interval.txt'])
  ITI = bml_annot_read(['annot/' SUBJECT '_inter_trial_interval.txt']);
  pre_ITI = bml_annot_rename(ITI(:,{'id','starts','ends','trial_id_post','session_id'}),'trial_id_post','trial_id');
  pre_ITI.starts_ITI = pre_ITI.starts;
%  trial = bml_annot_left_join(produced_triplet(:,{'id','starts','ends','session_id','trial_id'}),pre_ITI);
  
  ITI = ITI(~ismissing(ITI.duration),:);
  ITI = ITI(ITI.duration>0,:);
else
  ITI = produced_triplet(:,{'id','session_id','trial_id'});
  ITI.starts = produced_triplet.starts - 4;
  ITI.ends = produced_triplet.starts - 2;
  pre_ITI = ITI;
%  stimulus_triplet.ends = produced_triplet.starts - 0.5;
%  trial = produced_triplet;
%  trial.starts_ITI = produced_triplet.starts - 4;
end

pre_ITI.midpoint = (pre_ITI.starts + pre_ITI.ends)/2;
pre_ITI.starts = pre_ITI.midpoint - 0.25;
pre_ITI.ends = pre_ITI.midpoint + 0.25;
pre_ITI = bml_annot_table(pre_ITI);



% compiling epochs of intereset for xcorr analysis
selcols = {'id','starts','ends','duration','session_id','trial_id','type',...
  'syl_id','stim','stim_volume','duration_C','duration_V','rms_audio_p'};
produced_syllable.type(:) = {'prod'};
stimulus_syllable.type(:) = {'stim'};
pre_ITI.type(:) = {'ITI'};
epoch = bml_annot_rowbind(produced_syllable(:,selcols), stimulus_syllable, pre_ITI);
epoch.epoch_id = epoch.id;
%epoch = bml_annot_extend(epoch,0.15,0.15); %extending epoch to do xcorr more discriminative

cd(PATH_PROTOCOL)
  
%% loading continuous preprocessed data
load([PATH_CLEAN filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'],'D');

% if isfile([PATH_CLEAN filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'])
%     
% elseif isfile([PATH_CLEAN_v2 filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'])
%         load([PATH_CLEAN_v2 filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'],'D');
% else
%     warning("file does not exist")
% end
cfg=[];
cfg.remask_nan = true;
cfg.value = 0;
Df = bml_mask(cfg, D);

%% applying high pass filter at 70Hz
cfg=[];
cfg.bpfilter = 'yes';
cfg.bpfreq = [70,180];
cfg.bpfiltord = 5;
cfg.bpfilttype = 'but';
cfg.bpfiltdir = 'twopass';
cfg.bsfilter = 'yes';
cfg.bsfreq = [59 61; 118 122; 177 186; 236 244];
cfg.bsfiltdir = 'twopass';
cfg.bsfilttype = 'but';
cfg.bsfiltord = 5;
Df=ft_preprocessing(cfg,Df);

%combining artifacts with empty_electrode table
% cfg=[];
% cfg.groupby='label';
% artifact_empty = bml_annot_union(cfg, artifact, empty_electrode);
% 
% %masking artifacts and empty_electrodes with NaNs
% cfg=[];
% cfg.annot=artifact_empty;
% cfg.label_colname = 'label';
% cfg.complete_trial = false; %masks entire trials
% cfg.value=NaN;
% Df = bml_mask(cfg,Df);


%% creating table of channels per epochs
epoch_ch = table();
for i=1:length(Df.label)
  epoc_ch_i = epoch;
  epoc_ch_i.channel1(:) = Df.label(i);
  epoch_ch = [epoch_ch; epoc_ch_i];
end
epoch_ch.id=[];
epoch_ch = bml_annot_table(epoch_ch);
epoch_ch.channel2(:) = {'audio_p'};

%% running coherence with bml_annot_calculate2
cfg=[];
cfg.epoch=epoch_ch;
coherence_audio_p = bml_annot_calculate2(cfg,Df,{'coherence','norm1','norm2'}, @get_coherence);

%saving data
writetable(coherence_audio_p, [PATH_ANALYSIS '/data/' SUBJECT '_' CLEANMETHOD '_coherence_audio_p.txt'],'Delimiter','\t');

