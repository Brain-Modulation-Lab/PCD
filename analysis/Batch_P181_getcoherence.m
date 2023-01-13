ANALYSIS_PATH = 'Z:\Commits\Vibration_artifacts\denoising_method\coherence-analysis\AllSubjects';

SUBJECTS_TBL = readtable(fullfile(ANALYSIS_PATH,'Subjects_list.txt'));
PROTOCOL_TBL = readtable(fullfile(ANALYSIS_PATH,'Conditions_list.txt'));

PROTOCOL_FUNCTION = 'apply_audio_p_coherence_batch';

for subj_i = 1 : height(SUBJECTS_TBL)
    SUBJECT = SUBJECTS_TBL.Subject{subj_i};
    PROTOCOL_PATH = fullfile(ANALYSIS_PATH, SUBJECT);
    exe_daytime = datestr(now,'yyyymmdd_HHMM');
    diary([PROTOCOL_PATH filesep 'batch_' PROTOCOL_FUNCTION '_' exe_daytime '.log'])

% addpath(PROTOCOL_PATH);
% 
% cd(PROTOCOL_PATH)

     
fprintf('=== Running protocol %s ===\n',PROTOCOL_FUNCTION)

for i=1:height(PROTOCOL_TBL)
  CONDITION = PROTOCOL_TBL.Conditions{i};
  fprintf('Running protocol.')
  %running protocol
  try
    proto = str2func(PROTOCOL_FUNCTION);
    proto(SUBJECT,CONDITION);
    fprintf('OK\n')
    
  catch err
    fprintf('FAILED: %s\n',err.message)
  end
  cd(ANALYSIS_PATH)
end
diary('off')
end