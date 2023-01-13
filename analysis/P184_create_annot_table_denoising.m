% read coherence table and count number of electrode to remove

clearvars
close all
clc

%%
SUBJECT_LIST = readtable('X:\Commits\Vibration_artifacts\denoising_method\Denoising_allsubjects\Subjets-ALL.txt');
SUBJECT_LIST = SUBJECT_LIST.subject;
% idx_noclean=[29:length()];
SUBJECT_LIST(30:end) = [];

thr  = 3.08;
modalities = {"Raw", "Clean"};
methods = {"None", "PCD Audio"};

for s=1:length(SUBJECT_LIST)
    annot_vibration_artifact_electrodes = table();
     
    SUBJECT = SUBJECT_LIST{s};
    save_path = fullfile('X:\DBS', SUBJECT, '\Preprocessed Data\Sync\annot');
    save_path_clean = fullfile('X:\Commits\Vibration_artifacts\denoising_method\Data\Clean\Data', SUBJECT);
    save_path_raw = fullfile('X:\Commits\Vibration_artifacts\denoising_method\Data\Raw\Data', SUBJECT);


    elec_table.subj = string(SUBJECT);
    archive  = fullfile('X:\Commits\Vibration_artifacts\denoising_method\coherence-analysis\AllSubjects', SUBJECT, 'data\coherence_audio_p_by_electrode.tsv');

    CoTable = tdfread(archive,'\t');
    channel_names = cellstr(CoTable.electrode);
    
    
    session = [];
    label = [];
    coherence = [];
    method = [];
    artefactual = [];
    datatype = [];

    for m=1:length(modalities)
        for ss=1:max(CoTable.session_id)
            
            if any(ischar(CoTable.rNC))
                coherence_val = cellstr(CoTable.rNC);
                idx = find(CoTable.session_id ==ss & contains(string(CoTable.subject),modalities{m}));
                idx_dirty = find(CoTable.session_id ==ss & contains(string(CoTable.subject),modalities{m}) & str2double(cellstr(CoTable.rNC))>thr) - min(idx) + 1;
            else
                idx = find(CoTable.session_id ==ss & contains(string(CoTable.subject),modalities{m}));
                idx_dirty = find(CoTable.session_id ==ss & contains(string(CoTable.subject),modalities{m}) & CoTable.rNC>thr) - min(idx) + 1;
                coherence_val = CoTable.rNC;
            end
            binary = zeros(length(idx),1);
            binary(idx_dirty) = 1;


            session = [session; ss*ones(length(idx),1)];
            label = [label; channel_names(idx)];
            coherence = [coherence; coherence_val(idx) ];
            method = [method; repmat(methods{m}, length(idx), 1)];
            artefactual = [artefactual; binary];
            datatype = [datatype; repmat(modalities{m}, length(idx), 1)];
  
        end
    end
    %create table annot
    annot_vibration_artifact_electrodes.data = datatype;
    annot_vibration_artifact_electrodes.session = session;
    annot_vibration_artifact_electrodes.label = label;
    annot_vibration_artifact_electrodes.coherence_val = coherence;
    annot_vibration_artifact_electrodes.method = method;
    annot_vibration_artifact_electrodes.artifact = artefactual;

    writetable(annot_vibration_artifact_electrodes, fullfile(save_path, [SUBJECT, '_electrodes_vibration_artifact.txt']));
    
end


