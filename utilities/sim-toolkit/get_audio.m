function [audio, epochs] = get_audio(PATH_AUDIO)

% load audio
load(fullfile(PATH_AUDIO,'Audio.mat'),'z','t')
% load audio events
event_table = readtable(fullfile(PATH_AUDIO,'Audio_epochs.txt'));

audio = struct();
audio.trial = z;
audio.time = t;
audio.label = cell({'audio_p'});


t0 = event_table.t0;
epochs = event_table;
epochs.starts = t0 - 4.5;
epochs.ends = t0 + 4.5 - 1/1000;

epochs = bml_annot_table(epochs);





