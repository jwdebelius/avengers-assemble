clear

% ************************* Sample's PARAMETERS *********************
base_samples_dir = '../vaginal_samples/';
%sample_name = 'vaginal_samples';
primer_set_name = 'optivag';
kmer_len = 250;

% *********************** LOAD METHOD's PARAMETERS *******************
run('../Configs/params_script2')

% *************************** LOAD DB and TAXONOMY FILE *********************
% Load the taxonomy
if exist([uniS16_dir '/' AlgoConfig.taxaFileName],'file')
    Taxonomy = load([uniS16_dir '/' AlgoConfig.taxaFileName]);
else
    Taxonomy = [];
end


% Load the 16S sequences
if isempty(Taxonomy) && ~exist('Sequence_uni','var')
    uniS16_file = [uniS16_dir '/' db_filename '.fasta'];
    [Header_uni, Sequence_uni] = fastaread(uniS16_file);
else
    Header_uni = {};
    Sequence_uni = {};
end

% *************************** PROFILE SAMPLES  ***************************

samples = dir(base_samples_dir);
num_samples = length(samples);

% % *************************** PROFILES EACH SAMPLE *********************
% Skips the hidden files
for i = 1:num_samples
    if contains(samples(i).name, 'pool')
        samples(i).name
        SampleConfig = struct;
        SampleConfig.sample_name = samples(i).name;
        SampleConfig.sample_dir = samples(i).folder;
        SampleConfig.primers_seq = primers_seq;
    
        main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig, Header_uni,Sequence_uni,Taxonomy)
    end
end

% SampleConfig = struct;
% SampleConfig.sample_name = sample_name;
% SampleConfig.sample_dir = [base_samples_dir '/' sample_name];
% SampleConfig.primers_seq = primers_seq;
% 
% main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig,Header_uni,Sequence_uni,Taxonomy)
