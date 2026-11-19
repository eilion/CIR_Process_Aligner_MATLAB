function draw_Samples(inputFile,ALPHA,BETA0,RHO0,INT,SEED)

% inputFile: the name of folder in 'Inputs/'. For example, 'DATA_02282025'.

% ALPHA: the value of \alpha in the main manuscript. For example, '2.0'.

% BETA0: the initial value of \beta in the main manuscript.

% RHO0: the initial value of \rho in the main manuscript.

% INT: the length of unitless interval. For example, '0.5'.

% SEED: random seed. If you want to "randomize" the running, set SEED to be
% 'NaN'.



addpath('Codes/');
disp('-----------------------------------------------------------------');
disp('## This software is designed to learn the accumulation rate model (gamma process).');
disp(['#  The input file is ',inputFile,'.']);

if isnan(SEED) || strcmp(SEED,'NaN')
    SEED = 'shuffle';
end

% random_seed:
STREAM = RandStream('mrg32k3a','Seed',SEED);
RandStream.setGlobalStream(STREAM);

% Load setting:
[setting] = getSetting(inputFile,ALPHA,INT,0);

% Load data:
[data,param,SAMPLES,target,setting] = getData(inputFile,setting,[],BETA0,RHO0);

parfor ll = 1:length(SAMPLES)
    % random_seed:
    workerStream = RandStream.getGlobalStream();
    workerStream.Substream = 1 + ll;

    [data(ll),SAMPLES(ll)] = SVGD(data(ll),param,SAMPLES(ll),target,setting,5000,100);
    [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,1000,100);
    [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,1000,1000);
    disp(['   ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'% and sedRate = ',num2str(data(ll).R),'.']);
end

outputFile = saveResults(inputFile,data,param,SAMPLES,target,setting);

disp(['## Results are saved in the folder Outputs/',outputFile,'/.']);


disp('-----------------------------------------------------------------');


end