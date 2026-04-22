function draw_Samples(inputFile,ALPHA,BETA0,RHO0,INT)

% inputFile: the name of folder in 'Initializations/'. For example,
% 'DATA_02282025_K_4_INT_0.5'. Here, K = 2*ALPHA.
% ALPHA: the value of \alpha in the main manuscript. For example, '2.0'.
% BETA0: the initial value of \beta in the main manuscript.
% RHO0: the initial value of \rho in the main manuscript.
% INT: the length of unitless interval. For example, '0.5'.

addpath('Codes/');
disp('-----------------------------------------------------------------');
disp('## This software is designed to learn the accumulation rate model (gamma process).');
disp(['#  The input file is ',inputFile,'.']);

% Load setting:
[setting] = getSetting(inputFile,ALPHA,INT,0);

% Load data:
[data,param,SAMPLES,target,setting] = getData(inputFile,setting,[],BETA0,RHO0);

parfor ll = 1:length(SAMPLES)
    [data(ll),SAMPLES(ll)] = SVGD(data(ll),param,SAMPLES(ll),target,setting,5000,100);
    [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,1000,100);
    [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,1000,1000);
    disp(['   ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'% and sedRate = ',num2str(data(ll).R),'.']);
end

outputFile = saveResults(inputFile,data,param,SAMPLES,target,setting);

disp(['## Results are saved in the folder Outputs/',outputFile,'/.']);


disp('-----------------------------------------------------------------');


end