function draw_Samples(inputFile,ALPHA,BETA0,RHO0,INT)

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