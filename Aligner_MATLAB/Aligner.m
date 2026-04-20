function [outputFile] = Aligner(inputFile)

addpath('Codes/');
disp('-----------------------------------------------------------------');
disp('## This software is designed to infer age models from d18O and 14C proxies.');
disp(['#  The input file is ',inputFile,'.']);

% Load data:
[data,SAMPLES,target,setting] = getData(inputFile);
disp('#  Initial Parameters are:');
for ll = 1:length(data)
    disp(['   ',data(ll).name,': Average Sedrate = ',num2str(1./data(ll).R),', Scale = ',num2str(data(ll).scale),', Shift = ',num2str(data(ll).shift),'.']);
end

disp('## Accumulation rates and age samples are being drawn...');
if length(SAMPLES) > 1
    parfor ll = 1:length(SAMPLES)
        if strcmp(setting(ll).islearn_average_sedrate,'Yes')||strcmp(setting(ll).islearn_shift,'Yes')||strcmp(setting(ll).islearn_scale,'Yes')
            for r = 1:10
                [data(ll),SAMPLES(ll)] = SVGD(data(ll),SAMPLES(ll),target,setting(ll),3000,100);
                [data(ll),SAMPLES(ll)] = GPST(data(ll),SAMPLES(ll),target,setting(ll),300,100);
                [data(ll),SAMPLES(ll)] = learnParam(data(ll),SAMPLES(ll),target,setting(ll));
            end
            disp(['#  ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'%.']);
            disp(['   Average Sedrate = ',num2str(1./data(ll).R),', Scale = ',num2str(data(ll).scale),', Shift = ',num2str(data(ll).shift),'.']);
        end
    end
else
    if strcmp(setting.islearn_average_sedrate,'Yes')||strcmp(setting.islearn_shift,'Yes')||strcmp(setting.islearn_scale,'Yes')
        for r = 1:10
            [data,SAMPLES] = SVGD(data,SAMPLES,target,setting,3000,100);
            [data,SAMPLES] = GPST(data,SAMPLES,target,setting,300,100);
            [data,SAMPLES] = learnParam(data,SAMPLES,target,setting);
        end
        disp(['#  ',data.name,' has been processed with the average accRate = ',num2str(mean(SAMPLES.APP_RATE)),'%.']);
        disp(['   Average Sedrate = ',num2str(1./data.R),', Scale = ',num2str(data.scale),', Shift = ',num2str(data.shift),'.']);
    end
end


disp('## Average accumulation rates have been estimated. Samples are being drawn...');
if length(SAMPLES) > 1
    parfor ll = 1:length(SAMPLES)
        [data(ll),SAMPLES(ll)] = SVGD(data(ll),SAMPLES(ll),target,setting(ll),3000,100);
        [data(ll),SAMPLES(ll)] = GPST(data(ll),SAMPLES(ll),target,setting(ll),300,100);
        [data(ll),SAMPLES(ll)] = GPST(data(ll),SAMPLES(ll),target,setting(ll),200,setting(ll).nSamples);
        disp(['   ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'% and sedRate = ',num2str(1./data(ll).R),'.']);
    end
else
    [data,SAMPLES] = SVGD(data,SAMPLES,target,setting,3000,100);
    [data,SAMPLES] = GPST(data,SAMPLES,target,setting,300,100);
    [data,SAMPLES] = GPST(data,SAMPLES,target,setting,200,setting.nSamples);
    disp(['   ',data.name,' has been processed with the average accRate = ',num2str(mean(SAMPLES.APP_RATE)),'% and sedRate = ',num2str(1./data.R),'.']);
end


CI_C14 = getCI_C14(data,target,setting);


outputFile = saveResults(inputFile,data,SAMPLES,target,CI_C14,setting);

disp(['## Results are saved in the folder Outputs/',outputFile,'/.']);


disp('-----------------------------------------------------------------');


end