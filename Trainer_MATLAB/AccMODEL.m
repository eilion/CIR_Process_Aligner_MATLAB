function AccMODEL(inputFile,ALPHA,BETA0,RHO0,initial,INT,max_iters,iters_for_R)

addpath('Codes/');
disp('-----------------------------------------------------------------');
disp('## This software is designed to learn the accumulation rate model (gamma process).');
disp(['#  The input file is ',inputFile,'.']);

% Load setting:
[setting] = getSetting(inputFile,ALPHA,INT,max_iters);

% Load data:
[data,param,SAMPLES,target,setting] = getData(inputFile,setting,initial,BETA0,RHO0);

% Initialization:
if isempty(initial)
    disp('#  Accumulation rates are being initialized...');
    parfor ll = 1:length(SAMPLES)
        [data(ll),SAMPLES(ll)] = SVGD(data(ll),param,SAMPLES(ll),target,setting,5000,100);
        [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,1000,100);
        disp(['   ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'% and sedRate = ',num2str(data(ll).R),'.']);
    end
    path = ['Initializations/',inputFile,'_K_',num2str(2*ALPHA),'_INT_',num2str(INT),'.mat'];
    save(path,'data','param','target','setting','SAMPLES');
end

R0 = length(setting.BETA);

for r = 1:setting.max_iters
    if ~isinf(iters_for_R)
        if rem(R0+r,iters_for_R)==1
            disp('## Average accumulation rates are being updated...');
            [data,SAMPLES] = learnRR(data,SAMPLES,target,setting);

            disp('## Accumulation rates are being initialized...');
            parfor ll = 1:length(SAMPLES)
                [data(ll),SAMPLES(ll)] = SVGD(data(ll),param,SAMPLES(ll),target,setting,5000,100);
                [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,1000,100);
                disp(['   ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'% and sedRate = ',num2str(data(ll).R),'.']);
            end
            path = ['Initializations/',inputFile,'_K_',num2str(2*ALPHA),'_INT_',num2str(INT),'.mat'];
            save(path,'data','param','target','setting','SAMPLES');
        end
    end

    disp(['## Iteration ',num2str(r),'/',num2str(setting.max_iters),':']);
    disp('#  Samples are being drawn...');
    parfor ll = 1:length(SAMPLES)
        [data(ll),SAMPLES(ll)] = GPST(data(ll),param,SAMPLES(ll),target,setting,250,100);
        disp(['   ',data(ll).name,' has been processed with the average accRate = ',num2str(mean(SAMPLES(ll).APP_RATE)),'% and sedRate = ',num2str(data(ll).R),'.']);
    end
    disp('#  Parameters are being updated...');

    param = learnParam(data,param,SAMPLES,setting);
    disp(['   BETA = ',num2str(param.BETA),'.']);
    disp(['   RHO = ',num2str(param.RHO),'.']);

    setting.BETA = [setting.BETA;param.BETA];
    setting.RHO = [setting.RHO;param.RHO];

    if rem(r,10) == 0
        path = ['Initializations/',inputFile,'_K_',num2str(2*ALPHA),'_INT_',num2str(INT),'.mat'];
        save(path,'data','param','target','setting','SAMPLES');
    end
end

path = ['Initializations/',inputFile,'_K_',num2str(2*ALPHA),'_INT_',num2str(INT),'.mat'];
save(path,'data','param','target','setting','SAMPLES');


disp('-----------------------------------------------------------------');


end