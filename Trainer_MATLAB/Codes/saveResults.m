function [outputFile] = saveResults(inputFile,data,param,SAMPLES,target,setting)

path = ['Outputs/',inputFile];

if exist(path,'dir') ~= 7
    mkdir(path);
end


path = ['Outputs/',inputFile,'/KK_',num2str(setting.K),'_',num2str(setting.interval),'/Trial_001'];
if exist(path,'dir') == 7
    n = 1;
    DET = 1;
    while DET == 1
        n = n + 1;
        if n < 10
            path = ['Outputs/',inputFile,'/KK_',num2str(setting.K),'_',num2str(setting.interval),'/Trial_00',num2str(n)];
        elseif n < 100
            path = ['Outputs/',inputFile,'/KK_',num2str(setting.K),'_',num2str(setting.interval),'/Trial_0',num2str(n)];
        else
            path = ['Outputs/',inputFile,'/KK_',num2str(setting.K),'_',num2str(setting.interval),'/Trial_',num2str(n)];
        end
        
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end

mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'data','param','target','setting','SAMPLES');

outputFile = path(9:end);


end