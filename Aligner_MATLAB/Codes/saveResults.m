function [outputFile] = saveResults(inputFile,data,SAMPLES,target,C14_CI,setting)

path = ['Outputs/',inputFile];

if exist(path,'dir') ~= 7
    mkdir(path);
end


for ll = 1:length(data)
    data(ll).R = 1./data(ll).R;
end


path = ['Outputs/',inputFile,'/Trial_001'];
if exist(path,'dir') == 7
    n = 1;
    DET = 1;
    while DET == 1
        n = n + 1;
        if n < 10
            path = ['Outputs/',inputFile,'/Trial_00',num2str(n)];
        elseif n < 100
            path = ['Outputs/',inputFile,'/Trial_0',num2str(n)];
        else
            path = ['Outputs/',inputFile,'/Trial_',num2str(n)];
        end
        
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end

mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'data','target','setting','SAMPLES','C14_CI');

outputFile = path(9:end);


end