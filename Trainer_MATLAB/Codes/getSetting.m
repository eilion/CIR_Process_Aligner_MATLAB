function [setting] = getSetting(inputFile,alpha,INT,max_iters)

setting = struct('max_iters',cell(1,1));

path = 'Defaults/setting.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

setting.K = round(2*alpha);
setting.max_iters = max_iters;
setting.nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples:')==1});
setting.a = str2double(INFO{2}{strcmp(INFO{1},'a:')==1});
setting.b = str2double(INFO{2}{strcmp(INFO{1},'b:')==1});
setting.interval = INT;

path = ['Inputs/',inputFile,'/setting.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);

    if sum(strcmp(INFO{1},'nSamples:')==1) == 1
        setting.nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples:')==1});
    end

    if sum(strcmp(INFO{1},'a:')==1) == 1
        setting.a = str2double(INFO{2}{strcmp(INFO{1},'a:')==1});
    end

    if sum(strcmp(INFO{1},'b:')==1) == 1
        setting.b = str2double(INFO{2}{strcmp(INFO{1},'b:')==1});
    end
end

setting.BETA = [];
setting.RHO = [];


end
