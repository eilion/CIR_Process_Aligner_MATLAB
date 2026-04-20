function [data,SAMPLES,target,setting] = getData(inputFile)


path = ['Inputs/',inputFile,'/Records'];
list = dir(path);
list(strcmp({list.name},'.')) = [];
list(strcmp({list.name},'..')) = [];
list(strcmp({list.name},'.DS_Store')) = [];

L = length(list);

data = struct('name',cell(L,1),'depth',cell(L,1),'d18O',cell(L,1),'C14',cell(L,1),'ETC',cell(L,1),'scale',cell(L,1),'shift',cell(L,1),'R',cell(L,1),'MIN',cell(L,1),'MAX',cell(L,1));
setting = struct('name',cell(L,1),'nSamples',cell(L,1),'islearn_average_sedrate',cell(L,1),'islearn_shift',cell(L,1),'islearn_scale',cell(L,1),'alpha',cell(L,1),'beta',cell(L,1),'rho',cell(L,1),'a_d18O',cell(L,1),'b_d18O',cell(L,1),'a_14C',cell(L,1),'b_14C',cell(L,1));

for ll = 1:L
    setting(ll).name = list(ll).name;

    path = 'Defaults/setting_record.txt';
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);

    setting(ll).nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples:')==1});
    setting(ll).islearn_average_sedrate = INFO{2}{strcmp(INFO{1},'islearn_average_sedrate:')==1};
    setting(ll).islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
    setting(ll).islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};

    setting(ll).alpha = str2double(INFO{2}{strcmp(INFO{1},'alpha:')==1});
    setting(ll).beta = str2double(INFO{2}{strcmp(INFO{1},'beta:')==1});
    setting(ll).rho = str2double(INFO{2}{strcmp(INFO{1},'rho:')==1});
    setting(ll).a_d18O = str2double(INFO{2}{strcmp(INFO{1},'a_d18O:')==1});
    setting(ll).b_d18O = str2double(INFO{2}{strcmp(INFO{1},'b_d18O:')==1});
    setting(ll).a_14C = str2double(INFO{2}{strcmp(INFO{1},'a_14C:')==1});
    setting(ll).b_14C = str2double(INFO{2}{strcmp(INFO{1},'b_14C:')==1});

    setting(ll).start_depth = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
    setting(ll).min_start_age = str2double(INFO{2}{strcmp(INFO{1},'min_start_age:')==1});

    path = ['Inputs/',inputFile,'/Records/',list(ll).name,'/setting_record.txt'];
    if exist(path,'file') == 2
        fileID = fopen(path);
        INFO = textscan(fileID,'%s %s');
        fclose(fileID);

        if sum(strcmp(INFO{1},'nSamples:')==1) == 1
            setting(ll).nSamples = str2double(INFO{2}{strcmp(INFO{1},'nSamples:')==1});
        end

        if sum(strcmp(INFO{1},'islearn_average_sedrate:')==1) == 1
            setting(ll).islearn_average_sedrate = INFO{2}{strcmp(INFO{1},'islearn_average_sedrate:')==1};
        end

        if sum(strcmp(INFO{1},'islearn_shift:')==1) == 1
            setting(ll).islearn_shift = INFO{2}{strcmp(INFO{1},'islearn_shift:')==1};
        end

        if sum(strcmp(INFO{1},'islearn_scale:')==1) == 1
            setting(ll).islearn_scale = INFO{2}{strcmp(INFO{1},'islearn_scale:')==1};
        end

        if sum(strcmp(INFO{1},'alpha:')==1) == 1
            setting(ll).alpha = str2double(INFO{2}{strcmp(INFO{1},'alpha:')==1});
        end

        if sum(strcmp(INFO{1},'beta:')==1) == 1
            setting(ll).beta = str2double(INFO{2}{strcmp(INFO{1},'beta:')==1});
        end

        if sum(strcmp(INFO{1},'rho:')==1) == 1
            setting(ll).rho = str2double(INFO{2}{strcmp(INFO{1},'rho:')==1});
        end

        if sum(strcmp(INFO{1},'a_d18O:')==1) == 1
            setting(ll).a_d18O = str2double(INFO{2}{strcmp(INFO{1},'a_d18O:')==1});
        end

        if sum(strcmp(INFO{1},'b_d18O:')==1) == 1
            setting(ll).b_d18O = str2double(INFO{2}{strcmp(INFO{1},'b_d18O:')==1});
        end

        if sum(strcmp(INFO{1},'a_14C:')==1) == 1
            setting(ll).a_14C = str2double(INFO{2}{strcmp(INFO{1},'a_14C:')==1});
        end

        if sum(strcmp(INFO{1},'b_14C:')==1) == 1
            setting(ll).b_14C = str2double(INFO{2}{strcmp(INFO{1},'b_14C:')==1});
        end

        if sum(strcmp(INFO{1},'start_depth:')==1) == 1
            setting(ll).start_depth = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
        end

        if sum(strcmp(INFO{1},'min_start_age:')==1) == 1
            setting(ll).min_start_age = str2double(INFO{2}{strcmp(INFO{1},'min_start_age:')==1});
        end
    end


    data(ll).name = list(ll).name;

    path = ['Inputs/',inputFile,'/Records/',list(ll).name,'/record_info.txt'];
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);

    data(ll).R = str2double(INFO{2}{strcmp(INFO{1},'init_sedrate:')==1});

    if sum(strcmp(INFO{1},'num_Z:')==1) == 1
        data(ll).num_Z = str2double(INFO{2}{strcmp(INFO{1},'num_Z:')==1});
    else
        data(ll).num_Z = NaN;
    end

    data(ll).R = 1./data(ll).R;

    data(ll).scale = str2double(INFO{2}{strcmp(INFO{1},'init_scale:')==1});
    data(ll).shift = str2double(INFO{2}{strcmp(INFO{1},'init_shift:')==1});

    path = ['Inputs/',inputFile,'/Records/',list(ll).name,'/d18O_data.txt'];
    if exist(path,'file') == 2
        AA = readtable(path);
        DD_d18O = AA.depth;
        YY_d18O = AA.d18O;
    else
        DD_d18O = [];
        YY_d18O = [];
        data(ll).scale = NaN;
        data(ll).shift = NaN;
    end

    path = ['Inputs/',inputFile,'/Records/',list(ll).name,'/C14_data.txt'];
    if exist(path,'file') == 2
        AA = readtable(path);
        DD_C14 = AA.depth;
        YY_C14 = [AA.age,AA.error,AA.dR,AA.dSTD,AA.cc];
    else
        DD_C14 = [];
        YY_C14 = [];
    end

    path = ['Inputs/',inputFile,'/Records/',list(ll).name,'/additional_ages.txt'];
    if exist(path,'file') == 2
        AA = readtable(path);
        DD_add = AA.depth;
        YY_add_mean = AA.age;
        YY_add_stdv = AA.unct;
    else
        DD_add = [];
        YY_add_mean = [];
        YY_add_stdv = [];
    end


    DD = [DD_d18O;DD_C14;DD_add];
    DD = unique(DD);

    data(ll).d18O = [YY_d18O,zeros(size(YY_d18O,1),1)];
    for n = 1:size(YY_d18O,1)
        data(ll).d18O(n,end) = find(DD==DD_d18O(n));
    end

    data(ll).C14 = [YY_C14,zeros(size(YY_C14,1),1)];
    for n = 1:size(YY_C14,1)
        data(ll).C14(n,end) = find(DD==DD_C14(n));
    end

    data(ll).ETC = [YY_add_mean,YY_add_stdv,zeros(size(YY_add_mean,1),1)];
    for n = 1:size(YY_add_mean,1)
        data(ll).ETC(n,end) = find(DD==DD_add(n));
    end

    data(ll).depth = DD;
end






target = struct('stack',cell(1,1),'calibration',cell(1,1));

target.stack = struct('A',cell(1,1),'MU',cell(1,1),'PMU',cell(1,1),'SIG',cell(1,1),'PSIG',cell(1,1));

path = ['Inputs/',inputFile,'/stack.csv'];
AA = readtable(path);
target.stack.A = AA.age;
target.stack.MU = AA.mean;
target.stack.SIG = AA.sigma;
N = size(target.stack.A,1);
target.stack.PMU = zeros(N,1);
target.stack.PMU(1) = (target.stack.MU(2)-target.stack.MU(1))./(target.stack.A(2)-target.stack.A(1));
target.stack.PMU(N) = (target.stack.MU(N)-target.stack.MU(N-1))./(target.stack.A(N)-target.stack.A(N-1));
target.stack.PMU(2:N-1) = (target.stack.MU(3:N)-target.stack.MU(2:N-1))./(target.stack.A(3:N)-target.stack.A(2:N-1)) + (target.stack.MU(2:N-1)-target.stack.MU(1:N-2))./(target.stack.A(2:N-1)-target.stack.A(1:N-2));
target.stack.PMU(2:N-1) = target.stack.PMU(2:N-1)/2;
target.stack.PSIG = zeros(N,1);
target.stack.PSIG(1) = (target.stack.SIG(2)-target.stack.SIG(1))./(target.stack.A(2)-target.stack.A(1));
target.stack.PSIG(N) = (target.stack.SIG(N)-target.stack.SIG(N-1))./(target.stack.A(N)-target.stack.A(N-1));
target.stack.PSIG(2:N-1) = (target.stack.SIG(3:N)-target.stack.SIG(2:N-1))./(target.stack.A(3:N)-target.stack.A(2:N-1)) + (target.stack.SIG(2:N-1)-target.stack.SIG(1:N-2))./(target.stack.A(2:N-1)-target.stack.A(1:N-2));
target.stack.PSIG(2:N-1) = target.stack.PSIG(2:N-1)/2;


target.calibration = struct('A',cell(4,1),'MU',cell(4,1),'PMU',cell(4,1),'SIG',cell(4,1),'PSIG',cell(4,1));

path = 'Defaults/calibration_curves/IntCal20.csv';
AA = readtable(path);

target.calibration(1).A = AA.cal_bp/1000;
target.calibration(1).MU = AA.a14C_age/1000;
target.calibration(1).SIG = AA.error/1000;

N = size(target.calibration(1).A,1);
target.calibration(1).PMU = zeros(N,1);
target.calibration(1).PMU(1) = (target.calibration(1).MU(2)-target.calibration(1).MU(1))./(target.calibration(1).A(2)-target.calibration(1).A(1));
target.calibration(1).PMU(N) = (target.calibration(1).MU(N)-target.calibration(1).MU(N-1))./(target.calibration(1).A(N)-target.calibration(1).A(N-1));
target.calibration(1).PMU(2:N-1) = (target.calibration(1).MU(3:N)-target.calibration(1).MU(2:N-1))./(target.calibration(1).A(3:N)-target.calibration(1).A(2:N-1)) + (target.calibration(1).MU(2:N-1)-target.calibration(1).MU(1:N-2))./(target.calibration(1).A(2:N-1)-target.calibration(1).A(1:N-2));
target.calibration(1).PMU(2:N-1) = target.calibration(1).PMU(2:N-1)/2;
target.calibration(1).PSIG = zeros(N,1);
target.calibration(1).PSIG(1) = (target.calibration(1).SIG(2)-target.calibration(1).SIG(1))./(target.calibration(1).A(2)-target.calibration(1).A(1));
target.calibration(1).PSIG(N) = (target.calibration(1).SIG(N)-target.calibration(1).SIG(N-1))./(target.calibration(1).A(N)-target.calibration(1).A(N-1));
target.calibration(1).PSIG(2:N-1) = (target.calibration(1).SIG(3:N)-target.calibration(1).SIG(2:N-1))./(target.calibration(1).A(3:N)-target.calibration(1).A(2:N-1)) + (target.calibration(1).SIG(2:N-1)-target.calibration(1).SIG(1:N-2))./(target.calibration(1).A(2:N-1)-target.calibration(1).A(1:N-2));
target.calibration(1).PSIG(2:N-1) = target.calibration(1).PSIG(2:N-1)/2;


path = 'Defaults/calibration_curves/Marine20.csv';
AA = readtable(path);

target.calibration(2).A = AA.cal_bp/1000;
target.calibration(2).MU = AA.a14C_age/1000;
target.calibration(2).SIG = AA.error/1000;

N = size(target.calibration(2).A,1);
target.calibration(2).PMU = zeros(N,1);
target.calibration(2).PMU(1) = (target.calibration(2).MU(2)-target.calibration(2).MU(1))./(target.calibration(2).A(2)-target.calibration(2).A(1));
target.calibration(2).PMU(N) = (target.calibration(2).MU(N)-target.calibration(2).MU(N-1))./(target.calibration(2).A(N)-target.calibration(2).A(N-1));
target.calibration(2).PMU(2:N-1) = (target.calibration(2).MU(3:N)-target.calibration(2).MU(2:N-1))./(target.calibration(2).A(3:N)-target.calibration(2).A(2:N-1)) + (target.calibration(2).MU(2:N-1)-target.calibration(2).MU(1:N-2))./(target.calibration(2).A(2:N-1)-target.calibration(2).A(1:N-2));
target.calibration(2).PMU(2:N-1) = target.calibration(2).PMU(2:N-1)/2;
target.calibration(2).PSIG = zeros(N,1);
target.calibration(2).PSIG(1) = (target.calibration(2).SIG(2)-target.calibration(2).SIG(1))./(target.calibration(2).A(2)-target.calibration(2).A(1));
target.calibration(2).PSIG(N) = (target.calibration(2).SIG(N)-target.calibration(2).SIG(N-1))./(target.calibration(2).A(N)-target.calibration(2).A(N-1));
target.calibration(2).PSIG(2:N-1) = (target.calibration(2).SIG(3:N)-target.calibration(2).SIG(2:N-1))./(target.calibration(2).A(3:N)-target.calibration(2).A(2:N-1)) + (target.calibration(2).SIG(2:N-1)-target.calibration(2).SIG(1:N-2))./(target.calibration(2).A(2:N-1)-target.calibration(2).A(1:N-2));
target.calibration(2).PSIG(2:N-1) = target.calibration(2).PSIG(2:N-1)/2;


path = 'Defaults/calibration_curves/SHCal20.csv';
AA = readtable(path);

target.calibration(3).A = AA.cal_bp/1000;
target.calibration(3).MU = AA.a14C_age/1000;
target.calibration(3).SIG = AA.error/1000;

N = size(target.calibration(3).A,1);
target.calibration(3).PMU = zeros(N,1);
target.calibration(3).PMU(1) = (target.calibration(3).MU(2)-target.calibration(3).MU(1))./(target.calibration(3).A(2)-target.calibration(3).A(1));
target.calibration(3).PMU(N) = (target.calibration(3).MU(N)-target.calibration(3).MU(N-1))./(target.calibration(3).A(N)-target.calibration(3).A(N-1));
target.calibration(3).PMU(2:N-1) = (target.calibration(3).MU(3:N)-target.calibration(3).MU(2:N-1))./(target.calibration(3).A(3:N)-target.calibration(3).A(2:N-1)) + (target.calibration(3).MU(2:N-1)-target.calibration(3).MU(1:N-2))./(target.calibration(3).A(2:N-1)-target.calibration(3).A(1:N-2));
target.calibration(3).PMU(2:N-1) = target.calibration(3).PMU(2:N-1)/2;
target.calibration(3).PSIG = zeros(N,1);
target.calibration(3).PSIG(1) = (target.calibration(3).SIG(2)-target.calibration(3).SIG(1))./(target.calibration(3).A(2)-target.calibration(3).A(1));
target.calibration(3).PSIG(N) = (target.calibration(3).SIG(N)-target.calibration(3).SIG(N-1))./(target.calibration(3).A(N)-target.calibration(3).A(N-1));
target.calibration(3).PSIG(2:N-1) = (target.calibration(3).SIG(3:N)-target.calibration(3).SIG(2:N-1))./(target.calibration(3).A(3:N)-target.calibration(3).A(2:N-1)) + (target.calibration(3).SIG(2:N-1)-target.calibration(3).SIG(1:N-2))./(target.calibration(3).A(2:N-1)-target.calibration(3).A(1:N-2));
target.calibration(3).PSIG(2:N-1) = target.calibration(3).PSIG(2:N-1)/2;



path = ['Inputs/',inputFile,'/calibration_curve_14C.txt'];
if exist(path,'file') == 2
    AA = readtable(path);

    target.calibration(4).A = AA.cal_bp/1000;
    target.calibration(4).MU = AA.a14C_age/1000;
    target.calibration(4).SIG = AA.error/1000;


    N = size(target.calibration(4).A,1);
    target.calibration(4).PMU = zeros(N,1);
    target.calibration(4).PMU(1) = (target.calibration(4).MU(2)-target.calibration(4).MU(1))./(target.calibration(4).A(2)-target.calibration(4).A(1));
    target.calibration(4).PMU(N) = (target.calibration(4).MU(N)-target.calibration(4).MU(N-1))./(target.calibration(4).A(N)-target.calibration(4).A(N-1));
    target.calibration(4).PMU(2:N-1) = (target.calibration(4).MU(3:N)-target.calibration(4).MU(2:N-1))./(target.calibration(4).A(3:N)-target.calibration(4).A(2:N-1)) + (target.calibration(4).MU(2:N-1)-target.calibration(4).MU(1:N-2))./(target.calibration(4).A(2:N-1)-target.calibration(4).A(1:N-2));
    target.calibration(4).PMU(2:N-1) = target.calibration(4).PMU(2:N-1)/2;
    target.calibration(4).PSIG = zeros(N,1);
    target.calibration(4).PSIG(1) = (target.calibration(4).SIG(2)-target.calibration(4).SIG(1))./(target.calibration(4).A(2)-target.calibration(4).A(1));
    target.calibration(4).PSIG(N) = (target.calibration(4).SIG(N)-target.calibration(4).SIG(N-1))./(target.calibration(4).A(N)-target.calibration(4).A(N-1));
    target.calibration(4).PSIG(2:N-1) = (target.calibration(4).SIG(3:N)-target.calibration(4).SIG(2:N-1))./(target.calibration(4).A(3:N)-target.calibration(4).A(2:N-1)) + (target.calibration(4).SIG(2:N-1)-target.calibration(4).SIG(1:N-2))./(target.calibration(4).A(2:N-1)-target.calibration(4).A(1:N-2));
    target.calibration(4).PSIG(2:N-1) = target.calibration(4).PSIG(2:N-1)/2;
end


SAMPLES = struct('name',cell(L,1));
for ll = 1:L

    SAMPLES(ll).name = data(ll).name;

    if isnan(data(ll).num_Z)
        setting(ll).interval = (data(ll).depth(end)+1e-8-setting(ll).start_depth)./(length(data(ll).depth)-1);
        data(ll).num_Z = length(data(ll).depth);
    else
        setting(ll).interval = (data(ll).depth(end)+1e-8-setting(ll).start_depth)./(data(ll).num_Z-1);
    end
    SAMPLES(ll).D = linspace(setting(ll).start_depth+setting(ll).interval/2,setting(ll).start_depth+(data(ll).num_Z+0.5)*setting(ll).interval,data(ll).num_Z)';

    SAMPLES(ll).D0 = SAMPLES(ll).D - setting(ll).interval/2;
    SAMPLES(ll).D0 = SAMPLES(ll).D0(1:end-1);

    N0 = size(SAMPLES(ll).D,1);
    SAMPLES(ll).F = sqrt(data(ll).R).*ones(N0,100,2.*setting(ll).alpha)./sqrt(2.*setting(ll).alpha);
    SAMPLES(ll).F = SAMPLES(ll).F + normrnd(0,0.01,[N0,100,2.*setting(ll).alpha]);

    N = size(data(ll).depth,1);
    SAMPLES(ll).AGE = [];
    SAMPLES(ll).AGE_D = [];

    D1 = data(ll).depth;
    D1 = sort(D1,'ascend');

    AP = 0:setting(ll).interval/setting(ll).nSamples:(setting(ll).interval-setting(ll).interval/setting(ll).nSamples);
    SAMPLES(ll).D = SAMPLES(ll).D - AP;

    AP = 0:setting(ll).interval/100:(setting(ll).interval-setting(ll).interval/100);
    HH = zeros(N,N0,100);
    for k = 1:100
        for n = 1:N
            CNT = floor((D1(n)+AP(k))/setting(ll).interval);
            HH(n,1:CNT,k) = setting(ll).interval;
            HH(n,CNT+1,k) = D1(n) + AP(k) - CNT.*setting(ll).interval;
            HH(n,1,k) = HH(n,1,k) - AP(k);
        end
    end

    SAM_G = sum(SAMPLES(ll).F.^2,3);

    SAMPLES(ll).BIAS = ones(1,100) + normrnd(0,0.1,[1,100]);

    AA = SAMPLES(ll).BIAS.^2 - 0.1 + squeeze(sum(HH.*reshape(SAM_G,[1,N0,100]),2));

    SAMPLES(ll).AGE_D = [];
    SAMPLES(ll).AGE = AA;
    SAMPLES(ll).ACC_RATE = SAM_G;

    SAMPLES(ll).APP_RATE = [];
end


end