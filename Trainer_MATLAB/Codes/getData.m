function [data,param,SAMPLES,target,setting] = getData(inputFile,setting,initial,BETA,RHO)

path = ['Initializations/',initial,'.mat'];

if exist(path,'file') == 2

    max_iters = setting.max_iters;

    AA = load(path);
    data = AA.data;
    param = AA.param;
    setting = AA.setting;
    target = AA.target;
    SAMPLES = AA.SAMPLES;

    setting.max_iters = max_iters;

    L = length(SAMPLES);

    LENGTH = zeros(L,1);

    for ll = 1:L
        if size(SAMPLES(ll).F,2) < setting.nSamples
            RAND_SEED = ceil(size(SAMPLES(ll).F,2).*rand(setting.nSamples,1));

            SAMPLES(ll).F = SAMPLES(ll).F(:,RAND_SEED,:);
            SAMPLES(ll).BIAS = SAMPLES(ll).BIAS(:,RAND_SEED);
        elseif size(SAMPLES(ll).F,2) > setting.nSamples
            SAMPLES(ll).F = SAMPLES(ll).F(:,1:setting.nSamples,:);
            SAMPLES(ll).BIAS = SAMPLES(ll).BIAS(:,1:setting.nSamples);
        end
        
        LENGTH(ll) = max(data(ll).depth.*data(ll).R);
    end

    [~,ORDER] = sort(LENGTH,'descend');
    data = data(ORDER);
    SAMPLES = SAMPLES(ORDER);

else

    path = ['Inputs/',inputFile,'/average_sedrates.txt'];
    fileID = fopen(path,'r');
    AA = textscan(fileID,'%s %f');
    fclose(fileID);
    record_names = AA{1};
    ave_sedrates = AA{2};

    L = length(record_names);

    data = struct('name',cell(L,1),'depth',cell(L,1),'C14',cell(L,1),'R',cell(L,1));

    LENGTH = zeros(L,1);

    for ll = 1:L
        data(ll).name = record_names{ll};

        path = ['Inputs/',inputFile,'/records/',record_names{ll},'.txt'];
        AA = readtable(path);
        data(ll).depth = AA.depth;
        data(ll).C14 = [AA.age,AA.error,AA.dR,AA.dSTD,AA.cc];

        [~,order] = sort(data(ll).depth,'ascend');
        data(ll).depth = data(ll).depth(order);
        data(ll).C14 = data(ll).C14(order,:);

        data(ll).R = 1./ave_sedrates(ll);

        LENGTH(ll) = max(data(ll).depth.*data(ll).R);
    end
    [~,ORDER] = sort(LENGTH,'descend');
    data = data(ORDER);


    target = struct('A',cell(1,1),'MU',cell(1,1),'PMU',cell(1,1),'SIG',cell(1,1),'PSIG',cell(1,1));


    path = ['Inputs/',inputFile,'/calibration.txt'];
    AA = readtable(path);

    target.A = AA.Var1/1000;
    target.MU = AA.Var2/1000;
    target.SIG = AA.Var3/1000;
    

    N = size(target.A,1);
    target.PMU = zeros(N,1);
    target.PMU(1) = (target.MU(2)-target.MU(1))./(target.A(2)-target.A(1));
    target.PMU(N) = (target.MU(N)-target.MU(N-1))./(target.A(N)-target.A(N-1));
    target.PMU(2:N-1) = (target.MU(3:N)-target.MU(2:N-1))./(target.A(3:N)-target.A(2:N-1)) + (target.MU(2:N-1)-target.MU(1:N-2))./(target.A(2:N-1)-target.A(1:N-2));
    target.PMU(2:N-1) = target.PMU(2:N-1)/2;
    target.PSIG = zeros(N,1);
    target.PSIG(1) = (target.SIG(2)-target.SIG(1))./(target.A(2)-target.A(1));
    target.PSIG(N) = (target.SIG(N)-target.SIG(N-1))./(target.A(N)-target.A(N-1));
    target.PSIG(2:N-1) = (target.SIG(3:N)-target.SIG(2:N-1))./(target.A(3:N)-target.A(2:N-1)) + (target.SIG(2:N-1)-target.SIG(1:N-2))./(target.A(2:N-1)-target.A(1:N-2));
    target.PSIG(2:N-1) = target.PSIG(2:N-1)/2;

    MIN = 0;
    MAX = max(abs(target.A));


    param = struct();

    param.BETA = BETA;
    param.RHO = RHO;

    SAMPLES = struct('name',cell(L,1));
    for ll = 1:L
        SAMPLES(ll).name = data(ll).name;

        SAMPLES(ll).D = (MIN+setting.interval/2:setting.interval:MAX+1.5*setting.interval-1e-24)';

        ID = (SAMPLES(ll).D<max(data(ll).depth*data(ll).R)+1.5*setting.interval);
        SAMPLES(ll).D = SAMPLES(ll).D(ID,:);
        SAMPLES(ll).D0 = SAMPLES(ll).D - setting.interval/2;
        SAMPLES(ll).D0 = SAMPLES(ll).D0(1:end-1);

        N0 = size(SAMPLES(ll).D,1);
        SAMPLES(ll).F = ones(N0,setting.nSamples,setting.K)./sqrt(setting.K);
        SAMPLES(ll).F = SAMPLES(ll).F + normrnd(0,0.01,[N0,setting.nSamples,setting.K]);

        N = size(data(ll).depth,1);
        SAMPLES(ll).AGE = [];
        SAMPLES(ll).AGE_D = [];

        D1 = data(ll).depth.*data(ll).R;
        D1 = sort(D1,'ascend');

        AP = 0:setting.interval/setting.nSamples:(setting.interval-setting.interval/setting.nSamples);

        SAMPLES(ll).D = SAMPLES(ll).D - AP;

        HH = zeros(N,N0,setting.nSamples);
        for k = 1:setting.nSamples
            for n = 1:N
                CNT = floor((D1(n)+AP(k))/setting.interval);
                HH(n,1:CNT,k) = setting.interval;
                HH(n,CNT+1,k) = D1(n) + AP(k) - CNT.*setting.interval;
                HH(n,1,k) = HH(n,1,k) - AP(k);
            end
        end

        SAM_G = sum(SAMPLES(ll).F.^2,3);

        SAMPLES(ll).BIAS = ones(1,setting.nSamples) + normrnd(0,0.1,[1,setting.nSamples]);

        AA = SAMPLES(ll).BIAS.^2 - 0.1 + squeeze(sum(HH.*reshape(SAM_G,[1,N0,setting.nSamples]),2));

        SAMPLES(ll).AGE_D = [];
        SAMPLES(ll).AGE = AA;
        SAMPLES(ll).ACC_RATE = SAM_G;
        SAMPLES(ll).APP_RATE = [];
    end
end


end