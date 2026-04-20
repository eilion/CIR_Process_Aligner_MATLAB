function [data,SAMPLES] = learnRR(data,SAMPLES,target,setting)

Q = length(data);

for q = 1:Q
    data(q).R = median((SAMPLES(q).AGE(end,:)-SAMPLES(q).AGE(1,:))./(data(q).depth(end)-data(q).depth(1)));
end

MIN = 0;
MAX = max(abs(target.A));

for ll = 1:Q
    SAMPLES(ll).D = (MIN+setting.interval/2:setting.interval:MAX+setting.interval/2-1e-24)';

    ID = (SAMPLES(ll).D<max(data(ll).depth*data(ll).R)+setting.interval/2);
    SAMPLES(ll).D = SAMPLES(ll).D(ID,:);

    N0 = size(SAMPLES(ll).D,1);
    SAMPLES(ll).F = ones(N0,setting.nSamples,setting.K)./sqrt(setting.K);
    SAMPLES(ll).F = SAMPLES(ll).F + normrnd(0,0.01,[N0,setting.nSamples,setting.K]);

    SAMPLES(ll).BIAS = ones(1,setting.nSamples) + normrnd(0,0.1,[1,setting.nSamples]);
end


end