function [data,SAMPLES] = learnParam(data,SAMPLES,target,setting)

A2 = setting.a_d18O;
B2 = setting.b_d18O;

% average accumulation rates:

if strcmp(setting.islearn_average_sedrate,'Yes')
    data.R = median((SAMPLES.AGE(end,:)-SAMPLES.AGE(1,:))./(data.depth(end)-data.depth(1)));
end

% shift and scale parameters:
if strcmp(setting.islearn_shift,'Yes')||strcmp(setting.islearn_scale,'Yes')
    A = target.stack.A;
    MU = target.stack.MU;
    SIG = target.stack.SIG;

    beta1 = 0.9;
    beta2 = 0.999;
    epsilon = 1e-8;
    gamma = 1e-3;

    AGE = SAMPLES.AGE;
    YY = data.d18O;
    if ~isempty(YY)
        ID_Y = YY(:,end);
        YY = YY(:,1:end-1);
    else
        ID_Y = [];
    end

    if ~isempty(ID_Y)
        Mw = 0;
        Vw = 0;

        AGE = AGE(ID_Y,:);

        ID = (AGE>=min(A))&(AGE<=max(A));
        AGE = AGE(ID);
        YY = YY.*ones(1,size(ID,2));
        YY = YY(ID);

        mu = interp1(A,MU,AGE);
        sig = interp1(A,SIG,AGE);

        CC = data.scale;
        HH = data.shift;

        for r = 1:10000
            PDEV = zeros(2,1);

            if isinf(A2)&&isinf(B2)
                QQ = (YY-CC.*mu-HH)./(CC.^2.*sig.^2);
            else
                QQ = (2.*A2+1)*(YY-CC.*mu-HH)./(2.*B2.*CC.^2*sig.^2+(YY-CC.*mu-HH).^2);
            end
            PDEV(1) = sum(QQ.*(YY-HH)./CC-1./CC,"all")./size(AGE,2);
            PDEV(2) = sum(QQ,"all")./size(AGE,2);

            Mw = beta1*Mw - (1-beta1).*PDEV;
            Vw = beta2*Vw + (1-beta2).*PDEV.*PDEV;

            CC = CC - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon);
            HH = HH - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon);
        end

        if strcmp(setting.islearn_scale,'Yes')
            data.scale = CC;
        end

        if strcmp(setting.islearn_shift,'Yes')
            data.shift = HH;
        end
    end
end


end