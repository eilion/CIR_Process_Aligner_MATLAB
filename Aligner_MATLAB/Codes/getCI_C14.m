function [CI] = getCI_C14(data,target,setting)

cal_curve_0 = target.calibration;

L = length(data);
CI = struct('name',cell(L,1),'depth',cell(L,1),'lower',cell(L,1),'median',cell(L,1),'upper',cell(L,1),'samples',cell(L,1));

if L == 1
    a = 3;
    b = 4;

    CI.name = data.name;

    if ~isempty(data.C14)
        cal_curve = cal_curve_0(data.C14(1,end-1));

        AA = cal_curve.A;
        MU = cal_curve.MU;
        SIG = cal_curve.SIG;

        MU_1 = MU(1);
        MU_end = MU(end);

        M = size(data.C14,1);

        depth = data.depth(data.C14(:,end));

        lower = zeros(M,2);
        MEDIAN = zeros(M,2);
        upper = zeros(M,2);

        SAMPLES = zeros(M,1000);

        for m = 1:M
            Table = data.C14(m,1:end-1);

            rand_seed = log(rand(100,10000));

            c14_age = min(max(Table(1)-Table(3),MU_1),MU_end);
            x = interp1q(MU,AA,c14_age)*ones(1,10000);

            for iters = 1:100
                c14_x = min(max(x,AA(1)),AA(end));
                det_x = interp1q(AA,MU,c14_x')';
                err_x = interp1q(AA,SIG,c14_x')';
                RR_x = - (a+0.5)*log(2*b+(det_x+Table(3)-Table(1)).^2./(err_x.^2+Table(2)^2+Table(4)^2)) - 0.5*log(err_x.^2+Table(2)^2+Table(4)^2);

                z = normrnd(x,0.3);
                c14_z = min(max(z,AA(1)),AA(end));
                det_z = interp1q(AA,MU,c14_z')';
                err_z = interp1q(AA,SIG,c14_z')';
                RR_z = - (a+0.5)*log(2*b+(det_z+Table(3)-Table(1)).^2./(err_z.^2+Table(2)^2+Table(4)^2)) - 0.5*log(err_z.^2+Table(2)^2+Table(4)^2);

                RR_z(z<AA(1)|z>AA(end)) = -inf;

                index = (RR_z-RR_x>rand_seed(iters,:));
                x(index) = z(index);
            end

            x_samples = x;
            SAMPLES(m,:) = x_samples(1:1000);

            lower(m,1) = quantile(x_samples,0.025);
            lower(m,2) = quantile(x_samples,0.16);
            MEDIAN(m,1) = quantile(x_samples,0.5);
            MEDIAN(m,2) = mean(x_samples);
            upper(m,1) = quantile(x_samples,0.84);
            upper(m,2) = quantile(x_samples,0.975);
        end


        CI.depth = depth;
        CI.lower = lower;
        CI.median = MEDIAN;
        CI.upper = upper;
        CI.samples = SAMPLES;
    end
else
    parfor ll = 1:L
        a = 3;
        b = 4;

        CI(ll).name = data(ll).name;

        if ~isempty(data(ll).C14)
            cal_curve = cal_curve_0(data(ll).C14(1,end-1));

            AA = cal_curve.A;
            MU = cal_curve.MU;
            SIG = cal_curve.SIG;

            MU_1 = MU(1);
            MU_end = MU(end);

            M = size(data(ll).C14,1);

            depth = data(ll).depth(data(ll).C14(:,end));

            lower = zeros(M,2);
            MEDIAN = zeros(M,2);
            upper = zeros(M,2);

            SAMPLES = zeros(M,1000);

            for m = 1:M
                Table = data(ll).C14(m,1:end-1);

                rand_seed = log(rand(100,10000));

                c14_age = min(max(Table(1)-Table(3),MU_1),MU_end);
                x = interp1q(MU,AA,c14_age)*ones(1,10000);

                for iters = 1:100
                    c14_x = min(max(x,AA(1)),AA(end));
                    det_x = interp1q(AA,MU,c14_x')';
                    err_x = interp1q(AA,SIG,c14_x')';
                    RR_x = - (a+0.5)*log(2*b+(det_x+Table(3)-Table(1)).^2./(err_x.^2+Table(2)^2+Table(4)^2)) - 0.5*log(err_x.^2+Table(2)^2+Table(4)^2);

                    z = normrnd(x,0.3);
                    c14_z = min(max(z,AA(1)),AA(end));
                    det_z = interp1q(AA,MU,c14_z')';
                    err_z = interp1q(AA,SIG,c14_z')';
                    RR_z = - (a+0.5)*log(2*b+(det_z+Table(3)-Table(1)).^2./(err_z.^2+Table(2)^2+Table(4)^2)) - 0.5*log(err_z.^2+Table(2)^2+Table(4)^2);

                    RR_z(z<AA(1)|z>AA(end)) = -inf;

                    index = (RR_z-RR_x>rand_seed(iters,:));
                    x(index) = z(index);
                end

                x_samples = x;
                SAMPLES(m,:) = x_samples(1:1000);

                lower(m,1) = quantile(x_samples,0.025);
                lower(m,2) = quantile(x_samples,0.16);
                MEDIAN(m,1) = quantile(x_samples,0.5);
                MEDIAN(m,2) = mean(x_samples);
                upper(m,1) = quantile(x_samples,0.84);
                upper(m,2) = quantile(x_samples,0.975);
            end


            CI(ll).depth = depth;
            CI(ll).lower = lower;
            CI(ll).median = MEDIAN;
            CI(ll).upper = upper;
            CI(ll).samples = SAMPLES;
        end
    end
end

end