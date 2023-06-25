function flap_rate=FlapFromAcc(accx,accy,accz)
params=ReadNNParams();
[inc,raccx,raccy,raccz]=Acc1HzToFeature(accx,accy,accz);
accstat=GetAccStats(raccx,raccy,raccz);
naccstat=(accstat-params.ss_mean)./params.ss_std;
flap_rate=nan(1,size(accx,1));
flap_rate(inc)=max(0,max(0,naccstat*params.w1+params.w1_)*params.w2+params.w2_);
end

function accstat=GetAccStats(accx,accy,accz)
GetAccStat=@(acc) [mean(acc);var(acc,1);min(acc);max(acc)];
accstat=[GetAccStat(accx);GetAccStat(accy);GetAccStat(accz)]';
end

function [inc,raccx,raccy,raccz]=Acc1HzToFeature(accx,accy,accz)
inc=find(~isnan(accx(1:end-4)) & ~isnan(accy(1:end-4)) & ~isnan(accz(1:end-4)));
for i=1:4
    inc=inc(~isnan(accx(inc+i)) & ~isnan(accy(inc+i)) & ~isnan(accz(inc+i)));
end
raccx=[accx(inc),accx(inc+1),accx(inc+2),accx(inc+3),   accx(inc+4)]';
raccy=[accy(inc),accy(inc+1),accy(inc+2),accy(inc+3),accy(inc+4)]';
raccz=[accz(inc),accz(inc+1),accz(inc+2),accz(inc+3),accz(inc+4)]';
inc=inc+2;
end

function params=ReadNNParams()
a=xlsread('fit_params.csv');
params.w1=reshape(a(2,:),a(1,2),a(1,1))';
params.w1_=a(4,1:a(3,1));
params.w2=a(6,1:a(5,1))';
params.w2_=a(8,1);
params.ss_mean=a(9,1:12);
params.ss_std=a(10,1:12);
end