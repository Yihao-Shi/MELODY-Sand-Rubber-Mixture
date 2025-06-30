function [gamma_n,gammaDD_n,gammaDR_n,gammaRR_n,atDD,atDR,atRR,at]=Fabric_at(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac)
Total=[DD;DR;RR];
delta=[1,0;0,1];
gamma_n=zeros(2,2);gammaDD_n=zeros(2,2);gammaDR_n=zeros(2,2);gammaRR_n=zeros(2,2);
gamma_t=zeros(2,2);gammaDD_t=zeros(2,2);gammaDR_t=zeros(2,2);gammaRR_t=zeros(2,2);

for i=1:length(Total)
    temp_nn=kron(Total(i,1:2),Total(i,1:2).');
    temp_tn=kron([-Total(i,2),Total(i,1)],Total(i,1:2).');
    gamma_n=gamma_n+(Total(i,3)*temp_nn)/(1+sum(sum(devia_ac.*temp_nn)));
    gamma_t=gamma_t+(Total(i,4)*temp_tn)/(1+sum(sum(devia_ac.*temp_nn)));
end
for i=1:length(DD)
    temp_nn=kron(DD(i,1:2),DD(i,1:2).');
    temp_tn=kron([-DD(i,2),DD(i,1)],DD(i,1:2).');
    gammaDD_n=gammaDD_n+(DD(i,3)*temp_nn)/(1+sum(sum(devia_acDD.*temp_nn)));
    gammaDD_t=gammaDD_t+(DD(i,4)*temp_tn)/(1+sum(sum(devia_acDD.*temp_nn)));
end
for i=1:length(DR)
    temp_nn=kron(DR(i,1:2),DR(i,1:2).');
    temp_tn=kron([-DR(i,2),DR(i,1)],DR(i,1:2).');
    gammaDR_n=gammaDR_n+(DR(i,3)*temp_nn)/(1+sum(sum(devia_acDR.*temp_nn)));
    gammaDR_t=gammaDR_t+(DR(i,4)*temp_tn)/(1+sum(sum(devia_acDR.*temp_nn)));
end
for i=1:length(RR)
    temp_nn=kron(RR(i,1:2),RR(i,1:2).');
    temp_tn=kron([-RR(i,2),RR(i,1)],RR(i,1:2).');
    gammaRR_n=gammaRR_n+(RR(i,3)*temp_nn)/(1+sum(sum(devia_acRR.*temp_nn)));
    gammaRR_t=gammaRR_t+(RR(i,4)*temp_tn)/(1+sum(sum(devia_acRR.*temp_nn)));
end

gamma_n=gamma_n./length(Total);gamma_t=gamma_t./length(Total);
gammaDD_n=gammaDD_n./length(DD);gammaDD_t=gammaDD_t./length(DD);
gammaDR_n=gammaDR_n./length(DR);gammaDR_t=gammaDR_t./length(DR);
gammaRR_n=gammaRR_n./length(RR);gammaRR_t=gammaRR_t./length(RR);

para_t=0.5*(gamma_t(1,1)+gamma_t(2,2));
paraDD_t=0.5*(gammaDD_t(1,1)+gammaDD_t(2,2));
paraDR_t=0.5*(gammaDR_t(1,1)+gammaDR_t(2,2));
paraRR_t=0.5*(gammaRR_t(1,1)+gammaRR_t(2,2));

devia_at=8/3*(gamma_t-para_t*delta)/(gamma_n(1,1)+gamma_n(2,2));
devia_atDD=8/3*(gammaDD_t-paraDD_t*delta)/(gammaDD_n(1,1)+gammaDD_n(2,2));
devia_atDR=8/3*(gammaDR_t-paraDR_t*delta)/(gammaDR_n(1,1)+gammaDR_n(2,2));
devia_atRR=8/3*(gammaRR_t-paraRR_t*delta)/(gammaRR_n(1,1)+gammaRR_n(2,2));

at=sqrt(0.5*sum(sum(devia_at.*devia_at)));
atDD=sqrt(0.5*sum(sum(devia_atDD.*devia_atDD)));
atDR=sqrt(0.5*sum(sum(devia_atDR.*devia_atDR)));
atRR=sqrt(0.5*sum(sum(devia_atRR.*devia_atRR)));
end