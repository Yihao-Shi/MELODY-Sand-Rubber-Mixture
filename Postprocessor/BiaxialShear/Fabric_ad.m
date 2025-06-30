function [adDD,adDR,adRR,ad]=Fabric_ad(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac)
Total=[DD;DR;RR];
delta=[1,0;0,1];
gamma_d=zeros(2,2);gammaDD_d=zeros(2,2);gammaDR_d=zeros(2,2);gammaRR_d=zeros(2,2);

for i=1:length(Total)
    temp_nn=kron(Total(i,1:2),Total(i,1:2).');
    gamma_d=gamma_d+(Total(i,5)*temp_nn)/(1+sum(sum(devia_ac.*temp_nn)));
end
for i=1:length(DD)
    temp_nn=kron(DD(i,1:2),DD(i,1:2).');
    gammaDD_d=gammaDD_d+(DD(i,5)*temp_nn)/(1+sum(sum(devia_acDD.*temp_nn)));
end
for i=1:length(DR)
    temp_nn=kron(DR(i,1:2),DR(i,1:2).');
    gammaDR_d=gammaDR_d+(DR(i,5)*temp_nn)/(1+sum(sum(devia_acDR.*temp_nn)));
end
for i=1:length(RR)
    temp_nn=kron(RR(i,1:2),RR(i,1:2).');
    gammaRR_d=gammaRR_d+(RR(i,5)*temp_nn)/(1+sum(sum(devia_acRR.*temp_nn)));
end

gamma_d=gamma_d./length(Total);
gammaDD_d=gammaDD_d./length(DD);
gammaDR_d=gammaDR_d./length(DR);
gammaRR_d=gammaRR_d./length(RR);

para_d=0.5*(gamma_d(1,1)+gamma_d(2,2));
paraDD_d=0.5*(gammaDD_d(1,1)+gammaDD_d(2,2));
paraDR_d=0.5*(gammaDR_d(1,1)+gammaDR_d(2,2));
paraRR_d=0.5*(gammaRR_d(1,1)+gammaRR_d(2,2));

devia_ad=4*(gamma_d-para_d*delta)/(gamma_d(1,1)+gamma_d(2,2));
devia_adDD=4*(gammaDD_d-paraDD_d*delta)/(gammaDD_d(1,1)+gammaDD_d(2,2));
devia_adDR=4*(gammaDR_d-paraDR_d*delta)/(gammaDR_d(1,1)+gammaDR_d(2,2));
devia_adRR=4*(gammaRR_d-paraRR_d*delta)/(gammaRR_d(1,1)+gammaRR_d(2,2));

ad=sqrt(0.5*sum(sum(devia_ad.*devia_ad)));
adDD=sqrt(0.5*sum(sum(devia_adDD.*devia_adDD)));
adDR=sqrt(0.5*sum(sum(devia_adDR.*devia_adDR)));
adRR=sqrt(0.5*sum(sum(devia_adRR.*devia_adRR)));
end