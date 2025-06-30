function [anDD,anDR,anRR,an]=Fabric_an(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac)
Total=[DD;DR;RR];
delta=[1,0;0,1];
gamma_n=zeros(2,2);gammaDD_n=zeros(2,2);gammaDR_n=zeros(2,2);gammaRR_n=zeros(2,2);

for i=1:length(Total)
    temp_nn=kron(Total(i,1:2),Total(i,1:2).');
    gamma_n=gamma_n+(Total(i,3)*temp_nn)/(1+sum(sum(devia_ac.*temp_nn)));
end
for i=1:length(DD)
    temp_nn=kron(DD(i,1:2),DD(i,1:2).');
    gammaDD_n=gammaDD_n+(DD(i,3)*temp_nn)/(1+sum(sum(devia_acDD.*temp_nn)));
end
for i=1:length(DR)
    temp_nn=kron(DR(i,1:2),DR(i,1:2).');
    gammaDR_n=gammaDR_n+(DR(i,3)*temp_nn)/(1+sum(sum(devia_acDR.*temp_nn)));
end
for i=1:length(RR)
    temp_nn=kron(RR(i,1:2),RR(i,1:2).');
    gammaRR_n=gammaRR_n+(RR(i,3)*temp_nn)/(1+sum(sum(devia_acRR.*temp_nn)));
end

gamma_n=gamma_n./length(Total);
gammaDD_n=gammaDD_n./length(DD);
gammaDR_n=gammaDR_n./length(DR);
gammaRR_n=gammaRR_n./length(RR);

para_n=0.5*(gamma_n(1,1)+gamma_n(2,2));
paraDD_n=0.5*(gammaDD_n(1,1)+gammaDD_n(2,2));
paraDR_n=0.5*(gammaDR_n(1,1)+gammaDR_n(2,2));
paraRR_n=0.5*(gammaRR_n(1,1)+gammaRR_n(2,2));

devia_an=4*(gamma_n-para_n*delta)/(gamma_n(1,1)+gamma_n(2,2));
devia_anDD=4*(gammaDD_n-paraDD_n*delta)/(gammaDD_n(1,1)+gammaDD_n(2,2));
devia_anDR=4*(gammaDR_n-paraDR_n*delta)/(gammaDR_n(1,1)+gammaDR_n(2,2));
devia_anRR=4*(gammaRR_n-paraRR_n*delta)/(gammaRR_n(1,1)+gammaRR_n(2,2));

an=sqrt(0.5*sum(sum(devia_an.*devia_an)));
anDD=sqrt(0.5*sum(sum(devia_anDD.*devia_anDD)));
anDR=sqrt(0.5*sum(sum(devia_anDR.*devia_anDR)));
anRR=sqrt(0.5*sum(sum(devia_anRR.*devia_anRR)));
end