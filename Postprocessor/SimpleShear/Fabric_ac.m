function [acDD,acDR,acRR,ac,devia_acDD,devia_acDR,devia_acRR,devia_ac]=Fabric_ac(DD,DR,RR)
Total=[DD;DR;RR];
delta=[1,0;0,1];
fai=zeros(2,2);faiDD=zeros(2,2);faiDR=zeros(2,2);faiRR=zeros(2,2);

for i=1:length(Total)
    fai=fai+kron(Total(i,1:2),Total(i,1:2).');
end
for i=1:length(DD)
    faiDD=faiDD+kron(DD(i,1:2),DD(i,1:2).');
end
for i=1:length(DR)
    faiDR=faiDR+kron(DR(i,1:2),DR(i,1:2).');
end
for i=1:length(RR)
    faiRR=faiRR+kron(RR(i,1:2),RR(i,1:2).');
end

fai=fai./length(Total);
faiDD=faiDD./length(DD);
faiDR=faiDR./length(DR);
faiRR=faiRR./length(RR);

para=0.5*(fai(1,1)+fai(2,2));
paraDD=0.5*(faiDD(1,1)+faiDD(2,2));
paraDR=0.5*(faiDR(1,1)+faiDR(2,2));
paraRR=0.5*(faiRR(1,1)+faiRR(2,2));

devia_ac=4*(fai-para*delta);
devia_acDD=4*(faiDD-paraDD*delta);
devia_acDR=4*(faiDR-paraDR*delta);
devia_acRR=4*(faiRR-paraRR*delta);

ac=sqrt(0.5*sum(sum(devia_ac.*devia_ac)));
acDD=sqrt(0.5*sum(sum(devia_acDD.*devia_acDD)));
acDR=sqrt(0.5*sum(sum(devia_acDR.*devia_acDR)));
acRR=sqrt(0.5*sum(sum(devia_acRR.*devia_acRR)));
end