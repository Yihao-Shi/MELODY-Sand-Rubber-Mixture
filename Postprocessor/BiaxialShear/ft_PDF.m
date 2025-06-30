clear;
proportion='30';
timestep=36;
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
fn=fx.*xnorm+fy.*ynorm;
[masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);

fxPerSec=zeros(particle_num,max(countorNode));
fyPerSec=zeros(particle_num,max(countorNode));
xnormPerSec=zeros(particle_num,max(countorNode));
ynormPerSec=zeros(particle_num,max(countorNode));
for i=1:particle_num
    for j=1:countorNode(i)
        fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
        fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
        xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
        ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
    end
end
ftPerSec=sqrt(abs(fxPerSec.^2+fyPerSec.^2-fnPerSec.^2));

DDcount=0;DRcount=0;RRcount=0;
fc=zeros(size(ftPerSec,1)*size(ftPerSec,2),4);
for i=1:particle_num
    for j=1:countorNode(i)
        if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
            fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,5)=fxPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,6)=fyPerSec(i,j);
        end
    end
end
fc(fc(:,1)==0&fc(:,2)==0,:)=[];

temp=zeros(size(fc,1),11);count=1;
[slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
for i=1:length(slavePar)
    test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
    [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
    for j=1:length(masterPar)
        temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
        temp(count,3)=sum(test(mrows(j):mrows(j+1)-1,3));
        temp(count,4)=sum(test(mrows(j):mrows(j+1)-1,4));
        normalize=sqrt(temp(count,3)^2+temp(count,4)^2);
        temp(count,3)=temp(count,3)/normalize;
        temp(count,4)=temp(count,4)/normalize;
        temp(count,5)=sum(test(mrows(j):mrows(j+1)-1,5));
        temp(count,6)=sum(test(mrows(j):mrows(j+1)-1,6));
        count=count+1;
    end
end
temp(count-1:end,:)=[];
temp(:,7)=abs(-temp(:,3).*temp(:,6)+temp(:,4).*temp(:,5));

RR_force=[];DR_force=[];DD_force=[]);
for i=1:size(fc,1)
    if fc(i,1)<soft_num&&fc(i,2)<soft_num
        DD_force=[DD_force;temp(i,7)];
    elseif fc(i,1)<soft_num&&fc(i,2)>=soft_num||fc(i,1)>=soft_num&&fc(i,2)<soft_num
        DR_force=[DR_force;temp(i,7)];
    else
        RR_force=[RR_force;temp(i,7)];
    end
end

mforce0=899;
if timestep==36
    result=HistRate([DD_force;DR_force;RR_force]/mforce0);
    resultDD=HistRate(DD_force/mforce0);
    resultDR=HistRate(DR_force/mforce0);
    resultRR=HistRate(RR_force/mforce0);
end

figure,semilogy(result(:,1),result(:,2),'g');hold on;
semilogy(resultDD(:,1),resultDD(:,2),'r');
semilogy(resultDR(:,1),resultDR(:,2),'b');
semilogy(resultRR(:,1),resultRR(:,2),'black');