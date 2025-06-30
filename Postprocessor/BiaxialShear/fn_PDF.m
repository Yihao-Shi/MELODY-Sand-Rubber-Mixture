clear;
proportion='30';
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
DD=zeros(36,1);DR=zeros(36,1);RR=zeros(36,1);Aver=zeros(36,1);
strong_DD=zeros(36,1);strong_DR=zeros(36,1);strong_RR=zeros(36,1);
for timestep=1:36%1:size(master,1)
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
    
    DDcount=0;DRcount=0;RRcount=0;
    fc=zeros(size(fnPerSec,1)*size(fnPerSec,2),4);
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
    temp(:,7)=abs(temp(:,3).*temp(:,5)+temp(:,4).*temp(:,6));

    RR_force=[];DR_force=[];DD_force=[];
    for i=1:size(temp,1)
        if temp(i,1)<soft_num&&temp(i,2)<soft_num
            DD_force=[DD_force;temp(i,7)];
        elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
            DR_force=[DR_force;temp(i,7)];
        else
            RR_force=[RR_force;temp(i,7)];
        end
    end
    
    mForce=(sum(DD_force)+sum(DR_force)+sum(RR_force))/(length(DD_force)+length(DR_force)+length(RR_force));
    DDstrong_num=length(find(DD_force>mForce));DRstrong_num=length(find(DR_force>mForce));RRstrong_num=length(find(RR_force>mForce));
    strong_DD(timestep)=DDstrong_num/length(DD_force);DD(timestep)=mean(DD_force);
    strong_DR(timestep)=DRstrong_num/length(DR_force);DR(timestep)=mean(DR_force);
    strong_RR(timestep)=RRstrong_num/length(RR_force);RR(timestep)=mean(RR_force);
    Aver(timestep)=mForce;

    mforce0=12691;
    if timestep==36
        result=HistRate([DD_force;DR_force;RR_force]/mforce0);
        resultDD=HistRate(DD_force/mforce0);
        resultDR=HistRate(DR_force/mforce0);
        resultRR=HistRate(RR_force/mforce0);
    end
end
weyy=Axialstrain(proportion);weyy=weyy([1:10:350,350]);
figure,plot(weyy,strong_DD,'b',weyy,strong_DR,'g',weyy,strong_RR,'r');
figure,plot(weyy,DD,'b',weyy,DR,'g',weyy,RR,'r');
figure,semilogx(result(:,1),result(:,2),'g');hold on;
semilogx(resultDD(:,1),resultDD(:,2),'r');
semilogx(resultDR(:,1),resultDR(:,2),'b');
semilogx(resultRR(:,1),resultRR(:,2),'black');