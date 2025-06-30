function [DD_strain,DR_strain,RR_strain,E_strain,E_nstrain,E_tstrain]=E_cstrain(particle_num,soft_num,FRACTION)
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
gapns=dlmread([FRACTION 'gapn_info.asc'],' ');
gapts=dlmread([FRACTION 'gapt_info.asc'],' ');
DD_k=1000000000;
DR_k=1820000000;
RR_k=10000000000;

%% extract data
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);
for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   gapn(i,:)=gapns((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   gapt(i,:)=gapts((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
DD_strain=zeros(36,2);DR_strain=zeros(36,2);RR_strain=zeros(36,2);
E_strain=zeros(36,2);E_nstrain=zeros(36,2);E_tstrain=zeros(36,2);
for timestep=1:size(master,1)
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,~]=master_process(particle_num,countorNode,master,fn,timestep);
    
    lengthPerSec=zeros(particle_num,max(countorNode));
    gapnPerSec=zeros(particle_num,max(countorNode));
    gaptPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            gapnPerSec(i,j)=gapn(timestep,sum(countorNode(1:i-1))+j);
            gaptPerSec(i,j)=gapt(timestep,sum(countorNode(1:i-1))+j);
            lengthPerSec(i,j)=len(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    
    fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),8);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,3)=gapnPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,4)=gaptPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,5)=lengthPerSec(i,j);
            end
        end
    end
    fc(fc(:,3)>=0|fc(:,4)>=0,:)=[];
    
    for i=1:size(fc,1)
        if fc(i,1)<soft_num&&fc(i,2)<soft_num
            fc(i,6)=0.5*DD_k*(fc(i,3)^2+fc(i,4)^2)*fc(i,5);
            fc(i,7)=0.5*DD_k*fc(i,3)^2*fc(i,5);
            fc(i,8)=0.5*DD_k*fc(i,4)^2*fc(i,5);
        elseif fc(i,1)<soft_num&&fc(i,2)>=soft_num||fc(i,1)>=soft_num&&fc(i,2)<soft_num
            fc(i,6)=0.5*DR_k*(fc(i,3)^2+fc(i,4)^2)*fc(i,5);
            fc(i,7)=0.5*DR_k*fc(i,3)^2*fc(i,5);
            fc(i,8)=0.5*DR_k*fc(i,4)^2*fc(i,5);
        else
            fc(i,6)=0.5*RR_k*(fc(i,3)^2+fc(i,4)^2)*fc(i,5);
            fc(i,7)=0.5*RR_k*fc(i,3)^2*fc(i,5);
            fc(i,8)=0.5*RR_k*fc(i,4)^2*fc(i,5);
        end
    end
    
    temp=zeros(size(fc,1),5);count=1;
    [slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
    for i=1:length(slavePar)
        test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
        [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
        for j=1:length(masterPar)
            temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
            temp(count,3)=sum(test(mrows(j):mrows(j+1)-1,6));
            temp(count,4)=sum(test(mrows(j):mrows(j+1)-1,7));
            temp(count,5)=sum(test(mrows(j):mrows(j+1)-1,8));
            count=count+1;
        end
    end
    temp(temp(:,1)==0&temp(:,2)==0,:)=[];
    
    RR=0;DR=0;DD=0;
    RRnum=0;DRnum=0;DDnum=0;
    for i=1:size(temp,1)
        if temp(i,1)<soft_num&&temp(i,2)<soft_num
            DD=DD+temp(i,3);
            DDnum=DDnum+1;
        elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
            DR=DR+temp(i,3);
            DRnum=DRnum+1;
        else
            RR=RR+temp(i,3);
            RRnum=RRnum+1;
        end
    end
    DD_strain(timestep,1)=DD/2;DD_strain(timestep,2)=DD/2/DDnum;
    DR_strain(timestep,1)=DR/2;DR_strain(timestep,2)=DR/2/DRnum;
    RR_strain(timestep,1)=RR/2;RR_strain(timestep,2)=RR/2/RRnum;
    E_nstrain(timestep,1)=sum(temp(:,4));E_nstrain(timestep,2)=mean(temp(:,4));
    E_tstrain(timestep,1)=sum(temp(:,5));E_tstrain(timestep,2)=mean(temp(:,5));
    E_strain(timestep,1)=sum(temp(:,3));
    E_strain(timestep,2)=mean(temp(:,3));
end
DD_strain(:,1)=DD_strain(:,1)-DD_strain(1,1);DD_strain(:,2)=DD_strain(:,2)-DD_strain(1,2);
DR_strain(:,1)=DR_strain(:,1)-DR_strain(1,1);DR_strain(:,2)=DR_strain(:,2)-DR_strain(1,2);
RR_strain(:,1)=RR_strain(:,1)-RR_strain(1,1);RR_strain(:,2)=RR_strain(:,2)-RR_strain(1,2);
E_strain(:,1)=E_strain(:,1)-E_strain(1,1);E_strain(:,2)=E_strain(:,2)-E_strain(1,2);
end