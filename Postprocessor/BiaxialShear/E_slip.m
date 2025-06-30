function [DD_s,DR_s,RR_s,E_s,aniso]=E_slip(proportion,FRACTION)
disp=dlmread([FRACTION 'displacement_info.asc'],' ');
bdisp=zeros(36,11);
for i=1:35 
    bdisp(i+1,:)=sum(disp(1:i*10,:));
end
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
velxs=dlmread([FRACTION 'xvel_info.asc'],' ');
velys=dlmread([FRACTION 'yvel_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
gapts=dlmread([FRACTION 'gapt_info.asc'],' ');
DD_k=1000000000;
DR_k=1820000000;
RR_k=10000000000;

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   velx(i,:)=velxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   vely(i,:)=velys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   gapt(i,:)=gapts((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
energy=zeros(36,8);aniso=cell(1,1);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);
    
    velxPerSec=zeros(particle_num,max(countorNode));
    velyPerSec=zeros(particle_num,max(countorNode));
    fxPerSec=zeros(particle_num,max(countorNode));
    fyPerSec=zeros(particle_num,max(countorNode));
    xnormPerSec=zeros(particle_num,max(countorNode));
    ynormPerSec=zeros(particle_num,max(countorNode));
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    gaptPerSec=zeros(particle_num,max(countorNode));
    lengthPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            velxPerSec(i,j)=velx(timestep,sum(countorNode(1:i-1))+j);
            velyPerSec(i,j)=vely(timestep,sum(countorNode(1:i-1))+j);
            fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
            fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
            xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
            ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
            gaptPerSec(i,j)=gapt(timestep,sum(countorNode(1:i-1))+j);
            lengthPerSec(i,j)=len(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    ftPerSec=abs(-fxPerSec.*ynormPerSec+fyPerSec.*xnormPerSec);

    fc=zeros(size(ftPerSec,1)*size(ftPerSec,2),13);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)~=-1
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                if fnPerSec(i,j)~=0&&ftPerSec(i,j)~=0
                    fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,5)=fnPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,6)=ftPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,7)=locxPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,8)=locyPerSec(i,j);
                    if masterPerSec(i,j)<5025
                        fc(sum(countorNode(1:i-1))+j,9)=mean(locxPerSec(masterPerSec(i,j)+1,masterPerSec(masterPerSec(i,j)+1,:)==i-1));
                        fc(sum(countorNode(1:i-1))+j,10)=mean(locyPerSec(masterPerSec(i,j)+1,masterPerSec(masterPerSec(i,j)+1,:)==i-1));
                    else
                        fc(sum(countorNode(1:i-1))+j,9)=bdisp(timestep,2*(masterPerSec(i,j)-5025)+2);
                        fc(sum(countorNode(1:i-1))+j,10)=bdisp(timestep,2*(masterPerSec(i,j)-5025)+3);
                    end
                    fc(sum(countorNode(1:i-1))+j,11)=gaptPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,13)=lengthPerSec(i,j);
                end
            end
        end
    end
    
    for i=1:size(fc,1)
        if timestep>1
            if fc(i,6)/fc(i,5)>0.499&&fc(i,6)/fc(i,5)<0.501
                d=abs(-((fc(i,7)-fcIni(i,7))-(fc(i,9)-fcIni(i,9)))*fc(i,4) ...
                    +((fc(i,8)-fcIni(i,8))-(fc(i,10)-fcIni(i,10)))*fc(i,3));
%                 if d-abs(fc(i,11)-fcIni(i,11))>0
                if fc(i,1)<soft_num&&fc(i,2)<soft_num
                    k=DD_k;
                elseif fc(i,1)<soft_num&&fc(i,2)>=soft_num||fc(i,1)>=soft_num&&fc(i,2)<soft_num
                    k=DR_k;
                elseif fc(i,1)>=soft_num&&fc(i,2)>=soft_num
                    k=RR_k;
                end
                fc(i,12)=0.5*abs(fcIni(i,6)+fc(i,6))*(d-(fc(i,6)-fcIni(i,6))/k/fc(i,13));
%                 end
            end
        end 
    end
    fcIni=fc;
    
    fc(fc(:,3)==0&fc(:,4)==0,:)=[];
    temp=zeros(size(fc,1),6);count=1;
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
            temp(count,5)=sum(test(mrows(j):mrows(j+1)-1,12));
            temp(count,6)=mean(test(mrows(j):mrows(j+1)-1,7));
            temp(count,7)=mean(test(mrows(j):mrows(j+1)-1,8));
            count=count+1;
        end
    end
    temp(temp(:,1)==0&temp(:,2)==0,:)=[];
    
    RR_slip=[];DR_slip=[];DD_slip=[];
    for i=1:size(temp,1)
        if ismissing(temp(i,5))==0&&temp(i,5)>0
            if temp(i,1)<soft_num&&temp(i,2)<soft_num
                DD_slip=[DD_slip;temp(i,5)];
            elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
                DR_slip=[DR_slip;temp(i,5)];
            elseif temp(i,1)>=soft_num&&temp(i,2)>=soft_num
                RR_slip=[RR_slip;temp(i,5)];
            end
        end
    end
    temp(isnan(temp(:,5))|temp(:,5)==0,:)=[];
    aniso{timestep,1}=temp;
    
    if ~isempty(DD_slip)
        energy(timestep,1)=sum(DD_slip);energy(timestep,5)=mean(DD_slip);
    end
    if ~isempty(DR_slip)
        energy(timestep,2)=sum(DR_slip);energy(timestep,6)=mean(DR_slip);
    end
    if ~isempty(RR_slip)
        energy(timestep,3)=sum(RR_slip);energy(timestep,7)=mean(RR_slip);
    end
    energy(timestep,4)=energy(timestep,1)+energy(timestep,2)+energy(timestep,3);
    if ~isempty(temp(:,5))
        energy(timestep,8)=mean(temp(:,5));
    end
end
DD_s=zeros(36,1);DR_s=zeros(36,1);RR_s=zeros(36,1);E_s=zeros(36,1);
for i=1:36
    DD_s(i,1)=sum(energy(1:i,1));
    DR_s(i,1)=sum(energy(1:i,2));
    RR_s(i,1)=sum(energy(1:i,3));
    E_s(i,1)=sum(energy(1:i,4));
    DD_s(i,2)=sum(energy(1:i,5));
    DR_s(i,2)=sum(energy(1:i,6));
    RR_s(i,2)=sum(energy(1:i,7));
    E_s(i,2)=sum(energy(1:i,8));
end
end