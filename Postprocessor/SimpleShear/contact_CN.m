clear
proportion='40';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
master=master(:,3:end-1);
len=len(:,3:end-1);
fx=dlmread([FRACTION 'fx_info.asc'],' ');
fy=dlmread([FRACTION 'fy_info.asc'],' ');
xnorm=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorm=dlmread([FRACTION 'ynorm_info.asc'],' ');
xnorm=xnorm(:,3:end-1);
ynorm=ynorm(:,3:end-1);
fx=fx(:,3:end-1);
fy=fy(:,3:end-1);

%% extract data
nodes=64;
particle_num=size(master,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%%
RRZm=zeros(size(master,1),1);Z=zeros(size(master,1),1);
DRZm=zeros(size(master,1),1);Zm=zeros(size(master,1),1);
DDZm=zeros(size(master,1),1);
for timestep=1:size(master,1)
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnormPerSec]=master_process(particle_num,nodes,master,fn,timestep);
    contact_info=[];
    for i=1:size(masterPerSec,1)
        temp=masterPerSec(i,:);
        temp(find(temp(1,:)==-1))=[];
        end1=i-1;end2=unique(temp);
        for j=1:length(end2)
            contact_info=[contact_info;end1,end2(j),0];
        end
    end
    [~,idx]=sort(contact_info(:,1))
    contact_info=contact_info(idx,:)
    contact_info(find(contact_info(:,1)>=797|contact_info(:,2)>=797),:)=[];
    
    for i=1:size(masterPerSec,1)
       for j=1:size(masterPerSec,2)
           end1=i-1;
           end2=masterPerSec(i,j);
           for k=1:size(contact_info,1)
               if end1==contact_info(k,1)&end2==contact_info(k,2)
                   contact_info(k,3)=contact_info(k,3)+1;
                   break;
               end
           end
       end
    end
    
    num=zeros(particle_num,1);
    RR_num=zeros(rigid_num,1);DR_num=zeros(particle_num,1);DD_num=zeros(soft_num,1);
    for p=1:particle_num
        info=contact_info(find(contact_info(:,1)==p-1),:)
        for i=1:size(info,1)
            num(p)=num(p)+sign(info(i,3));
            if p-1<soft_num&info(i,2)-1<soft_num
                DD_num(p)=DD_num(p)+sign(info(i,3));
            elseif (p-1<soft_num&info(i,2)-1>=soft_num)|(p-1>=soft_num&info(i,2)-1<soft_num)
                DR_num(p)=DR_num(p)+sign(info(i,3));
            elseif p-1>=soft_num&info(i,2)-1>=soft_num
                RR_num(p-soft_num)=RR_num(p-soft_num)+sign(info(i,3));
            end
        end
    end
    Dnum=num(1:soft_num);Rnum=num(soft_num+1:particle_num);

    n0=length(find(num==0));n1=length(find(num==1));
    RR0=length(find(num(1:soft_num)==0));RR1=length(find(num(1:soft_num)==1));
    DD0=length(find(num(soft_num+1:end)==0));DD1=length(find(num(soft_num+1:end)==1));
    RRZm(timestep)=(sum(RR_num)-RR1)/(rigid_num-RR1-RR0);
    DRZm(timestep)=(sum(DR_num)-n1)/(particle_num-n1-n0);
    DDZm(timestep)=(sum(DD_num)-DD1)/(soft_num-DD1-DD0);
    Z(timestep)=sum(num)/particle_num;
    Zm(timestep)=(sum(num)-n1)/(particle_num-n1-n0);
end
[ShearStrain]=Shearstrain(proportion);
figure,plot(ShearStrain,DDZm,'r',ShearStrain,DRZm,'g',ShearStrain,RRZm,'b');
figure,plot(ShearStrain,Z,'r',ShearStrain,Zm,'b');