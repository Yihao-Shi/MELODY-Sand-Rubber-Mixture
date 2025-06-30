clear;
proportion='20';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];

particle_num=5025;
total_num=5029;
rigid_num=ceil((1-str2double(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
dt=1;

[D_strain,D_damp,E_contact,E_boundary,E_k]=E_work(FRACTION);
con=E_contact+D_strain+D_damp-E_k+E_boundary;
contact_work=dlmread([FRACTION 'contactWork_info.asc'],' ');contact_work=contact_work-contact_work(1,:);
damping_work=dlmread([FRACTION 'dampingWork_info.asc'],' ');damping_work=damping_work-damping_work(1,:);
internal_work=dlmread([FRACTION 'InternalWork_info.asc'],' ');internal_work=internal_work-internal_work(1,:);
kinetic=dlmread([FRACTION 'kinetic_info.asc'],' ');kinetic=kinetic-kinetic(1,:);

end_ID=5027;
adamping=sum(damping_work(:,3:end_ID),2);
acontact=sum(contact_work(:,[3:4414,4416:end_ID]),2);
ainternal=sum(internal_work(:,3:end_ID),2);
akinetic=sum(kinetic(:,[3:4414,4416:end_ID]),2);
cons=(adamping+ainternal-akinetic+acontact)./E_boundary;
val=sum(contact_work(:,5028:5031),2)-sum(kinetic(:,5028:5031),2)+E_boundary;

o=zeros(36,5025);
for i=3:5027
    a=contact_work(:,i);b=kinetic(:,i)-damping_work(:,i)-internal_work(:,i);
    o(:,i-2)=(a-b)./a;
end

cw_work=sum(contact_work(:,5028:5031),2);
kp_work=sum(kinetic(:,3:5027),2);
kw_work=sum(kinetic(:,5028:5031),2);
rest_work=E_boundary+E_contact;

source=load(['source' proportion '.mat']);
velxs=dlmread([FRACTION 'xvel_info.asc'],' ');
velys=dlmread([FRACTION 'yvel_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
xpos=dlmread([FRACTION 'posx_info.asc'],' ');
ypos=dlmread([FRACTION 'posy_info.asc'],' ');
rpos=dlmread([FRACTION 'posr_info.asc'],' ');
for i=size(rpos,1):-1:1
    rpos(i,:)=rpos(i,:)-rpos(1,:);
end
DD_k=1000000000;
DR_k=1820000000;
RR_k=10000000000;
DW_k=2000000000;
RW_k=20000000000;
DD_rho=1100e5;RR_rho=2650e5;

% extract data
[~,countorNode]=NodeInfo(FRACTION,total_num,soft_num);

parfor i=1:36
   velx(i,:)=velxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   vely(i,:)=velys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
end

%
energy=zeros(36,14);contact_energy=zeros(36,5029);
energy(:,14)=rest_work;

inixposPerSec=zeros(particle_num);
iniyposPerSec=zeros(particle_num);
inirposPerSec=zeros(particle_num);
inilocxPerSec=zeros(particle_num,max(countorNode));
inilocyPerSec=zeros(particle_num,max(countorNode));
for i=1:total_num
    inixposPerSec(i)=xpos(1,i+2);iniyposPerSec(i)=ypos(1,i+2);inirposPerSec(i)=rpos(1,i+2);
    for j=1:countorNode(i)
        inilocxPerSec(i,j)=locx(1,sum(countorNode(1:i-1))+j);
        inilocyPerSec(i,j)=locy(1,sum(countorNode(1:i-1))+j);
    end
end

for timestep=1:36
    xposPerSec=zeros(particle_num);
    yposPerSec=zeros(particle_num);
    rposPerSec=zeros(particle_num);
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    velxPerSec=zeros(particle_num,max(countorNode));
    velyPerSec=zeros(particle_num,max(countorNode));
    dispxPerSec=zeros(particle_num,max(countorNode));
    dispyPerSec=zeros(particle_num,max(countorNode));
    disprPerSec=zeros(particle_num,max(countorNode));
    for i=1:total_num
        xposPerSec(i)=xpos(timestep,i+2);yposPerSec(i)=ypos(timestep,i+2);rposPerSec(i)=rpos(timestep,i+2);
        for j=1:countorNode(i)
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
            velxPerSec(i,j)=velx(timestep,sum(countorNode(1:i-1))+j);
            velyPerSec(i,j)=vely(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    
    for i=1:soft_num
        kinematic=source.kinematics{timestep,i};
        for j=1:countorNode(i)
            dispxPerSec(i,j)=kinematic(j,2);
            dispyPerSec(i,j)=kinematic(j,3);
            disprPerSec(i,j)=0;
        end
    end

    for i=soft_num+1:total_num
        kinematic=source.kinematics{timestep,i};
        for j=1:countorNode(i)
            dispxPerSec(i,j)=kinematic(1,1);
            dispyPerSec(i,j)=kinematic(1,2);
            disprPerSec(i,j)=kinematic(1,3);
        end
    end
    
    for i=soft_num+1:total_num
        dx=inilocxPerSec(i,:)-inixposPerSec(i);
        dy=inilocyPerSec(i,:)-iniyposPerSec(i);
        dr=disprPerSec(i,:)-inirposPerSec(i);
        dispxPerSec(i,:)=dispxPerSec(i,:)-dx+dx.*cos(dr)-dy.*sin(dr);
        dispyPerSec(i,:)=dispyPerSec(i,:)-dy+dx.*sin(dr)+dy.*cos(dr);
    end

    fc=zeros(total_num*167,50);
    count=zeros(total_num,1);
    for i=1:total_num
        contact=source.contacts{timestep,i};
        for j=1:size(contact,1)
            slave_node=contact(j,2);master_segement=contact(j,7);
            if contact(j,8)<0&&(i-1<particle_num||contact(j,4)<particle_num)&&i-1~=contact(j,4)&&slave_node<countorNode(i)&&i-1~=4412&&contact(j,4)~=4412
                count(i)=count(i)+1;
                cur_node=(i-1)*167+count(i);

                fc(cur_node,1)=i-1;fc(cur_node,2)=contact(j,4);                % particle ID
                fc(cur_node,3)=slave_node;fc(cur_node,4)=master_segement;
                fc(cur_node,5)=contact(j,8);fc(cur_node,6)=contact(j,9);       % gap

                if fc(cur_node,1)<soft_num&&fc(cur_node,2)<soft_num
                    k=DD_k;
                    eta=0;
                    miu=0.5;
                elseif fc(cur_node,1)<soft_num&&fc(cur_node,2)>=soft_num&&fc(cur_node,2)<particle_num ...
                       ||fc(cur_node,1)>=soft_num&&fc(cur_node,1)<particle_num&&fc(cur_node,2)<soft_num
                    k=DR_k;
                    eta=0;
                    miu=0.5;
                elseif fc(cur_node,1)>=soft_num&&fc(cur_node,1)<particle_num&&fc(cur_node,2)>=soft_num&&fc(cur_node,2)<particle_num
                    k=RR_k;
                    eta=1;
                    miu=0.5;
                elseif fc(cur_node,1)>=soft_num&&fc(cur_node,1)<particle_num&&fc(cur_node,2)>=particle_num|| ...
                       fc(cur_node,1)>=particle_num&&fc(cur_node,2)>=soft_num&&fc(cur_node,2)<particle_num
                    k=RW_k;
                    eta=1;
                    miu=0;
                elseif fc(cur_node,1)<soft_num&&fc(cur_node,2)>=particle_num ...
                       ||fc(cur_node,1)>=particle_num&&fc(cur_node,2)<soft_num
                    k=DW_k;
                    eta=0;
                    miu=0;
                else
                    k=0;
                    eta=0;
                    miu=0;
                end
                fc(cur_node,7)=k;
                fc(cur_node,51)=eta;
                fc(cur_node,52)=miu;
                
                fc(cur_node,8)=contact(j,11);fc(cur_node,9)=contact(j,12);                        % norm vector
                fc(cur_node,10)=contact(j,13);fc(cur_node,11)=contact(j,14);                      % tangent vector
                fc(cur_node,12)=contact(j,15);                                                    % length
                fc(cur_node,13)=contact(j,16);fc(cur_node,14)=contact(j,17);                      % contact force
                fc(cur_node,15)=fc(cur_node,13)*fc(cur_node,8)+fc(cur_node,14)*fc(cur_node,9);    % n-fc
                fc(cur_node,16)=fc(cur_node,13)*fc(cur_node,10)+fc(cur_node,14)*fc(cur_node,11);  % t-fc
                if i>soft_num
                    dx=locxPerSec(i,contact(j,18)+1)-xposPerSec(i);
                    dy=locyPerSec(i,contact(j,18)+1)-yposPerSec(i);
                    fc(cur_node,17)=-fc(cur_node,13)*dy+fc(cur_node,14)*dx;                       % r-fc
                end

                fc(cur_node,18)=contact(j,18);fc(cur_node,19)=contact(j,19);                      % master node ID/shape function
                fc(cur_node,20)=contact(j,21);fc(cur_node,21)=contact(j,22);
                fc(cur_node,22)=contact(j,24);fc(cur_node,23)=contact(j,25);
                fc(cur_node,24)=contact(j,27);fc(cur_node,25)=contact(j,28);
                fc(cur_node,26)=contact(j,30);                                                    % effective mass

                fc(cur_node,27)=velxPerSec(i,slave_node+1);                                       % x-velocity
                fc(cur_node,28)=velxPerSec(fc(cur_node,2)+1,fc(cur_node,18)+1);
                fc(cur_node,29)=velxPerSec(fc(cur_node,2)+1,fc(cur_node,20)+1);
                fc(cur_node,30)=velxPerSec(fc(cur_node,2)+1,fc(cur_node,22)+1);
                fc(cur_node,31)=velxPerSec(fc(cur_node,2)+1,fc(cur_node,24)+1);

                fc(cur_node,32)=velyPerSec(i,slave_node+1);                                       % y-velocity
                fc(cur_node,33)=velyPerSec(fc(cur_node,2)+1,fc(cur_node,18)+1);
                fc(cur_node,34)=velyPerSec(fc(cur_node,2)+1,fc(cur_node,20)+1);
                fc(cur_node,35)=velyPerSec(fc(cur_node,2)+1,fc(cur_node,22)+1);
                fc(cur_node,36)=velyPerSec(fc(cur_node,2)+1,fc(cur_node,24)+1);

                fc(cur_node,37)=fc(cur_node,27)-(fc(cur_node,19)*fc(cur_node,28) +...             % x-rel_vel
                                                 fc(cur_node,21)*fc(cur_node,29) +...
                                                 fc(cur_node,23)*fc(cur_node,30) +...
                                                 fc(cur_node,25)*fc(cur_node,31));
                fc(cur_node,38)=fc(cur_node,32)-(fc(cur_node,19)*fc(cur_node,33) +...             % y-rel_vel
                                                 fc(cur_node,21)*fc(cur_node,34) +...
                                                 fc(cur_node,23)*fc(cur_node,35) +...
                                                 fc(cur_node,25)*fc(cur_node,36));
                fc(cur_node,39)=fc(cur_node,37)*fc(cur_node,8)+fc(cur_node,38)*fc(cur_node,9);    % n-rel_vel
                fc(cur_node,40)=fc(cur_node,37)*fc(cur_node,10)+fc(cur_node,38)*fc(cur_node,11);  % t-rel_vel
                
                fc(cur_node,53)=dispxPerSec(i,slave_node+1);                                          % x-displacement
                master_dispx1=dispxPerSec(fc(cur_node,2)+1,fc(cur_node,18)+1);
                master_dispx2=dispxPerSec(fc(cur_node,2)+1,fc(cur_node,20)+1);
                master_dispx3=dispxPerSec(fc(cur_node,2)+1,fc(cur_node,22)+1);
                master_dispx4=dispxPerSec(fc(cur_node,2)+1,fc(cur_node,24)+1);

                fc(cur_node,54)=dispyPerSec(i,slave_node+1);                                          % y-displacement
                master_dispy1=dispyPerSec(fc(cur_node,2)+1,fc(cur_node,18)+1);
                master_dispy2=dispyPerSec(fc(cur_node,2)+1,fc(cur_node,20)+1);
                master_dispy3=dispyPerSec(fc(cur_node,2)+1,fc(cur_node,22)+1);
                master_dispy4=dispyPerSec(fc(cur_node,2)+1,fc(cur_node,24)+1);
                               
                rel_dispx=fc(cur_node,53)-(fc(cur_node,19)*master_dispx1 +...                       % x-rel_disp
                                       fc(cur_node,21)*master_dispx2 +...
                                       fc(cur_node,23)*master_dispx3 +...
                                       fc(cur_node,25)*master_dispx4);
                rel_dispy=fc(cur_node,54)-(fc(cur_node,19)*master_dispy1 +...                       % y-rel_disp
                                       fc(cur_node,21)*master_dispy2 +...
                                       fc(cur_node,23)*master_dispy3 +...
                                       fc(cur_node,25)*master_dispy4);
                fc(cur_node,41)=rel_dispx*fc(cur_node,8)+rel_dispy*fc(cur_node,9);              % n-rel_disp
                fc(cur_node,42)=rel_dispx*fc(cur_node,10)+rel_dispy*fc(cur_node,11);            % t-rel_disp

                fc(cur_node,46)=-2*fc(cur_node,51)*fc(cur_node,39)*sqrt(fc(cur_node,26)*k*fc(cur_node,12));
                fc(cur_node,47)=-2*fc(cur_node,51)*fc(cur_node,40)*sqrt(fc(cur_node,26)*k*fc(cur_node,12));
                
                if fc(cur_node,15)>0
                    fn=abs(fc(cur_node,15)-fc(cur_node,46));
                    ft=abs(fc(cur_node,16)-fc(cur_node,47));
                    fc(cur_node,43)=0.5*(fn^2+ft^2)/fc(cur_node,7)/fc(cur_node,12);    % dE(^e_c)
                    fc(cur_node,44)=0.5*fn^2/fc(cur_node,7)/fc(cur_node,12);
                    fc(cur_node,45)=0.5*ft^2/fc(cur_node,7)/fc(cur_node,12);
                    fc(i,49)=fc(i,39)*fc(i,46);
                else
                    fc(cur_node,46)=0;fc(cur_node,47)=0;
                end
            end
        end
    end
    
    if timestep>1
        for i=1:size(fc,1)
            if fc(i,15)>0&&fc(i,16)>0
                row=fc(i,1)*167+1;
                for j=row:row+count(fc(i,1)+1,1)
                    if fcIni(row,1)==fc(i,1)&&fcIni(row,2)==fc(i,2)&&fcIni(row,3)==fc(i,3)&&fcIni(row,4)==fc(i,4)
                        row=j;break;
                    end
                end
                
                if j~=row
                    fc(i,55)=(fc(i,16)-fc(i,47))/fc(i,7)/fc(i,12);
                    if abs(-fc(i,7)*fc(i,12)*fc(i,40)+fc(i,47))>fc(i,52)*abs(fc(i,15))
                        p_delta=abs(fc(i,40)-fc(i,55));
                        fc(i,55)=0;
                    else
                        p_delta=0;
                    end
                else
                    fc(i,55)=(fc(i,16)-fc(i,47)-fcIni(row,16))/fc(i,7)/fc(i,12);
                    if abs(-fc(i,7)*fc(i,12)*(fc(i,40)+fcIni(row,6))+fc(i,47))>fc(i,52)*abs(fc(i,15))
                        p_delta=abs(fc(i,40)-fc(i,55));
                        fc(i,55)=0;
                    else
                        p_delta=0;
                    end
                end
                fc(i,50)=abs(fc(i,16))*p_delta;
                fc(i,49)=fc(i,49)+fc(i,55)*fc(i,47);
            end
        end 
    end
    fcIni=fc;
    fc(fc(:,15)==0&fc(:,16)==0,:)=[];
    
%     temp=zeros(size(fc,1),3);count=1;
%     [slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
%     for i=1:length(slavePar)
%         test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
%         [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
%         for j=1:length(masterPar)
%             temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
%             temp(count,3)=sum(test(mrows(j):mrows(j+1)-1,50));
%             count=count+1;
%         end
%     end
%     temp(temp(:,1)==0&temp(:,2)==0,:)=[];
%     
%     RR_slip=[];DR_slip=[];DD_slip=[];
%     for i=1:size(temp,1)
%         if ismissing(temp(i,3))==0&&temp(i,3)>0
%             if temp(i,1)<soft_num&&temp(i,2)<soft_num
%                 DD_slip=[DD_slip;temp(i,3)];
%             elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
%                 DR_slip=[DR_slip;temp(i,3)];
%             elseif temp(i,1)>=soft_num&&temp(i,2)>=soft_num
%                 RR_slip=[RR_slip;temp(i,3)];
%             end
%         end
%     end
%     temp(isnan(temp(:,3))|temp(:,3)==0,:)=[];
    
    if timestep>1
%         if ~isempty(DD_slip)
%             energy(timestep,1)=energy(timestep-1,1)+mean(DD_slip);
%         end
%         if ~isempty(DR_slip)
%             energy(timestep,2)=energy(timestep-1,2)+mean(DR_slip);
%         end
%         if ~isempty(RR_slip)
%             energy(timestep,3)=energy(timestep-1,3)+mean(RR_slip);
%         end
%         if ~isempty(temp(:,3))
%             energy(timestep,4)=energy(timestep-1,4)+sum(temp(:,3));
%         end
        energy(timestep,8)=energy(timestep-1,8)+sum(fc(:,49));
    end
    energy(timestep,5)=sum(fc(:,43));energy(timestep,6)=sum(fc(:,44));energy(timestep,7)=sum(fc(:,45));
    energy(timestep,9)=D_damp(timestep);
    energy(timestep,10)=D_strain(timestep);
    energy(timestep,11)=E_k(timestep);
    energy(timestep,12)=E_boundary(timestep);
end
energy(:,5)=energy(:,5)-energy(1,5);energy(:,6)=energy(:,6)-energy(1,6);energy(:,7)=energy(:,7)-energy(1,7);
energy(:,5)=0.5*energy(:,5);energy(:,6)=0.5*energy(:,6);energy(:,7)=0.5*energy(:,7);
energy(:,4)=E_boundary-acontact-energy(:,8)-energy(:,5);