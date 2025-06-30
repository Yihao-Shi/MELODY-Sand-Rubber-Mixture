clear
proportion='0';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];

particle_num=5025;
total_num=5029;
rigid_num=ceil((1-str2double(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

contacts=cell(1,1);kinematics=cell(1,1);forces=cell(1,1);
for timestep=26:26
    filenum=num2str(timestep-1);
    if timestep<11
        data=fopen([FRACTION 'DYNAMIC_0000' filenum '.asc']);
    else
        data=fopen([FRACTION 'DYNAMIC_000' filenum '.asc']);
    end
    
    scanned_particle=0;
    while ~feof(data)
        tline=fgets(data);
%         if strfind(tline,'KINEMATICS')~=0
%             if scanned_particle<soft_num
%                 particle_num=str2double(fgets(data));
%                 data_row=[];
%                 for i=1:particle_num
%                     str_row=fgets(data);
%                     data_row=[data_row;str2num(char(regexp(str_row,'\s+','split')))'];
%                 end
%                 data_row=data_row(:,[1,4:7]);
%                 scanned_particle=scanned_particle+1;
%                 kinematics{timestep, scanned_particle}=data_row;
%             else
%                 data_row=[];tline=fgets(data);
%                 for i=1:2
%                     str_row=fgets(data);
%                     data_row=[data_row;str2num(char(regexp(str_row,'\s+','split')))'];
%                 end
%                 scanned_particle=scanned_particle+1;
%                 kinematics{timestep, scanned_particle}=data_row;
%             end
%         end
        
%         if strfind(tline,'FORCES')~=0
%             if scanned_particle<soft_num
%                 particle_num=str2double(fgets(data));
%                 data_row=[];
%                 for i=1:particle_num
%                     str_row=fgets(data);
%                     data_row=[data_row;str2num(char(regexp(str_row,'\s+','split')))'];
%                 end
%                 data_row=data_row(:,[1,4:5]);
%                 forces{timestep, scanned_particle}=data_row;
%             else
%                 data_row=[];
%                 str_row=fgets(data);
%                 data_row=[data_row;str2num(char(regexp(str_row,'\s+','split')))'];
%                 forces{timestep, scanned_particle}=data_row;
%             end
%         end
            
        if strfind(tline,'CONTACTS')~=0
            contact_num=str2double(fgets(data));
            
            data_row=[];
            for i=1:contact_num
                str_row=fgets(data);
                data_row=[data_row;str2num(char(regexp(str_row,'\s+','split')))'];
            end
            data_row(data_row(:,8)>0,:)=[];
            scanned_particle=scanned_particle+1;
            contacts{timestep, scanned_particle}=data_row;
        end
    end
    fclose(data);
end
for timestep=36:-1:1
    for i=1:total_num
        kinematics{timestep,i}=kinematics{timestep,i}-kinematics{1,i};
    end
end