clear
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];

particle_num=5025;
total_num=5029;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
%%
[softFieldNode,countorNode]=NodeInfo(FRACTION,total_num,soft_num);
%%
Spies=cell(1,1);

% Spies{1,1}={'kinetic_info'};
% Spies{1,2}={int2str(total_num)};
% Spies{1,3}={'1'};
% for i=1:total_num
%     Spies{1,4}{i,1}={'Energy'};
%     Spies{1,4}{i,2}={'Kinetic'};
%     Spies{1,4}{i,3}={int2str(i-1)};
% end
% 
% Spies{2,1}={'internalWork_info'};
% Spies{2,2}={int2str(total_num)};
% Spies{2,3}={'1'};
% for i=1:total_num
%     Spies{2,4}{i,1}={'Work'};
%     Spies{2,4}{i,2}={'Internal'};
%     Spies{2,4}{i,3}={int2str(i-1)};
% end
% 
% Spies{3,1}={'contactWork_info'};
% Spies{3,2}={int2str(total_num)};
% Spies{3,3}={'1'};
% for i=1:total_num
%     Spies{3,4}{i,1}={'Work'};
%     Spies{3,4}{i,2}={'Contact'};
%     Spies{3,4}{i,3}={int2str(i-1)};
% end
% 
% Spies{4,1}={'damping_info'};
% Spies{4,2}={int2str(total_num)};
% Spies{4,3}={'1'};
% for i=1:total_num
%     Spies{4,4}{i,1}={'Work'};
%     Spies{4,4}{i,2}={'Damping'};
%     Spies{4,4}{i,3}={int2str(i-1)};
% end

Spies{1,1}={'mass_info'};
Spies{1,2}={int2str(particle_num)};
Spies{1,3}={'1'};
for i=1:particle_num
    for j=1:countorNode(i,1)
        Spies{1,4}{sum(countorNode(1:i-1))+j,1}={'MassScaling'};
        Spies{1,4}{sum(countorNode(1:i-1))+j,2}={'Mass'};
        Spies{1,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
        Spies{1,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
    end
end

% Spies{2,1}={'posy_info'};
% Spies{2,2}={int2str(total_num)};
% Spies{2,3}={'1'};
% for i=1:total_num
%     Spies{2,4}{i,1}={'Position'};
%     Spies{2,4}{i,2}={'Y'};
%     Spies{2,4}{i,3}={int2str(i-1)};
%     Spies{2,4}{i,4}={int2str(-1)};
% end

% 
% Spies{3,1}={'locx_info'};
% Spies{3,2}={int2str(sum(countorNode))};
% Spies{3,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{3,4}{sum(countorNode(1:i-1))+j,1}={'Position'};
%         Spies{3,4}{sum(countorNode(1:i-1))+j,2}={'X'};
%         Spies{3,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{3,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{4,1}={'locy_info'};
% Spies{4,2}={int2str(sum(countorNode))};
% Spies{4,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{4,4}{sum(countorNode(1:i-1))+j,1}={'Position'};
%         Spies{4,4}{sum(countorNode(1:i-1))+j,2}={'Y'};
%         Spies{4,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{4,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{5,1}={'xvel_info'};
% Spies{5,2}={int2str(sum(countorNode))};
% Spies{5,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{5,4}{sum(countorNode(1:i-1))+j,1}={'Velocity'};
%         Spies{5,4}{sum(countorNode(1:i-1))+j,2}={'X'};
%         Spies{5,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{5,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{6,1}={'yvel_info'};
% Spies{6,2}={int2str(sum(countorNode))};
% Spies{6,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{6,4}{sum(countorNode(1:i-1))+j,1}={'Velocity'};
%         Spies{6,4}{sum(countorNode(1:i-1))+j,2}={'Y'};
%         Spies{6,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{6,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{7,1}={'rvel_info'};
% Spies{7,2}={int2str(sum(countorNode))};
% Spies{7,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{7,4}{sum(countorNode(1:i-1))+j,1}={'Velocity'};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,2}={'R'};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{8,1}={'dispx_info'};
% Spies{8,2}={int2str(sum(countorNode))};
% Spies{8,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{8,4}{sum(countorNode(1:i-1))+j,1}={'Displacement'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,2}={'X'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{9,1}={'dispy_info'};
% Spies{9,2}={int2str(sum(countorNode))};
% Spies{9,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{9,4}{sum(countorNode(1:i-1))+j,1}={'Displacement'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,2}={'Y'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{10,1}={'dispr_info'};
% Spies{10,2}={int2str(sum(countorNode))};
% Spies{10,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{10,4}{sum(countorNode(1:i-1))+j,1}={'Displacement'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,2}={'R'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end

% Spies{8,1}={'fx_info'};
% Spies{8,2}={int2str(sum(countorNode))};
% Spies{8,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{8,4}{sum(countorNode(1:i-1))+j,1}={'Force'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,2}={'X'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,3}={'Contact'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,4}={int2str(i-1)};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{9,1}={'fy_info'};
% Spies{9,2}={int2str(sum(countorNode))};
% Spies{9,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{9,4}{sum(countorNode(1:i-1))+j,1}={'Force'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,2}={'Y'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,3}={'Contact'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,4}={int2str(i-1)};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{10,1}={'fr_info'};
% Spies{10,2}={int2str(sum(countorNode))};
% Spies{10,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{10,4}{sum(countorNode(1:i-1))+j,1}={'Force'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,2}={'R'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,3}={'Contact'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,4}={int2str(i-1)};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end

% Spies{1,1}={'gapn_info'};
% Spies{1,2}={int2str(sum(countorNode))};
% Spies{1,3}={'1'};
% for i=1:particle_num
%     for j=1:countorNode(i,1)
%         Spies{1,4}{sum(countorNode(1:i-1))+j,1}={'Contact'};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,2}={'Gapn'};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,4}={int2str(0)};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{2,1}={'slave_info'};
% Spies{2,2}={int2str(sum(countorNode))};
% Spies{2,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{2,4}{sum(countorNode(1:i-1))+j,1}={'Contact'};
%         Spies{2,4}{sum(countorNode(1:i-1))+j,2}={'Slave'};
%         Spies{2,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{2,4}{sum(countorNode(1:i-1))+j,4}={int2str(0)};
%         Spies{2,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{7,1}={'fx_info'};
% Spies{7,2}={int2str(sum(countorNode))};
% Spies{7,3}={'1'};
% for i=1:particle_num
%     for j=1:countorNode(i,1)
%         Spies{7,4}{sum(countorNode(1:i-1))+j,1}={'Contact'};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,2}={'Fx'};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,4}={int2str(0)};
%         Spies{7,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{8,1}={'fy_info'};
% Spies{8,2}={int2str(sum(countorNode))};
% Spies{8,3}={'1'};
% for i=1:particle_num
%     for j=1:countorNode(i,1)
%         Spies{8,4}{sum(countorNode(1:i-1))+j,1}={'Contact'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,2}={'Fy'};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,4}={int2str(0)};
%         Spies{8,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{9,1}={'master_info'};
% Spies{9,2}={int2str(sum(countorNode))};
% Spies{9,3}={'1'};
% for i=1:particle_num
%     for j=1:countorNode(i,1)
%         Spies{9,4}{sum(countorNode(1:i-1))+j,1}={'Contact'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,2}={'Master'};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,4}={int2str(0)};
%         Spies{9,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{10,1}={'length_info'};
% Spies{10,2}={int2str(sum(countorNode))};
% Spies{10,3}={'1'};
% for i=1:particle_num
%     for j=1:countorNode(i,1)
%         Spies{10,4}{sum(countorNode(1:i-1))+j,1}={'Contact'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,2}={'Length'};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,4}={int2str(0)};
%         Spies{10,4}{sum(countorNode(1:i-1))+j,5}={int2str(j-1)};
%     end
% end
% 
% Spies{1,1}={'mass_info'};
% Spies{1,2}={int2str(sum(countorNode))};
% Spies{1,3}={'1'};
% for i=1:total_num
%     for j=1:countorNode(i,1)
%         Spies{1,4}{sum(countorNode(1:i-1))+j,1}={'MassScaling'};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,2}={'Mass'};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,3}={int2str(i-1)};
%         Spies{1,4}{sum(countorNode(1:i-1))+j,4}={int2str(j-1)};
%     end
% end
% 
% Spies{10,1}={'jacobian_info'};
% Spies{10,2}={int2str(sum(softFieldNode))};
% Spies{10,3}={'1'};
% for i=1:soft_num
%     for j=1:softFieldNode(i,1)
%         Spies{10,4}{sum(softFieldNode(1:i-1))+j,1}={'Jacobian'};
%         Spies{10,4}{sum(softFieldNode(1:i-1))+j,2}={int2str(i-1)};
%         Spies{10,4}{sum(softFieldNode(1:i-1))+j,3}={int2str(j-1)};
%     end
% end
%%
fid_data = fopen(['Biaxial_Spies_Mass' proportion '.txt'],'w');
fprintf(fid_data,'%0.16g\n',size(Spies,1));
for n=1:size(Spies,1)
    fprintf(fid_data,'%s\n',' ');
    fprintf(fid_data,'%s',cell2mat(Spies{n,1}));
    fprintf(fid_data,'%s',' ');
    fprintf(fid_data,'%s',cell2mat(Spies{n,2}));
    fprintf(fid_data,'%s',' ');
    fprintf(fid_data,'%s\n',cell2mat(Spies{n,3}));
    for i=1:size(Spies{n,4},1)
%         if n<2
%         fprintf(fid_data,'%s %s %s\n',cell2mat(Spies{n,4}{i,1}),cell2mat(Spies{n,4}{i,2}),cell2mat(Spies{n,4}{i,3}));
%         elseif n<8
        fprintf(fid_data,'%s %s %s %s\n',cell2mat(Spies{n,4}{i,1}),cell2mat(Spies{n,4}{i,2}),cell2mat(Spies{n,4}{i,3}),cell2mat(Spies{n,4}{i,4}));
%         else
%             fprintf(fid_data,'%s %s %s %s %s\n',cell2mat(Spies{n,4}{i,1}),cell2mat(Spies{n,4}{i,2}),cell2mat(Spies{n,4}{i,3}),cell2mat(Spies{n,4}{i,4}),cell2mat(Spies{n,4}{i,5}));
%         end
    end
end
fclose(fid_data);