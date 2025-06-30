% 1. Simulation name
Simulation_Name='Example';


% 2. Initilization
PFCball;
Contours=cell(1,1);
Distributions=cell(1,1);
Interpolations=cell(1,1);
Integrations=cell(1,1);
Bodies_Materials=cell(1,1);
Imposed_Pressures=cell(1,1);
Imposed_Velocities=cell(1,1);
Alid=cell(1,1);
Deactivations=cell(1,1);
To_Plot=zeros(44,1);
row_d=1;                                      %count the number of deformable particles
row_r=1;                                      %count the number of rigid particles
xdelta=(x_max-x_min)/10;
ydelta=(y_max-y_min)/20;


% 3. Bodies
for i=1:size(deformable_pos,1)
    Contours{i,1}={'Closed',cat(2,deformable_pos(row_d,1)+deformable_pos(row_d,3)*cos([0:1:360]'/180*pi),deformable_pos(row_d,2)+deformable_pos(row_d,3)*sin([0:1:360]'/180*pi)),'Linear'};
    Distributions{i,1}={'Unstructured'};
    Distributions{i,2}=deformable_pos(row_r,3)/10;
    Distributions{i,3}=8;
    Interpolations{i,1}='MLS';
    Interpolations{i,2}=10;
    Integrations{i,1}='Gauss';
    Integrations{i,2}=3;
    Detections(i,1)=deformable_pos(row_r,3)/5;
    Detections(i,2)=deformable_pos(row_r,3)/5;
    Bodies_Materials{i,1}='Material1';
    Imposed_Pressures{i,1}={[],'None',[],'None';...
                           };
    Imposed_Velocities{i,1}={[],'None',[],[],'None',[];...
                             };
    Initial_Velocities{i,1}=[0,0];
    Mesh_Ratios(i,1:2)=[3,1];
    Status{i,1}='active';
    Alid{i,1}=[];
    Alid{i,2}=[];
    Alid{i,3}=[];
    row_d=row_d+1;
end

for j=size(deformable_pos,1)+1:size(deformable_pos,1)+size(rigid_pos,1)
    Contours{j,1}={'Closed',cat(2,rigid_pos(row_r,1)+rigid_pos(row_r,3)*cos([0:1:360]'/180*pi),rigid_pos(row_r,2)+rigid_pos(row_r,3)*sin([0:1:360]'/180*pi)),'Linear'};
    Distributions{j,1}={'Rigid'};
    Distributions{j,2}=rigid_pos(row_r,3)/10;
    Distributions{j,3}=8;
    Interpolations{j,1}='None';
    Interpolations{j,2}=10;
    Integrations{j,1}='Gauss';
    Integrations{j,2}=3;
    Detections(j,1)=rigid_pos(row_r,3)/5;
    Detections(j,2)=rigid_pos(row_r,3)/5;
    Bodies_Materials{j,1}='Material2';
    Imposed_Pressures{j,1}={[],'None',[],'None';...
                            };
    Imposed_Velocities{j,1}={[],'None',[],[],'None',[];...
                             };
    Initial_Velocities{j,1}=[0,0];
    Mesh_Ratios(j,1:2)=[1,1];
    Status{j,1}='active';
    Alid{j,1}=[];
    Alid{j,2}=[];
    Alid{j,3}=[];
    row_r=row_r+1;
end

contour_row=size(deformable_pos,1)+size(rigid_pos,1)+1;
Contours{contour_row,1}={'Closed',[x_min_expand,y_min-xdelta;x_max_expand,y_min-xdelta;x_max_expand,y_min;x_min_expand,y_min],'Linear'};
Distributions{contour_row,1}={'Rigid'};
Distributions{contour_row,2}=xdelta/5;
Distributions{contour_row,3}=8;
Interpolations{contour_row,1}='None';
Interpolations{contour_row,2}=10;
Integrations{contour_row,1}='Gauss';
Integrations{contour_row,2}=3;
Detections(contour_row,1)=rhi/5;
Detections(contour_row,2)=rhi/5;
Bodies_Materials{contour_row,1}='Material3';
Imposed_Pressures{contour_row,1}={[],'None',[],'None';...
                        };
Imposed_Velocities{contour_row,1}={[0,0;1e6,0],'Soft',[2e10,0],[0,0;1e6,0],'Soft',[2e10,0];...
                         };
Initial_Velocities{contour_row,1}=[0,0];
Mesh_Ratios(contour_row,1:2)=[4,1];
Status{contour_row,1}='active';
Alid{contour_row,1}=[];
Alid{contour_row,2}=[];
Alid{contour_row,3}=[];

contour_row=contour_row+1;
Contours{contour_row,1}={'Closed',[x_min_expand,y_max;x_max_expand,y_max;x_max_expand,y_max+xdelta;x_min_expand,y_max+xdelta],'Linear'};
Distributions{contour_row,1}={'Rigid'};
Distributions{contour_row,2}=xdelta/5;
Distributions{contour_row,3}=8;
Interpolations{contour_row,1}='None';
Interpolations{contour_row,2}=10;
Integrations{contour_row,1}='Gauss';
Integrations{contour_row,2}=3;
Detections(contour_row,1)=rhi/5;
Detections(contour_row,2)=rhi/5;
Bodies_Materials{contour_row,1}='Material3';
Imposed_Pressures{contour_row,1}={[],'None',[],'None';...
                        };
Imposed_Velocities{contour_row,1}={[0,0;1e6,0],'Soft',[2e10,0],[0,0;1e6,0],'Soft',[2e10,0];...
                         };
Initial_Velocities{contour_row,1}=[0,0];
Mesh_Ratios(contour_row,1:2)=[1,4];
Status{contour_row,1}='active';
Alid{contour_row,1}=[];
Alid{contour_row,2}=[];
Alid{contour_row,3}=[];

contour_row=contour_row+1;
Contours{contour_row,1}={'Closed',[x_min-ydelta,y_min_expand;x_min,y_min_expand;x_min,y_max_expand;x_min-ydelta,y_max_expand],'Linear'};
Distributions{contour_row,1}={'Rigid'};
Distributions{contour_row,2}=ydelta/2.5;
Distributions{contour_row,3}=8;
Interpolations{contour_row,1}='None';
Interpolations{contour_row,2}=10;
Integrations{contour_row,1}='Gauss';
Integrations{contour_row,2}=3;
Detections(contour_row,1)=rhi/5;
Detections(contour_row,2)=rhi/5;
Bodies_Materials{contour_row,1}='Material3';
Imposed_Pressures{contour_row,1}={[0,0;5,1e3],'Oriented',[],'None';...
                        };
Imposed_Velocities{contour_row,1}={[],'None',[],[0,0;1e6,0],'Soft',[2e10,0];...
                         };
Initial_Velocities{contour_row,1}=[0,0];
Mesh_Ratios(contour_row,1:2)=[1,4];
Status{contour_row,1}='active';
Alid{contour_row,1}=[];
Alid{contour_row,2}=[];
Alid{contour_row,3}=[];

contour_row=contour_row+1;
Contours{contour_row,1}={'Closed',[x_max,y_min_expand;x_max+ydelta,y_min_expand;x_max+ydelta,y_max_expand;x_max,y_max_expand],'Linear'};
Distributions{contour_row,1}={'Rigid'};
Distributions{contour_row,2}=ydelta/2.5;
Distributions{contour_row,3}=8;
Interpolations{contour_row,1}='None';
Interpolations{contour_row,2}=10;
Integrations{contour_row,1}='Gauss';
Integrations{contour_row,2}=3;
Detections(contour_row,1)=rhi/5;
Detections(contour_row,2)=rhi/5;
Bodies_Materials{contour_row,1}='Material3';
Imposed_Pressures{contour_row,1}={[0,0;5,-1e3],'Oriented',[],'None';...
                        };
Imposed_Velocities{contour_row,1}={[],'None',[],[0,0;1e6,0],'Soft',[2e10,0];...
                         };
Initial_Velocities{contour_row,1}=[0,0];
Mesh_Ratios(contour_row,1:2)=[1,4];
Status{contour_row,1}='active';
Alid{contour_row,1}=[];
Alid{contour_row,2}=[];
Alid{contour_row,3}=[];


% 4. Materials and contact laws
Materials={'Material1','NeoHookean',[1100e5,0.7,0.07,7e5,0.495,0];...
           'Material2','ElasticLinear',[2650e5,0.7,0.07,7e8,0.3,0];...
           'Material3','ElasticLinear',[2650e5,0.7,0.07,8e8,0.3,0];...
           };
Contact_Laws={'Material1','Material1','MohrCoulomb','Evolutive',[1e9,1e9,0.1,0,0];...
              'Material1','Material2','MohrCoulomb','Evolutive',[1.82e9,1.82e9,0.1,0,0];...
              'Material2','Material2','DampedMohrCoulomb','Evolutive',[1e10,1e10,0.1,0,0,1.0];...
              'Material1','Material3','MohrCoulomb','Evolutive',[2e9,2e9,0.1,0,0];...
              'Material2','Material3','DampedMohrCoulomb','Evolutive',[2e10,2e10,0.1,0,0,1.0];...
              };


% 5. General boundary conditions
Periodic_Boundaries=[-1e6,1e6];
Gravity=[0;0];


% 6. Text outputs
Monitorings=cell(1,1);
Spies=cell(1,1);


% 7. Graphic outputs
To_Plot(1)=0;   % Body Index
To_Plot(2)=1;   % Initial Position
To_Plot(3)=0;   % Current Position
To_Plot(4)=1;   % Displacement
To_Plot(5)=1;   % Velocity
To_Plot(6)=0;   % Acceleration
To_Plot(7)=0;   % Force
To_Plot(8)=0;   % Internal Force
To_Plot(9)=1;   % Contact Force
To_Plot(10)=0;  % Body Force
To_Plot(11)=0;  % Dirichlet Force
To_Plot(12)=0;  % Neumann Force
To_Plot(13)=0;  % Damping Force
To_Plot(14)=0;  % Alid Force
To_Plot(15)=0;  % Jacobian
To_Plot(16)=1;  % Cauchy XX Stress
To_Plot(17)=1;  % Cauchy YY Stress
To_Plot(18)=1;  % Cauchy XY Stress
To_Plot(19)=0;  % Cauchy ZZ Stress
To_Plot(20)=0;  % Tresca Stress
To_Plot(21)=1;  % Von Mises Stress
To_Plot(22)=0;  % Major Principal Stress
To_Plot(23)=0;  % Intermediate Principal Stress
To_Plot(24)=0;  % Minor Principal Stress
To_Plot(25)=0;  % Spherical Stress
To_Plot(26)=0;  % Green-Lagrange XX strain
To_Plot(27)=0;  % Green-Lagrange YY strain
To_Plot(28)=0;  % Green-Lagrange XY strain
To_Plot(29)=0;  % Norm of the Green-Lagrange strain tensor
To_Plot(30)=0;  % Body Damage
To_Plot(31)=0;  % Body Relative Damage
To_Plot(32)=0;  % Normalized Displacement Error
To_Plot(33)=0;  % Displacement Error
To_Plot(34)=0;  % Internal Work
To_Plot(35)=0;  % Contact Work
To_Plot(36)=0;  % Body Work
To_Plot(37)=0;  % Dirichlet Work
To_Plot(38)=0;  % Neumann Work
To_Plot(39)=0;  % Damping Work
To_Plot(40)=0;  % Alid Work
To_Plot(41)=0;  % Temperature
To_Plot(42)=0;  % Scaling Parameter
To_Plot(43)=0;  % Active Contacts
To_Plot(44)=0;  % Contacting Bodies
To_Plot(45)=0;  % Internal Variable 0
To_Plot(46)=0;  % Internal Variable 1
To_Plot(47)=0;  % Internal Variable 2
To_Plot(48)=0;  % Internal Variable 3
To_Plot(49)=0;  % Internal Variable 4
To_Plot(50)=0;  % Internal Variable 5
To_Plot(51)=0;  % Internal Variable 6
To_Plot(52)=0;  % Internal Variable 7
To_Plot(53)=0;  % Internal Variable 8
To_Plot(54)=0;  % Internal Variable 9
Chains_Parameters=[1000,0];
Fields_Parameters=[x_min,x_max,y_min,y_max,1e-2,1e-1];



% 8. Numerical parameters
Scheme='Adaptive_Euler';
Scheme_Parameters=[1e-4,0.1,10];
Contact_Updating_Period=0.001;
Time_Stepping_Parameters=[0,1e-4,15];
Save_Periods=[1,1];


% 9. Flags
Activate_Plot=1;
Initialize_CZM=0;