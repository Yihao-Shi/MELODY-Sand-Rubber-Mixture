%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                       %%
%%                      Packing 2D                       %%
%%                                                       %%
%%                     Main Program                      %%
%%                Version 1.0 ; July 2012                %%
%%                                                       %%
%%                Author: Guilhem Mollon                 %%
%%               Supervisor: Jidong Zhao                 %%
%%                                                       %%
%%          Realized at Hong Kong University             %%
%%              of Science and Technology                %%
%%                     Year 2012                         %%
%%                                                       %%
%%            Please read enclosed .txt file             %%
%%                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% A. INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Created data files %%%%
FileName='Packing2D_Example.mat';                   %Name of the created MATLAB data file
PFCFileCreation=1;                                  %0: no PFC2D file created ; 1: creation of a PFC2D file
PFCFileName='Packing2D_Example.dat';                %Name of the created PFC2D file

%%%% Constrained Voronoi Tessellation parameters %%%%
Nparticles=500;                                     %Number of generated points in the domain
TargetSurfaceCOV=0.8;                               %Target COV of the grains surfaces
TargetMainOrientation=0;                            %Target angle of main orientation
TargetAnisotropy=1;                                 %Target ratio of anisotropy
DomainPoints=[0,0;1,0;1,1;0,1;0,0];                 %Packing domain. Syntax: [x1,y1;x2,y2;   ;xn,yn;x1,y1] (must be closed, with points order anticlockwise)
NiterMax=1000000;                                   %Maximum number of iterations for Constrained Voronoi Tessellation
TargetError=0.002;                                  %Target error for Constrained Voronoi Tessellation

%%%% Cell Filling Parameters %%%%
TargetSolidFraction=0.65;                           %Target solid fraction
NvarOptim=6;                                        %Number of optimization variables for the filling algorithm (optimized for n=2 to n=NvarOptim+1)

%%%% Fourier Spectrum Properties %%%%
COVSpectrum=0;                                      %COV of the Fourier descriptors
TypeCOVSpectrum=1;                                  %0:Modes vary individually ; 1:Spectrum varies altogether
TypeSpectrum=0;                                     %0:Taylored spectrum ; 1:Existing spectrum

%%%% If TypeSpectrum is 0 %%%%
DescriptorD2=0.1;                                   %Fourier descriptor n=2
DescriptorD3=0.1;                                   %Fourier descriptor n=3
SpectrumDecay1=-2;                                  %Exponential decay from n=3 to n=7
DescriptorD8=0.015;                                 %Fourier descriptor n=8
SpectrumDecay2=-2;                                  %Exponential decay from n>8

%%%% If TypeSpectrum is 1 %%%%
FileSpectrum='Spectrum_Toyoura.mat';                %Name of a Spectrum file

%%%% ODECS Algorithm Parameters %%%%
Dmax=0.02;                                          %Maximum distance for which a contour point is considered as "covered" by a given disc
Rmin=0.02;                                          %Minimum radius of a given disc
Pmax=0.05;                                          %Accepted proportion of "uncovered" contour points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% B. COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%     PRESS F5    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[VoronoiCells,VoronoiVertices,ConvergenceHistory,ODECS,SampleProperties,SolidFraction,CellsSurfaces,CellsOrientations,SampleD50,SampleCu,PDFOrientation,PDFSurface,PDFElongation,PDFRoundness,PDFCircularity,PDFRegularity]=Secondary_Program(FileName,PFCFileCreation,PFCFileName,Nparticles,TargetSurfaceCOV,TargetMainOrientation,TargetAnisotropy,DomainPoints,NiterMax,TargetError,TargetSolidFraction,NvarOptim,COVSpectrum,TypeCOVSpectrum,TypeSpectrum,FileSpectrum,DescriptorD2,DescriptorD3,SpectrumDecay1,DescriptorD8,SpectrumDecay2,Dmax,Rmin,Pmax);
save(FileName)