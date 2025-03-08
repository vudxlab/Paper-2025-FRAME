% StaBIL manual
% Example 1.2: static analysis of a 3D frame
% Units: m, kN
% Developed in MATLAB R2022a
%
% Author and programmer: Dx Laboratory
%
%         Adress: University of Transport and Communications
%
%
%   Main paper:
%   DOI:
%__________________________________________________________________________

clear
close all
clc

% Types={EltTypID EltName}
Types=  {1        'beam';
         2        'truss'};

% Sections
bCol=0.2;  % Column section width
hCol=0.2;  % Column section height
bBeam=0.2; % Beam section width
hBeam=0.3; % Beam section height
Atruss=0.001;

% Sections=[SecID A ky kz Ixx Iyy Izz yt yb zt zb]
Sections=  [1  hCol*bCol   Inf Inf 0.196*bCol^3*hCol   bCol^3*hCol/12 ...
               hCol^3*bCol/12  hCol/2  hCol/2  bCol/2  bCol/2;   % Columns
            2  hBeam*bBeam Inf Inf 0.249*bBeam^3*hBeam bBeam^3*hBeam/12 ...
               hBeam^3*bBeam/12 hBeam/2 hBeam/2 bBeam/2 bBeam/2; % Beams
            3  Atruss NaN NaN NaN NaN NaN NaN NaN NaN NaN];

% Materials=[MatID E      nu      rho];
Materials=  [
1     35e9  0.2     2500   1.5e11;   %concrete
2     35e9  0.2     2500   1.5e11
];  %steel

L=5;
H=3.5;
B=4;

% Nodes=[NodID X  Y  Z]
Nodes=  [1     0  0  0;
         2     L  0  0]      
Nodes=reprow(Nodes,1:2,2,[2 0 0 H])
Nodes=[Nodes; 
         10    2  0  2]             % reference node
Nodes=reprow(Nodes,1:7,2,[100 0 B 0])

figure
plotnodes(Nodes);

% Elements=[EltID TypID SecID MatID n1 n2 n3]
Elements=[ 1      1     1     1     1  3  10;
           2      1     1     1     2  4  10];
Elements=reprow(Elements,1:2,1,[2 0 0 0 2 2 0])
Elements=[ Elements;
           5      1     2     1     3  4  10; 
           6      1     2     1     5  6  10];   
Elements=reprow(Elements,1:6,2,[100 0 0 0 100 100 100]) 
Elements=[ Elements;
          301     2     2     2     3  103 NaN; 
          302     2     2     2     4  104 NaN];     
Elements=reprow(Elements,19:20,1,[2 0 0 0 2 2 0]) 
Elements=reprow(Elements,19:22,1,[4 0 0 0 100 100 0]) 
% Elements=[ Elements;
%           309     2     3     2      4 106 NaN; 
%           310     2     3     2      6 104 NaN];           

hold('on');
plotelem(Nodes,Elements,Types);
title('Nodes and elements');

% Plot elements in different colors in order to check the section definitions
figure
plotelem(Nodes,Elements(find(Elements(:,3)==1),:),Types,'r','Numbering','on','LineWidth',2);
hold('on');
plotelem(Nodes,Elements(find(Elements(:,3)==2),:),Types,'b','Numbering','on','LineWidth',2);
plotelem(Nodes,Elements(find(Elements(:,3)==3),:),Types,'g','Numbering','on','LineWidth',2);
title('Phần tử: Mặt cắt')

% Degrees of freedom
DOF=getdof(Elements,Types);

% Boundary conditions: hinges 
seldof=[  1.01;   1.02;   1.03;  1.04;   1.05;   1.06;  
          2.01;   2.02;   2.03;  2.04;   2.05;   2.06; 
        101.01; 101.02; 101.03;  101.04;   101.05;   101.06; 
        102.01; 102.02; 102.03;  102.04;   102.05;   102.06; 
        201.01; 201.02; 201.03;  201.04;   201.05;   201.06; 
        202.01; 202.02; 202.03;  202.04;   202.05;   202.06;
        ];

DOF=removedof(DOF,seldof);

% Assembly of stiffness matrix K 
[K,M] =asmkm(Nodes,Elements,Types,Sections,Materials,DOF);


% Eigenvalue problem
nMode=10;
[phi,omega]=eigfem(K,M,nMode);

% Display eigenfrequenties
disp('Lowest eigenfrequencies [Hz]');
f0=omega/2/pi;
disp(f0);

% Plot eigenmodes

% for i = 1:10
%     figure;
%     animdisp(Nodes,Elements,Types,DOF,phi(:,i),'DispMax','off','LineWidth',2)
%     title({['Hình thái dao động ' num2str(i)], ['f= ' num2str(f0(i)) ' Hz']})
% end




% Loads

% % Own weight
% DLoadsOwn=accel([0 0 9.81],Elements,Types,Sections,Materials);
% 
% % Wind load
% 
% % DLoads=[EltID n1globalX n1globalY n1globalZ ...]
% DLoadsWind =[1 0 0    0 0 1500 0;
%              2 0 0    0 0 1500 0;
%              3 0 1500 0 0 1500 0;
%              4 0 1500 0 0 1500 0];
%              
% DLoads=multdloads(DLoadsOwn,DLoadsWind);
% 
% P=elemloads(DLoads,Nodes,Elements,Types,DOF);
% 
% % Solve K * U = P 
% U=K\P;

% figure
% plotdisp(Nodes,Elements,Types,DOF,U(:,1),DLoads(:,:,1),Sections,Materials)
% title('Displacements: own weight')

% figure
% plotdisp(Nodes,Elements,Types,DOF,U(:,2),DLoads(:,:,2),Sections,Materials)
% title('Displacements: wind')
% 
% % Compute forces
% [ForcesLCS,ForcesGCS]=elemforces(Nodes,Elements,Types,Sections,Materials,DOF,U,DLoads);
% 
% % Compute reaction forces for load case 1
% Freac=reaction(Elements,ForcesGCS(:,:,1),[1.03; 2.03; 101.03; 102.03; 201.03; 202.03])
% 
% % Plot element forces for load case 1
% figure
% plotforc('norm',Nodes,Elements,Types,ForcesLCS(:,:,1),DLoads(:,:,1))
% title('Normal forces: Own weight')
% figure
% plotforc('sheary',Nodes,Elements,Types,ForcesLCS(:,:,1),DLoads(:,:,1))
% title('Shear forces along y: Own weight')
% figure
% plotforc('shearz',Nodes,Elements,Types,ForcesLCS(:,:,1),DLoads(:,:,1))
% title('Shear forces along z: Own weight')
% figure
% plotforc('momx',Nodes,Elements,Types,ForcesLCS(:,:,1),DLoads(:,:,1))
% title('Torsional moments: Own weight')
% figure
% plotforc('momy',Nodes,Elements,Types,ForcesLCS(:,:,1),DLoads(:,:,1))
% title('Bending moments around y: Own weight')
% figure
% plotforc('momz',Nodes,Elements,Types,ForcesLCS(:,:,1),DLoads(:,:,1))
% title('Bending moments around z: Own weight')
% 
% % Plot element forces for load case 2
% figure
% plotforc('norm',Nodes,Elements,Types,ForcesLCS(:,:,2),DLoads(:,:,2))
% title('Normal forces: Wind')
% figure
% plotforc('sheary',Nodes,Elements,Types,ForcesLCS(:,:,2),DLoads(:,:,2))
% title('Shear forces along y: Wind')
% figure
% plotforc('shearz',Nodes,Elements,Types,ForcesLCS(:,:,2),DLoads(:,:,2))
% title('Shear forces along z: Wind')
% figure
% plotforc('momx',Nodes,Elements,Types,ForcesLCS(:,:,2),DLoads(:,:,2))
% title('Torsional moments: Wind')
% figure
% plotforc('momy',Nodes,Elements,Types,ForcesLCS(:,:,2),DLoads(:,:,2))
% title('Bending moments around y : Wind')
% figure
% plotforc('momz',Nodes,Elements,Types,ForcesLCS(:,:,2),DLoads(:,:,2))
% title('Bending moments around z: Wind')
% 
% % Load combinations
% 
% % Safety factors
% gamma_own=1.35;
% gamma_wind=1.5;
% 
% % Combination factors
% psi_wind=1;
% 
% % Load combination (Ultimate Limit State, ULS)
% U_ULS=gamma_own*U(:,1)+gamma_wind*psi_wind*U(:,2);
% Forces_ULS=gamma_own*ForcesLCS(:,:,1)+gamma_wind*psi_wind*ForcesLCS(:,:,2);
% DLoads_ULS(:,1)=DLoads(:,1,1)
% DLoads_ULS(:,2:7)=gamma_own*DLoads(:,2:7,1)+gamma_wind*psi_wind*DLoads(:,2:7,2);
% 
% figure
% plotdisp(Nodes,Elements,Types,DOF,U_ULS,DLoads_ULS,Sections,Materials)
% 
% printdisp(Nodes,DOF,U_ULS);
% 
% printforc(Elements,Forces_ULS);
% 
% % Plot element forces
% figure
% plotforc('norm',Nodes,Elements,Types,Forces_ULS,DLoads_ULS)
% title('Normal forces: ULS')
% figure
% plotforc('sheary',Nodes,Elements,Types,Forces_ULS,DLoads_ULS)
% title('Shear forces along y: ULS')
% figure
% plotforc('shearz',Nodes,Elements,Types,Forces_ULS,DLoads_ULS)
% title('Shear forces along z: ULS')
% figure
% plotforc('momx',Nodes,Elements,Types,Forces_ULS,DLoads_ULS)
% title('Torsional moments: ULS')
% figure
% plotforc('momy',Nodes,Elements,Types,Forces_ULS,DLoads_ULS)
% title('Bending moments around y: ULS')
% figure
% plotforc('momz',Nodes,Elements,Types,Forces_ULS,DLoads_ULS)
% title('Bending moments around z: ULS')
% 
% % Plot stresses
% figure
% plotstress('snorm',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Normal stresses due to normal forces')
% figure
% plotstress('smomyt',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Normal stresses due to bending moments around y: top')
% figure
% plotstress('smomyb',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Normal stresses due to bending moments around y: bottom')
% figure
% plotstress('smomzt',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Normal stresses due to bending moments around z: top')
% figure
% plotstress('smomzb',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Normal stresses due to bending moments around z: bottom')
% figure
% plotstress('smax',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Maximal normal stresses')
% figure
% plotstress('smin',Nodes,Elements,Types,Sections,Forces_ULS,DLoads_ULS)
% title('Minimal normal stresses')
