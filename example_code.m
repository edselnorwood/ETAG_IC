% Alexander Moreno
% 02-20-2020
%
% Description: 
% A brief overview of the functions used to create the multi-coil 
% (elliptical/circular) wire antenna and B-Fields. This is missing some 
% functions. The code is a Matlab scripts 
%
% START: MAIN FUNC
%% Construct Wire Antenna
ra=0.3; % [m] radius_a (think of it as the "width" if it was a rect)
ri=0.3; % [m] radius_i (think of it as the "length" if it was a rect)
phi=2;  % [deg] sets the "pitch"
N=2;    % number of turns going along the z-direction
O=1;    % orientable (clock-wise or counter clock-wise)
wT=0.1; % [m] wire thickness
h=(1.1)*(2*wT*N); % height of the multi-coiled wire antenna
Nxy=1;  % number of turns along the xy-plane
% function to construct the multi-coil (elliptical/circular) wire antenna
% creates 4 points around the center, to create the wire thickness
[xS,yS,zS] = constrWireAnt(h,ra,ri,phi,N,O,wT,Nxy);
% note: xS,yS, and zS are all 1D array [1 x arraySize]
%% Calc B-Fields
I = 1;      % current [A]
Nx = 50;    % total number of segments along x-axis
Ny = 50;    % total number of segments along y-axis
Nz = 50;    % total number of segments along z-axis
Ns = [Nx,Ny,Nz];
% lower bounds for cuboid to calc the b-fields 
xminb=-(h+ra); yminb=-(h+ra); zminb=-(h+ra);  
% upper bounds for cuboid to calc the b-fields 
xmaxb=h+ra;    ymaxb=h+ra;    zmaxb= h+ra;
bBox = [xminb,yminb,zminb; xmaxb,ymaxb,zmaxb];
% Calc the b-fields for the given points in space defined above 
[X,Y,Z,BX,BY,BZ] = CalcB_WireAnt(I,xS,yS,zS,bBox,Ns);
% X,Y,Z,BX,BY,BZ dim size are Ny by Nx by Nz
normB=sqrt(BX.^2+BY.^2+BZ.^2);
nBX = BX./normB;
nBY = BY./normB;
nBZ = BZ./normB;
%% Plot
% multi-coil (elliptical/circular) wire antenna
hold all;
figure(1)
% only plotting the first set to make it visually easier to see
% this plot is not ploting the addtional points to "add" thickness
% this plot is not plotting the multi-coils (N or Nxy)
H=plot3(xS(1:361),zS(1:361),yS(1:361),'-');
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 
set(H,'linewidth',5); set(H,'color','r');
title('Antenna Structure');
grid on; axis equal; %axis tight;
%% B-Fields Quiver
% finding the cut-plane that i want to plot
% this is easy to find because total segments in each axis are the same
nn=ceil(Nx)/2; nn=25; 
BXn=BX./normB; BYn=BY./normB; BZn=BZ./normB;
% plots the arrows
quiver3(X(nn,:,:),Z(nn,:,:),Y(nn,:,:),BXn(nn,:,:),BZn(nn,:,:),BYn(nn,:,:),'w'); % matlab
view(0,90)
% B-Fields
% squeeze function reduces a matrix(multi-dim array) by 1dim
X2 = squeeze(X(nn,:,:));
Y2 = squeeze(Y(nn,:,:));
Z2 = squeeze(Z(nn,:,:));
B2 = squeeze(normB(nn,:,:));
% plots the contour/magnitude  
[M0,c]=contourf(X2,Z2,B2);
contourcbar;
view(0,90); grid on; axis tight;
xlabel('x[m]','FontWeight','bold','FontSize', 24); 
ylabel('z[m]','FontWeight','bold','FontSize', 24);
zlabel('z[m]','FontWeight','bold','FontSize', 24);
title('Coiled Wire Antenna:B-Fields (Model)','FontSize', 16);
view(0,90)
ylim([-0.5 0.9]); xlim([-0.9 0.9]);

% b-fields
% re-arranges the multi-dim matrix(array) into 1D array
% so i can export as a CSV file 
X  = reshape(X, [1,Nx*Ny*Nz]);
Y  = reshape(Y, [1,Nx*Ny*Nz]);
Z  = reshape(Z, [1,Nx*Ny*Nz]);
BX = reshape(BX,[1,Nx*Ny*Nz]);
BY = reshape(BY,[1,Nx*Ny*Nz]);
BZ = reshape(BZ,[1,Nx*Ny*Nz]);

% END: MAIN FUNC
%%
% fuction that creates the wire antenna 
function [xS0,yS0,zS0] = constrWireAnt(h,ra,ri,phi,N,O,wT,Nxy)
    helixSTEP = phi*(pi/180);
    start=0; fin = N*(2*pi) + helixSTEP/2;
    cst_xxx = start:helixSTEP:fin;
    size(cst_xxx)
    xS0=[];yS0=[];zS0=[];
    if(O==1) % clock wise
        for nx=1:Nxy
            % checking: 1st iteration 
            t   = wT/2;
            txy = (3/2)*wT*(nx-1);
            % flips back and forth of the winding (CW or CCW)
            % trying to simulate the actually winding 
            if(mod(nx,2)~=0) 
                %+z
                xS0 = [xS0,(ra+txy).*sin(cst_xxx)];
                yS0 = [yS0,(ri+txy).*cos(cst_xxx)];
                zS0 = [zS0,((h+t)*cst_xxx)./(2*pi*N)];
                %-z
                xS0 = [xS0,(ra+txy).*sin(cst_xxx)];
                yS0 = [yS0,(ri+txy).*cos(cst_xxx)];
                zS0 = [zS0,((h-t)*cst_xxx)./(2*pi*N)];
                %+xy
                xS0 = [xS0,(ra+t+txy).*sin(cst_xxx)];
                yS0 = [yS0,(ri+t+txy).*cos(cst_xxx)];
                zS0 = [zS0,((h)*cst_xxx)./(2*pi*N)];
                %-xy
                xS0 = [xS0,(ra-t+txy).*sin(cst_xxx)];
                yS0 = [yS0,(ri-t+txy).*cos(cst_xxx)];
                zS0 = [zS0,((h)*cst_xxx)./(2*pi*N)];
            else
                %+z
                xS0 = [xS0,-(ra+txy).*sin(cst_xxx)];
                yS0 = [yS0, (ri+txy).*cos(cst_xxx)];
                zS0 = [zS0, ((h+t)*cst_xxx)./(2*pi*N)];
                %-z
                xS0 = [xS0,-(ra+txy).*sin(cst_xxx)];
                yS0 = [yS0, (ri+txy).*cos(cst_xxx)];
                zS0 = [zS0, ((h-t)*cst_xxx)./(2*pi*N)];
                %+xy
                xS0 = [xS0,-(ra+t+txy).*sin(cst_xxx)];
                yS0 = [yS0, (ri+t+txy).*cos(cst_xxx)];
                zS0 = [zS0, ((h)*cst_xxx)./(2*pi*N)];
                %-xy
                xS0 = [xS0,-(ra-t+txy).*sin(cst_xxx)];
                yS0 = [yS0, (ri-t+txy).*cos(cst_xxx)];
                zS0 = [zS0, ((h)*cst_xxx)./(2*pi*N)];
                
            end % END: IF
        end %END: FOR
        
        else % counter clock wise
            for nx=1:Nxy
                t   = wT/2;
                txy = (3/2)*wT*(nx-1);
                
                if(mod(nx,2)~=0)
                    %+z
                    xS0 = [xS0,-(ra+txy).*sin(cst_xxx)];
                    yS0 = [yS0,(ri+txy).*cos(cst_xxx)];
                    zS0 = [zS0,((h+t)*cst_xxx)./(2*pi*N)];
                    %-z
                    xS0 = [xS0,-(ra+txy).*sin(cst_xxx)];
                    yS0 = [yS0, (ri+txy).*cos(cst_xxx)];
                    zS0 = [zS0, ((h-t)*cst_xxx)./(2*pi*N)];
                    %+xy
                    xS0 = [xS0,-(ra+t+txy).*sin(cst_xxx)];
                    yS0 = [yS0, (ri+t+txy).*cos(cst_xxx)];
                    zS0 = [zS0, ((h)*cst_xxx)./(2*pi*N)];
                    %-xy
                    xS0 = [xS0,-(ra-t+txy).*sin(cst_xxx)];
                    yS0 = [yS0, (ri-t+txy).*cos(cst_xxx)];
                    zS0 = [zS0, ((h)*cst_xxx)./(2*pi*N)];        
                    
                else
                    %+z
                    xS0 = [xS0,(ra+txy).*sin(cst_xxx)];
                    yS0 = [yS0,(ri+txy).*cos(cst_xxx)];
                    zS0 = [zS0,((h+t)*cst_xxx)./(2*pi*N)];
                    %-z
                    xS0 = [xS0,(ra+txy).*sin(cst_xxx)];
                    yS0 = [yS0,(ri+txy).*cos(cst_xxx)];
                    zS0 = [zS0,((h-t)*cst_xxx)./(2*pi*N)];
                    %+xy
                    xS0 = [xS0,(ra+t+txy).*sin(cst_xxx)];
                    yS0 = [yS0,(ri+t+txy).*cos(cst_xxx)];
                    zS0 = [zS0,((h)*cst_xxx)./(2*pi*N)]; 
                    %-xy
                    xS0 = [xS0,(ra-t+txy).*sin(cst_xxx)];
                    yS0 = [yS0,(ri-t+txy).*cos(cst_xxx)];
                    zS0 = [zS0,((h)*cst_xxx)./(2*pi*N)];                 
                end % END: IF
            end % END: FOR        
    end % END: 
    
end % end of constrWireAnt_10_25_2018