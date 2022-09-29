
% 2-dimensional analytical solution - Transient (using the analytical solution for
% instantaneous point source integrated with respect to time)

function NumVSAnalytical
clc
close all
clear ll

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 3;
ny = 100;
media ='SW'; %GW or SW

% Analytical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C0 = 150;     % concentration of the injected water - units: ppm
Q = 86.4;      % injection rate - units: 1m3/day
n = 1;      % porosity 
Dx = 0.01;      % longiitudinal dispersivity - y-direction m2/s
Dy = 0.01;      % longiitudinal dispersivity - y-direction m2/s
r = 0;         % reaction rate
Xc = 2;        % injection location loation
Yc = 2;        % injection location loation
tsim = 200;      % time after injection - units: days
Zgw = 70;      % aquifer thinkness
V = -0.5;        % GW seepage flow velocity (m/day) (average intersticial velocity)

% Numerical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_simNUMERICAL = 2000; % time step to load from the model results (to compare the analytical solution against)
address_sim = strcat('D:\FCL_SEC\Models\S2R2\Playground\Visualization_VisIT\1_S2R2_direct_results\',media,'\',int2str(t_simNUMERICAL),'_',media,'.csv');
    
% Storing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storage.C0 = C0;
storage.Q = Q;
storage.n = n;
storage.Dx = Dx;
storage.Dy = Dy;
storage.r = r;
storage.Xc = Xc;
storage.Yc = Yc;
storage.tsim = tsim;
storage.Zgw = Zgw;
storage.V = V;
storage.nx = nx;
storage.ny = ny;
storage.address_sim = address_sim;
save InputData storage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Analytical_calc
%Numerical_calc
Compare_plot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Analytical_calc

load InputData storage
C0=storage.C0;
Q=storage.Q;
n=storage.n;
Dx=storage.Dx;
Dy=storage.Dy;
r=storage.r;
Xc=storage.Xc;
Yc=storage.Yc;
tsim=storage.tsim;
Zgw=storage.Zgw;
V=storage.V;
nx=storage.nx;
ny=storage.ny;

Cfinal=zeros(nx,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');
for x=1:nx
    waitbar(x/nx)
    for y=1:ny
        
    % calculation of A
    Qprime = Q / Zgw; %fluid injection rate per unit thinkness of aquifer
    
    A = C0*Qprime / ...
        ( 4*n*pi()*(Dx*Dy)^(0.5) );

    % calculation of B
    B = exp(...
            V*(x-Xc)/ ...
            (2*Dx));

    % Integral to solve
    fun = @(t) exp( ...
                    -( V^2 / (4*Dx) + r) * t ...
                    -( (x-Xc)^2 ./ (4*Dx*t) ) ...
                    -( (y-Yc)^2 ./ (4*Dy*t))) ;
                
    D = quad(fun,0,tsim);

    % Final calculation
    C = A*B*D;
    
    % Saving results
    Cfinal(x,y) = C;
   
    end
end
close(h)

disp('Analytical solution calculation: OK')

storage.C_analytical = Cfinal';
save AnalytSOL storage;

%figure
%contour(Cfinal','ShowText','on');
%colormap
%colorbar('location','Eastoutside')
%grid on
%xlabel('Distance along the x-direction')
%ylabel('Distance along the y-direction')
%legend('Analytical Solution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Numerical_calc

load InputData storage
nx=storage.nx;
ny=storage.ny;
address_sim=storage.address_sim;

A=importdata(address_sim);
Adata=A.data;

% Transforming into a matricial form (MODFLOW grid convention)
c_num=zeros(ny,nx);
h=waitbar(0,'NUMERICAL Solution: Importing and transforming...');
for j=1:ny
     waitbar(j/ny)
    for i=1:nx
        l=find(Adata(:,2)==i & Adata(:,3)==j);
        c_num(j,i)=Adata(l,4);
    end
end
close(h)

disp('NUMERICAL solution loading and transformation: OK')

storage.C_numerical=c_num;
save NumericalSOL storage;

%figure
%contour(c_num,'ShowText','on');
%colormap
%colorbar('location','Eastoutside')
%grid on
%xlabel('Distance along the x-direction')
%ylabel('Distance along the y-direction')
%legend('Numerical Solution')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Compare_plot

load InputData storage
InputData = storage;

load AnalytSOL storage
C_analytical = storage.C_analytical;

%load NumericalSOL storage
%C_numerical = storage.C_numerical;

colormap(gray)
cmap = colormap
cmap = flipud(cmap)

cm = [jet(64);gray(64)];

figure
subplot(3,3,[1,2,4,5])
[~, hh]=contour(C_analytical);
set(hh,'LineColor','k','ShowText','on')

hold on

%{
pcolor(C_numerical)%,'ShowText','on');
shading flat
colormap(cmap)
grid on
%}

colorbar('location','East')
xlabel('Distance along the x-direction')
ylabel('Distance along the y-direction')
legend('Analytical','Numerical')

% Longitudinal profile
subplot(3,3,[3,6])
plot(C_analytical(:,InputData.Yc),1:1:InputData.ny,'k','linewidth',1)
ylim([1 InputData.ny])
xlabel('Concentration')
ylabel('Distance along the y-direction')
legend('Analytical Solution')
set (gca,'Xdir','reverse')
grid on

% Transverse profile
subplot(3,3,[7,8])
plot(1:1:InputData.nx, C_analytical(round(InputData.nx/2),:),'k','linewidth',1)
xlim([1 InputData.nx])
legend('Analytical Solution')
ylabel('Concentration')
xlabel('Distance along the x-direction')
grid on

disp('ANALYTICAL and NUMERICAL solutions compared: OK')
disp('Finished')