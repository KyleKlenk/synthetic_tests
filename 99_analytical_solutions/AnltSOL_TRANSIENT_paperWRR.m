
% 2-dimensional analytical solution - Transient (using the analytical solution for
% instantaneous point source integrated with respect to time)

function NumVSAnalytical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 12;
ny = 600;
media ='SW'; %GW or SW

% Analytical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mLayerDepth_summa_m         = 0.006;                        % m
hru_area_m2                 = 32700;                        % m2

%M       = 350;        % load; g
mLayerVolFracWat_summa_m3 = 39.376154283055975;
%UpLayVol = mLayerVolFracWat_summa_m3;
C0 = 2;
n       = 1;       	% porosity 
Dx      = 0.0006;       	% longiitudinal dispersivity - y-direction m2/s
Dy      = 0.0001;      	% longiitudinal dispersivity - y-direction m2/s
r       = 0;       	% reaction rate
%r       = 0.01/(60*60*24);       	% reaction rate
Xc      = 175;     	% injection location loation
Yc      = 1;    	% injection location loation
%tsim    = 24*24*60;    	% time after injection - units: days
Zgw     = mLayerDepth_summa_m;        % aquifer thinkness
specificYield = 0.22;
iLayerLiqFluxSoil_summa_m_s = 1.0359700931174889e-05 / 0.25 ;       % m/s
%iLayerLiqFluxSoil_summa_m3  = 0.30488599840555514;         % m3/s

%data_step                   = 900;                          % Seconds

V = iLayerLiqFluxSoil_summa_m_s;	% flow velocity (m/day) (average intersticial velocity)
%V = 0.6;
%V = -iLayerLiqFluxSoil_summa_m3/data_step * (60*60*24)/hru_area_m2;	% flow velocity (m/day) (average intersticial velocity)

tsim = [...
        ...0,...
        %60 * 60,...
        %60 * 60 * 2,...
        %60 * 60 * 5,...
        %60 * 60 * 10,...
        ...60 * 60 * 24,...
        %60 * 60 * 24 * 3,...
        %60 * 60 * 24 * 8,...
        60 * 60 * 24 * 15,...
        60 * 60 * 24 * 70,...
        60 * 60 * 24 * 120,...
        ];

    
V = iLayerLiqFluxSoil_summa_m_s
V = 5.1476e-08 / 1000;
%Dy = 0.06
C0 = 2;
%tsim = [2.5, 5, 10, 15, 20];

VolLayer = 39.376154283074612;   % m3
wflux_s2r = 0.30488599840048847; % m3/s

V = 6 / (VolLayer / (wflux_s2r / 900));
V = iLayerLiqFluxSoil_summa_m_s
L = V * 60 * 60 * 24 * 20 / 1000;

% Numerical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_simNUMERICAL = 2000; % time step to load from the model results (to compare the analytical solution against)
address_sim = strcat('D:\FCL_SEC\Models\S2R2\Playground\Visualization_VisIT\1_S2R2_direct_results\',media,'\',int2str(t_simNUMERICAL),'_',media,'.csv');
    
% Storing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storage.C0 = C0;
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
n_ti = numel(tsim);
for ti = 1:n_ti
    
    tsim_i = tsim(ti);
    
    Analytical_calc_1D_contS(tsim_i)
    %Analytical_calc_1D_2(tsim_i)
    %Analytical_calc_1D_3(tsim_i)
    %Analytical_calc_1D_justAdv(tsim_i, mLayerDepth_summa_m)
    %Analytical_calc_2D
    %Numerical_calc
    Compare_plot(ti, n_ti, tsim_i,mLayerDepth_summa_m)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1D, continuous source, Eq. 60 of 
% USGS Chapter B7
% ANALYTICAL SOLUTIONS FOR ONE-, TWO-, AND THREE-DIMENSIONAL SOLUTE
%TRANSPORT IN GROUND-WATER SYSTEMS WITH UNIFORM FLOW
function Analytical_calc_1D_contS(tsim_i)

load InputData storage
C0=storage.C0;
Dy=storage.Dy;
r=storage.r;
V=storage.V;
ny=storage.ny;

Cfinal=zeros(1,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');

for y=1:ny
    waitbar(y/ny)

    % calculation of A

    %A = (M / UpLayVol) /...
    %    ( 4*n*pi()*(Dx*Dy)^(0.5) );
    A =  C0 / 2;

    % calculation of B
    U = (V^2 + 4*r*Dy)^0.5;
    
    B1 = exp(...
            y * (V - U)/ ...
            (2 * Dy)...
            );
        
    B2 = erfc(...
            (y - (U * tsim_i) ) / ...
            (2 * (Dy * tsim_i)^0.5));
    
    C1 = exp(...
            y * (V + U)/ ...
            (2 * Dy));
    
    C2 = erfc(...
            (y + (U * tsim_i) ) / ...
            (2 * (Dy * tsim_i)^0.5));

    % Final calculation
    C = A * (B1 * B2 + C1 * C2);

    % Saving results
    Cfinal(1,y) = C;

end

close(h)

disp('Analytical solution calculation: OK')

storage.C_analytical = Cfinal';
save AnalytSOL storage;

function Analytical_calc_1D_2(tsim_i)

load InputData storage
C0=storage.C0;
Dy=storage.Dy;
r=storage.r;
V=storage.V;
ny=storage.ny;

Cfinal=zeros(1,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');

for y=1:ny
    waitbar(y/ny)

    % calculation of A

    %A = (M / UpLayVol) /...
    %    ( 4*n*pi()*(Dx*Dy)^(0.5) );
    
    U = (V^2 + 4*r*Dy)^0.5;
    
    A =  C0 * V^2 / (4 * r * Dy);

    % calculation of B
    
    
    B1 = 2 * exp((y * V / Dy) - r * tsim_i);
        
    B2 = erfc(...
            (y + V * tsim_i) / ...
            (2 * (Dy * tsim_i)^0.5));
    
   C1 = U / V - 1; 
        
   C2 = exp(...
            y * (V - U)/ ...
            (2 * Dy));
    
    C3 = erfc(...
            (y - U * tsim_i) / ...
            (2 * (Dy * tsim_i)^0.5));
        
    D1 = U / V + 1;
    
    D2 = exp(...
            y * (V + U)/ ...
            (2 * Dy));
    
    D3 = erfc(...
            (y + U * tsim_i) / ...
            (2 * (Dy * tsim_i)^0.5));

    % Final calculation
    C = A * (B1 * B2 + C1 * C2 * C3 + D1 * D2 * D3);

    % Saving results
    Cfinal(1,y) = C;

end

close(h)

disp('Analytical solution calculation: OK')

storage.C_analytical = Cfinal';
save AnalytSOL storage;

function Analytical_calc_1D_3(tsim_i)

load InputData storage
C0=storage.C0;
Dy=storage.Dy;
r=storage.r;
V=storage.V;
ny=storage.ny;

Cfinal=zeros(1,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');

for y=1:ny
    waitbar(y/ny)

    % calculation of A

    %A = (M / UpLayVol) /...
    %    ( 4*n*pi()*(Dx*Dy)^(0.5) );
    

    
    A =  C0;

    % calculation of B
    
    
    B1 = 0.5 *  erfc(...
                (y - V * tsim_i) / ...
                (2 * (Dy * tsim_i)^0.5));
    
   C1 = (V^2 * tsim_i / (pi() * Dy))^0.5;
   
   C2 = exp(-(...
            (y - V * tsim_i)^2)/ ...
            (4 * Dy * tsim_i));
        
   D1  = 0.5 * (1 + (V * y) / Dy + (V^2 * tsim_i) / Dy);
   
   D2  = exp(V * y / Dy);
    
    D3 = erfc(...
            (y + V * tsim_i) / ...
            (2 * (Dy * tsim_i)^0.5));
        

    % Final calculation
    C = A * (B1 +  C1 * C2 - D1 * D2 * D3);

    % Saving results
    Cfinal(1,y) = C;

end

close(h)

disp('Analytical solution calculation: OK')

storage.C_analytical = Cfinal';
save AnalytSOL storage;

function Analytical_calc_1D_justAdv(tsim_i, mLayerDepth_summa_m)

load InputData storage
C0=storage.C0;
Dy=storage.Dy;
r=storage.r;
V=storage.V;
ny=storage.ny;

Cfinal=zeros(1,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');

for y=1:ny
    waitbar(y/ny)

    % calculation of A

    %A = (M / UpLayVol) /...
    %    ( 4*n*pi()*(Dx*Dy)^(0.5) );
  

    % Final calculation
    C = C0 * (y * mLayerDepth_summa_m - V * tsim_i);

    % Saving results
    Cfinal(1,y) = C;

end

close(h)

disp('Analytical solution calculation: OK')

storage.C_analytical = Cfinal';
save AnalytSOL storage;

function Analytical_calc_2D

load InputData storage
M=storage.M;
UpLayVol=storage.UpLayVol;
n=storage.n;
Dx=storage.Dx;
Dy=storage.Dy;
r=storage.r;
Xc=storage.Xc;
Yc=storage.Yc;
tsim=storage.tsim;
V=storage.V;
nx=storage.nx;
ny=storage.ny;

Cfinal=zeros(nx,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');
for x=1:nx
    waitbar(x/nx)
    for y=1:ny
        
    % calculation of A
    
    A = (M / UpLayVol) /...
        ( 4*n*pi()*(Dx*Dy)^(0.5) );

    % calculation of B
    B = exp(...
            V*(y-Yc)/ ...
            (2*Dy));

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

function Compare_plot(ti, n_ti, tsim_i,mLayerDepth_summa_m)

load InputData storage
InputData = storage;

load AnalytSOL storage
C_analytical = storage.C_analytical;

%load NumericalSOL storage
%C_numerical = storage.C_numerical;

colormap(gray)
cmap = colormap
cmap = flipud(cmap);

cm = [jet(64);gray(64)];

nsubplot_x = ceil((n_ti)^0.5);
nsubplot_y = ceil(n_ti / nsubplot_x);

if ti==1
     figure
end
%subplot(3,3,[1,2,4,5])
%[~, hh]=contour(C_analytical);
%set(gca, 'Ydir', 'reverse'); 
%set(hh,'LineColor','k','ShowText','on')
%grid on
%hold on

%{
pcolor(C_numerical)%,'ShowText','on');
shading flat
colormap(cmap)
grid on
%}

%colorbar('location','East')
%xlabel('Distance along the x-direction')
%ylabel('Distance along the y-direction')
%legend('Analytical','Numerical')

% Longitudinal profile
subplot(nsubplot_x,nsubplot_y,ti)
plot(C_analytical,(0:1:InputData.ny-1)/1000,'k--','linewidth',2)
ylim([0 InputData.ny/1000])
%xlim([0 2])
set(gca, 'Ydir', 'reverse'); 
xlabel('Concentration')
ylabel('Distance along the y-direction')
title(['Time since load start: ', num2str(tsim_i/(24*60*60)), ' days'])
legend('Analytical Solution')
set (gca,'Xdir','reverse')
grid on

% Transverse profile
%subplot(3,3,[7,8])
%plot(0:1:InputData.nx-1, C_analytical(round(InputData.Yc),:),'k','linewidth',1)
%xlim([1 InputData.nx])
%legend('Analytical Solution')
%ylabel('Concentration')
%xlabel('Distance along the x-direction')
%grid on

disp('ANALYTICAL and NUMERICAL solutions compared: OK')
disp('Finished')