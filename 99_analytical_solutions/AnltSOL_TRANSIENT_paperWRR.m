
% 2-dimensional analytical solution - Transient (using the analytical solution for
% instantaneous point source integrated with respect to time)

function AnltSOL_TRANSIENT_paperWRR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test = 8; % 2, 4, 6, 8, 

% General %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny = 600;
C0 = 2;
mLayerDepth_summa_m         = 0.006;                        % m
Dy      = 0.0001;      	% longiitudinal dispersivity - y-direction m2/s
specificYield = 0.25;
iLayerLiqFluxSoil_summa_m_s = 1.0359700931174889e-05 / specificYield;       % m/s
%iLayerLiqFluxSoil_summa_m3  = 0.30488599840555514;         % m3/s
V = iLayerLiqFluxSoil_summa_m_s;
%L = V * 60 * 60 * 24 * 20 / 1000;

tsim = [...
        ...0,...
        60 * 60 * 24 * 15,...
        60 * 60 * 24 * 70,...
        60 * 60 * 24 * 120,...
        ];


% Storing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storage.C0 = C0;
storage.Dy = Dy;
storage.tsim = tsim;
storage.V = V;
storage.ny = ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_ti = numel(tsim);
for ti = 1:n_ti
    
    tsim_i = tsim(ti);
    
    if test == 4 || test == 8 
        
        % Set reation rate
        if test == 4; r = 0; end                % conservative
        if test == 8; r = 0.01/(60*60*24); end  % linear decay
        storage.r = r;
        save InputData storage
        
        % Call analytical solver
        Analytical_calc_1D_contS(tsim_i)
        
    elseif test ==2 || test == 6
        
    end
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
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Compare_plot(ti, n_ti, tsim_i,mLayerDepth_summa_m)

load InputData storage
InputData = storage;

load AnalytSOL storage
C_analytical = storage.C_analytical;


colormap(gray)
cmap = colormap
cmap = flipud(cmap);

cm = [jet(64);gray(64)];

nsubplot_x = ceil((n_ti)^0.5);
nsubplot_y = ceil(n_ti / nsubplot_x);

if ti==1
     figure
end

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


disp('ANALYTICAL and NUMERICAL solutions compared: OK')
disp('Finished')