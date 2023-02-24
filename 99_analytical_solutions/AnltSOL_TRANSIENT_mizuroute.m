
% 2-dimensional analytical solution - Transient (using the analytical solution for
% instantaneous point source integrated with respect to time)

function AnltSOL_TRANSIENT_mizuroute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test = 2; % 2_nrTrans_instS_PorMedia
% test = 4; % 4_nrTrans_contS_PorMedia
% test = 6; % 6_nrTrans_instS_PorMedia_linDecay
% test = 8; % 8_nrTrans_contS_PorMedia_linDecay
 
test = 8; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General mdoel setup %%%%%%%%%%%%%%%%%%%%%
% Don't change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dist_between_nodes = cumsum([80.348, 113.724, 134.015]); % km
timeDays_travel_between_nodes = [120, 250, 370]; % days

% Preliminary calcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Time steps
if test == 4 || test == 8 
    %x_request = [6430, 6753, 6893];
    %x_dist_request = [1, 400, 500];
    %x_dist_request = cumsum([80.348, 113.724, 134.015]); % km
    %x_dist_request = [150]; %km
    start_day =  '2017-08-28 12:15:00';
    tsim_nDays_sinceStart = 1100; % time reachis nodes: 120, 250 and 370;
    
    if test == 4
        timeDays_travel_between_nodes = timeDays_travel_between_nodes - 20; % days
    end
    
    V_kmday = dist_between_nodes ./ (timeDays_travel_between_nodes);
    Dy = mean(V_kmday)*0.5;

elseif test == 2 || test == 6 
    tsim = [...
            60 * 60 * 24 * 55,...
            60 * 60 * 24 * 85,...
            60 * 60 * 24 * 120,...
            ];
end

% Storing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storage.Dy = Dy;
storage.tsim_nDays_sinceStart = tsim_nDays_sinceStart;
storage.V_kmday = V_kmday;
%storage.ny = ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_x_request = numel(dist_between_nodes);
for x_reqi = 1:n_x_request
    
    
    if test == 4 || test == 8 
        
        % Input concentration (continuous)
        C0 = 10;
         
        % Set reation rate
        if test == 4; r = 0; end                % conservative
        if test == 8; r = 0.01; end  % linear decay
        storage.C0 = C0;
        storage.r = r;
        save InputData storage;
        
        % Call analytical solver
        Analytical_calc_1D_contS(dist_between_nodes,x_reqi,tsim_nDays_sinceStart)
        
    elseif test == 2 || test == 6
        
        % Input mass (instantaneous)
        M = 350;
        
        % Mass of the receiving node
        mLayerVolFracLiq = [0.2007, 0.2014, 0.2020, 0.2027, 0.2034, 0.2041, 0.2048, 0.2055, 0.2063, 0.2070, 0.2077, 0.2085, 0.2093, 0.2101, 0.2109, 0.2117, 0.2125, 0.2133, 0.2141, 0.2150, 0.2159, 0.2167, 0.2176, 0.2185, 0.2194, 0.2204, 0.2213, 0.2223, 0.2232, 0.2242, 0.2252, 0.2262, 0.2272, 0.2283, 0.2293, 0.2304, 0.2315, 0.2326, 0.2337, 0.2348, 0.2360, 0.2372, 0.2384, 0.2396, 0.2408, 0.2420, 0.2433, 0.2446, 0.2459, 0.2472, 0.2485, 0.2499, 0.2513, 0.2527, 0.2541, 0.2555, 0.2570, 0.2584, 0.2599, 0.2615, 0.2630, 0.2646, 0.2662, 0.2678, 0.2694, 0.2711, 0.2727, 0.2744, 0.2761, 0.2779, 0.2796, 0.2814, 0.2832, 0.2850, 0.2869, 0.2888, 0.2906, 0.2925, 0.2945, 0.2964, 0.2984, 0.3003, 0.3023, 0.3043, 0.3063, 0.3083, 0.3104, 0.3124, 0.3145, 0.3165, 0.3186, 0.3206, 0.3227, 0.3247, 0.3268, 0.3288, 0.3308, 0.3328, 0.3348, 0.3368];
        hru_area_m2      = 32700;
        
        mLayerVolFracWat_summa_m3 = mLayerVolFracLiq * hru_area_m2 * mLayerDepth_summa_m;
        mLayerVolFracWat_summa_mm_m2 = mLayerVolFracLiq * hru_area_m2 * mLayerDepth_summa_m * 1000;
        UpLayVol_mm_m2 = mLayerVolFracWat_summa_m3 * 1000; % to mm*m2
        
        iLayerLiqFluxSoil_summa_mm_s = iLayerLiqFluxSoil_summa_m_s * 1000;
        
        % Set reaction rate
        if test == 2; r = 0; end                % conservative
        if test == 6; r = 0.01/(60*60*24); end  % linear decay
        
        storage.M = M;
        storage.UpLayVol_mm_m2 = UpLayVol_mm_m2;
        storage.mLayerVolFracWat_summa_mm_m2 = mLayerVolFracWat_summa_mm_m2;
        storage.mLayerVolFracLiq = mLayerVolFracLiq;
        storage.r = r;
        %storage.mLayerDepth_summa_mm = mLayerDepth_summa_mm;
        storage.iLayerLiqFluxSoil_summa_mm_s = iLayerLiqFluxSoil_summa_mm_s;
        storage.hru_area_m2 =  hru_area_m2;
        save InputData storage;
        
        % Call analytical solver
        Analytical_calc_1D_instS(tsim_i)
        
    end
    %Analytical_calc_1D_2(tsim_i)
    %Analytical_calc_1D_3(tsim_i)
    %Analytical_calc_1D_justAdv(tsim_i, mLayerDepth_summa_m)
    %Analytical_calc_2D
    %Numerical_calc
    Compare_plot(dist_between_nodes,x_reqi,...
                 start_day)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1D, continuous source, Eq. 60 of 
% USGS Chapter B7
% ANALYTICAL SOLUTIONS FOR ONE-, TWO-, AND THREE-DIMENSIONAL SOLUTE
%TRANSPORT IN GROUND-WATER SYSTEMS WITH UNIFORM FLOW
function Analytical_calc_1D_contS(...
    dist_between_nodes,...
    x_reqi,...
    tsim_nDays_sinceStart)

load InputData storage
C0=storage.C0;
Dy=storage.Dy;
r=storage.r;
V_kmday=storage.V_kmday;

Cfinal=zeros(1,tsim_nDays_sinceStart);

x_request_i = dist_between_nodes(x_reqi);

V = V_kmday(x_reqi);
V = mean(V_kmday);

h=waitbar(0,'ANALYTICAL solution: calculating...');

for t_i=1:tsim_nDays_sinceStart
    waitbar(t_i/tsim_nDays_sinceStart)
    
    % Days since start (in seconds)
    tsim_i = t_i;

    % calculation of A

    %A = (M / UpLayVol) /...
    %    ( 4*n*pi()*(Dx*Dy)^(0.5) );
    A =  C0 / 2;

    % calculation of B
    U = (V^2 + 4*r*Dy)^0.5;
    
    B1 = exp(...
            x_request_i * (V - U)/ ...
            (2 * Dy)...
            );
        
    B2 = erfc(...
            (x_request_i - (U * tsim_i) ) / ...
            (2 * (Dy * tsim_i)^0.5));
    
    C1 = exp(...
            x_request_i * (V + U)/ ...
            (2 * Dy));
    
    C2 = erfc(...
            (x_request_i + (U * tsim_i) ) / ...
            (2 * (Dy * tsim_i)^0.5));

    % Final calculation
    C = A * (B1 * B2 + C1 * C2);

    % Saving results
    Cfinal(1,t_i) = C;

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


% 1D, instantaneous source, Eq. 76 of 
% USGS Chapter B7
% ANALYTICAL SOLUTIONS FOR ONE-, TWO-, AND THREE-DIMENSIONAL SOLUTE
%TRANSPORT IN GROUND-WATER SYSTEMS WITH UNIFORM FLOW
function Analytical_calc_1D_instS(tsim_i)

load InputData storage
M = storage.M;
UpLayVol_mm_m2 = storage.UpLayVol_mm_m2;
mLayerVolFracLiq = storage.mLayerVolFracLiq;
Dy = storage.Dy;
r = storage.r;
V_kmday = storage.V_kmday;
ny = storage.ny;
hru_area_m2 = storage.hru_area_m2;
mLayerDepth_summa_mm = storage.mLayerDepth_summa_mm;
iLayerLiqFluxSoil_summa_mm_s = storage.iLayerLiqFluxSoil_summa_mm_s;

% Loading point
Yc=1;
% porosity
n = 1; % not needed because it's accounted for in the intersticil flow velocity

% Mass injection unit of aquifer thickness (and area because it's 1D)
% M (g/m3), Vol = m3/s, M_unitVol = 
UpLayVol_m3 = UpLayVol_mm_m2 / 1000;
M_unitVol = M / UpLayVol_m3(1);
M_unitVol = M * (iLayerLiqFluxSoil_summa_mm_s);

M_unitVol = M / ((mLayerDepth_summa_mm / 1000) * mLayerVolFracLiq(1) * hru_area_m2);
M_unitVol = M_unitVol * ( 4 * n * pi() * Dy );

M_unitVol = 1.8/1000 * ( 4 * n * pi() * Dy ) / 1.75; %=1.2925e-06%


M_unitVol = (2 * iLayerLiqFluxSoil_summa_mm_s) ;

M_unitVol = (0.006 * iLayerLiqFluxSoil_summa_mm_s / 6);

Cinit = (2 * iLayerLiqFluxSoil_summa_mm_s * 900 * hru_area_m2/1000) ...
        / (iLayerLiqFluxSoil_summa_mm_s * 900 * hru_area_m2/1000 + UpLayVol_m3(1)/6);
M_unitVol = Cinit * V;
M_unitVol = M_unitVol / 3.2%* 1.8e6 
%M_unitVol = 1.2925e-06;

Cfinal=zeros(1,ny);

h=waitbar(0,'ANALYTICAL solution: calculating...');

for y=1:ny
    
    waitbar(y/ny)

    % calculation of A

    A = (M_unitVol) /...
        ( 4 * n * pi() * Dy);

    % calculation of B
    B = exp(...
            V*(y-Yc)/ ...
            (2*Dy));

    D = exp( ...
            -( V^2 / (4*Dy) + r) * tsim_i ...
            -( (y-Yc)^2 ./ (4 * Dy * tsim_i))) ;

    % Final calculation
    C = A * B * D;

    % Saving results
    Cfinal(1,y) = C;
   
end
close(h)

% Get the effect of dilution because of different layer volumes
% 1) Using cubic spline to extrapolate for all sub-layers because the
% analytical solution has 1 mm of resolution, while the model has 6 mm
mLayerVolFracLiq_interp = spline(...
    1:6:600,...
    mLayerVolFracLiq,...
    1:1:600);
%figure;plot(mLayerVolFracLiq_interp^0.5,[1:1:600]/1000); set(gca,"Ydir","reverse"); set(gca,"Xdir","reverse");

Cfinal_corr = Cfinal * mLayerVolFracLiq_interp(1) ./ mLayerVolFracLiq_interp;

disp('Analytical solution calculation: OK')

storage.C_analytical = Cfinal_corr';
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

function Compare_plot(dist_between_nodes,x_reqi,...
                 start_day)

load InputData storage;
InputData = storage;

load AnalytSOL storage;
C_analytical = storage.C_analytical;


colormap(gray);

n_x_request = numel(dist_between_nodes);

nsubplot_x = ceil((n_x_request)^0.5);
nsubplot_y = ceil(n_x_request / nsubplot_x);

if x_reqi==1
     figure
end

% Longitudinal profile
subplot(nsubplot_x,nsubplot_y,x_reqi)
plot(C_analytical,'k--','linewidth',2)
%ylim([0 InputData.ny/1000])
%xlim([0 2])
%set(gca, 'Ydir', 'reverse'); 
ylabel('Concentration')
xlabel(['Time since load start: ', start_day, ' days'])
title(['River reach: ', num2str(x_reqi)])
legend('Analytical Solution')
%set (gca,'Xdir','reverse')
grid on
hold on


disp('ANALYTICAL and NUMERICAL solutions compared: OK')
disp('Finished')