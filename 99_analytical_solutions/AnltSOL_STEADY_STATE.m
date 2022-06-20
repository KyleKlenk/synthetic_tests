
% 2-dimensional analytical solution - STEADY-STATE (using the analytical solution for
% instantaneous point source integrated with respect to time)

clc
close all
clear all

% input data
C0 = 200;     % concentration of the injected water - units: ppm
Q = 1;         % injection rate - units: 1m3/day
n = 0.1;       % porosity
Dx = 12.5;      % longiitudinal dispersivity - y-direction 
Dy = 12.5;      % longiitudinal dispersivity - y-direction 
r = 0;         % reaction rate
Xc = 100;        % injection location location
Yc = 200;        % injection location location
Zgw = 10;      % aquifer thinkness
V = 0.3;        % GW seepage flow velocity (m/day) (average intersticial velocity)

nx = 400;
ny = 400;

Cfinal=zeros(nx,ny);

for x=1:nx
    
    for y=1:ny
        
    % calculation of A
    Qprime = Q/Zgw; %fluid injection rate per unit thinkness of aquifer
    
    A = C0*Qprime/(2*n*pi()*(Dx*Dy)^(0.5));

    % calculation of B
    B = exp(V*(x-Xc)/(2*Dx));

    % Calculation of D - Bessel function
    B1 = V^2/(4*Dx)+r;
    B21 = (x-Xc)^2./Dx;
    B22 = (y-Yc)^2./Dy;
    B2 = (B21+B22);
    BesselArg = (B1*B2).^0.5;

    D= besselk(0, BesselArg);

    % Final calculation
    Cfinal(x,y) = A * B * D;
    
    end
end

% Plotting
figure

colormap(gray)
cmap = colormap
cmap = flipud(cmap)

% Contour
subplot(3,3,[1,2,4,5])
contour(Cfinal,'ShowText','on');
xlabel('Distance along the x-direction')
ylabel('Distance along the y-direction')
colormap
colorbar('location','East')
grid on
legend('Analytical Solution')

% Longitudinal profile
subplot(3,3,[3,6])
plot(Cfinal(:,Yc),1:1:ny,'k','linewidth',1)
ylim([1 ny])
xlabel('Concentration')
ylabel('Distance along the y-direction')
legend('Analytical Solution')
set (gca,'Xdir','reverse')
grid on

% Transverse profile
subplot(3,3,[7,8])
plot(1:1:nx, Cfinal(Xc,:),'k','linewidth',1)
xlim([1 nx])
legend('Analytical Solution')
ylabel('Concentration')
xlabel('Distance along the x-direction')
grid on