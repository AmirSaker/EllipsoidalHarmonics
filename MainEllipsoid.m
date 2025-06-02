% Main file for Ellipsoidal Harmonic Surface Reconstruction
% Copyright (c) 2025, SAKER AMIR
% Based on: 
% - Lamé function implementations (computeLameK/L/M/N.m)
% - Coordinate transformations (cartesian2Ellipsoidal.m, ellipsoidal2Cart.m)
% - Compute Surface Basis and reconstruction

clear; clc; close all;
%% Add necessary paths
addpath(genpath('C:\Users\PC\Desktop\UL\M2\Mon année\TIR\CODE\Ellipsoidal-harmonics-decomposition-main\code'));
addpath(genpath('C:\Users\PC\Desktop\UL\M2\Mon année\TIR\CODE\spheroidal-harmonics-decomposition-main\stlTools'));
addpath('C:\Users\PC\Desktop\UL\M2\Mon année\TIR\CODE\spheroidal-harmonics-decomposition-main\Input_geometry');
addpath('C:\Users\PC\Desktop\UL\M2\Mon année\TIR\CODE\spheroidal-harmonics-decomposition-main\Input_geometry\data');
addpath('C:\Users\PC\Desktop\UL\M2\Mon année\TIR\CODE\spheroidal-harmonics-decomposition-main\Input_geometry\data');
addpath('C:\Users\PC\Desktop\UL\M2\Mon année\TIR\CODE\spheroidal-harmonics-decomposition-main\Input_geometry\test_sable');

%% Input surface
filename = 'B1.stl'; % Example file
[~, baseName, ext] = fileparts(filename);

% Load mesh
if strcmpi(ext, '.mat')
    load(filename);
elseif strcmpi(ext, '.stl')
    [v, f] = stlRead(filename);
end

% Visualization
plot_mesh(v, f);
view([-70 10]);
title('Input Surface');
rec_refinement = 5;
%% Surface registration and ellipsoid fitting
[v_reg, centroid, Rots] = register_surface(v, false);

% Fit general ellipsoid (using your ft_ellipsoid.m)
[a1, a2, a3] = fit_ellipsoid(v_reg,false); 

% Calculate focal parameters
h1 = sqrt(a2^2 - a3^2);
h2 = sqrt(a1^2 - a3^2);
h3 = sqrt(a1^2 - a2^2);

x=v_reg(:,1);
y=v_reg(:,2);
z=v_reg(:,3);

% Visualize fitting
scatter3(v_reg(:,1), v_reg(:,2), v_reg(:,3), 'k.');
hold on;
[xe,ye,ze] = ellipsoid(0,0,0,a1,a2,a3,50);
surf(xe, ye, ze, 'FaceAlpha', 0.3);
axis equal; title('Ellipsoid Fitting');

%% Surface Parameterization 
map = ellipsoidal_conformal_map(v_reg, f, a1, a2, a3);
plot_mesh(map, f);
view([-70 10]);
title('Ellipsoidal conformal parameterization');

angle_distortion = angle_distortion(v_reg, f, map);
area_distortion = area_distortion(v_reg, f, map);

figure;
histogram(angle_distortion,-180:1:180);
xlim([-180 180])
title('Angle Distortion');
xlabel('Angle difference (degree)')
ylabel('Number of angles')
set(gca,'FontSize',12);

figure;
histogram(area_distortion,-5:0.1:5);
xlim([-5 5])
title('Area Distortion');
xlabel('log(final area/initial area)')
ylabel('Number of faces')
set(gca,'FontSize',12);

%% Compute Ellipsoidal Harmonic Basis
max_degree = 1; % degree
% [rho, mu, nu] = cartesian2Ellipsoidal(a1, a2, a3, x, y, z);
% Build basis matrix
D = ComputeBasis(map, a1, a2, a3, max_degree);

%% Harmonic Decomposition
% Solve for coefficients (separately for x,y,z)
% coeffs = zeros(size(D,2), 3);
% coeffs = D \ v_reg;
%% Harmonic Decomposition 2
%%% Solve for coefficients (separately for x,y,z)
coeffs = pinv(D) * v_reg;  % More stable numericaly
%% Shape Descriptors Analysis
% Dl = zeros(3, max_degree+1);
% for n = 0:max_degree
%     for m = 1:2*n+1
%         idx = n^2 + n + m + 1;
%         Dl(:,n+1) = Dl(:,n+1) + abs(coeffs(idx,:)).^2;
%     end
%     Dl(:,n+1) = sqrt(Dl(:,n+1));
% end
% 
% % Plot descriptors
% figure;
% loglog(0:max_degree, Dl(1,:), 'r', 0:max_degree, Dl(2,:), 'g', 0:max_degree, Dl(3,:), 'b');
% legend('X','Y','Z'); title('Shape Descriptors');

%% Surface Reconstruction
[rec_v, rec_f] = icosahedron_sphere(rec_refinement, false); % Base mesh

% Scale to ellipsoid
rec_v(:,1) = rec_v(:,1) * a1;
rec_v(:,2) = rec_v(:,2) * a2;
rec_v(:,3) = rec_v(:,3) * a3;

% Compute basis for reconstruction points
%D_rec = computeSurfaceBasis(rec_v, a1, a2, a3, max_degree);
D_rec = ComputeBasis(rec_v, a1, a2, a3, max_degree);
% Reconstruct
rec_v = D_rec * coeffs;
rec_v = real(D_rec * coeffs(1:(max_degree+1)^2, :));

% Transform back to original space
rec_v = (Rots' * rec_v')';
for i=1:3  
    rec_v(:,i) = rec_v(:,i)+ centroid(i);
end

paraview_patch(rec_v, rec_f)
%% Output and Visualization
stlWrite(fullfile('output', ['reconstructed_', baseName, '.stl']), rec_f, rec_v);
figure;
plot_mesh(rec_v, rec_f);
title('Reconstructed Surface');