% ME 614 CFD -- Course project
% Solve the 2d incompressible flow around a square
% --------------------------------------------------------------
% governing equations:
% u_t + u u_x + v u_y = -p_x + (u_xx + u_yy) / Re
% v_t + u v_x + v v_y = -p_y + (v_xx + v_yy) / Re
% u_x + v_y = 0
% --------------------------------------------------------------
% boundary equations:
% u, v = 0 on the upper and lower plate
% u = 4 y / H * (1 - y / H), v = 0, p_x = 0 on the inlet
% u_x = v_x = p = 0 on the outlet
% u = v = p_n = 0 on the square
% --------------------------------------------------------------

%% ------------------------------ geometry -------------------------------
clc; clear all; close all;

nd1 = 20;     % number of cells along the horizontal boundary of the square 
nd2 =20;    % number of cells along the vertical boundary of the square 
nx1 = 40;    % number of cells in the inlet region
nx2 = 400;    % number of cells in the outlet region
ny1 = 38;    % number of cells between bottom of square and the lower plate
ny2 = 52;    % number of cells between top of square and the top plate

nx = nx1 + nd1 + nx2;    % number of cells along the x-direction
ny = ny1 + nd2 + ny2;  % number of cells along the y-direction

dx = 1 / nd1;    dy = 1 / nd2;    
Re = 100;     % Reynolds number 
CFL = 0.8;  % CFL number
dt = CFL / (1/dx);
max_iter = 2500;

%% ----------------------------- initialization ----------------------------
ubar = zeros(1, ny+2);  % store the value of U in the inlet such that ubar_j = U_{0,j}
for j = 2 : (ny+1)
    ubar(j) = (2*j-3) / ny * (2 - (2*j-3) / ny);
end

ubar(1) = - ubar(2);
ubar(ny+2) = -ubar(ny+1);

U = zeros(nx+1, ny+2);  
V = zeros(nx+2, ny+1);  
P = zeros(nx+2, ny+2);
U(1, :) = ubar;

coe_p = zeros((nx+2) * (ny+2), (nx+2) * (ny+2));    % store the coefficient matrix of the possion equation of p
rhs_p = zeros((nx+2) * (ny+2), 1);      % coe_p * P = rhs_p

coe_u = zeros((nx+1) * (ny+2), (nx+1) * (ny+2));    % store the coefficient matrix of the possion equation of u
rhs_u = zeros((nx+1) * (ny+2), 1);

coe_v = zeros((nx+2) * (ny+1), (nx+2) * (ny+1));    % store the coefficient matrix of the possion equation of v
rhs_v = zeros((nx+2) * (ny+1), 1);

U_1 = zeros(nx+1, ny+2);    % U*
U_2 = zeros(nx+1, ny+2);    % U**

V_1 = zeros(nx+2, ny+1);    % V*
V_2 = zeros(nx+2, ny+1);    % V**

new_P = zeros(nx+2, ny+2);
new_U = zeros(nx+1, ny+2);
new_V = zeros(nx+2, ny+1);

x_ord = zeros(nx+1, ny+1); % for streamlines and contour plots
y_ord = zeros(nx+1, ny+1);
u_ord = zeros(nx+1, ny+1);
v_ord = zeros(nx+1, ny+1);
vort_ord = zeros(nx+1, ny+1);

temp_t = []; temp_v = []; temp_v2 = []; temp_v3 = []; % for analysis of Strouhal number

%% ----------------------------- start calculation ----------------------------
for iter = 1 : max_iter
    display(iter);
    % get U_{i,j}^*
    for i = 2 : nx
        for j = 2 : (ny+1)
            if i > nx1 && i < (nx1 + nd1 + 2) && j > (ny1 + 1) && j < (ny1 + nd2 + 2)
                ;
            else
                U_1(i, j) = dt / dx * max((U(i-1,j) + U(i,j)) / 2, 0) * U(i-1,j) + ...
                    dt / dx * max(-(U(i,j) + U(i+1,j)) / 2, 0) * U(i+1,j) + ...
                    dt / dy * max((V(i, j-1) + V(i+1, j-1)) / 2, 0) * U(i, j-1) + ...
                    dt / dy * max(-(V(i,j) + V(i+1, j)) / 2, 0) * U(i, j+1) + ...
                    (1 - dt / dx * (max((U(i,j) + U(i+1,j)) / 2, 0) + max(-(U(i-1,j) + U(i,j)) / 2, 0)) - ...
                    dt / dy * (max(-(V(i, j-1) + V(i+1, j-1)) / 2, 0) + max((V(i,j) + V(i+1, j)) / 2, 0))) * U(i,j);
            end
        end
    end
    
    U_1(1, :) = ubar;   % left boundary
    U_1(nx+1, 2:(ny+1)) = U_1(nx, 2:(ny+1)); % right boundary
    U_1(:, ny+2) = -U_1(:, ny+1);   % top boundary
    U_1(:, 1) = -U_1(:, 2); % bottom boundary
    
    % left and right boundary of square is zero velocity, so do not change
    U_1((nx1+2):(nx1+nd1), ny1+2) = - U_1((nx1+2):(nx1+nd1), ny1+1);    % bottom boundary
    U_1((nx1+2):(nx1+nd1), ny1+nd2+1) = - U_1((nx1+2):(nx1+nd1), ny1+nd2+2);    % top boundary
    
    % get U**
    for i = 2 : nx
        for j = 2 : (ny+1)
            if i > nx1 && i < (nx1 + nd1 + 2) && j > (ny1 + 1) && j < (ny1 + nd2 + 2)
                ;
            else
                m = i + (nx+1)*(j-1);
                coe_u(m, m) = 1 + 2 / Re * dt / dx^2 + 2 / Re * dt / dy^2;
                coe_u(m, m-1) = -dt / Re / dx^2;
                coe_u(m, m+1) = -dt / Re / dx^2;
                coe_u(m, m-nx-1) = -dt / Re / dy^2;
                coe_u(m, m+nx+1) = -dt / Re / dy^2;
                rhs_u(m) = U_1(i, j);
            end
        end
    end
    
    for j = 2 : (ny+1)  % left boundary
        m = 1 + (nx+1)*(j-1);   coe_u(m, m) = 1;    rhs_u(m) = ubar(j); 
    end
    
    for j = 2 : (ny+1)  % right boundary
        m = nx + 1 + (nx+1)*(j-1);  coe_u(m, m) = 1;    coe_u(m, m-1) = -1;
    end
    
    for i = 1 : (nx+1)  % top boundary
        m = i + (ny + 1) * (nx + 1);    coe_u(m, m) = 1; coe_u(m, m - nx - 1) = 1;
    end
    
    for i = 1 : (nx+1)  % bottom boundary
        m = i;  coe_u(m, m) = 1; coe_u(m, m + nx + 1) = 1;
    end
    
    for i = (nx1+1) : (nx1+nd1+1)
        for j = (ny1 + 2) : (ny1 + nd2 + 1)    
            m = i + (j-1) * (nx + 1);   coe_u(m, m) = 1;
            if j == (ny1+2) % bottom boundary of square
                  coe_u(m, m - nx - 1) = 1;
            elseif j == (ny1 + nd2 + 1) % top boundary of square
                  coe_u(m, m + nx + 1) = 1;
            end
        end
    end
    
    coe_u = sparse(coe_u);
    temp = coe_u \ rhs_u;
    
    for i = 1 : (nx + 1)
        for j = 1 : (ny + 2)
            U_2(i, j) = temp(i + (j-1) * (nx + 1));
        end
    end
    
    % get V_{i,j}^*
    for i = 2 : (nx+1)
        for j = 2 : ny
            if i > (nx1 + 1) && i < (nx1 + nd1 +2) && j > ny1 && j < (ny1 + nd2 + 2)
                ;
            else
                V_1(i, j) = dt / dx * max((U(i-1, j) + U(i-1, j+1)) / 2, 0) * V(i-1, j) + ...
                    dt / dx * max(-(U(i,j) + U(i, j+1)) / 2, 0) * V(i+1, j) + ...
                    dt / dy * max(-(V(i, j) + V(i, j+1)) / 2, 0) * V(i, j+1) + ...
                    dt / dy * max((V(i, j-1) + V(i,j)) / 2, 0) * V(i, j-1) + ...
                    (1 - dt / dx * (max((U(i,j) + U(i, j+1)) / 2, 0) + max(-(U(i-1, j) + U(i-1, j+1)) / 2, 0)) - ...
                    dt / dy * (max((V(i, j) + V(i, j+1)) / 2, 0) + max(-(V(i, j-1) + V(i,j)) / 2, 0))) * V(i,j);
            end
        end
    end
    
    V_1(1, :) = -V_1(2, :); % left boundary
    V_1(nx+2, :) = V_1(nx+1, :);    % right boundary
    V_1(nx1 + 2, (ny1+1):(ny1+nd2+1)) = -V_1(nx1 + 1, (ny1+1):(ny1+nd2+1));     % left boundary of square
    V_1(nx1 + nd1 + 1, (ny1+1):(ny1+nd2+1)) = -V_1(nx1 + nd1 + 2, (ny1+1):(ny1+nd2+1));     % right boundary of square
    
    % get V_{i,j}**
    for i = 2 : (nx+1)
        for j = 2 : ny
            if i > (nx1 + 1) && i < (nx1 + nd1 +2) && j > ny1 && j < (ny1 + nd2 + 2)
                ;
            else
                m = i + (j-1) * (nx + 2);
                coe_v(m, m-1) = -dt / dx^2 / Re;
                coe_v(m, m+1) = -dt / dx^2 / Re;
                coe_v(m, m-nx-2) = -dt / dy^2 / Re;
                coe_v(m, m+nx+2) = -dt / dy^2 / Re;
                coe_v(m, m) = 1 + 2 * dt / dx^2 / Re + 2 * dt / dy^2 / Re;
                rhs_v(m) = V_1(i, j);   
            end
        end
    end
    
    
    for j = 1 : (ny + 1)
        m = 1 + (j-1) * (nx + 2);   coe_v(m, m) = 1;    coe_v(m, m+1) = 1;  % left boundary
        m = nx + 2 + (j-1) * (nx+2);    coe_v(m, m) = 1;    coe_v(m, m-1) = -1; % right boundary
    end
    
    for i = 2 : (nx+1)
        m = i + ny * (nx + 2);  coe_v(m, m) = 1;    % top boundary
        m = i;  coe_v(m, m) = 1;    % bottom boundary
    end
    
    for i = (nx1 + 2) : (nx1 + nd1 + 1)
        for j = (ny1 + 1) : (ny1 + nd2 + 1)
            m = i + (j-1) * (nx + 2);   coe_v(m, m) = 1; 
            if i == (nx1 + 2)
                coe_v(m, m-1) = 1;  % left boundary of square
            elseif i == (nx1 + nd1 + 1)
                coe_v(m, m+1) = 1;  % right boundary of square
            end
        end
    end
    
    coe_v = sparse(coe_v);
    temp = coe_v \ rhs_v;
    for i = 1 : (nx + 2)
        for j = 1 : (ny + 1)
            V_2(i, j) = temp(i + (j-1) * (nx + 2));
        end
    end
    
    % get p_{n+1}
    for i = 2 : (nx + 1)
        for j = 2 : (ny + 1)
            if i > (nx1+1) && i < (nx1 + nd1 + 2) && j > (ny1 + 1) && j < (ny1 + nd2 + 2)
                ;
            else
                m = i + (j-1) * (nx + 2);
                coe_p(m, m-1) = dt / dx^2;   coe_p(m, m+1) = dt / dx^2;
                coe_p(m, m - nx - 2) = dt / dy^2;    coe_p(m, m + nx + 2) = dt / dy^2;
                coe_p(m, m) = - 2 * dt / dx^2 - 2 * dt / dy^2;
                rhs_p(m) = (U_2(i, j) - U_2(i-1, j)) / dx + (V_2(i, j) - V_2(i, j-1)) / dy; 
            end
        end
    end

    for j = 1 : (ny + 2)
        m = 1 + (j - 1) * (nx + 2);     coe_p(m, m) = 1;    coe_p(m, m+1) = -1;     % left boundary
        m = nx + 2 + (j-1) * (nx + 2);  coe_p(m, m) = 1;    % right boundary
    end
    
    for i = 2 : (nx+1)
        m = i;  coe_p(m, m) = 1;    coe_p(m, m + nx + 2) = -1;  % bottom boundary
        m = i + (ny + 1) * (nx + 2); coe_p(m, m) = 1;   coe_p(m, m - nx - 2) = -1;  % top boundary
    end
    
    for i = (nx1 + 2) : (nx1 + nd1 + 1)
        for j = (ny1 + 2) : (ny1 + nd2 + 1)
            m = i + (j-1) * (nx + 2);   coe_p(m, m) = 1;
            if i == (nx1 + 2)
                coe_p(m, m-1) = -1; % left boundary of square
            elseif i == (nx1 + nd1 + 1)
                coe_p(m, m+1) = -1; % right boundary of square
            elseif j == (ny1 + 2)
                coe_p(m, m - nx - 2) = -1;  % bottom boundary of square
            elseif j == (ny1 + nd2 + 1)
                coe_p(m, m + nx + 2) = -1;  % top boundary of square
            end
        end
    end
    
    coe_p = sparse(coe_p);
    temp = coe_p \ rhs_p;
    
    for i = 1 : (nx + 2)
        for j = 1 : (ny + 2)
            new_P(i, j) = temp(i + (j-1) * (nx +2));
        end
    end
    
    % get U^{n+1}
    for i = 2 : nx
        for j = 2 : (ny+1)
            if i > nx1 && i < (nx1 + nd1 + 2) && j > (ny1 + 1) && j < (ny1 + nd2 + 2)
                ;
            else
                new_U(i, j) = U_2(i, j) + dt / dx * (new_P(i,j) - new_P(i+1, j));
            end
        end
    end
    
    new_U(1, :) = ubar;   % left boundary
    new_U(nx+1, 2:(ny+1)) = new_U(nx, 2:(ny+1)); % right boundary
    new_U(:, ny+2) = -new_U(:, ny+1);   % top boundary
    new_U(:, 1) = -new_U(:, 2); % bottom boundary
    
    % left and right boundary of square is zero velocity, so do not change
    new_U((nx1+2):(nx1+nd1), ny1+2) = - new_U((nx1+2):(nx1+nd1), ny1+1);    % bottom boundary
    new_U((nx1+2):(nx1+nd1), ny1+nd2+1) = - new_U((nx1+2):(nx1+nd1), ny1+nd2+2);    % top boundary
    
    % get V^{n+1}
    for i = 2 : (nx+1)
        for j = 2 : ny
            if i > (nx1 + 1) && i < (nx1 + nd1 +2) && j > ny1 && j < (ny1 + nd2 + 2)
                ;
            else
                new_V(i, j) = V_2(i, j) + dt / dy * (new_P(i, j) - new_P(i, j+1));
            end
        end
    end
    
    new_V(1, :) = -new_V(2, :); % left boundary
    new_V(nx+2, :) = new_V(nx+1, :);    % right boundary
    new_V(nx1 + 2, (ny1+1):(ny1+nd2+1)) = -new_V(nx1 + 1, (ny1+1):(ny1+nd2+1));     % left boundary of square
    new_V(nx1 + nd1 + 1, (ny1+1):(ny1+nd2+1)) = -new_V(nx1 + nd1 + 2, (ny1+1):(ny1+nd2+1));     % right boundary of square
    
    U = new_U;  V = new_V;  P = new_P;
    
    %%  -------------------------- post process -----------------------
    temp_t = [temp_t, iter * dt];
    temp_v = [temp_v, V(nx1 + nd1*2 + 1, ny1 + nd2 + 1)];
    temp_v2 = [temp_v2, V(nx1 + nd1 * 4 + 1, ny1 + nd2 +1 )];
    temp_v3 = [temp_v3, V(nx1 + nd1 * 11 + 1, ny1 + nd2 +1 )];
    
    if mod(iter, 100) == 0 
        for i = 1 : (nx + 1)
            for j = 1 : (ny + 1)
                x_ord (i, j) = (i-1) * dx;  y_ord(i, j) = (j-1) * dy;
                if i == 1 || i == (nx+1) || j == 1 || j == (ny+1)
                    u_ord(i, j) = NaN; v_ord(i, j) = NaN; vort_ord(i, j) = NaN; 
                elseif i > nx1 && i < (nx1+nd1+2) && j > ny1 && j < (ny1 + nd2 + 2)
                    u_ord(i, j) = NaN; v_ord(i, j) = NaN; vort_ord(i, j) = NaN; 
                else
                    u_ord(i, j) = (U(i, j) + U(i, j+1)) / 2;
                    v_ord(i, j) = (V(i, j) + V(i+1, j)) / 2;
                    vort_ord(i, j) = (V(i+1, j) - V(i, j)) / dx - (U(i, j+1) - U(i, j)) / dy;
                end
            end
        end
       
        figure(1);
        gap = 1;
        contourf(x_ord(1:gap:end, 1:gap:end), y_ord(1:gap:end, 1:gap:end), u_ord(1:gap:end, 1:gap:end), 'LineColor', 'none');
        colorbar;
        caxis([-0.3,1])
        axis equal
        title([sprintf('t = %0.1f',iter*dt), ',  Re = ', num2str(Re), sprintf(',  Ly1/D =%0.1f',ny1/nd2)], 'fontsize', 40, 'fontweight', 'bold');
        xlabel('x', 'fontsize', 40)
        ylabel('y', 'fontsize', 40)
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times','fontsize',25)
        set(0, 'DefaultLineLineWidth', 3);
        drawnow
        
    end                 
    
end

%% export time and v for analysis of Strouhal number
% write data of time and v_velocity to txt for comparison
% fid = fopen('time_coarse.txt','wt');
% fprintf(fid, '%.8f\n', temp_t);
% fclose(fid);
% 
% fid = fopen('v_velocity_coarse.txt','wt');
% fprintf(fid, '%.8f\n', temp_v);
% fclose(fid);

%% plot streamlines
% gap = 3;        
% h1 = quiver(x_ord(1:gap:end, 1:gap:end), y_ord(1:gap:end, 1:gap:end), u_ord(1:gap:end, 1:gap:end), v_ord(1:gap:end, 1:gap:end));
% set(h1, 'AutoScale', 'on', 'AutoScaleFactor', 0.2);
% axis equal
% title([sprintf('t=%0.1f',iter*dt), ', Re = ', num2str(Re), sprintf(',  Ly1/D =%0.1f',ny1/nd2)], 'fontsize', 40, 'fontweight', 'bold');
% xlabel('x', 'fontsize', 40)
% ylabel('y', 'fontsize', 40)
% xlim([0,14])
% ylim([0,5])
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',25)
% set(0, 'DefaultLineLineWidth', 2);
% % 
% starty = 0:0.2:5;
% startx = zeros(size(starty)) + 0.05;
% 
% streamline(x_ord', y_ord', u_ord', v_ord', startx, starty)
