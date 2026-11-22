% first clean figures and memory from all previous calculations
clc;
clear;
%clf;

N = 20000; % time

Delta = 0.01; % step size
D_rot = 0.01; % rotational diffusion
R = 200; % number of rods (displayed as red dots)
length_of_rod = 0.005;

global flow_rate;
flow_rate = 1.0;
frame_frequency = 100;


% ellipsoidal parameters
global a_ellipse;
global b_ellipse;
a_ellipse = 2.0;
b_ellipse = 1.0;

hist_resolution = int16(sqrt(R));
n_time = zeros(R);

velocity_matrix = zeros(R);

x = zeros(R,N);
y = zeros(R,N);
z = zeros(R,N);

rod_col_size = N/frame_frequency;
x_rod = zeros(R,rod_col_size);
y_rod = zeros(R,rod_col_size);
z_rod = zeros(R,rod_col_size);


theta = zeros(R,N);
phi = zeros(R,N);

for r=1:R
    % initial conditions
    frame_indicator=0;

    r_x0 = a_ellipse*(2*rand()-1);
    r_y0 = b_ellipse*(2*rand()-1);

    while (r_x0^2/(a_ellipse^2) + r_y0^2/(b_ellipse^2) >= 1) 
        r_x0 = a_ellipse * (2*rand()-1);
        r_y0 = b_ellipse * (2*rand()-1);
    end

    x0 = r_x0;
    y0 = r_y0;
    z0 = -0.01;

    
    phi0 = 2*pi*rand(); % angle on the xy-plane
    theta0 = (pi/2) *rand(); % angle on the xz and yz plane
    
    % EULERS METHOD 
    x(r,1) = x0 + Delta*x_deriv(phi0, theta0);
    y(r,1) = y0 + Delta*y_deriv(phi0, theta0);
    z(r,1) = z0 + Delta*z_deriv(theta0, x0, y0);
    theta(r,1) = theta0 + Delta*theta_deriv(theta0, phi0, x0, y0);
    phi(r,1) = phi0 + Delta*0;
    
    n_time(r) = int16(N * rand()) + 1;

    for n=2:N
        % Euler's method
        r1 = rand();
        r2 = rand();

    
        x(r,n) = x(r,n-1) + Delta*x_deriv(phi(r,n-1), theta(r,n-1));
        y(r,n) = y(r,n-1) + Delta*y_deriv(phi(r,n-1), theta(r,n-1));
        z(r,n) = z(r,n-1) + Delta*z_deriv(theta(r,n-1), x(r,n-1), y(r,n-1));

        theta(r,n) = theta(r,n-1) +...
            Delta*...
            theta_deriv(theta(r,n-1),phi(r,n-1), x(r,n-1), y(r,n-1))+...
            sqrt(2*D_rot*Delta)*(2*r1-1);

        phi(r,n) = phi(r,n-1) + Delta*0 +...
            sqrt(2*D_rot*Delta)*(2*r2-1)/sin(theta(r,n-1));

        if n > n_time(r)
            % adjustment of trajectory due to the wall x^2+2*y^2=1
            if ((x(r,n))^2/(a_ellipse^2) + (y(r,n))^2/(b_ellipse^2) > 1)
                d = sqrt((x(r,n))^2/(a_ellipse^2) +...
                    (y(r,n))^2/(b_ellipse^2));
                x(r,n) = x(r,n)/d;
                y(r,n) = y(r,n)/d;
            end
    
            if (n ~= 1 && mod(n,frame_frequency)==1)
                frame_indicator = frame_indicator+1;
                x_rod(r,frame_indicator) = x(r,n);
                y_rod(r,frame_indicator) = y(r,n);
                z_rod(r,frame_indicator) = z(r,n);
            end
    
            if(n == N)
                x_rod(r,end) = x(r,n);
                y_rod(r,end) = y(r,n);
                z_rod(r,end) = z(r,n);
            end
        end

    end


end


%Calcuate the velocities of each rod
for k=1:R
    total_distance = z_rod(k,end);
    velocity_matrix(k) = calc_velocity(total_distance, N);
end

avg_velocity = calc_avg_velocity(R,velocity_matrix);

maxz=max(z_rod,[],"all");
h_z = 1;
z_resolution = int16(maxz/h_z)+1;

for it=2:N/frame_frequency
    
    for k=1:z_resolution
        z_pdf(k) = 0;  
    end

    for r=1:R
        % calculating density, 
        for z_body = (z_rod(r,it)-length_of_rod):h_z:(z_rod(r,it)+length_of_rod) 
            if (z_rod(r,it)>0)
                z_index = int16(z_body/h_z)+1;
                if (z_index<z_resolution)&&(z_index>0)
                    z_pdf(z_index) = z_pdf(z_index) + 1.0;
                end
            end
        end
    end
    

    figure(1)
    clf()
    z_hist = h_z/2.0:h_z:(z_resolution*h_z);
    plot(z_hist,z_pdf,'black'); hold on;
    xlabel("z");
    axis([0 maxz 0 max(z_pdf)+10]);
    ch=sprintf("Time %f",it*frame_frequency*Delta_t); 
    title(ch);
end









% FUNCTION DEFINITIONS

function [] = mishas_histogram(z_rod,frame_frequency,Delta_t,N)

    maxz=max(z_rod,[],"all");
    h_z = 1;
    z_resolution = int16(maxz/h_z)+1;
    
    for it=2:N/frame_frequency
        
        for k=1:z_resolution
            z_pdf(k) = 0;  
        end
    
        for r=1:R
            % calculating density, 
            for z_body = (z_rod(r,it)-length_of_rod):h_z:(z_rod(r,it)+length_of_rod) 
                if (z_rod(r,it)>0)
                    z_index = int16(z_body/h_z)+1;
                    if (z_index<z_resolution)&&(z_index>0)
                        z_pdf(z_index) = z_pdf(z_index) + 1.0;
                    end
                end
            end
        end
        
    
        figure(1)
        clf()
        z_hist = h_z/2.0:h_z:(z_resolution*h_z);
        plot(z_hist,z_pdf,'black'); hold on;
        xlabel("z");
        axis([0 maxz 0 max(z_pdf)+10]);
        ch=sprintf("Time %f",it*frame_frequency*Delta_t); 
        title(ch);
    end
    

end

function [] = Graphical_Simulation()

%     % plotting wall
%     d_phi_wall = 0.05;
%     phi_wall = 0:d_phi_wall:2*pi;
%     
%     for i = 1:length(phi_wall)
%         x_wall(i) = a_ellipse*cos(phi_wall(i));
%         y_wall(i) = b_ellipse*sin(phi_wall(i));
%     end
%     
%     minz=min(z_rod,[],"all"); % minimum over all elements of z_rod
%     maxz=max(z_rod,[],"all"); % maximum over all elements of z_rod
%     z_wall = [minz maxz];
%     xz_wall = [-a_ellipse a_ellipse];

    % for i=1:frame_indicator
    %     f1=figure(1);
    % 
    %     scrsz = get(groot,'ScreenSize');
    %     maxscrsz=min(scrsz(3),scrsz(4));
    %     set(f1,'Position',[scrsz(3)/3 0 maxscrsz maxscrsz],'Color','w')
    % 
    %     clf;
    %     subplot(1,3,1);
    %     plot(x_wall,y_wall); hold on;
    %     plot(x_rod(:,i),y_rod(:,i),'red.');hold on;
    %     grid on
    %     daspect([1 1 1]);
    %     
    %     subplot(1,3,2); 
    %     plot([-a_ellipse -a_ellipse],z_wall); hold on;
    %     plot([a_ellipse a_ellipse],z_wall); hold on;
    %     plot(x_rod(:,i),z_rod(:,i),'red.');hold on;
    %     grid on
    %     daspect([1 1 1]);
    %     
    %     subplot(1,3,3); 
    %     plot([-b_ellipse -b_ellipse],z_wall); hold on;
    %     plot([b_ellipse b_ellipse],z_wall); hold on;
    %     plot(y_rod(:,i),z_rod(:,i),'red.');hold on;
    %     grid on
    %     daspect([1 1 1]);
        
    %     ch=sprintf("%d.png",i);
    %     saveas(gcf,ch);

end

function a = x_deriv(phi, theta)

    a = cos(phi)*sin(theta);

end

function b = y_deriv(phi, theta)
    b = sin(phi)*sin(theta);
end

function c = z_deriv(theta, x, y)
    global flow_rate;
    global a_ellipse;
    global b_ellipse;
    c = cos(theta) + flow_rate*(x^2/(a_ellipse^2) + y^2/(b_ellipse^2) - 1);
end

function q = theta_deriv(theta, phi, x, y)
    global flow_rate;
    global a_ellipse;
    global b_ellipse;
    q = -(sin(theta))^2 * (cos(phi)*(2*x/a_ellipse^2) + sin(phi)*(2*y/b_ellipse^2));
    q = flow_rate * q;
end




% Velocity functions

function velocity = calc_velocity(total_distance,total_time)
    % caluclate the velocity of one rod, given the total_distance the rod
    % was able to swim and the duration of time it took the rod to swim
    % that distance.
    velocity = total_distance/total_time;
end

function avg_velocity = calc_avg_velocity(number_of_rods,velocity_arr)
    % calculates the average velocity of a simulation, passing in an array
    % of already calcuated velocities of the rods, function will get a sum
    % of the velocity values and divide it by the given number of rods.
    
    v_sum = 0;
    for i=1:length(velocity_arr)
        v_sum = v_sum + velocity_arr(i);
    end

    avg_velocity = v_sum/number_of_rods;
end


% Distribution function
function [] = xyz_distribution(x_rod,y_rod,z_rod,flowRate)

    figure('Name',num2str(flowRate));
    
    subplot(1, 3, 1);
    histogram(x_rod(:), 'Normalization','count');
    title('X Position Distribution');
    xlabel('x');
    ylabel('counted');
    
    subplot(1, 3, 2);
    histogram(y_rod(:), 'Normalization', 'count');
    title('Y Position Distribution');
    xlabel('y');
    ylabel('counted');
    
    subplot(1, 3, 3);
    histogram(z_rod(:), 'Normalization', 'count');
    title('Z Position Distribution');
    xlabel('z');
    ylabel('counted');


end

