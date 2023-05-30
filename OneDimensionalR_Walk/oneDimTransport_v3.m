% first clean figures and memory from all previous calculations
clc;
clf;
clear;


R = 10000;
Nt = 20000;

Delta_t = 0.01; 
Ts = Nt*Delta_t / R;
v_prop = 1.0; 
t_downstream = 5;
frame_frequency = 100;
length_of_rod = 0.005;


for r=1:R 
    z(r) = 0; 
    direction_of_swimming(r) = 0; 
    time_remaining(r) = 0; 
    time_start(r) = r * Ts;
end

index_times = 0;
for r=1:R 
    for nt=1:Nt 
        z_old(r) = z(r); 
        
        if (mod(nt-1,frame_frequency)==0)
            snap_num = int16((nt-1)/frame_frequency)+1;
            z_rod(r,snap_num) = z_old(r);
        end
        % while current_time is greater than the time the rod started
        % swimming.

        if (nt*Delta_t >= time_start(r))
            % First jump - upstream swimming
            if (direction_of_swimming(r) == 0)
                % Metropolis algorithm 
                accepted = 0;
                while (accepted == 0) 
                    t_upstream = 1000.0*rand();
                    checker = rand();
                    if checker < power_law(t_upstream)
                        accepted = 1;
                        index_times = index_times + 1;
                        times_upstream(index_times) = t_upstream;
                    end
                end
                time_remaining(r) = t_upstream;
                direction_of_swimming(r) = 1;
            end
            
            % If we are currently swimming upstream
            if (direction_of_swimming(r) == 1)
                time_remaining(r) = time_remaining(r) - Delta_t;
                z(r)=z_old(r)+v_prop*Delta_t;
                if time_remaining(r) < 0
                    % swimming direction changed to down stream
                    direction_of_swimming(r) = -1; 
                    time_remaining(r) = t_downstream;
                end
            end
            
            % If we are currently swimming downstream
            if (direction_of_swimming(r) == -1)
                time_remaining(r) = time_remaining(r) - Delta_t;
                z(r) = z_old(r) - v_prop*Delta_t;

                if time_remaining(r) < 0
                    % Metropolis algorithm 
                    accepted = 0;
                    while (accepted == 0) 
                        t_upstream = 1000.0*rand();
                        checker = rand();
                        if checker < power_law(t_upstream)
                            accepted = 1;
                        end
                    end
                    time_remaining(r) = t_upstream;
                    direction_of_swimming(r) = 1;    
                end
            end
            
            
        end
        
    end
    
end



maxz=max(z_rod,[],"all");

h_z = 5;
z_resolution = int16(maxz/h_z)+1;
z_resolution = double(z_resolution);

% for k=1:z_resolution
%     z_hist2(k) = (k-0.5)*h_z;  
% end


for it=2:Nt/frame_frequency
    
    for k=1:z_resolution
        z_pdf(k) = 0;  
    end

    for r=1:R
        for z_body = (z_rod(r,it)-length_of_rod):h_z:(z_rod(r,it)+length_of_rod) 
            if (z_rod(r,it) > 0)
                z_index = int16(z_body/h_z)+1;

                if (z_index< z_resolution)&&(z_index>0)
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
%     ch=sprintf("%d.png",it);
%     saveas(gcf,ch);
   
    
end


function z = psi(t)
    z = exp(-t);

end

function z = power_law(t)
    g = 1.2; % gamma
    u = (1 + t).^(g + 1);
    z = 1/u;
end