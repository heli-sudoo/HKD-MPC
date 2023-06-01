addpath("PreProcessedData");
addpath("PostProcessedData/");

% load("Inplace-trot.mat");

%%
tau_sz = size(body_states, 1);
status_durations = zeros(tau_sz, 4);
GRFs = zeros(tau_sz, 12);
torques = zeros(tau_sz, 12);

% Calculate contact status durations for each leg 
dt = t(2) - t(1);
for leg = 1:4
    status_durations(:, leg) = Induce_status_duration_per_leg(contacts(:,leg), dt);
end

% Calculate the GRF reference
mass = 9; g = 10;
for k = 1:tau_sz
    F = [0, 0, mass * g / sum(contacts(k, :))];
    for leg = 1: 4
        if contacts(k, leg)
            GRFs(k, 3*(leg-1)+1:3*leg) = F;
        end
    end    
end

%% Write to file
% write contact information to csv file
fname = "PostProcessedData/"+"quad_reference.csv";
fid = fopen(fname, 'w');
fprintf(fid, 'dt\n');
fprintf(fid, '%4.3f\n', dt);
for i = 1:tau_sz
    fprintf(fid, 'body_state \n');
    fprintf_array(fid, body_states(i, :), '%6.3f ');       
    
    fprintf(fid, 'qJ\n');
    fprintf_array(fid, qJs(i, :), '%6.3f ');

    fprintf(fid, 'foot_placements\n');
    fprintf_array(fid, foot_placements(i, :), '%6.3f ');

    fprintf(fid, 'grf\n');
    fprintf_array(fid, GRFs(i, :), '%6.3f ');

    fprintf(fid, 'torque\n');
    fprintf_array(fid, torques(i, :), '%6.3f ');   

    fprintf(fid, 'contact\n');
    fprintf_array(fid, contacts(i, :), '%d ');

    fprintf(fid, 'status_dur\n');
    fprintf_array(fid, status_durations(i, :), '%6.3f ');
end
fclose(fid);




%% Help functions
function status_durs = Induce_status_duration_per_leg(contacts_leg, dt)
tau_sz = length(contacts_leg);
status_durs = zeros(tau_sz, 1);
status_dur = 0;
status_start = 1;
contact_prev = contacts_leg(1);

for k = 2:tau_sz
    status_dur = status_dur + dt;
    contact_cur = contacts_leg(k);    
    if (contact_cur ~= contact_prev) 
        status_end = k - 1;
        status_durs(status_start:status_end) = status_dur;
        status_start = k;
        status_dur = 0;
        contact_prev = contact_cur;
    end
    
    if k == tau_sz
        status_end = k;
        status_durs(status_start:status_end) = status_dur;
    end        
end
end

function fprintf_array(fid, a, format)
for j = 1:length(a)
    fprintf(fid, format, a(j));
end
fprintf(fid, '\n');
end
