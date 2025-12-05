function vmi_sim_visualize_particle_paths(input_crate)
desolve_struct(input_crate);

colors = lines;




v = sqrt(vx_components.^2 + vy_components.^2 +vz_components.^2);
R = sqrt(x_components.^2 + y_components.^2 +z_components.^2);

%



b_v = v(1:num_molecules,end)>18 | v(num_molecules+1:end, end)>18;
b_v = true(size(b_v));

%b_v = v(1:num_molecules,end)<2 | v(num_molecules+1:end, end)<2;


b_v = R(1:num_molecules,end) < droplet_radii(1:num_molecules)   | R(num_molecules+1:end,end) < droplet_radii(num_molecules+1:end);
b_v = R(1:num_molecules,end) < droplet_radii(1:num_molecules)   & R(num_molecules+1:end,end) < droplet_radii(num_molecules+1:end);

b_v = true(size(b_v));

indices = 1:num_molecules;

try
id_select = randsample(indices(b_v), 10);
catch ME
disp(ME);

end

path_colors = [];

figure
for id = id_select
line = plot(time, v(id,:));
hold on
plot(time, v(id+num_molecules,:), 'color', line.Color, 'LineStyle','--');

path_colors = [path_colors; line.Color];

end

xlabel('t / ps');
ylabel(['v / ', char(0197), '/ps']);



counter = 1;

for id = id_select
figure


r0 = [x_components(id,1); y_components(id,1); z_components(id,1)];
r =  [x_components(id,:); y_components(id,:); z_components(id,:)];


%r = r- r0;

line = plot3(r(1,:), r(2,:), r(3,:),'color',   [path_colors(counter,:),1]);

hold on
scatter3(r(1,1), r(2,1), r(3,1), 4,'MarkerEdgeColor',path_colors(counter,:) ,'MarkerEdgeAlpha',1);


test = 1;

id = id + num_molecules;



r0 = [x_components(id,1); y_components(id,1); z_components(id,1)];
r =  [x_components(id,:); y_components(id,:); z_components(id,:)];



line =  plot3(r(1,:), r(2,:), r(3,:),'color',  [path_colors(counter,:),1]);

hold on
scatter3(r(1,1), r(2,1), r(3,1), 4,'MarkerEdgeColor',path_colors(counter,:),'MarkerEdgeAlpha',1);


counter = counter + 1;


xlim([-80,80]);
ylim([-80,80]);
zlim([-80,80]);
xlabel(['x / ', char(0197)])
ylabel(['y / ', char(0197)])
zlabel(['z / ', char(0197)])


pbaspect([1,1,1]);

        [X,Y,Z] = sphere;
        X = X*mean(droplet_radii(id));
        Y = Y*mean(droplet_radii(id));
        Z = Z*mean(droplet_radii(id));

        s=  surface(X,Y,Z, 'FaceAlpha',0.1, 'FaceLighting','gouraud', 'EdgeColor',[0.1,0.3,0.9], 'EdgeAlpha',0.2, 'EdgeLighting','gouraud');
        s.CData = Z*0;
        colormap('winter')

end

% figure
% 
% v = sqrt(vx_components.^2 + vy_components.^2 +vz_components.^2);
% 
% plot(time, v);

end