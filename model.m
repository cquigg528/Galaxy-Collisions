% Core masses
mc = [4 3];

% Stars per core
Nstars = 3500;

% Core positions
core1 = [70 -90 -40];
core2 = [40 35 20];
gr_0(1, :) = core1;
gr_0(2, :) = core2;


% Core initial velocities (vx, vy, vz) and star rotation control:
% Counter-clockwise if last elements are ~=0, else clockwise
vcore1 = [0.05, 0, 0, 1];
vcore2 = [-0.05, -2, 5, 1];
gv_0(1, :) = vcore1;
gv_0(2, :) = vcore2;

% Run code
[t, r] = galaxy2(1600.0, 7, mc, Nstars, gr_0, gv_0);


% Run simulation
% === Graphics ===
plotenable = 1;
pausesecs = 0.0025;

starsize = 1;
starcolor = 'c';
starmarker = 'o'; 
coresize = 12;
corecolor = 'y';

starcolor2 = 'm';
corecolor2 = 'b';

avienable = 1;
avifilename = 'galaxy.avi';
aviframerate = 25;


if avienable
    aviobj = VideoWriter(avifilename);
    aviobj.Quality = 100;
    open(aviobj);
end

for i=1 : length(r(1,1,:))
    if plotenable
        clf;
        hold on;
        axis square;
        box on;
        set(gca,'Color','k');
%         view(-29, 5);    % change camera (azimuth, elevation)
%         view(0, 80);
        view(-81, 9);
        xlabel("x");
        ylabel('y');
        zlabel('z');
        grid on;
        xlim([-200 200]);
        ylim([-200 200]);
        zlim([-150 150]);
        plot3(squeeze(r(2:Nstars+1, 1, i)), squeeze(r(2:Nstars + 1, 2, i)), ...
            squeeze(r(2:Nstars+1, 3, i)), "Marker", starmarker, ... 
            "Markersize", starsize, "MarkerEdgeColor", starcolor2, ... 
            "MarkerFaceColor", starcolor2, "Linestyle", 'none');
        plot3(squeeze(r(1, 1, i)), squeeze(r(1, 2, i)), ...
            squeeze(r(1, 3, i)), "Marker", starmarker, "Markersize", ...  
            coresize, "MarkerEdgeColor", corecolor2, "MarkerFaceColor", ...
            corecolor2);
        plot3(squeeze(r(Nstars + 3:end, 1, i)), squeeze(r(Nstars + 3:end, 2, i)), ...
            squeeze(r(Nstars + 3:end, 3, i)), "Marker", starmarker, ... 
            "Markersize", starsize, "MarkerEdgeColor", starcolor, ... 
            "MarkerFaceColor", starcolor, "Linestyle", 'none');
        plot3(squeeze(r(Nstars + 2, 1, i)), squeeze(r(Nstars + 2, 2, i)), ...
            squeeze(r(Nstars + 2,3, i)), "Marker", starmarker, "Markersize", ...  
            coresize, "MarkerEdgeColor", corecolor, "MarkerFaceColor", ...
            corecolor);
        drawnow;
        if avienable
            if i == 1
                framecount = 3 * aviframerate;
            else
                framecount = 1;
            end
            for iframe = 1 : framecount
                writeVideo(aviobj, getframe(gcf));
            end
        end
        
        pause(pausesecs);
    end
end

        
if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end
