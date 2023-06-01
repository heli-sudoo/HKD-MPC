function drawFloor()

% Create the ground. It reaches from -8 to +8
    v = []; f = []; c = [];
    % Number of tiles on each side of the origin
    nTilesX = 8;  
    nTilesY = 8;  
    % Size of the tiles
    sTiles = 1;  
    % Create vertices:
    for i = -sTiles*nTilesX:sTiles:sTiles*nTilesX
        for j = -sTiles*nTilesY:sTiles:sTiles*nTilesY
            v = [v;[i,j,0]];                         
        end
    end
    % Connect them and color them as checkerboard:
    for i = 1:2*nTilesY
        for j = 1:2*nTilesX
            f = [f;[i,i+1,i+2+2*nTilesY,i+1+2*nTilesY]+(j-1)*(2*nTilesY+1)];	% Connect vertices
            if mod(j+i,2)==0	% Color faces alternating
                c = [c;[1,1,1]];
            else
                c = [c;[.8,.8,.8]];
            end
        end
    end

p = patch('faces', f, 'vertices', v, 'FaceVertexCData', c, 'FaceColor', 'flat');
%set(p, 'FaceAlpha',0.7);    % Set to partly transparent
set(p, 'FaceLighting','phong'); % Set the renderer
set(p, 'FaceColor','flat');     % Set the face coloring
set(p, 'EdgeColor','none');     % Don't show edges    
set(p,'AmbientStrength',.6);
end