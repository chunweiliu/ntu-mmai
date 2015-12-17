function [num_shots, shot_frames] = shot_detect( index, filename )
%The function will read a mpeg vedio and detect where the shot is.
%   Load the vedio as an object and load all vedio frames.
vedio_object = mmreader( filename );
vedio_frames = read( vedio_object );
num_shots = 0;
shot_frames = [];

if ( index == 1 ) 
%Process the pixel compression throught each frame.
%   For all frames, find the difference between the neighbor frame
%       video_frames( row, col, color_channel, frame_index )
%       frame( row, col, color_channel )

    frame = vedio_frames(:, :, 1, 1);
    [height width] = size(frame);
    xy = height * width * 256 ;
    threshold_h = 0.06;
    
    num_frames = get( vedio_object, 'NumberOfFrames' );
    for i = 2: num_frames
        prevf = vedio_frames(:, :, :, i-1);
        frame = vedio_frames(:, :, :, i);
        % RGB 3 chanel
        pixel_dif = 0;
        for j = 1: 3
            pixel_dif = pixel_dif + abs( sum( sum( frame(:, :, j) - prevf(:, :, j) ) ) );
        end
        
        % Simple pixel difference algorithm
        s(i) = pixel_dif / xy;
        shot = s(i) > threshold_h; 
        
        if (shot)
            num_shots = num_shots + 1;
            shot_frames(num_shots) = i;
        end
    end
    % drow the difference couve
    x = linspace(1, num_frames-1, num_frames-1);
    plot(x, s(x));
    
elseif ( index == 2 )
%The Color historgram mehtod
%   use RGB historgram
%   hist(frame, nbins): draw the histrogram
%   histc(frame, nbins): count the histrogram

    frame = vedio_frames(:, :, 1, 1);
    xy = 256 * 256;
    threshold_h = 0.2;
    bins = [64, 128, 192, 256];
    
    num_frames = get( vedio_object, 'NumberOfFrames' );
    for i = 2: num_frames
        prevf = vedio_frames(:, :, :, i-1);
        frame = vedio_frames(:, :, :, i);
        hist_dif = 0;
        for j = 1: 3
            
            % L1 distance
            hist_dif = hist_dif + abs( sum( sum( histc( frame(:, :, j), bins ) - histc( prevf(:, :, j), bins ) ) ) );
        end
        
        % Simple histogram difference algorithm
        s(i) = hist_dif / xy;
        shot = s(i) > threshold_h;
        
        if(shot)
            num_shots = num_shots + 1;
            shot_frames(num_shots) = i;
        end
    end
    x = linspace(1, num_frames-1, num_frames-1);
    plot(x, s(x));

elseif( index == 3 )
%REGION HISTOGRAMS
    frame = vedio_frames(:, :, 1, 1);
    [height width] = size(frame);
    xy = height * width;
    threshold_h = 0.4;
    bins = [64, 128, 192, 256];
    bx = [1, floor(width/4), floor(width/2), floor(width*3/4), width];
    by = [1, floor(height/4), floor(height/2), floor(height*3/4), height];
    
    num_frames = get( vedio_object, 'NumberOfFrames' );
    for i = 2: num_frames
        prevf = vedio_frames(:, :, :, i-1);
        frame = vedio_frames(:, :, :, i);
        hist_dif = 0;
        for j = 1: 3
            
            % consider 16 region
            for  m = 1: 4
                for n = 1: 4
                    
                    % L1 distance
                    
                    %hist_dif = hist_dif + abs( sum( sum( histc( frame(1:240, 1:320, j), bins ) - histc( prevf(1: 240, 1:320, j), bins ) ) ) );
                    hist_dif = hist_dif + abs( sum( sum( ...
                        histc( frame(by(m):by(m+1), bx(n):bx(n+1), j), bins ) - histc( prevf(by(m): by(m+1), bx(n):bx(n+1), j), bins ) ) ) );      
                end
            end
        end
        
        % Simple histogram difference algorithm
        s(i) = hist_dif / xy;
        shot = s(i) > threshold_h;
        
        if(shot)
            num_shots = num_shots + 1;
            shot_frames(num_shots) = i;
        end
    end

    x = linspace(1, num_frames-1, num_frames-1);
    plot(x, s(x));
end

%Print the result
fprintf('The file "%s" has "%d" shots.\n', filename, num_shots);
fprintf('The shots are detectived below\n');
for i = 1: num_shots,
fprintf('%d ', shot_frames(i));
end
fprintf('\n');
