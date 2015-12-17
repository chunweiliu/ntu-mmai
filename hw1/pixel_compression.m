function [num_shots, shot_frames] = shot_detection( index, filename )
%The function will read a mpeg vedio and detect where the shot is.
%   Load the vedio as an object and load all vedio frames.
vedio_object = mmreader( filename );
vedio_frames = read( vedio_object );
num_shots = 0;

if ( index == 1 ) 
%Process the pixel compression throught each frame.
%   video_frames( x, y, color_channel, frame_index )
%   frame( x, y, color_channel )
    num_frames = get( vedio_object, 'NumberOfFrames' );
    frame = vedio_frames(:, :, :, 1);
    [height width] = size(frame);
    xy = height * width;
    threshold1 = 20;
    shot = 0;
    for i = 2: num_frames
        prevf = vedio_frames(:, :, :, i-1);
        frame = vedio_frames(:, :, :, i);
        % RGB 3 chanel
        pixel_dif = 0;
        for j = 1: 3
            pixel_dif = pixel_dif + abs(sum( sum( frame(:, :, j) - prevf(:, :, j) ) ));
        end
        
        s(i) = pixel_dif / xy;
        shot = s(i) > threshold1; 
        
        if (shot)
            num_shots = num_shots + 1;
            %shot_frames(num_shots) = frame;
        end
    end
    % drow the difference couve
    x = linspace(1, num_frames-1, num_frames-1);
    plot(x, s(x));
end