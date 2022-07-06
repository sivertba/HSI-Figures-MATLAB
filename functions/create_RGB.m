function [x_img] = create_RGB(cube,R,G,B,FRA, OS)
% Create RGB Image
x_img = cube(:,:,[R G B]);


% The variable FRA roughly controls the contrast.
% Increase the value of FRA to make the image brighter.
switch nargin
    case 4
        FRA = 0.02;
        OS = 0.5;
    case 5
        OS = 0.5;
end

% scale the radiances from 0 to 1
RA = (1.0*x_img(:,:,1))/max(max(x_img(:,:,1)));
GA = (1.0*x_img(:,:,2))/max(max(x_img(:,:,2)));
BA = (1.0*x_img(:,:,3))/max(max(x_img(:,:,3)));

RA = OS+(FRA/std(std(RA))*(RA - mean(mean(RA))));
GA = OS+(FRA/std(std(GA))*(GA - mean(mean(GA))));
BA = OS+(FRA/std(std(BA))*(BA - mean(mean(BA))));

RA(RA>1) = 1; RA(RA<0) = 0;
GA(GA>1) = 1; GA(GA<0) = 0;
BA(BA>1) = 1; BA(BA<0) = 0;

x_img(:,:,1) = RA;
x_img(:,:,2) = GA;
x_img(:,:,3) = BA;

end