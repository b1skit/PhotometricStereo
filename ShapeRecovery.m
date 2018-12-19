%% CMPT 412 Assignment 03
% Photometric stereo recovery of shape & relative distance
% By Adam Badke
% 301310785

%% Clear the workspace:
clearvars

%% Globals constants
global TABLE_DIM;
global RATIO_MIN;   
global RATIO_MAX;

TABLE_DIM = 512;    % Square dimension of the lookup table
RATIO_MIN = 0.001;  % Minimum E_i/E_j ratio value to map
RATIO_MAX = 185;    % Maximum E_i/E_j ratio value to map

QVR_DNSTY = 10;     % Spacing to sample when plotting quiver plots

%% README: DATA SET SELECTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program batch processes images stored in an array in sets of 3.
% For best results, synthetic and real images are processed with different
% light measurements. The program must be configured in advance to select
% which set of images you wish to process.
%
% To select the image set & corresponding light configurations you wish to
% process, comment/uncomment the relevant sections below.
%
% Uncommenting the first section will process the set of synthetic images.
%
% Uncommenting the second section will process the set of real images.
%
% Only 1 section should be uncommented at a time.
%
% The .tif images must be lit in the order of left, center, right:
% image1 = image lit from the left
% image2 = image lit from the center
% image3 = image lit from the right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% SYNTHETIC IMAGE SET PROCESSING
% % Uncomment this block to process the set of synthetic images
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Extracted measurements from calibration sphere:
% lightPositions(:,:,1) = [226,160, 0]; 
% lightPositions(:,:,2) = [320, 125, 0];
% lightPositions(:,:,3) = [415, 160, 0];
% 
% % Center position
% CENTER_X = 320;
% CENTER_Y = 240;
% CAL_RADIUS = 153; % Measured radius of the calibration sphere
% 
% % Synthetic images:
% tifFiles = {
%     'cone1.tif';
%     'cone2.tif';
%     'cone3.tif';
%     
%     'hexagon1.tif';
%     'hexagon2.tif';
%     'hexagon3.tif';
%     
%     'octogon1.tif';
%     'octogon2.tif';
%     'octogon3.tif';
% 
%     'sphere1.tif';
%     'sphere2.tif';
%     'sphere3.tif';
%     
%     'torus1.tif';
%     'torus2.tif';
%     'torus3.tif';    
% };


% REAL IMAGE SET PROCESSING
% Uncomment this block to process the set of real images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracted measurements from calibration sphere:
lightPositions(:,:,1) = [290,290, 0];
lightPositions(:,:,2) = [318, 261, 0];
lightPositions(:,:,3) = [357, 289, 0];

% Center position
CENTER_X = 318; 
CENTER_Y = 296;
CAL_RADIUS = 153; % Measured radius of the calibration sphere

% Real images:
tifFiles = {
    'cone2-lamp1.tif';
    'cone2-lamp3.tif'; % Switched order
    'cone2-lamp2.tif'; % Switched order
    
    'cone-lamp1.tif';
    'cone-lamp3.tif'; % Switched order
    'cone-lamp2.tif'; % Switched order

    'cylinder-lamp1.tif';
    'cylinder-lamp3.tif'; % Switched order
    'cylinder-lamp2.tif'; % Switched order
    
    'ellipsoid-lamp1.tif';
    'ellipsoid-lamp3.tif'; % Switched order
    'ellipsoid-lamp2.tif'; % Switched order
    
    'hex1-lamp1.tif';
    'hex1-lamp3.tif'; % Switched order
    'hex1-lamp2.tif'; % Switched order
    
    'hex2-lamp1.tif';
    'hex2-lamp3.tif'; % Switched order
    'hex2-lamp2.tif'; % Switched order    
    
    'sphere-lamp1.tif';
    'sphere-lamp3.tif'; % Switched order
    'sphere-lamp2.tif'; % Switched order
};


%% END OF IMAGE SELECTION BLOCKS

%% Compute the Z element of the light direction vector:
lightDirections = zeros(1, 3, 3); % [x,y,z]x3
for i = 1:3
    lightPositions(1,3,i) = sqrt(CAL_RADIUS^2 - (lightPositions(1,1,i) ...
        - CENTER_X)^2 - (lightPositions(1,2,i) - CENTER_Y)^2);
    
    % Calculate the direction of the light from the object (light pop
    % NOTE: The TOP of the image is 0,0
    lightDirections(:,:,i) = ( lightPositions(:,:,i) ...
        - [CENTER_X, CENTER_Y, 0] )/norm(lightPositions(:,:,i) ...
        - [CENTER_X, CENTER_Y, 0]); 
end

% Print the light positions and directions:
lightPositions
lightDirections

%% Build a lookup table:
gradientLookup = zeros(TABLE_DIM, TABLE_DIM, 2);
for row = 1: TABLE_DIM
    for col = 1: TABLE_DIM

        E1oE2 = indexToRatio(row);
        E2oE3 = indexToRatio(col);
        
        ps1 = lightDirections(1,1,1);
        ps2 = lightDirections(1,1,2);
        ps3 = lightDirections(1,1,3);
        
        qs1 = lightDirections(1,2,1);
        qs2 = lightDirections(1,2,2);
        qs3 = lightDirections(1,2,3);
        
        % Derived formulas for p, q:
        p = (qs1*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            qs2*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            E2oE3*qs1*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E2oE3*qs3*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E1oE2*E2oE3*qs2*(ps1^2 + qs1^2 + 1)^(1/2) - ...
            E1oE2*E2oE3*qs3*(ps1^2 + qs1^2 + 1)^(1/2))/...
            (ps1*qs2*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            ps2*qs1*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            E2oE3*ps1*qs3*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E2oE3*ps3*qs1*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E1oE2*E2oE3*ps2*qs3*(ps1^2 + qs1^2 + 1)^(1/2) - ...
            E1oE2*E2oE3*ps3*qs2*(ps1^2 + qs1^2 + 1)^(1/2));
        
        q = -(ps1*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            ps2*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            E2oE3*ps1*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E2oE3*ps3*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E1oE2*E2oE3*ps2*(ps1^2 + qs1^2 + 1)^(1/2) - ...
            E1oE2*E2oE3*ps3*(ps1^2 + qs1^2 + 1)^(1/2))/...
            (ps1*qs2*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            ps2*qs1*(ps3^2 + qs3^2 + 1)^(1/2) - ...
            E2oE3*ps1*qs3*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E2oE3*ps3*qs1*(ps2^2 + qs2^2 + 1)^(1/2) + ...
            E1oE2*E2oE3*ps2*qs3*(ps1^2 + qs1^2 + 1)^(1/2) - ...
            E1oE2*E2oE3*ps3*qs2*(ps1^2 + qs1^2 + 1)^(1/2));
        
        gradientLookup(row, col, 1) = -p; % Negate due to the flipped axis
        gradientLookup(row, col, 2) = q;
    end
end


%% Loop through each image, processing 1 at a time
index = 1;
while index <= size(tifFiles,1)
    % Load the 3 consecutive images:
    image1 = double(imread(char(tifFiles(index))));
    index = index + 1;
    image2 = double(imread(char(tifFiles(index))));
    index = index + 1;
    image3 = double(imread(char(tifFiles(index))));
    index = index + 1;
    
    % Convert to greyscale, [0,1] in a single matrix for easier handling:
    clear combinedImage;
    combinedImg(:,:,1) = (image1(:,:,1) + ...
        image1(:,:,2) + image1(:,:,3))/(3 * 255);
    combinedImg(:,:,2) = (image2(:,:,1) + ...
        image2(:,:,2) + image2(:,:,3))/(3 * 255);
    combinedImg(:,:,3) = (image3(:,:,1) + ...
        image3(:,:,2) + image3(:,:,3))/(3 * 255);

    
    %% Assemble a gradient and normal map:
    normalMap = zeros(size(combinedImg,1), size(combinedImg,2), 3);
    gradientMap = zeros(size(combinedImg,1), size(combinedImg,2), 2);
    for y = 1:size(combinedImg,1)
        for x = 1:size(combinedImg,2)

            E1 = combinedImg(y, x, 1);
            E2 = combinedImg(y, x, 2);
            E3 = combinedImg(y, x, 3);

            if E1 ~= 0 && E3 ~= 0 && E3 ~=0
                p = gradientLookup(ratioToIndex(E1/E2), ...
                    ratioToIndex(E2/E3), 1);
                q = gradientLookup(ratioToIndex(E1/E2), ...
                    ratioToIndex(E2/E3), 2);

                gradientMap(y, x, 1) = p;
                gradientMap(y, x, 2) = q;

                theNormal = [-p, -q, 1];
                theNormal = theNormal/norm(theNormal);

                % Shift normal map R & G to [0,1], & B to [0.5,1]
                normalMap(y, x, 1) = (theNormal(1,1) + 1) / 2;
                normalMap(y, x, 2) = (theNormal(1,2) + 1) / 2;
                normalMap(y, x, 3) = (theNormal(1,3) / 2) + 0.5;
            end
        end
    end
    
    %% Calculate the depth map:
    depthMap = zeros(size(combinedImg,1), size(combinedImg,2));
    % Integrate across the rows:
    for y = 1:size(combinedImg,1)
        if y > 1
           depthMap(y, x) = depthMap(y - 1, x) + gradientMap(y, x, 2);
        end
        for x = 2:size(combinedImg,2)
            depthMap(y, x) = depthMap(y, x - 1) + gradientMap(y, x, 1);
        end
    end

    % Integrate down the columns, and average the results:
    for x = 1:size(combinedImg,2)
    	for y = 1:size(combinedImg,1)
            if x > 1 || y > 1
                depthMap(y, x) = (depthMap(y, x) + gradientMap(y, x, 2))/2;
            end
    	end
    end

    % Normalize to [0,1]:
    depthMap = depthMap - min(min(depthMap));
    depthMap = depthMap / max(max(depthMap));

    %% Display the results:
    % Display the quiver plot:
    [xGrid,yGrid] = meshgrid(1:size(gradientMap,2), 1:size(gradientMap,1));
    figure;
    quiver(xGrid(1:QVR_DNSTY:end,1:QVR_DNSTY:end), ...
        yGrid(1:QVR_DNSTY:end,1:QVR_DNSTY:end), ...
        gradientMap(1:QVR_DNSTY:end,1:QVR_DNSTY:end,1), ...
        gradientMap(1:QVR_DNSTY:end,1:QVR_DNSTY:end,2), 2);
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    title('Quiver plot of sampled gradient (2x scale)');

    % Display a comparison of the source, resulting normal map & depth map
    figure;
	title('Source image, normal map, & depth map results');
    montage({imread(char(tifFiles(index - 2))), normalMap, depthMap},...
        'Size', [3, 1]);

end

%% Functions:

% Convert a E_i/E_j ratio to a lookup table query index
% ratioValue: Either E1/E2, or E2/E3
% Returns a lookup table query index in [1,TABLE_DIM]
function index = ratioToIndex(ratioValue)
    global TABLE_DIM;
    global RATIO_MIN;
    global RATIO_MAX;
    
    if ratioValue >= RATIO_MIN && ratioValue <= RATIO_MAX
        index = log(ratioValue);
        index = index - log(RATIO_MIN);
        index = index / (log(RATIO_MAX) - log(RATIO_MIN));

        index = index * (TABLE_DIM - 1);
        index = index + 1;

        index = round(index);
    else
        % Handle out of bounds requests:
        if ratioValue < RATIO_MIN 
            index = 1;
        elseif ratioValue > RATIO_MAX
            index = TABLE_DIM;
        end
    end
end


% Convert a lookup table query index to a E_i/E_j ratio
% index: [1, TABLE_DIM]
% Returns an Ei/Ej ratio from a lookup table query index
function ratio = indexToRatio(index)
    global TABLE_DIM;
    global RATIO_MIN;
    global RATIO_MAX;
    
    if index >= 1 && index <= TABLE_DIM
        
        index = index - 1;
        index = index / (TABLE_DIM - 1);
        
        index = index * (log(RATIO_MAX) - log(RATIO_MIN));
        index = index + log(RATIO_MIN);
        index = exp(index);
        
        ratio = index;
        
    else
        % Handle out of bounds requests:
        if index < 1
           ratio = RATIO_MIN; 
        end
        if index > TABLE_DIM
            ratio = RATIO_MAX;
        end
    end
end


