
baseWD = '/rigel/tzsts/projects/mhk2154/';

load( strcat(baseWD, 'data/chbt.mat') );

data         = chbt;
[d1, d2, d3] = size(data);
samplingRate = 10.40;
sTime = (1:d3)/samplingRate;

%NOTE: removing the amplitude A
maxFrame    = reshape(prctile(reshape(data, d1*d2, d3).', 95).', d1, d2);
maxFrameInv = boolean(maxFrame>eps).*(maxFrame.^-1);
data        = data .* repmat(maxFrameInv, 1, 1, d3);


%NOTE: setting up the output director
dirname      = strcat( baseWD, num2str(d1), '-chbt_wp_newmask');
bfilename    = strcat( dirname, '/pasf' );
mkdir( dirname );


fprintf('The input data size is [%d,%d,%d]. \n', d1, d2, d3);
fprintf('The output file names start with %s. \n', bfilename);
fprintf('The number of computational threads is: %d. \n', maxNumCompThreads);


%NOTE: loading the spatial mask to mask out artifacts
%spatialMask = double(boolean(imread('mask.bmp')));
load('mask.mat');


pasfOpts = struct('nTopEVs', 3, ...
	'saveData', 2, 'boolUseSavedData', 1, ...
	'SampleTimes', sTime, ...
	'cmethod', 'weightedPhase', ... 
	'outlierThreshold', 0.2, ...
	'boolDemeanInput', 2, ...
	'spans', 15, ...
	'errorRate', 0.1, 'bfilename', bfilename);


tic;
[Components, Clusters, ClusterInfo, SDFInfo] = pasf(data, 5, pasfOpts, spatialMask);
toc;

%load( strcat(bfilename, '_output.mat') );

plot_PASF_outputs(Components, Clusters, ClusterInfo, SDFInfo, false, bfilename);

[d1, d2, d3, d4] = size(Components);

parfor i=1:d4,
       cName = strcat(bfilename, '_cmp', num2str(i));
       if i == d4-1,
	       cName = strcat(bfilename, '_reminder');
       elseif i == d4,
	       cName = strcat(bfilename, '_signal');
       end

       saveVideo(Components(:, :, :, i), 1:d3,  cName);
end

disp('done.');

