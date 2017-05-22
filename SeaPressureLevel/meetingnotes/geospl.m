%This will generate images and video in subfolder called 'images'!
%
function sign = geospl(nc_name) % Designed for data downloaded from NCEP.
	lat = ncread(nc_name, 'lat'); % Read latitdue, longitude and pressure matrix from NC file.
	lon = ncread(nc_name, 'lon');
	SeaPressure = ncread(nc_name, 'slp');
    SeaPressure = permute(SeaPressure, [2 1 3]);

	latlim = double([min(lat) max(lat)]); % Range of latitude and longitude.
	lonlim = double([min(lon) max(lon)]);
	Dimension = size(SeaPressure); % Dimension of dataset latitude*longitude*pages.
	num_obs = Dimension(3); % Number of observations.
    Latitude = repmat(lat, 1, Dimension(2));
    Longitude = repmat(transpose(lon), Dimension(1), 1);
    mat_name = 'SeaLevelPressure_Info.mat';
    save(mat_name, 'Dimension', 'Latitude', 'Longitude', 'SeaPressure'); % Save data as mat file.
    
    
    
	%imageNames = cell(num_obs,1);
	%imDirName = strcat(pwd,'/images/');
	%if exist(imDirName) ~= 7,  
	%	mkdir( imDirName );  
	%end

	crange = prctile(SeaPressure(:) , [1 99]);

	outputVideo = VideoWriter(strcat(pwd,'/SeaPressure.avi'));
	outputVideo.FrameRate = 5; % Set FPS manually.
	open(outputVideo);

	for i = 1:num_obs,

		disp( strcat( 'generating frame # ', num2str(i) ) );

		exp = SeaPressure(:,:,i); % Need transpose for lat*lon.
		% FIXME Check the following to see if it does what I assumed it does.
		R = georasterref('RasterSize', size(exp), 'Latlim', latlim, 'Lonlim', lonlim); % Set up raster reference from dataset.
		figure('Visible','off','Color','k'); 
		worldmap(latlim, lonlim); % Set up world map.
		load geoid
		geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');
		geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15]); % Add land areas and coastlines.
		geoshow(exp, R, 'DisplayType', 'texturemap');
		mesh(peaks);
		caxis(crange);
		colormap(parula);

		fig = gcf;
		fig.Color = 'white';
		F = getframe(gcf);
		writeVideo(outputVideo, F);
	end

	close(outputVideo);
	sign = 'Done';

end
