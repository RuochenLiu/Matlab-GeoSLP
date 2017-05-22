function [Components, Clusters, ClusterInfo, SDFInfo] = pasf(Z, nSpatialComponents, options, spatialMask)
% PASF Phased-Aligned Spectral Filtering for Decomposing Spatiotemporal Dynamics.  
% PASF seeks to extract dynamics from spatiotemporal signal via data
% assimilation and modeling.  PASF assumes that the observed spatiotemporal
% data represent superimposed lower-rank smooth oscillations and movements from
% a generative dynamic system, mixed with higher-rank random noises. Separating
% the signals from noises is essential to visualize, model and understand these
% lower-rank dynamic systems. PASF uses a filtering framework for identifying
% lower-rank dynamics and its components that are embedded in a high
% dimensional spatiotemporal system. The approach is based on the structural
% decomposition and phase-aligned construction in the frequency domain.  
% 
% INPUT 'Z' is the spatiotemporal signal of dimension d1 x d2 x T where the
% 		third dimension represents the time samples.
% 	'nSpatialComponents' is the number of dynamic spatial components. PASF
% 		decomposes Z into 'nSpatialComponents' spatial components where
% 		each represents a separate dynamics.
% 	'options' is a struct with optional input parameters. Parameters are:
%		'boolDemeanInput' is a logical variable that indicates if the input
%			should be demeaned before running the PASF. The default is 1.
%		'nTopEVs' is the number of top eigenvectors PASF picks at each
%			frequency from the estimated spectral density function. The
%			default value is nSpatialComponents.
%		'EigenThresholdMethod' is the method that is used to compute
%			the hard threshold for dynamic eigenvalues to select dynamic
%			eigenvectors corresponding to the signals. The default is ErrRate which 
%			uses errorRate variable. The options are [ErrRate|MaxGap]
%		'errorRate' is the maximum error rate. The function computes
%			the eigenvalue threshold based on this rate. It should be
%			determine based on SNR.
%		'taperLen' is the taper length that is used for estimating the spectral
%			density function. The default value is zero.
%		'kmethod' is the type of smoothing kernel that is used to estimate
%			the spectral density function. The options are [box|gaussian|triangle].
%			The default is 'triangle', a.k.a Daniell.
%		'spans' is the bandwidth of the smoothing kernel that is used
%			to estimate the spectral density function. The default value is
%			floor(log(size(Z,3))).
%		'detrend' is a variable that indicates if the trend of the signal
%			should be removed for estimating the spectral density
%			function. 0 does not detrend. 1 subtracts the mean. 2
%			subtracts the linear trend. The default is 2.
%		'kmeansReplicates' is number of kmeans replications. The default is 128.
%		'outlierThreshold' is number between 0 and 2 that indicates
%			what distnace from cluster centroid in kmeans should be treated
%			as outliers and be removed from the cluster. The default is 2
%			and that means accept all. Note that the similarity is
%			1-correlation and all the distances are between 0 and 2.
%		'saveData' is a variable that indicates if the
%			PASF's intermediate and output results should be saved on the hard-drive.
%			The default is 1.
%		'boolUseSavedData" is a logical variable that indicates if
%			saved intermediate results should be used. The looks for a file
%			with name (bfilename)_intermediate_data.met.
%		'bfilename' is the base file name that function uses to save data if
%			indicated. The default is 'pasf_out'.
%		'boolParfor' is a logical variable indicates that if the parfor
%			should be used to accelerate the computations or not. The
%			default is true.
%		'boolQuietly' is a logical variable that indicates if the progress 
%			method should not be reported. The default is 0. 
% 		
% OUTPUT 'Components' is the input spatiotemporal signal that is decomposed
% 	into nSpatialComponents+2 components. The last two components are noise and
% 	the input signal after demeaning if specified, respectively. Components is an
% 	array of size d1 x d2 x T x (nSpatialComponents+2). 
%
% 	'Clusters' is an array of size nTopEVs x d3 where nTopEVs is the number
% 	of selected eigenvectors and d3 is the number of time samples.
% 	Clusters(k,f) shows the cluster label of the kth eigenvector of the
% 	spectral density function at frequency f. If the cluster number is zero
% 	that means the corresponding eigenvector didn't get clustered.
%
% 	'ClusterInfo' is a cell type object contains different information about
% 	the clustering performance.
%
% 	'SDFInfo' is a cell type object contains eigenvectors, eigevalues, phases
% 	of the spectral density function at all frequency. 
%
% Details	
% 	For more details about the PASF refer to https://arxiv.org/abs/1604.04899
% 
% Examples:
% 
% simOpts = struct('noiseVar', 1, 'noiseCorrCoeff', 0);
% data    = rotating_energy_sources(500, simOpts);
% Z       = pasf(data, 2);
% 
% %% Now plot the output 
% set(gca, 'XTickLabel','')
% 
% [d1, d2, d3, d4] = size(Z);
% minVec = prctile(reshape(Z,d1*d2*d3,d4), 1);
% maxVec = prctile(reshape(Z,d1*d2*d3,d4), 99);
% 
% for i = 1:d3,
% 	for j = 1:d4,
% 		subplot(2, 2, j);
% 		imagesc( Z(: ,: ,i ,j), [minVec(j), maxVec(j)] );
% 		colorbar;
% 
% 		if j == d4,
% 			title('Input Signal');
% 		elseif j == d4-1,
% 			title('Error');
% 		else
% 			title( strcat('Signal #', num2str(j)) );
% 		end
% 	end
% 
% 	axesHandles = findobj(get(figureHandle, 'Children'), 'flat', 'Type', 'axes');
% 	axis(axesHandles,'square')
% 	pause(0.02);
% end % 
%
%% AUTHOR    : Mohammad Khabbazian
%% $DATE     : 24-Apr-2017 $
%% $Revision : 1.00 $
%% DEVELOPED : R2016a
%% FILENAME  : pasf.m
%

	if nargin == 1
		nSpatialComponents = 1;
	end

	%%TODO separate pasf options from mvspec options
	%% setting up the default values
	pasfOpts = struct('nTopEVs', nSpatialComponents, ...
		'EigenThresholdMethod', 'ErrRate', ...
		'Normalize', 0, ...
		'SampleTimes', 1:size(Z, 3), ...
		'errorRate', 0.01, ...
		'cmethod', 'phase', ...
		'taperLen', 0, 'kmethod', 'triangle', 'spans', floor(log(size(Z,3))), ...
		'boolDemeanInput', 1, 'kmeansReplicates', 128, 'outlierThreshold', 2, ...
		'detrend', 2, 'boolZeroPad', 0, 'saveData', 1, 'boolUseSavedData', 0, ...
		'boolParfor', 1, ...
		'boolQuietly', 0, 'bfilename', 'pasf_out');


	if nargin > 2,
		fields = fieldnames(options);
		for f=1:length(fields),
			pasfOpts.(fields{f}) = options.(fields{f});
		end
	end

	if ~pasfOpts.boolQuietly,
		disp(pasfOpts);
	end


	[d1, d2, d3] = size(Z);

	%% masks
	if nargin < 4,
		% zero (false) means mask
		spatialMask = ones(d1, d2);
	end


	%%TODO just for test add the to the function input arguments
	%frequencies = pasfOpts.SampleTimes;
	%freqMask = ones(pasfOpts.nTopEVs, length(frequencies));
	%freqMask(:, frequencies > 10 ) = 0;

	if size(spatialMask, 1) ~= size(Z, 1) || size(spatialMask, 2) ~= size(Z, 2),
		error('spatialMask size is different from the input spatiotemproal signal.');
	end	
	%Z = Z .* repmat( spatialMask, [1, 1, size(Z, 3)] );

	Z = permute(reshape(Z, d1*d2, d3), [2, 1]); %NOTE: Now each column is a time series.

	%%TODO: Add noSpatialMask to avoid the extra computations.
	Z = Z .* repmat(reshape(spatialMask, [1, d1*d2]), [size(Z, 1), 1]);

	if pasfOpts.boolDemeanInput == 1,
		Z = detrend(Z, 'constant');
		if ~pasfOpts.boolQuietly, disp('demeaned the input data.');  end
	elseif pasfOpts.boolDemeanInput == 2,
		Z = detrend(Z);
		if ~pasfOpts.boolQuietly, disp('detrended the input data.'); end
	end


	if ~pasfOpts.boolQuietly,
		disp('done with the initial steps.');
	end


	%NOTE: Estimate the SDF and compute the top eigenvectors/values of the spectral density matrices.
	nTimes   = size(Z, 1);
	nSpatial = size(Z, 2);
	nTopEVs  = pasfOpts.nTopEVs;  

	if ~pasfOpts.boolUseSavedData,

		[PVec_FT_Cxx, PPhase_FT_Cxx, Eigenvalues, NF, Traces] = mvspec(Z, pasfOpts);

		if ~pasfOpts.boolQuietly,
			disp('done with estimating the spectral density matrices and their eigen decomposition.');
		end

		if pasfOpts.saveData > 0,
			filename = strcat(pasfOpts.bfilename, '_intermediate_data.mat');
			save(filename, 'PVec_FT_Cxx', 'PPhase_FT_Cxx', 'Eigenvalues', 'NF', 'Traces', '-v7.3');
			disp( strcat('done with saving the intermediate data in "', filename, '".') );
		end
	end

	if pasfOpts.boolUseSavedData,
		filename = strcat(pasfOpts.bfilename, '_intermediate_data.mat');
		load(filename, 'PVec_FT_Cxx', 'PPhase_FT_Cxx', 'Eigenvalues', 'NF', 'Traces');
		disp( strcat('done with loading the data from "', filename, '".') );
	end



	%NOTE: Decide which eigenvectors we should take into the account bases on their
	%eigenvalues.
	sortedEVs   = sort(Eigenvalues(:), 'descend');

	if strcmp( pasfOpts.EigenThresholdMethod, 'MaxGap' ), 
		diffED = sortedEVs(1:end-1) - sortedEVs(2:end);
		[maxV, maxIdx] = max( diffED );
		%%NOTE: to make sure we pick at least nSpatialComponents eigenvectors.
		maxIdx = max( maxIdx, min(2*nSpatialComponents, length(sortedEVs)-1) );
		eigenValueThreshold = sortedEVs( maxIdx + 1 );
		if ~pasfOpts.boolQuietly,
			fprintf('Max gap between sorted eigenvalues is %d. \n', maxV);
			fprintf('%d eigenvectors will be used to build the filter. \n', maxIdx);
		end
	else
		indicesMask = ( cumsum(sortedEVs) ./ sum(Traces) ) > (1 - pasfOpts.errorRate);
		if( sum(indicesMask) > 0 )
			eigenValueThreshold = max( sortedEVs(indicesMask) );
		else
			eigenValueThreshold = 0;
			warning('set the the eigenvalue threshold to zero.');
		end
	end

	if ~pasfOpts.boolQuietly,
		fprintf('Eigenvalue threshold is %d. \n', eigenValueThreshold);
	end

	scoreMask = Eigenvalues > eigenValueThreshold;



	%NOTE: unwrap the phases to have a smooth phase vector in 2D.
	if exist('unwrap2D') ~= 3,  mex unwrap2D.c;  end

	PPhase_FT_Cxx = reshape(PPhase_FT_Cxx, nSpatial, nTopEVs * NF )';
	SM            = boolean( spatialMask(:) );

	%NOTE: When boolParfor is true, it unwraps all the phase vectors but when false it only unwraps
	%those that are tagged for clustering. The first one is useful when one likes to plot and
	%investigate the phases in all the frequency but can be very slow.
	if pasfOpts.boolParfor == 1,
		%if isempty(gcp('nocreate')),  parpool(maxNumCompThreads);  end
		parfor i = 1:size(PPhase_FT_Cxx, 1),
			PPhase_FT_Cxx(i, :) = unwrap2D_phase_vector(PPhase_FT_Cxx(i, :), SM, d1, d2);
		end
	else
		for i = find( scoreMask(:)' ),
			PPhase_FT_Cxx(i, :) = unwrap2D_phase_vector(PPhase_FT_Cxx(i, :), SM, d1, d2);
		end
	end

	if ~pasfOpts.boolQuietly,  disp('done with the phase unwrapping.');  end

	if ~pasfOpts.boolQuietly,
		fprintf('done with the unwrapping.');
	end



	%%NOTE: kmeans/linkage accepts the point cloud as a matrix whose rows are the data point.
	Clusts = zeros(size(scoreMask)); 
	half   = 2 + floor(NF/2);
	scoreMask(:, half:end) = 0;

	if strcmp( pasfOpts.cmethod, 'phase'),
		datapoints = PPhase_FT_Cxx;
	elseif strcmp( pasfOpts.cmethod, 'weightedPhase'),
		datapoints = sqrt(reshape(abs(PVec_FT_Cxx), nSpatial, nTopEVs * NF )') .* PPhase_FT_Cxx;
	elseif strcmp( pasfOpts.cmethod, 'maskedPhase')
		datapoints = boolean(reshape(abs(PVec_FT_Cxx), nSpatial, nTopEVs * NF )'>1/sqrt(nSpatial)) .* PPhase_FT_Cxx;
	elseif strcmp( pasfOpts.cmethod, 'modulus'),
		datapoints = reshape(abs(PVec_FT_Cxx), nSpatial, nTopEVs * NF )';
	else
		error('undefined cmethod!');
	end

	%freqMask = boolean(freqMask .* scoreMask);
	freqMask = scoreMask;
	[Clusts, ClusterInfo] = robust_cluster(datapoints, ...
		Eigenvalues, ...
		PPhase_FT_Cxx, ...
		reshape(abs(PVec_FT_Cxx), nSpatial, nTopEVs * NF )',... 
		freqMask, spatialMask, nSpatialComponents, pasfOpts);

	if mod(NF, 2) == 0, %%NOTE: We ran clustering for the first half of frequencies.
		Clusts(:, half:end) = Clusts(:, half-2:-1:2);
	else
		Clusts(:, half:end) = Clusts(:, half-1:-1:2);
	end
	Clusts = Clusts(:); 
	ClusterInfo.CEnergy = 2*ClusterInfo.CEnergy;

	if ~pasfOpts.boolQuietly,
		disp('done with the clustering.');
	end


	%NOTE: Construct the filters and then apply it to the input spatio-temporal signal.
	C_omega = PVec_FT_Cxx;
	B_omega = zeros(nTopEVs, nSpatial, size(PVec_FT_Cxx, 3) );
	for i=1:size(PVec_FT_Cxx, 3),
		B_omega(:, :, i) = PVec_FT_Cxx(:, :, i)';
	end

	Components = filter_spatiotemporal_signal(Z, C_omega, B_omega, Clusts, nSpatialComponents);

	if ~pasfOpts.boolQuietly,
		disp('done with the filtering/separating dynamic components.');
	end


	%%preparing the outputs
	Components = permute(Components, [2, 1, 4, 3]);
	Components = reshape(Components, d1, d2, d3, nSpatialComponents);
	Z          = reshape(permute(Z, [2, 1]), d1, d2, d3);
	Err        = Z - sum(Components, 4);
	Components(:, :, :, nSpatialComponents + 1) = Err; 
	Components(:, :, :, nSpatialComponents + 2) = Z;

	Clusters     = reshape(Clusts, nTopEVs, NF);
	%NOTE: PVec_FT_Cxx \in C^{S x nTopEVs x T};
	SDFInfo.Eigenvectors = permute(reshape(abs(PVec_FT_Cxx), d1, d2, nTopEVs, NF), [1, 2, 4, 3]);
	SDFInfo.Eigenvalues = Eigenvalues;
	SDFInfo.Traces      = Traces;

	SDFInfo.PhaseArray = zeros(d1, d2, NF, nTopEVs);
	for j = 1:NF, 
		for k = 1:nTopEVs,
			SDFInfo.PhaseArray(:, :, j, k) = reshape(PPhase_FT_Cxx( (j-1)*nTopEVs+k, :), d1, d2);
		end
	end


	if ~pasfOpts.boolQuietly,
		disp( strcat('Error rate is: ', num2str(sum(Err(:).^2)/sum(Z(:).^2)), '.') );
	end


	if pasfOpts.saveData > 1,
		filename = strcat(pasfOpts.bfilename, '_output.mat');
		save(filename, 'Components', 'Clusters', 'SDFInfo', 'pasfOpts', 'ClusterInfo', '-v7.3');

		if pasfOpts.saveData > 2,
			Eigenvalues  = SDFInfo.Eigenvalues;
			Eigenvectors = SDFInfo.Eigenvectors;
			save(strcat(pasfOpts.bfilename, '_Components.mat'), 'Components',                  '-v7.3');
			save(strcat(pasfOpts.bfilename, '_Eigenpairs.mat'), 'Eigenvalues', 'Eigenvectors', '-v7.3');
			save(strcat(pasfOpts.bfilename, '_Clusters.mat'  ), 'Clusters', 'ClusterInfo',     '-v7.3');
			save(strcat(pasfOpts.bfilename, '_Filters.mat'  ),  'C_omega', 'B_omega',          '-v7.3');
		end
	end


end %end of pasf



function Components = filter_spatiotemporal_signal(Z, C_omega, B_omega, Clusts, nSpatialComponents)

	[nTopEVs, nSpatial, NF] = size(B_omega);
	nTimeSamples = size(Z, 1);

	Components = zeros(nTimeSamples, nSpatial, nSpatialComponents);

	for cNum = 1:nSpatialComponents,
		boolClust = (Clusts == cNum);
		mask = repmat( reshape(boolClust , nTopEVs, 1, nTimeSamples), 1, nSpatial, 1);

		X = zeros(NF, nTopEVs);
		tZ = zeros(nTimeSamples, 1);

		for i=1:size(B_omega,  1),
			for j=1:size(B_omega, 2),
				tZ(1:nTimeSamples) = squeeze(Z(:,j));
				X(:, i) = X(:, i) + squeeze(B_omega(i, j, :).*mask(i, j, :)) .* fft(tZ);
			end
		end

		mask  = permute(mask, [2, 1, 3]);
		Z_hat = zeros(NF, nSpatial);

		for i=1:size(C_omega, 1),
			for j=1:size(C_omega, 2),
				Z_hat(:, i) = Z_hat(:, i) + ifft(squeeze(C_omega(i, j, :) .* mask(i, j, :)) .* squeeze(X(:, j)) );
			end
		end
		Components(:, :, cNum) = real(Z_hat(1:nTimeSamples, :)); 
	end

end



function unwrappedPhaseVec = unwrap2D_phase_vector(phaseVec, mask, d1, d2)
	phaseVec = reshape(phaseVec, d1, d2);
	SM = boolean(mask);
	notSM = double(~SM); %NOTE: The concept of mask changes here cause unwrap2D(.) uses it differently.
	unwrappedPhaseVec = reshape( unwrap2D( phaseVec, reshape(notSM, d1, d2), 128), 1, d1*d2);
	%cosmic changes for ploting
	unwrappedPhaseVec(SM) = unwrappedPhaseVec(SM) - mean( unwrappedPhaseVec(SM) ); 	
end



function [Clusts, CI] = robust_cluster(vectors, weights, phase, modulus, ... 
		mask, spatialMask, nClusters, pasfOpts)

	%TODO merge it with the prev function
	Clusts  = zeros(size(mask)); 
	SM      = boolean(spatialMask(:));
	mask    = mask(:);

	CI.Tree = linkage(vectors, 'ward', 'corr');
	Memberships = mask(mask);
	kmeansOpt = statset('UseParallel', pasfOpts.boolParfor);

	for outlierThreshold = 1:-0.05:min(pasfOpts.outlierThreshold, 1),
		mask(mask)   = (Memberships > 0); 
		norm_vectors = normalize_rows( vectors(mask, SM) );

		[Memberships, ctrs, ~, cdist] = kmeans( norm_vectors, nClusters, ... 
			'Options', kmeansOpt, ... 
			'Distance', 'correlation', 'Replicates', 128);
		for c = 1:nClusters,
			outlierMask = (Memberships == c) & (cdist(:, c) > outlierThreshold);
			Memberships( outlierMask ) = 0; 
		end
	end
	CEnergy = zeros(nClusters, 1);


	%%sorting based on the energy level
	for c = 1:nClusters,
		CEnergy(c) = sum( weights(Memberships == c) );
	end
	CMap = rename_cluster_map(Memberships, CEnergy); 
	Memberships(Memberships>0) = CMap(Memberships(Memberships>0));



	modulus = modulus(:, SM);
	phase   = phase(:, SM);

	ctrs      = zeros(size(ctrs));
	modctrs   = zeros(size(ctrs));
	phasectrs = zeros(size(ctrs));

	totalEnergy = sum(weights(:));
	for c = 1:nClusters,

		vecs = repmat(weights(Memberships == c), 1, size(norm_vectors,2)) .* norm_vectors(Memberships == c, :);
		ctrs(c, :) = mean( vecs, 1);

		vecs = repmat(weights(Memberships == c), 1, size(modulus, 2)) .* modulus(Memberships == c, :);
		modctrs(c, :) = mean( vecs, 1);

		vecs = repmat(weights(Memberships == c), 1, size(phase, 2)) .* phase(Memberships == c, :);
		phasectrs(c, :) = mean( vecs, 1);

		CEnergy(c) = sum( weights(Memberships == c) )/totalEnergy;
	end


	CI.D                       = cdist;
	CI.D(:, CMap(1:nClusters)) = cdist(:, 1:nClusters);

	CI.Centers        = zeros( nClusters, length(SM) );
	CI.Centers(:, SM) = ctrs;

	CI.ModCenters        = zeros( nClusters, length(SM) );
	CI.ModCenters(:, SM) = modctrs;


	CI.PhaseCenters        = zeros( nClusters, length(SM) );
	CI.PhaseCenters(:, SM) = phasectrs;

	CI.Mem            = round(Memberships);
	CI.CSizes         = histc(Memberships, 1:nClusters);
	CI.CEnergy        = CEnergy;
	Clusts(mask)      = CI.Mem; 
end



function mat = normalize_rows( mat )
	nRows = size(mat, 2);
	mat = mat -  repmat( mean(mat, 2), 1, nRows ); 
	mat = mat ./ repmat( sqrt(sum(mat.^2, 2)), 1, nRows );
end


function CMap = rename_cluster_map(mem, weightVec)
	%% Reindex the cluster names in decreasing order of their sizes.
	nC      = max(mem);
	[~,J]   = sort(weightVec, 'descend');
	J(J)    = 1:nC; %build the map to rename the cluster names. 
	CMap    = J;
end


%%It estimates the multivariate spectral density function.
function [PVec_FT_Cxx, PPhase_FT_Cxx, Eigenvalues, N, Traces] =  mvspec(X, opts)

	N0 = size(X,1); 
	N  = N0;
	ncols  = size(X,2);
	xfreq  = 1;
	nTopEVs = opts.nTopEVs;

	if(opts.detrend == 1 )
		X = detrend(X, 'constant');
	elseif( opts.detrend == 2 )
		X = detrend(X);
	end

	if opts.taperLen > 0,
		w = taper(N, opts.taperLen)';
		X = X .* repmat(w, 1, ncols);
	end


	if opts.boolZeroPad ==  1, 
		X     = [X; zeros(N-1, ncols)];
		N     = size(X,1); 
		ncols = size(X,2);
	end

	xfft = fft(X);
	clearvars X; %we don't need it anymore.

	if ~opts.boolQuietly,
		disp('done with xfft computation.');
	end


	spans = 2*floor(opts.spans/2)+1; %%just to make the smoothing easier
	switch opts.kmethod,
		case 'box'
			g = ones(spans,1);
		case 'gaussian'
			%NOTE: e.g. gausswin(7)=[0.0439, 0.2494, 0.7066, 1.00, 0.7066, 0.2494, 0.0439]'
			g = gausswin(spans);  
		case 'triangle'
			%NOTE: e.g.   triang(7)=[0.25, 0.50, 0.75, 1.00, 0.75, 0.5, 0.25]' 
			g = triang(spans);   
		otherwise
			error('kmethod is undefined [box|gaussian|triangle]');
	end
	g = g/sum(g);

	PVec_FT_Cxx   = zeros(ncols, nTopEVs, N); 
	PPhase_FT_Cxx = zeros(ncols, nTopEVs, N); 
	Eigenvalues   = zeros(nTopEVs, N);
	Traces        = zeros(N, 1);

	half = ceil(spans/2);
	nFreq = floor(N/2)+1;
	for i=2:nFreq,   %NOTE: Estimate the spectral density function by smoothing the periodogram.

		M = zeros(ncols, ncols);

		%NOTE: Since all the coefficients, g_j,  are positive M is positive semi-definite.
		for j=1:spans, %with the assumption that span is an even number
			M = M + g(j) * get_periodogram(xfft, N, ncols, N0*xfreq, i, j-half);
		end
		Traces(i)     = trace(M);
		Traces(N-i+2) = Traces(i);

		[V, D] = svds(M, nTopEVs, 'L'); 
		D      = diag(D);

		PVec_FT_Cxx(:, :, i)     = V; 
		PVec_FT_Cxx(:, :, N-i+2) = conj(V); 

		Eigenvalues(:, i)      = D;  
		Eigenvalues(:, N-i+2)  = D;  

		PPhase_FT_Cxx(:, :, i)     = angle(V);
		PPhase_FT_Cxx(:, :, N-i+2) = angle(V); %%TODO I guess it should be -angle(V)!?

		if (~opts.boolQuietly) && (floor((i-1)*25/nFreq) < floor(i*25/nFreq)),
			fprintf('%d%c ', floor(i*100/nFreq), '%');
			if i == nFreq, fprintf('\n'); end
		end

	end

end 


function pgram = get_periodogram(xfft, N, ncols, const, i, offset)
	idx = i + offset;
	if( idx < 2 ) %ignoring zero freq;
		%%TODO
		idx = N + idx - 1; %-1 to ignore zero freq 
		%idx = 1-(idx - 1); %-1 to ignore zero freq 
	end
	if( idx > N )
		idx = idx - N;
	end
	pgram = zeros(ncols, ncols);
	%pgram = xfft(idx, :).' * conj( xfft(idx, :) )/(const);
	pgram = kron(xfft(idx, :).', conj(xfft(idx, :)))/(const); % eqv to the above line but faster!
end


function w = taper(n, p) 
	if( p > 0.5 || p < 0 )
		error('p must be between 0 and 0.5');
	end

	m = floor(n * p);
	w = 0.5 * (1 - cos(pi * (1:2:2*m-1)/(2*m)) );
	w = [w, ones(1, n - 2*m), flip(w)];
end

