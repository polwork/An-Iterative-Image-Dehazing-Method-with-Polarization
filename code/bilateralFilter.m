
%

function output = bilateralFilter( data, edge, edgeMin, edgeMax, sigmaSpatial, sigmaRange, ...
    samplingSpatial, samplingRange )

if( ndims( data ) > 2 ),
    error( 'data must be a greyscale image with size [ height, width ]' );
end

if( ~isa( data, 'double' ) ),
    error( 'data must be of class "double"' );
end

if ~exist( 'edge', 'var' ),
    edge = data;
elseif isempty( edge ),
    edge = data;
end

if( ndims( edge ) > 2 ),
    error( 'edge must be a greyscale image with size [ height, width ]' );
end

if( ~isa( edge, 'double' ) ),
    error( 'edge must be of class "double"' );
end

inputHeight = size( data, 1 );
inputWidth = size( data, 2 );

if ~exist( 'edgeMin', 'var' ),
    edgeMin = min( edge( : ) );
    warning( 'edgeMin not set!  Defaulting to: %f\n', edgeMin );
end

if ~exist( 'edgeMax', 'var' ),
    edgeMax = max( edge( : ) );
    warning( 'edgeMax not set!  Defaulting to: %f\n', edgeMax );
end

edgeDelta = edgeMax - edgeMin;

if ~exist( 'sigmaSpatial', 'var' ),
    sigmaSpatial = min( inputWidth, inputHeight ) / 14;
    fprintf( 'Using default sigmaSpatial of: %f\n', sigmaSpatial );
end

if ~exist( 'sigmaRange', 'var' ),
    sigmaRange = 0.006 * edgeDelta;
    fprintf( 'Using default sigmaRange of: %f\n', sigmaRange );
end

if ~exist( 'samplingSpatial', 'var' ),
    samplingSpatial = sigmaSpatial;
end

if ~exist( 'samplingRange', 'var' ),
    samplingRange = sigmaRange;
end

if size( data ) ~= size( edge ),
    error( 'data and edge must be of the same size' );
end

% parameters
derivedSigmaSpatial = sigmaSpatial / samplingSpatial;
derivedSigmaRange = sigmaRange / samplingRange;

paddingXY = floor( 2 * derivedSigmaSpatial ) + 1;
paddingZ = floor( 2 * derivedSigmaRange ) + 1;

% allocate 3D grid
downsampledWidth = floor( ( inputWidth - 1 ) / samplingSpatial ) + 1 + 2 * paddingXY;
downsampledHeight = floor( ( inputHeight - 1 ) / samplingSpatial ) + 1 + 2 * paddingXY;
downsampledDepth = floor( edgeDelta / samplingRange ) + 1 + 2 * paddingZ;

gridData = zeros( downsampledHeight, downsampledWidth, downsampledDepth );
gridWeights = zeros( downsampledHeight, downsampledWidth, downsampledDepth );

% compute downsampled indices
[ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 );

% ii =
% 0 0 0 0 0
% 1 1 1 1 1
% 2 2 2 2 2

% jj =
% 0 1 2 3 4
% 0 1 2 3 4
% 0 1 2 3 4

% so when iterating over ii( k ), jj( k )
% get: ( 0, 0 ), ( 1, 0 ), ( 2, 0 ), ... (down columns first)

di = round( ii / samplingSpatial ) + paddingXY + 1;
dj = round( jj / samplingSpatial ) + paddingXY + 1;
dz = round( ( edge - edgeMin ) / samplingRange ) + paddingZ + 1;

% perform scatter (there's probably a faster way than this)
% normally would do downsampledWeights( di, dj, dk ) = 1, but we have to
% perform a summation to do box downsampling
for k = 1 : numel( dz ),
       
    dataZ = data( k ); % traverses the image column wise, same as di( k )
    if ~isnan( dataZ  ),
        
        dik = di( k );
        djk = dj( k );
        dzk = dz( k );

        gridData( dik, djk, dzk ) = gridData( dik, djk, dzk ) + dataZ;
        gridWeights( dik, djk, dzk ) = gridWeights( dik, djk, dzk ) + 1;
        
    end
end

% make gaussian kernel
kernelWidth = 2 * derivedSigmaSpatial + 1;
kernelHeight = kernelWidth;
kernelDepth = 2 * derivedSigmaRange + 1;

halfKernelWidth = floor( kernelWidth / 2 );
halfKernelHeight = floor( kernelHeight / 2 );
halfKernelDepth = floor( kernelDepth / 2 );

[gridX, gridY, gridZ] = meshgrid( 0 : kernelWidth - 1, 0 : kernelHeight - 1, 0 : kernelDepth - 1 );
gridX = gridX - halfKernelWidth;
gridY = gridY - halfKernelHeight;
gridZ = gridZ - halfKernelDepth;
gridRSquared = ( gridX .* gridX + gridY .* gridY ) / ( derivedSigmaSpatial * derivedSigmaSpatial ) + ( gridZ .* gridZ ) / ( derivedSigmaRange * derivedSigmaRange );
kernel = exp( -0.5 * gridRSquared );

% convolve
blurredGridData = convn( gridData, kernel, 'same' );
blurredGridWeights = convn( gridWeights, kernel, 'same' );

% divide
blurredGridWeights( blurredGridWeights == 0 ) = -2; % avoid divide by 0, won't read there anyway
normalizedBlurredGrid = blurredGridData ./ blurredGridWeights;
normalizedBlurredGrid( blurredGridWeights < -1 ) = 0; % put 0s where it's undefined

% for debugging
% blurredGridWeights( blurredGridWeights < -1 ) = 0; % put zeros back

% upsample
[ jj, ii ] = meshgrid( 0 : inputWidth - 1, 0 : inputHeight - 1 ); % meshgrid does x, then y, so output arguments need to be reversed
% no rounding
di = ( ii / samplingSpatial ) + paddingXY + 1;
dj = ( jj / samplingSpatial ) + paddingXY + 1;
dz = ( edge - edgeMin ) / samplingRange + paddingZ + 1;

% interpn takes rows, then cols, etc
% i.e. size(v,1), then size(v,2), ...
output = interpn( normalizedBlurredGrid, di, dj, dz );
