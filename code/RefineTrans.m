
function trans = RefineTrans(Img, trans, lamda, sigma)
% check the input is a three-channel image
check3Dim(Img);

if ~exist( 'lamda', 'var' ),
    lamda = 2;
    warning( 'lamda not set!  Defaulting to: %f\n', lamda );
end

if ~exist( 'sigma', 'var' ),
    sigma = 0.5;
    warning( 'sigma not set!  Defaulting to: %f\n', sigma );
end
% estimating the transmission function
[nRows, nCols] = size(trans);

% first-order differential filter
% also we can construct high-order filters 
d{1} = [-1, -1, -1; -1, 8, -1; -1, -1, -1];

% normalization
num_filters = length(d); 
for k = 1 : num_filters
    d{k} = d{k} / norm(d{k}(:));
end

% calculating weighting function
for k = 1 : num_filters
    WFun{k} = WeighComputing(Img, d{k}, sigma);
end

% FFT for constant
Tf = fft2(trans);
DS = 0; 
for k = 1 : num_filters
    D{k} = psf2otf(d{k}, [nRows, nCols]);
    DS = DS + abs(D{k}).^2;
end

% cyclic looping
beta = 1; beta_rate = 2 * sqrt(2);
beta_max = 2^5;
Outiter = 0;

while beta < beta_max
    % updating parameters    
    gamma = lamda / beta;
    % show the results
    Outiter = Outiter + 1; 
    % fprintf('Outer iteration %d; beta %.3g\n', Outiter, beta);    
    figure(100), imshow(trans, []); title(num2str(Outiter)); pause(0.05);     

    % fixing t, solving u
    DU = 0;
    for k = 1 : num_filters
        dt{k} = imfilter(trans, d{k}, 'circular');
        u{k} = max(abs(dt{k}) - WFun{k} / beta / num_filters, 0) .* sign(dt{k}); 
        DU = DU + fft2(imfilter(u{k}, flipud(fliplr(d{k})), 'circular')); 
    end
    % fixing u, solving t;    
    trans = RemoveOuts(abs(ifft2((gamma * Tf + DU) ./ ( gamma + DS))));  
    % increasing beta
    beta = beta * beta_rate;
end

close(100);


