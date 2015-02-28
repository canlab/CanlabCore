function K = emp_feat_map_segments( kern, dat1, dat2, ind1, ind2, kerparam)

% K = emp_feat_map_segments( kern, dat1, dat2, ind1, ind2, kerparam)
% 
% The same as 'emp_feat_map' but can additionally handle segments of
% features seperately (e.g. data that contains several pictures in one long
% vector). 
%
% The kernel parameters contain the data- and the kernel-objects
% needed and the number of segments (all of same size -> can be improved).
%
% Kernelparameters (with defaults):
% kerparam = { prototypes, kernel_base, kernel_2nd, nofSegs }
% 
% prototypes                -- data object that contains the prototypes
% kernel_base               -- kernel object to compute feature vectors
%                              (needs not to be a proper kernel, btw)
% kernel_2nd = 'linear'     -- kernel object applied to these feature vectors
% nofSegs = 1               -- no. of segments in case features are composed of segments

prototypes = kerparam{ 1};
kernel_base = kerparam{ 2};

l = length( kerparam);
switch l
    case 4
        kernel_2nd = kerparam{ 3};
        nofSegs = kerparam{ 4};
    case 3
        switch class( kerparam{ 3})
            case 'kernel'
                kernel_2nd = kerparam{ 3};
                nofSegs = 1;
            case 'double'
                kernel_2nd = kernel;
                nofSegs = kerparam{ 3};
            otherwise
                error( 'wrong type of kernel parameter');
        end
    case 2
        kernel_2nd = kernel;
        nofSegs = 1;
    otherwise
        error( 'not enough kernel parameters');
end


[sz1, fsz1] = get_dim( dat1);
[sz2, fsz2] = get_dim( dat2);
[szprot, fszprot] = get_dim( prototypes);

% test if features of data and prototypes agree
if ( fsz1 ~= fszprot) | ( fsz1 ~= fsz2) 
    error( '!! no. of features for data and prototypes does not agree !!'); end

% compute length of segments from total length and nofSegs
if mod( fsz1, nofSegs) == 0 
    lengthOfSegs = fsz1 / nofSegs;
else error( 'no. of features is not a multiple of no. of segments'); end


for i = 1:nofSegs
    findx = ( i - 1)*lengthOfSegs + 1 : i*lengthOfSegs;
    feat_vec1{ i} = calc( kernel_base, get( prototypes, 1:szprot, findx), get( dat1, ind1, findx));
    feat_vec2{ i} = calc( kernel_base, get( prototypes, 1:szprot, findx), get( dat2, ind2, findx));
end

feat_vec1 = [ feat_vec1{ :}];
feat_vec2 = [ feat_vec2{ :}];

%keyboard;
K = calc( kernel_2nd, data( feat_vec1), data( feat_vec2));
