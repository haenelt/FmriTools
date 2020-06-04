% SPMT contrast transformation

% This script transforms t-maps to p, -log(p), r (correlation coefficient),
% d-maps (effect size) or z-maps. 
%
% The following formulas are used:
% (1) r = sign(t) / sqrt( df / t*t + 1 )
% (2) d = 2*r / sqrt( 1 - sqrt(r) )
% (3) p = 1 - spm_Tcdf
% (4) -log10(1-p) = -log( 1 - spm_Tcdf )
% For the last case of log transformation this means that a p-value of 0.01
% is transformed to a value of 2.
%
% Examples:
% p-value	-log10(1-P)
% 0.1		1
% 0.05		1.3
% 0.01		2
% 0.001		3
% 0.0001	4
%
% All maps can be thresholded using height and extent threshold (cluster
% size) and corrections for multiple comparisons based on family-wise error
% (FWE) or false discovery rate (FDR). If no correction for multiple
% comparisons is applied, p for conjunctions is the p of the conjunction
% SPM.
%
% All spmT files in the SPM.mat directory will be transformed. Outputs get 
% a basename following the convention: Type_Contrast_Pheight_Pextent_K_Neg.
%
%   Type:      P    - p-value
%              logP - log p-value
%              R    - correlation coefficient
%              D    - effect size
%              T    - t-value
%              Z    - z-value
%   Contrast:  name used in the contrast manager
%   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")
%              pFWE - p-value with FWE correction in %
%              pFDR - p-value with FDR correction in %
%   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")
%              pkFWE - extent p-value with FWE correction in %
%   K:         extent threshold in voxels
%   Neg:       image also shows thresholded inverse effects (e.g. neg. 
%              values)
%
% The script was taken from cg_spmT2x.m v1.25 by Christian Gaser and
% adapted to my own purposes.

% created by Daniel Haenelt
% Date created: 17-04-2020
% Last modified: 26-05-2020

% input
spm_mat_in = '/data/pt_01880/Experiment3_Stripes/p2/disparity/GE_EPI2/contrast_all/SPM.mat'; % SPM.mat
path_output = '/data/pt_01880/Experiment3_Stripes/p2/disparity/results/'; % path where output is written
name_output = 'all'; % basename of output contrasts
name_sess = 'GE_EPI2'; % name of session (if multiple sessions exist)

% add spm to path
pathSPM = '/data/pt_01880/source/spm12';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathSPM);

% load spm.mat
load(spm_mat_in);

% input path
path_spm_mat = fileparts(spm_mat_in);

% select transformation
str = '1-p|-log(1-p)|correlation coefficient cc|effect size d|z-score|apply thresholds without conversion';
sel = spm_input('Convert t value to?',1,'m',str,1:6,1);

% select threshold
str = 'yes|no';
threshold = spm_input('apply height threshold','+1','b',str,[],1);

if strcmp(threshold,'yes')

    % select height threshold adjustment
    str = 'FWE|FDR|none';
    adjustment = spm_input('p value adjustment','+1','b',str,[],1);
    
    switch adjustment
        case 'FWE' % family-wise false positive rate
            u0  = spm_input('p value (family-wise error)','+1','r',0.05,1,[0,1]);
        case 'FDR' % false discovery rate
            u0  = spm_input('p value (false discovery rate)','+1','r',0.05,1,[0,1]);
        otherwise  % no adjustment
            u0  = spm_input('threshold (T or p value)','+1','r',0.001,1);
    end

    pk = spm_input('cluster size (k or p value)','+1','r',0,1,[0,Inf]);
    if (pk < 1) && (pk > 0)
        str = 'uncorrected|FWE corrected';
        extent_FWE = spm_input('p value (extent)','+1','b',str,[0 1],1);
    end

    if sel > 2
        str = 'yes|no';
        neg_results = spm_input('Show also negative values','+1','b',str,[1 0],1);
    else
        neg_results = 0;
    end

end

% get number calculated contrast files
spmT = dir(fullfile(path_spm_mat,'/spmT*.nii'));
nspmT = length(spmT);

% check if correct files are listed
for i = 1:nspmT
    if ~strcmp(spmT(i).name(1:6),'spmT_0')
        error('Only spmT_0* files can be used');
    end
end

for Ic = 1:nspmT
    
    % load volume
    Vspm = spm_vol(fullfile(path_spm_mat,spmT(Ic).name));

    % get volume information
    xCon = SPM.xCon;                    % contrast definitions
    df   = [xCon(Ic).eidf SPM.xX.erdf]; % degrees of freedom
    STAT = xCon(Ic).STAT;               % statistic indicator
    R    = SPM.xVol.R;			        % vector of resel counts
    S    = SPM.xVol.S;			        % volume (in vox)
    XYZ  = SPM.xVol.XYZ;		        % in-mask XYZ coordinates
    FWHM = SPM.xVol.FWHM;               % smoothness of components (in vox)
    v2r  = 1/prod(FWHM(~isinf(FWHM)));  % voxels to resels
    
    % get data array
    Z = spm_get_data(Vspm,XYZ);
    
    % apply threshold
    if strcmp(threshold,'yes')
    
        switch adjustment
            case 'FWE' % family-wise false positive rate
                u  = spm_uc(u0,df,STAT,R,1,S);
            case 'FDR' % false discovery rate
                u  = spm_uc_FDR(u0,df,STAT,1,Vspm,0);
            otherwise  % no adjustment
                if u0 <= 1
                    u = spm_u(u0,df,STAT);
                else
                    u = u0;
                end
        end

        % calculate height threshold filtering
        if neg_results
            Q = find((Z > u) | (Z < -u));
        else
            Q = find(Z > u);
        end

        % apply height threshold
        Z = Z(:,Q);
        XYZ = XYZ(:,Q);
        if isempty(Q)
            fprintf('No voxels survive height threshold u=%0.2g\n',u);
            continue
        end 

        % extent threshold
        if (pk < 1) && (pk > 0)
            if extent_FWE
                Pk = 1;
                k = 0;
                while (Pk >= pk && k<S)
                    k = k + 1;
                    [Pk, Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
            else
                Pn = 1;
                k = 0;
                while (Pn >= pk && k<S)
                    k = k + 1;
                    [Pk, Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
            end
        else
            k = pk;
        end

        % calculate extent threshold filtering
        A = spm_clusters(XYZ);
        Q = [];
        for i = 1:max(A)
            j = find(A == i);
            if length(j) >= k
                Q = [Q j]; 
            end
        end

        % eliminate voxels
        Z     = Z(:,Q);
        XYZ   = XYZ(:,Q);
        if isempty(Q)
            fprintf('No voxels survived extent threshold k=%0.2g\n',k);
            continue
        end
    end

    % transform data
    switch sel
        case 1
            t2x = 1-spm_Tcdf(Z,df(2));
        case 2
            t2x = -log10(max(eps,1-spm_Tcdf(Z,df(2))));
        case 3
            t2x = sign(Z).*(1./((df(2)./((Z.*Z)+eps))+1)).^0.5;
        case 4
            tmp = (df(2)./((Z.*Z)+eps))+1;
            t2x = 2./((1-(1./tmp)).*tmp).^0.5;
        case 5
            t2x = spm_t2z(Z, df(2));
        case 6
            t2x = Z;
    end
        
    % output basename
    switch sel
        case 1
            t2x_name = 'P';
        case 2
            t2x_name = 'logP';
        case 3
            t2x_name = 'R';
        case 4
            t2x_name = 'D';
        case 5
            t2x_name = 'Z';
        case 6
            t2x_name = 'T';
    end

    if isempty(name_output) && isempty(name_sess)
        basename = xCon(Ic).name;
    elseif isempty(name_output) && ~isempty(name_sess)
        basename =  [xCon(Ic).name '_' name_sess];
    elseif ~isempty(name_output) && isempty(name_sess)
        basename = [name_output '_' xCon(Ic).name];
    else
        basename = [name_output '_' xCon(Ic).name '_' name_sess];
    end
        
    if strcmp(threshold,'yes')
        switch adjustment
            case 'FWE'
                p_height_str = '_pFWE';
            case 'FDR'
                p_height_str = '_pFDR';
            otherwise
                p_height_str = '_p';
        end
            
        if (pk < 1) && (pk > 0)
            if extent_FWE
                p_extent_str = ['_pkFWE' num2str(pk*100)];
            else
                p_extent_str = ['_pk' num2str(pk*100)];
            end
        else
            p_extent_str = '';
        end
        
        k_extent_str = ['_k' num2str(k)];
        
        if neg_results
            neg_str = '_bi';
        else
            neg_str = '';
        end

        name = [t2x_name '_' basename p_height_str num2str(u0*100) p_extent_str k_extent_str neg_str];
    else
        name = [t2x_name '_' basename];
    end
    fprintf('Save %s\n', name);        

    % output path
    path_res = fullfile(path_output,t2x_name,'native');
    if ~exist(path_res, 'dir')
        mkdir(path_res)
    end
        
    % output filename
    out = deblank(fullfile(path_res,[name '.nii']));

    % reconstruct (filtered) image from XYZ & Z pointlist
    Y = zeros(Vspm.dim(1:3));
    OFF = XYZ(1,:) + Vspm.dim(1)*(XYZ(2,:)-1 + Vspm.dim(2)*(XYZ(3,:)-1));
    Y(OFF) = t2x;

    VO = Vspm;
    VO.fname = out;
    VO.dt = [spm_type('float32') spm_platform('bigend')];
    spm_write_vol(VO,Y);
    
end
