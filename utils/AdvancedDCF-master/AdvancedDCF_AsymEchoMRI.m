
function [DCF_AsymEcho] = AdvancedDCF_AsymEchoMRI( AngleInfo_input, Option ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by KyungPyo Hong 
% %  For asymmetric echo MRI: e.g., ultrashort echo time (UTE)
% %
% %   - AngleInfo_input (2-dimension) = A x B
% %        A: angle information
% %        B: repeats (or time frames) 
% %   - Option.Nsamples: the number of sampled points in a ray (even number)
% %   - Option.AnglePrecision: round-up for matching up with the precision of AngleInfo_input
% %                            (0 = full precision)
% %   - Option.WeightedContrast: weight-adjust for various image contrasts 
% %                              (Song HK, et al. Magn Reson Med. 2000 Dec;44(6):825-32.) 
% %   - Option.AsymEcho: define asymmetric echo
% %       Fully sampled:   |*************$************|   ($ = the center of frequency encoding)
% %       Asym. echo:      |---------****$************|   (- = empty) 
% %         -> Option.AsymEcho = the length of star marks (*) including $ mark
% %                              (NM/2 <= Option.AsymEcho <= NM)
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NM = Option.Nsamples;
    [ NR  Nrep ] = size( AngleInfo_input );

    if ~isfield( Option, 'AsymEcho' ) 
        Option.AsymEcho = NM;
    end
   
    if ~isfield( Option, 'AnglePrecision' ) 
        Option.AnglePrecision = 0;
    end
    AnglePrecision = 10^(Option.AnglePrecision);

    if ~isfield( Option, 'WeightedContrast' ) || length( Option.WeightedContrast(:) ) == 1
        WeightedContrast = ones( NM, NR, Nrep );
    else
        WeightedContrast = Option.WeightedContrast;
    end
    
    if ~isfield( Option, 'Display' )
        Option.Display = 0;
    end
    
    if Option.Display == 1
        ray_info = -0.5:1/NM:0.5-1/NM;
        figure(101)
        for ii = 1:Nrep
            for jj = 1:NR
                GridInfo = ray_info*exp( 1i*AngleInfo_input(jj,ii) );
                if Option.AsymEcho < NM
                    GridInfo( 1:NM-Option.AsymEcho ) = 0;
                end
                plot( real( GridInfo ), imag( GridInfo ) )
                hold on
            end
            axis square; 
        end
    end

    % prepare angle info
    AngleInfo = mod( [AngleInfo_input;   AngleInfo_input+pi], 2*pi );  
    AngleInfo = permute( repmat( AngleInfo, [ 1  1  NM] ), [3 1 2] );
    if Option.AsymEcho < NM
        AngleInfo( 1:NM-Option.AsymEcho, 1:NR, : ) = -1;
    end
    if AnglePrecision > 1
        AngleInfo = round( AngleInfo * AnglePrecision )/AnglePrecision;
    end
    
    DCF_AsymEcho = zeros( NM, NR, Nrep ); 
    for TimeFrame = 1:Nrep
        fprintf( 'TimeFrame = %3d\n', TimeFrame );
        [tmp_ NR_ind ] =  find( WeightedContrast(NM/2+1,:,TimeFrame) > 0 ); 
        NumOfNR = length( NR_ind );
        DCF_AsymEcho( NM/2+1, NR_ind, TimeFrame ) = 1/NumOfNR;   % at radius = 0
        for Radius = 1:NM/2-1
% % %             fprintf( 'Radius = %3d\n', Radius )
            NumOfGrids = (2*Radius+1)^2 - (2*Radius-1)^2;
            AngleCrit = 2*pi/NumOfGrids;
            if AnglePrecision > 1
                AngleCrit = round( AngleCrit*AnglePrecision )/AnglePrecision;
            end
            
            tmp_WC = WeightedContrast( [NM/2+1-Radius  NM/2+1+Radius], :, TimeFrame )';
            [ WC_ind  tmp_ ] = find( tmp_WC(:) > 0 );
            DegreeOfCorrelations = zeros( 1, NR*2 );
            tmp_calc = [];
            for ray_loc = WC_ind'
                Correlations = zeros( 1, length(WC_ind)*2 );
                ang1 = AngleInfo( NM/2+1-Radius, ray_loc, TimeFrame );
                ang2 = AngleInfo( NM/2+1-Radius, WC_ind, TimeFrame );
                if ang1 == -1
                    Correlations = 0;
                else                    
                    ang2 = ang2( find( ang2 >= 0 ) );
                    ang3 = max( ang1, ang2 );
                    ang4 = min( ang1, ang2 );
                    tmp_calc(:,1) = abs( ang3 - ang4 );
                    tmp_calc(:,2) = abs( ang3 - (ang4+2*pi) );
                    if AnglePrecision > 1
                        tmp_calc = round( tmp_calc*AnglePrecision )/AnglePrecision;
                    end
                    tmp_sel = sort( tmp_calc, 2, 'ascend' );
                    AngularDistance = tmp_sel( :, 1 );
                    [tmp_ind1 tmp_] = find( AngularDistance >= AngleCrit );
                    [tmp_ind2 tmp_] = find( AngularDistance < AngleCrit );
                    Correlations( tmp_ind1 ) = 0;
                    Correlations( tmp_ind2 ) = (1 - (AngularDistance(tmp_ind2)/AngleCrit));
                end
                DegreeOfCorrelations(ray_loc) = DegreeOfCorrelations(ray_loc) + sum(Correlations);
            end
            DCF_AsymEcho( NM/2+1-Radius, 1:NR, TimeFrame ) = 1./DegreeOfCorrelations(1:NR);
            DCF_AsymEcho( NM/2+1+Radius, 1:NR, TimeFrame ) = 1./DegreeOfCorrelations(NR+1:end);
        end
    end
    DCF_AsymEcho( find( DCF_AsymEcho == inf ) ) = 0;
    DCF_AsymEcho(1,:,:) = DCF_AsymEcho(2,:,:);
end