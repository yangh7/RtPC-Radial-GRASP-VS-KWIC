
function [DCF] = AdvancedDCF_2DRadial( AngleInfo_input, Option ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by KyungPyo Hong
% % Version 1.0:  February, 2019.
% % Version 1.1:  added a weighted-contrast function (KWIC)
% %
% %   - AngleInfo_input:  [0, 2*pi) 
% %       : A x B = rays x repeats (or time frames)
% %   - Option.Nsamples: the number of sampled points in a ray (even number)
% %   - Option.AnglePrecision: round-up for matching up with the precision of AngleInfo_input
% %                            (0 = full precision)
% %   - Option.WeightedContrast: weight-adjust for various image contrasts 
% %                              (Song HK, et al. Magn Reson Med. 2000 Dec;44(6):825-32.) 
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unit metric
    %   # unit distance
    %       - a discrete length
    %          ex) dr (Nyquist distance=1/FOV) = dL (arc length) 
    %   # for 2D radial (arc length)
    %       - Unit metric = unit distance = arc length (dL)
    %   # for 3D radial (surface area)
    %       - Unit metric = unit area = (unit distance)^2
    UnitDistance = 1;
    UnitMetric = UnitDistance;

    
    NM = Option.Nsamples;
    [ Nray  Nrep ] = size( AngleInfo_input );
    
    if ~isfield( Option, 'AnglePrecision' ) 
        Option.AnglePrecision = 0;
    end
    CalcPrecision = 10^(Option.AnglePrecision);
    
    if ~isfield( Option, 'WeightedContrast' ) || length( Option.WeightedContrast(:) ) == 1
        WeightedContrast = ones( NM, Nray, Nrep );
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
            for jj = 1:Nray
                GridInfo = ray_info*exp( 1i*AngleInfo_input(jj,ii) );
                plot( real( GridInfo ), imag( GridInfo ) )
                hold on
            end
            axis square; 
        end
    end

    % prepare angle info
    AngleInfo = mod( [AngleInfo_input;   AngleInfo_input+pi], 2*pi );  
    if CalcPrecision > 1
        AngleInfo = round( AngleInfo * CalcPrecision )/CalcPrecision;
    end
    

    
    DCF = zeros( NM, Nray, Nrep ); 
    for TimeFrame = 1:Nrep
        fprintf( 'TimeFrame = %3d\n', TimeFrame );
        [tmp_ NR_ind ] =  find( WeightedContrast(NM/2+1,:,TimeFrame) > 0 ); 
        NumOfNR = length( NR_ind );
        DCF( NM/2+1, NR_ind, TimeFrame ) = 1/NumOfNR;   % at radius = 0
        tmp_DCF = zeros( NM/2-1, Nray*2 );
        for radius = 1:NM/2-1
            fprintf( 'radius = %3d\n', radius )
            
            tmp_WC = WeightedContrast( [NM/2+1-radius  NM/2+1+radius], :, TimeFrame )';
            [ WC_ind  tmp_ ] = find( tmp_WC(:) > 0 );
            SumOfOverlaps = zeros( 1, Nray*2 );
            tmp_calc = [];
            for ray_loc = WC_ind'
                % % dL = Unit Metric
                % %    = radius * |theta2 - theta1|
                                
                ang1 = AngleInfo( ray_loc, TimeFrame );
                ang2 = AngleInfo( WC_ind, TimeFrame );
                ang3 = max( ang1, ang2 );
                ang4 = min( ang1, ang2 );
                tmp_calc(:,1) = abs( ang3 - ang4 );
                tmp_calc(:,2) = abs( ang3 - (ang4+2*pi) );
                tmp_sel = sort( tmp_calc, 2, 'ascend' );
                if CalcPrecision > 1
                    tmp_sel = round( tmp_sel*CalcPrecision )/CalcPrecision;
                end
                Overlap = (UnitMetric - radius*tmp_sel(:,1))/UnitMetric;
                [tmp_ind, ~] = find( Overlap <= 0 );
                Overlap( tmp_ind ) = 0;
                SumOfOverlaps(ray_loc) = SumOfOverlaps(ray_loc) + sum(Overlap);
            end
            tmp_DCF( radius, : ) = 1./SumOfOverlaps;
        end
        DCF( NM/2:-1:2, :, TimeFrame ) = tmp_DCF(:,1:Nray);
        DCF( NM/2+2:end, :, TimeFrame ) = tmp_DCF(:,Nray+1:end);
    end
    DCF( find( DCF == inf ) ) = 0;
    DCF(1,:,:) = DCF(2,:,:);
end