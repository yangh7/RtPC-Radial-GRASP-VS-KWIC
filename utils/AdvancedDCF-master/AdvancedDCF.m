
function [DCF] = AdvancedDCF( AngleInfo_input, Option ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by KyungPyo Hong
% % Version 1.0: February, 2019.
% % Version 1.1: added a function of weighted contrast (KWIC)
% %
% %   - AngleInfo_input (2-dimension) = A x B
% %        A: angle information
% %        B: repeats (or time frames) 
% %   - Option.Nsamples: the number of sampled points in a ray (even number)
% %   - Option.AnglePrecision: round-up for matching up with the precision of AngleInfo_input
% %                            (0 = full precision)
% %   - Option.WeightedContrast: weight-adjust for various image contrasts 
% %                              (Song HK, et al. Magn Reson Med. 2000 Dec;44(6):825-32.) 
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unit metric
    %   # unit distance
    %       - a minimal discrete distance between two neighboring grids
    %   # for 2D radial (angular distance)
    %       - Unit metric = unit distance
    %   # for 3D radial (surface area)
    %       - Unit metric = unit area = (unit distance)^2
    UnitDistance = 1;
    
    if ~isfield( Option, 'Dimension' )
        Option.Dimension = '2D'
    end
    switch Option.Dimension
        case '3D'
            UnitMetric = UnitDistance^2;
            Dimension = 3;
        otherwise  % for 2D
            UnitMetric = UnitDistance;
            Dimension = 2;
    end
    
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

    if Dimension == 3 
        ray_info = -0.5:1/NM:0.5-1/NM;
        ray_grid = zeros( NM, Nray, Nrep, 3 );
        for ii = 1:Nrep
            for jj = 1:Nray
                ray_grid(:,jj,ii,1) = ray_info.*sin(AngleInfo_input(jj,2,ii)).*cos(AngleInfo_input(jj,1,ii)); 
                ray_grid(:,jj,ii,2) = ray_info.*sin(AngleInfo_input(jj,2,ii)).*sin(AngleInfo_input(jj,1,ii)); 
                ray_grid(:,jj,ii,3) = ray_info.*cos(AngleInfo_input(jj,2,ii)); 
            end
        end
    end    

    if Dimension == 2
        % prepare angle info
        AngleInfo = mod( [AngleInfo_input;   AngleInfo_input+pi], 2*pi );  
        if CalcPrecision > 1
            AngleInfo = round( AngleInfo * CalcPrecision )/CalcPrecision;
        end
    end    
    
    DCF = zeros( NM, Nray, Nrep ); 
    for TimeFrame = 1:Nrep
        fprintf( 'TimeFrame = %3d\n', TimeFrame );
        [tmp_ NR_ind ] =  find( WeightedContrast(NM/2+1,:,TimeFrame) > 0 ); 
        NumOfNR = length( NR_ind );
        DCF( NM/2+1, NR_ind, TimeFrame ) = 1/NumOfNR;   % at radius = 0
        
        if Dimension == 3
            VectorPos = [squeeze(ray_grid(NM/2+2,:,TimeFrame,1))' squeeze(ray_grid(NM/2+2,:,TimeFrame,2))' squeeze(ray_grid(NM/2+2,:,TimeFrame,3))' ];
            VectorNeg = [squeeze(ray_grid(NM/2,:,TimeFrame,1))' squeeze(ray_grid(NM/2,:,TimeFrame,2))' squeeze(ray_grid(NM/2,:,TimeFrame,3))' ];
            TwoVectors = [VectorPos; VectorNeg];
        end
        tmp_DCF = zeros( NM/2-1, Nray*2 );
        parfor radius = 1:NM/2-1
            fprintf( 'radius = %3d\n', radius )
            
            tmp_WC = WeightedContrast( [NM/2+1-radius  NM/2+1+radius], :, TimeFrame )';
            [ WC_ind  tmp_ ] = find( tmp_WC(:) > 0 );
            
            SumOfCorrelations = zeros( 1, Nray*2 );
            tmp_calc = [];
            Correlation = [];
            for ray_loc = WC_ind'
                if Dimension == 2
                    % % dL = Unit Metric
                    % %    = radius * |theta2 - theta1|
                    ang1 = AngleInfo( ray_loc, TimeFrame );
                    ang2 = AngleInfo( WC_ind, TimeFrame );
                    ang3 = max( ang1, ang2 );
                    ang4 = min( ang1, ang2 );
                    tmp_calc(:,1) = abs( ang3 - ang4 );
                    tmp_calc(:,2) = abs( ang3 - (ang4+2*pi) );
                    tmp_sel = sort( tmp_calc, 2, 'ascend' );
                    Correlation  =  (UnitMetric - radius*tmp_sel(:,1))/UnitMetric;
                elseif Dimension == 3
                    % % correlation for 3D radial: based on the surface area of sphere 
                    % %  - dA = Unit Metric = unit area = (unit distance)^2
                    % %       = radius^2 * (theta2-theta1) * |cos(phi1)-cos(phi2)|
                    % %       = 2*pi*(radius^2)*(1-cos(AngularDistanceBtwTwoPoints)) at the pole (r,0,0)

                    % % calculate an angular distance between two points on the surface at a given radius
                    dotAB = sum( TwoVectors(ray_loc,:).*TwoVectors(WC_ind,:), 2 );
                    tmp_AngDist = atan2( sqrt( sum( cross( ones(length(WC_ind),3).*TwoVectors(ray_loc,:), TwoVectors(WC_ind,:) ).^2,2)), dotAB );
                    AngularDistance = abs( tmp_AngDist );
                    Correlation = (UnitMetric - (2*pi*radius^2)*(1-cos(AngularDistance)))/UnitMetric;
                end
                if CalcPrecision > 1
                    Correlation = round( Correlation*CalcPrecision )/CalcPrecision;
                end
                [tmp_ind, ~] = find( Correlation <= 0 );
                Correlation( tmp_ind ) = 0;
                SumOfCorrelations(ray_loc) = SumOfCorrelations(ray_loc) + sum(Correlation);
            end
            tmp_DCF( radius, : ) = 1./SumOfCorrelations;
        end
        DCF( NM/2:-1:2, :, TimeFrame ) = tmp_DCF(:,1:Nray);
        DCF( NM/2+2:end, :, TimeFrame ) = tmp_DCF(:,Nray+1:end);
    end
    DCF( find( DCF == inf ) ) = 0;
    DCF(1,:,:) = DCF(2,:,:);
end