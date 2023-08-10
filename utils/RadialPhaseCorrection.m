function Results = RadialPhaseCorrection( PhaseCorrData, kSpaceData, option)
    [NM, NR, NI, NC] = size( kSpaceData );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate phase Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp( option.PhaseCorrection, 'CalTC' )
        
        kdata_data = kSpaceData;
        kdata_phase = PhaseCorrData;

        [NM1, NR1, NI1, NC1] = size( kdata_phase );
        count = 1;
        for ii = 1:NR1
            for jj = 1:NI1
                kdata_phase_reorder(:,count,:) = kdata_phase(:,ii,jj,:);
                count = count + 1;
            end
        end

        for jj = 1:NC1
            for kk = 1:2:size(kdata_phase_reorder,2)
                L1 = kdata_phase_reorder(:,kk,jj);
                L2 = kdata_phase_reorder(:,kk+1,jj);
                L2_1 = L2(NM:-1:1);
                L2_2 = circshift(L2_1,1,1);
                L1_fft = fft( abs( L1 ) );
                L2_fft = conj(fft(abs(L2_2)));

                if jj == 1 && kk == 1
                    %Huili
%                     figure(101)
%                     subplot(2,2,1); plot( abs(L1) )
%                     subplot(2,2,2); plot( abs(L2_1) )
                end        
                G = (L1_fft.*L2_fft);
                FT_G = fftshift( abs( fft(G) ) );
                % calculating support
                if 1
                    mask_sup = Supp( (FT_G), .1 );  % Hassan's 
                else
                    mask_sup = Supp( (FT_G), .9 );
% %                 mask_sup = Supp( (FT_G), .5 );
                end
                c_g_sup = angle( fftshift(G) ).*mask_sup;

                if jj == 1 && kk == 1
                    %Huili
%                     subplot(2,2,3); plot( c_g_sup )
                end

                % fitting by y
                index = find( mask_sup == 1 );
    % % %             data = [ ((index)),((c_g_sup(index))) ];
                data = [ index, c_g_sup(index) ];
                P = polyfit( data(:,1), data(:,2), 1 );
                phase_shift( floor(kk/2)+1, jj ) = P(1);
                yfit = polyval( P, data(:,1) );
                
                y = data(:,2);
                yresid = yfit - data(:,2);
                SSresid = sum( yresid.^2 );
                SStotal = (length(y)-1) * var(y);
                rsq = 1 - SSresid/SStotal;
                phase_shift_y_res( floor(kk/2)+1, jj ) = rsq;
                shift_info( floor(kk/2)+1, jj ) = -P(1)*NM/(4*pi);
                L1_norm( floor(kk/2)+1, jj ) = norm( L1, 2 );
            end
        end
        coil_all = sum( L1_norm, 2 );
        weights = L1_norm./repmat( coil_all, [1,NC] );
        shift_new = sum( shift_info.*weights, 2 );

        shift_actual = shift_new;
        fun = @(x,xdata)x(1)*sin(xdata).^2+x(2)*cos(xdata).^2+x(3)*sin(xdata).*cos(xdata);
        shift_actual = shift_actual(:,1);

        [NM2, NR2, NC2] = size( kdata_phase_reorder );
        theta = 0:NR2/2-1;
    %     theta_all = zeros( length(theta)*2, 1 );
    %     theta_all(1:2:NR2) = theta*traj_corr_angle;
    %     theta_all(2:2:NR2) = theta*traj_corr_angle + pi;
        theta_all = option.PC_angles;

        shift_actual_all = zeros( length(theta)*2,1);
        shift_actual_all(1:2:NR2) = shift_actual;
        shift_actual_all(2:2:NR2) = shift_actual;
        theta_all = mod( theta_all, 2*pi );

        x0 = [0,0,0];
        xdata = theta_all(:);
        ydata = shift_actual_all;

        x = lsqcurvefit( fun, x0, xdata, ydata );
        [S,I] = sort(xdata);

        %Huili
%         test = x(1).*sin(xdata(I)).^2+x(2).*cos(xdata(I)).^2+x(3).*sin(xdata(I)).*cos(xdata(I));
%         subplot(2,2,4); plot(xdata(I),test,'LineWidth',3,'color','k')
%         hold on
%         plot(xdata(I),ydata(I),'LineWidth',2,'color','r')
        Results.PC_coeff = x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif strcmp( option.PhaseCorrection, 'UseTC' )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct phase errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [NM,NR,NI,NC] = size( kSpaceData );
        kdata_ref_cor = kSpaceData;
        x = option.PC_coeff;

        counter_main =1:1:NR;
        % nt_act = lnt;
        nt_act = NI;

        % need angle here
        theta_radian = option.Angles;
        traj_corr_values = x(1).*sin(theta_radian).^2+x(2).*cos(theta_radian).^2+x(3).*sin(theta_radian).*cos(theta_radian);
        traj_corr_values_reshape = reshape( traj_corr_values, [1,NR,NI]);
        traj_corr_values_reshape_rep = repmat(traj_corr_values_reshape, [NM,1,1,NC]);

    %     lin = transpose( -NM/2:NM/2-1 );
        lin = transpose( -(NM/2-1):NM/2 );
        lin_rep = repmat( lin, [1,NR,NI,NC]);
        shift_max = traj_corr_values_reshape_rep./NM;
        traj_corr_values_rep = exp( 2*1i*shift_max*pi.*lin_rep);
        inv_ft = fft( kdata_ref_cor, [], 1 );
        inv_ft_aft = inv_ft.*fftshift(traj_corr_values_rep,1);
        Results.kSpace = ifft(inv_ft_aft);
% % %             figure;     imagesc( abs(kdata_ref_traj(:,:,1,1,1)) )
%             figure;     imagesc( abs(ifft(inv_ft_aft(:,:,1,1,1))) )
    end

end
