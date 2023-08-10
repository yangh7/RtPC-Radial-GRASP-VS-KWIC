function ress = mtimes(a,bb)

% % %  if a.adjoint,
% % %      % Multicoil non-Cartesian k-space to Cartesian image domain
% % %      % nufft for each coil and time point
% % %      b = bb(:,:,:);
% % %      ress(:,:,:) = a.st'*b;
% % %  else
% % %      % Cartesian image to multicoil non-Cartesian k-space 
% % %        ress =  a.st*bb;
% % %  end
            
 if a.adjoint,
     % Multicoil non-Cartesian k-space to Cartesian image domain
     % nufft for each coil and time point
     for tt = 1:size(bb,4)
        b = reshape( bb(:,:,:,tt), [a.dataSize(1)*a.dataSize(2)  a.imSize(3)] );
        ress(:,:,tt) = a.st(:,:,tt)'*b;
     end
 else
     % Cartesian image to multicoil non-Cartesian k-space 
     ress = zeros( a.dataSize(1), a.dataSize(2), a.imSize(3), a.dataSize(3));
     for tt = 1:size(bb,3)
         tmp_ress = a.st(:,:,tt)*bb(:,:,tt);
         ress(:,:,:,tt) = reshape( tmp_ress, [a.dataSize(1), a.dataSize(2), size(tmp_ress,2)] );
     end
 end
