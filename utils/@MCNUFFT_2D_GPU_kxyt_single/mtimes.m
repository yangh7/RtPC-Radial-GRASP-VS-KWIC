function ress = mtimes(a,bb)

 if a.adjoint,
     % Multicoil non-Cartesian k-space to Cartesian image domain
     % nufft for each coil and time point
        for ii = 1:size(bb,3)
            b = col( bb(:,:,ii) );
            ress(:,:,ii) = a.st'*b;
        end
 else
     % Cartesian image to multicoil non-Cartesian k-space 
        tmp_ress = a.st*bb;
        ress = reshape( tmp_res, [ a.dataSize(1), a.dataSize(2), size(tmp_ress,2) ] );
 end
            
