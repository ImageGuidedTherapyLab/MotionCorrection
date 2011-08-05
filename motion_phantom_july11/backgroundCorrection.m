% Script for iterative solving of laplace equation
% load data
 
 
% Get polygon ROI containing hot data
% (can be replaced with any method for excluding heating, noise, or
% artifact.)
imagesc(abs(dat0(:,:,40)))
roi0=roipoly;
dat3=dat0;
% Phase One: Laplace Integration across entire image:
%This iteration makes a low resolution background estimate including
%heating zone.  Essentially can be replaced with a low pass (Wiener) if
%desired.  Note this blurs heating and penalty is larger ROI.
% However, gives smother/better data for BC calcs (noise sucks).
% iterate through slices (right now testing 2)
for zz=[2 40]
for niter=1:3
    datlap=0;
    for xx=[-1 0 1],
        for yy=[-1 0 1],
            if(xx~=0 || yy~=0),
                datlap=datlap+circshift(dat3(:,:,zz),[xx yy]);
            end
        end
    end
    % add mean of neighbors back in
    dat3(:,:,zz)=datlap/8;
end
end
% Alternate approach ... just remove comment (long computation time)
%dat3=wiener3D(dat0,[3 3]);
%Phase 2: Interative Laplace to estimate background in zone of heating
% remove the hot roi from the data
% consider replacing 0's with estimate of solution
% as either mean of pre-heating region or FFT technique
dat3=roi3D(dat3,1-roi0);
% iterate through slices (right now testing 2)
for zz=[2 40]
% iterate average values until converge
% currently just being lazy and hardcoding
for niter=1:100
    datlap=0;
    for xx=[-1 0 1],
        for yy=[-1 0 1],
            if(xx~=0 || yy~=0),
                datlap=datlap+circshift(dat3(:,:,zz),[xx yy]);
            end
        end
    end
    % add computed grid back in while enforcing Dirichlet BC's
    dat3(:,:,zz)=roi3d(dat3(:,:,zz),1-roi0)+roi3d(datlap/8,roi0);
end
end
% checkin' stuff out
% technically this can be converted to temperature, but subtraction likely
% better.
imagesc(cat(2,angle(dat0(:,:,40)),angle(dat3(:,:,40)),angle(conj(dat3(:,:,40)).*dat0(:,:,40))))
axis image, axis off, caxis([-pi pi]), colorbar
datcorr=conj(dat3).*dat0;
imagesc(cat(2,angle(conj(dat0(:,:,2)).*dat0(:,:,40)),angle(conj(dat3(:,:,40)).*dat0(:,:,40)),angle(conj(datcorr(:,:,2)).*datcorr(:,:,40))))
axis image, axis off, caxis([-pi pi]), colorbar
