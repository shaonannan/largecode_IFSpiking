%% orientation hypercolumn and orientation-cluster
clear all
% nxhyp = 32; nyhyp = 6;  nhyp = nxhyp*nyhyp;
nxhyp = 20; nyhyp = 6;  nhyp = nxhyp*nyhyp;
nxc   = 24; nyc   = 24;
% nxhyp = 8; nyhyp = 4;  nhyp = nxhyp*nyhyp;
% nxc   = 20; nyc   = 20;
nmax  = nxc*nyc*nhyp;
nori  = 12; oriint = pi/nori;
nphase = 4;

theta = zeros(1,nmax); clustorien = zeros(1,nmax);clusthyp = zeros(1,nmax);

% the 1st orientation hypercolumn
nxh = nxc/2; nyh = nyc/2;
for i = 1:1:nxc
    for j = 1:1:nyc
        ind = (j-1)*nxc*nxhyp+i;
        theta(ind) = atan2((nxh-i),(j-nyh))/2.0;
        if theta(ind) < 0
            theta(ind)= theta(ind)+pi;
        end
        clustorien(ind) = floor(theta(ind)/oriint);
        if clustorien(ind) == nori
            clustorien(ind) = nori - 1;
        end
    end
end

for i = 1:1:nxc
    for j = nyc+1:1:2*nyc
        ind = (j-1)*nxc*nxhyp+i;
        yref= 2*nyc-j+1;
        indref = (yref-1)*nxc*nxhyp + i;
        theta(ind) = theta(indref);
        clustorien(ind) = clustorien(indref);
    end
end

for i = nxc+1:1:2*nxc
    for j = 1:1:2*nyc
        ind = (j-1)*nxc*nxhyp+i;
        xref= 2*nxc-i+1;
        indref = (j-1)*nxc*nxhyp + xref;
        theta(ind) = theta(indref);
        clustorien(ind) = clustorien(indref);
    end
end

template_theta = zeros(2*nxc,2*nyc);   template_clustorien = zeros(2*nxc,2*nyc);
for i = 1:1:2*nxc
    tmp = reshape(theta,nxc*nxhyp,nyc*nyhyp);
    template_theta(i,:) = tmp(i,1:2*nyc); 
    tmp = reshape(clustorien,nxc*nxhyp,nyc*nyhyp);
    template_clustorien(i,:) = tmp(i,1:2*nyc); 
end

for i = 1:1:nxc*nxhyp
    for j = 1:1:nyc*nyhyp
        xhyp = ceil(i/nxc);
        yhyp = ceil(j/nyc);
        ind  = (j-1)*nxc*nxhyp + i;
        clusthyp(ind) = (yhyp-1)*nxhyp + xhyp;
    end
end

        

theta = reshape(repmat(template_theta,nxhyp/2,nyhyp/2),1,nmax);
clustorien = reshape(repmat(template_clustorien,nxhyp/2,nyhyp/2),1,nmax);
phase = floor(nphase*rand(1,nmax));
size(theta)
% figure(1);imagesc(reshape(theta,nxc*nxhyp,nyc*nyhyp));axis off;axis equal;
% figure(2);imagesc(reshape(clustorien,nxc*nxhyp,nyc*nyhyp));axis off;axis equal;
% figure(3);imagesc(reshape(clusthyp,nxc*nxhyp,nyc*nyhyp));axis off;axis equal;
% % figure(3);imagesc(template_theta);axis off;axis equal;
% figure(4);imagesc(template_clustorien);axis off;axis equal;
save('clustorien.txt','clustorien','-ASCII');
save('clusthyp.txt','clusthyp','-ASCII');
save('theta.txt','theta','-ASCII');
save('phase.txt','phase','-ASCII');

figure(12); figlarge = (reshape(theta,nxc*nxhyp,nyc*nyhyp));
clustlarge = (reshape(clustorien,nxc*nxhyp,nyc*nyhyp));
subplot(1,2,1);imagesc(figlarge(1:30,1:20));
axis off;axis equal;
subplot(1,2,2);imagesc(clustlarge(1:30,1:20));
axis off;axis equal;





