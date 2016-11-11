function nPb= nrPb_fn(N,wp,bp)

% Function to return number of passbands in root-flipping algorithm
w = (-N/2+1:N/2)/N*2*pi;
idxPass = [];
for ii = 1:length(bp)
    idxPass = [idxPass find(w >= (bp(ii)-wp) & w <= (bp(ii)+wp))];
end
nPb=length(idxPass);