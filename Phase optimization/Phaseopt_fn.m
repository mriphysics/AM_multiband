function [ rfmb ] = Phaseopt_fn(rfsb,mb,tb,bs,AM_only)
%PHASEOPT_FN Summary of this function goes here
%   Detailed explanation goes here

N = length(rfsb);
t = 0:1/(N-1):1;

spos = (1:mb)-(mb+1)/2; 
phi_sel = 2*pi*tb*bs*t(:)*spos;

% Load in phase-offsets
if AM_only
    load('bmax_conj.mat')
else
    load('bmax_wong.mat')
end

phi_sol_PO = cell2mat(pstore(mb));
phi_sol_PO = angle(exp(1i*phi_sol_PO))+pi;

rfmb = sum( repmat(rfsb(:),[1 mb]).* exp(1i*phi_sel + repmat(1i*phi_sol_PO,[N 1])) ,2);

if AM_only
    rfmb = real(rfmb);
end

end

