function fy = gaussian(x,mu,var) 
fy = 1./sqrt(2*pi*var).*exp(-(x-mu).^2./2./var);
end
