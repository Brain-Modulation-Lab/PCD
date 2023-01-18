function y = my_gaussian(x,mu,fwhm,intensity)
y = intensity*exp(-((x-mu).^2)/(2*(fwhm)/2.355)^2);