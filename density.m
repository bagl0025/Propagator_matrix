clear all
model = 'karato.txt'
data = load (model);
n=size(data,1)
a=6371.0;

%load model values into arrays
for i = 1:1:n
    B(i) = data(i,2);
    z(i) = 6371.0 - data(i,1);
end

for i=1:1:n
    r=z(i);
	x = r / a;

%The following is based on PREM at 1 sec
% starting at the 400 km disc.

	if (r >= 5971.0 && r <= 6151.0)
		rho(i) = 7.1089 - 3.8045 * x;
    elseif (r >= 6151.0 && r <= 6291.0)
		rho(i) = 2.6910 + 0.6924 * x;
    elseif (r >= 6291.0 && r <= 6346.6)
		rho(i) = 2.6910 + 0.6924 * x;
    elseif (r >= 6346.6 && r <= 6356.0) 
    	rho(i) = 2.9;
    elseif (r >= 6356.0 && r <= 6368.0)
		rho(i) = 2.6;
    elseif (r >= 6368.0 && r <= 6371.0)
		rho(i) = 1.02;
    end
end

    output = fopen ('newmodel.txt','w');

	for i=1:1:n
        fprintf(output,'%.4f\t%.5f\t%.5f\n',6371.0-z(i),B(i),rho(i));
    end
    
	fclose (output);
