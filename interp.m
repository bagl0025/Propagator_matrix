clear all

%load values (z, B, rho) from model
model = 'karato.txt'
data = load (model);
n=size(data,1);

%load model values into arrays
for i = 1:1:n
    z(i) = data(i,1);
    B(i) = data(i,2);
    rho(i) = data(i,3);
end

mindepth = round(z(1));
maxdepth = round(z(i));
j=4;
for i=mindepth+1:1:maxdepth-1
    if (i == z(j))
        i=i+1;
        j=j+2;
    end
    xi(i)=i;
end

newB=interp1(z,B,xi,'linear');
newrho=interp1(z,rho,xi,'linear');

output = fopen('tmp.txt','w');
for i=1:1:size(xi,2)
    fprintf(output,'%.2f\t%f\t%f\n',xi(i), newB(i), newrho(i));
end
fclose (output);