function [ mag,phase ] = readdoris(path,nbc,nbl)
fid = fopen(path,'r');
if fid~=-1
    cint = (fread(fid,[nbc*2 nbl],'float32'))';
    fclose(fid);
    realpart = cint(:,1:2:nbc*2);
    cplxpart = cint(:,2:2:nbc*2);
    cint = realpart + sqrt(-1).*cplxpart;
    phase = angle(cint);
    mag=sqrt(realpart.^2+cplxpart.^2);
else
    mag=[];phase=[];
end 