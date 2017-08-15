function x = correct_cos_stim_by_same_fun(T,x,Fs,Fc,odd)

inter=[-0.1,2];
n=2;
segs=inter(1):((inter(2)-inter(1))/(n-1)):inter(2);

for i=1:n-1
    ind=find(T>segs(i) & T<=segs(i+1));
    seg_x=x(ind);
    xc = xcorr(x,seg_x);
    xc = xc(8193:end);
    p=[1:(min(ind)-length(ind)),(max(ind)+10*length(ind)):(length(xc)-length(ind))];
    [m,j] = max(xc(p));
    %ind=ind;
    x(ind)=x((p(j)+1):(p(j)+length(ind))); 
end


