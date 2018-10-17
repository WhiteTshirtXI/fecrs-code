function s = calc_arclength(x)

num_segments = (length(x)-2)/2;

xmid = (x(:,1:2:end-3) + x(:,3:2:end-1))/2;

for i=2:num_segments
   s(i) =norm(xmid(:,i)-xmid(:,i-1)); 
end

s = cumsum(s);

s = s(end);

end