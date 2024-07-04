h1=1;ACC1=h1;
ACC2 = 100000; acc = 0.01;
tmin=12.38;
i=0;
while true
i=i+1;
h2 = (ACC1 + ACC2) / 2;
h12 = h2 - h1; %tao2-tao1
if h2-tmin < 0
    ACC1 = h2;
    if ACC2 - ACC1 <= acc
        break;
    end
else 
    ACC2 = h2;

 end

end
