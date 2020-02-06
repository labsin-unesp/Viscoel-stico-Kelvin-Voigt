
for i=1:100
    rms_vec(i) = rms(aceleracao(14,25*i:25*(i+1)));
end

plot(rms_vec)