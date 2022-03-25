function fWN = genFreq(t) %t in ps;
dnu = 1/(max(t)-min(t));                % freq in THz
numax = (1/abs(t(2)-t(1)))/2;           % max freq in THz
fTHz = [-numax:dnu:-dnu 0:dnu:numax];   % freqs in THz
fWN = fTHz*33.356;                      %frequency in wavenumbers