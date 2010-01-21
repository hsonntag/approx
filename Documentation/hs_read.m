clc;
%[header] = argos_header_read('D:/hauptseminar/data/q00001708a.dat');

[data] = argos_data_read_start('D:/hauptseminar/data/q00001708.bin', 1, 1024);
for i=2:64
data = [data, argos_data_read_start('D:/hauptseminar/data/q00001708.bin', i*1024+1, 1024)];
end
data = transpose(data);
str = sprintf('D:/q00001708a.bin');
fd = fopen(str, 'wb');
fwrite(fd, data, 'int32');
fclose(fd);