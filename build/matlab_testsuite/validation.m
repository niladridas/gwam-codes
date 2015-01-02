% Validation test suite
clc
clear all
OUT_valid = load('/home/prabhupad/Dropbox/writings/Mycodes/build/out.txt');
test_data = load('/home/prabhupad/Dropbox/writings/Mycodes/build/data/validation.txt');
figure(1)
plot(OUT_valid(:,1));
hold on ;
plot(test_data(:,4));
figure(2)
plot(OUT_valid(:,2));
hold on ;
plot(test_data(:,5));
figure(3)
plot(OUT_valid(:,3));
hold on ;
plot(test_data(:,6));
