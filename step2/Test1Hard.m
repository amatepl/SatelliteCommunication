clear all;close all;clc;
% The parity-check matrix is of size N-K x N.
% The bit stream is devided into blocks of K bits (message vector).
% N : length of the codeword.
H =[1 1 0 1 1 0 0 1 0 0;
    0 1 1 0 1 1 1 0 0 0;
    0 0 0 1 0 0 0 1 1 1;
    1 1 0 0 0 1 1 0 1 0;
    0 0 1 0 0 1 0 1 0 1];

H1 = H;
% Transforming H into [IP']
H([3 5],:)=H([5 3],:);
H([4 5],:)=H([5 4],:);
H([1 5],:)=H([5 1],:);
H(5,:) = mod(H(5,:) + H(1,:),2);
H(1,:) = mod(H(1,:) + H(2,:),2);
H(1,:) = mod(H(1,:) + H(3,:),2);
H(2,:) = mod(H(2,:) + H(3,:),2);
H(5,:) = mod(H(5,:) + H(4,:),2);
H(1,:) = mod(H(1,:) + H(5,:),2);
H(2,:) = mod(H(2,:) + H(5,:),2);

% Parity array portion
P = H(:,6:end);

%Generator matrix (K x N-K)
G = [P' eye(5)];

bit_stream = randi([0,1],1,50);

% Dividing the bit stream into blocks of bits.
blocks = reshape(bit_stream,[],5);

encoded_bits = mod(blocks*G,2);

f1 = mod(sum(H(1,:)),2);
f2 = mod(sum(H(2,:)),2);
f3 = mod(sum(H(3,:)),2);
f4 = mod(sum(H(4,:)),2);
f5 = mod(sum(H(5,:)),2);

f_nodes = [f1 f2 f3 f4 f5];

% The ones in every column of H idicate the f node which is related to the
% given c. 
% Ex: the first column is : 1
%                           0
%                           0
%                           0
%                           1
%                           0
% So the f nodes that can check c1 are f1 and f4.
% Once we know this information we check all the c nodes that are connected
% to the given f nodes so we can deduce the value of the c node we want to  
% verify.