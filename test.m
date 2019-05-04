clear all; close all; clc;

A = rand([1,6]);
B = kron(A,ones([3,1])).';

C = (1:10).';
%C = kron(C,ones([6,1]));

%C(:,1:3)