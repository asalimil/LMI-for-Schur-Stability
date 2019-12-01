%% LMI for Schur Stabilization
clc; clear; close all;

A = [0  0   1  0
     0  0   0  1
    -0.5  0  -2.5  1
     0 -0.5   1 -2];

B = [0 0
     0 0
     1 0
     0 1];

gamma = sdpvar(1);
K = sdpvar(size(B,2),size(A,1));

C = [];
eta = 0.00001;
M = [-gamma*eye(size(A,1))    (A + B*K)
      (A + B*K)'      -gamma*eye(size(A,1))];
F = [C, M <= eta*eye(size(M))];
optimize(F,gamma);
gamma = value(gamma)
K = value(K)
