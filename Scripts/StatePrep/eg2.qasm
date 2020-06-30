//QASM netlist compiled by ../../qc
//Author: Nandakishore Santhi <nsanthi@lanl.gov>
//Date: 06-12-2017 10:52:58

OPENQASM 2.0;

include "qelib1.inc";

qreg q[2];
creg c[2];

//Quantum netlist implementing unitary operator 'eg2.state' using simple gates:
x q[0];
x q[1];
cx q[1], q[0];
u2(3.14159, -3.14159) q[0];
cx q[1], q[0];
h q[0];
x q[1];

measure q -> c;
