//QASM netlist compiled by ../../qc
//Author: Nandakishore Santhi <nsanthi@lanl.gov>
//Date: 06-12-2017 10:55:27

OPENQASM 2.0;

include "qelib1.inc";

qreg q[2];
creg c[2];

//Quantum netlist implementing unitary operator 'eg1.state' using simple gates:
x q[1];
cx q[1], q[0];
u2(3.14159, -3.14159) q[0];
cx q[1], q[0];
u3(1.5708, 0, 0) q[0];
x q[1];
cx q[0], q[1];
u3(0.785398, 3.14159, -3.14159) q[1];
cx q[0], q[1];
u3(0.785398, 0, 0) q[1];
x q[0];
cx q[1], q[0];
u2(3.14159, -3.14159) q[0];
cx q[1], q[0];
h q[0];

measure q -> c;
