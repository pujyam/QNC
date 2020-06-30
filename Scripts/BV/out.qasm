//QASM netlist compiled by ../../qc
//Author: Nandakishore Santhi <nsanthi@lanl.gov>
//Date: 30-06-2020 02:53:22

OPENQASM 2.0;

include "qelib1.inc";

//QASM file mainInline0.qasm linked below:
qreg q[5];
creg c[2];

//State preparation
x q[2];
h q[0];
h q[1];
h q[2];

//Quantum netlist implementing permutation operator 'P.op' using simple gates:
x q[1];
u3(0, 0.785398, 0.785398) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, 3.14159, -3.14159) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, -1.5708, 0) q[2];
u3(0, 0.392699, 0.392699) q[0];
cx q[1], q[0];
u3(0, -0.785398, -0.785398) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, 3.14159, -3.14159) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, 1.5708, 0) q[2];
u3(0, -0.392699, -0.392699) q[0];
cx q[1], q[0];
u3(0, 0.785398, 0.785398) q[2];
h q[2];
h q[1];
cx q[2], q[1];
h q[2];
h q[1];
u3(0.785398, 3.14159, -3.14159) q[2];
h q[2];
h q[1];
cx q[2], q[1];
h q[2];
h q[1];
u3(0.785398, -1.5708, 0) q[2];
u3(0, 0.392699, 0.392699) q[1];
x q[1];
u3(0, 0.785398, 0.785398) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, 3.14159, -3.14159) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, -1.5708, 0) q[2];
u3(0, 0.392699, 0.392699) q[0];
cx q[1], q[0];
u3(0, -0.785398, -0.785398) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, 3.14159, -3.14159) q[2];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
u3(0.785398, 1.5708, 0) q[2];
u3(0, -0.392699, -0.392699) q[0];
cx q[1], q[0];
u3(0, 0.785398, 0.785398) q[2];
h q[2];
h q[1];
cx q[2], q[1];
h q[2];
h q[1];
u3(0.785398, 3.14159, -3.14159) q[2];
h q[2];
h q[1];
cx q[2], q[1];
h q[2];
h q[1];
u3(0.785398, -1.5708, 0) q[2];
u3(0, 0.392699, 0.392699) q[1];

//QASM file mainInline1.qasm linked below:
h q[0];
h q[1];

//Measurement
barrier q;
measure q[0] -> c[0];
measure q[1] -> c[1];
