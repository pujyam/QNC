//QASM netlist compiled by ../../qc
//Author: Nandakishore Santhi <nsanthi@lanl.gov>
//Date: 08-12-2017 12:57:38

OPENQASM 2.0;

include "qelib1.inc";

//Quantum netlist implementing unitary operator 'F2.op' using simple gates:
x q[0];
x q[1];
cx q[0], q[1];
u3(0, -0.392699, -0.392699) q[1];
cx q[0], q[1];
u3(3.14159, -1.9635, 1.9635) q[1];
u3(0, -1.1781, -1.1781) q[0];
cx q[1], q[0];
u3(0.785398, 3.14159, -3.14159) q[0];
cx q[1], q[0];
u3(2.35619, -1.5708, -3.14159) q[0];
cx q[0], q[1];
u3(0.955317, -3.14159, 0) q[1];
cx q[0], q[1];
u3(0.955317, -3.14159, 0) q[1];
u3(3.14159, -1.1781, 1.1781) q[0];
cx q[1], q[0];
u3(0.659058, -3.14159, -1.89255) q[0];
cx q[1], q[0];
u3(2.48253, 2.03444, 3.14159) q[0];
u3(3.14159, 0, 0) q[1];
cx q[1], q[0];
u3(1.0472, 3.14159, -3.14159) q[0];
cx q[1], q[0];
u3(1.0472, 0, 0) q[0];
x q[1];
cx q[0], q[1];
u3(0.955317, 3.14159, -3.14159) q[1];
cx q[0], q[1];
u3(0.955317, 0, 0) q[1];
x q[0];
cx q[1], q[0];
u3(0.785398, 3.14159, -3.14159) q[0];
cx q[1], q[0];
u3(2.35619, 0, 3.14159) q[0];
