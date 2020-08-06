qreg q[4];
creg c[3];

//State preparation
x q[3];
h q[0];
h q[1];
h q[2];

operator q[0],q[1],q[2],q[3];

h q[0];
h q[1];
h q[2];

//Measurement
barrier q;
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
