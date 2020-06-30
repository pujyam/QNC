qreg q[3];
creg c[2];

//State preparation
x q[2];
h q[0];
h q[1];

operator q[0],q[1],q[2];

h q[0];
h q[1];

//Measurement
barrier q;
measure q[0] -> c[0];
measure q[1] -> c[1];
