h q[0];
h q[1];

//Measurement
barrier q;
measure q[0] -> c[0];
measure q[1] -> c[1];
