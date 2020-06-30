#!/usr/bin/env python3.5
##############################################################################
# © Copyright 2017-. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
#
# All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
#
# This is open source software; you can redistribute it and/or modify it under the terms of the BSD 3-clause License. If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL. Full text of the BSD 3-clause License can be found in the License file in the main development branch of the repository.
#
##############################################################################
# BSD 3-clause license:
# Copyright 2017- Triad National Security, LLC
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
##############################################################################
# Author: Nandakishore Santhi
# Date: 28 November, 2017
# Copyright: Open source, must acknowledge original author
# Purpose: Compile a Quantum Computer Algorithm (as a unitary operator) into a QASM code that can run directly
# on IBM-QX physical machines or QISKit simulator.
# LACC#: LANL LACC# C18030 - QNC: Quantum Netlist Compiler
#
##############################################################################
#
#Usage:
# ./apiRemoteExperiment.py ../Scripts/BV/out.qasm
#
import pprint, sys
pp = pprint.PrettyPrinter(indent=2)

QX_TOKEN = "" #TODO: Add your IBM web token here

from IBMQuantumExperience import IBMQuantumExperience

api = IBMQuantumExperience(QX_TOKEN, config = {"url":"https://quantumexperience.ng.bluemix.net/api"}, verify=True)

print("\nMy Credits:\n")
pp.pprint(api.get_my_credits())

print("\nAvailable Backends:\n")
pp.pprint(api.available_backends())

backend = 'ibmqx_qasm_simulator'
#backend = 'ibmqx4'
#backend = 'ibmqx5'

print("\nBackend Status (" + backend + "):\n")
pp.pprint(api.backend_status(backend))

print("\nBackend Calibration (" + backend + "):\n")
pp.pprint(api.backend_calibration(backend))

print("\nBackend Parameters (" + backend + "):\n")
pp.pprint(api.backend_parameters(backend))

#Single experiment
with open(sys.argv[1], 'r') as myfile: qasm = myfile.read()

print("\nRun experiment on (" + backend + "):\n")
print("\nQASM:\n" + qasm)
pp.pprint(api.run_experiment(qasm, backend, shots=1024, name="Example 0", timeout=60))

#pp.pprint(api.get_execution("e0207dfd02770b59b3e08a0cd71fe23e"))
#pp.pprint(api.get_result_from_execution("e0207dfd02770b59b3e08a0cd71fe23e"))
