/*
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
*/
#include <iomanip>
#include "qutils.h"
#include "qgate.h"
#include "qasm.h"

void linkMainFiles(const vector<string>& mainFilename, ofstream& outFile, uint64_t start, uint64_t stop) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    ifstream mainFile;
    for (int64_t i=start; i<stop; i++) {
        mainFile.open(mainFilename[i]);
        if (mainFile.is_open()) {
            cout << "Linking QASM file: " + mainFilename[i] + "\n";

            outFile << "\n//QASM file " + mainFilename[i] + " linked below:\n";
            string line;
            while (getline(mainFile, line)) outFile << line << '\n';
            mainFile.close();
        }
        else { cerr << "Could not open file " + mainFilename[i] + " for linking\n"; exit(1); }
    }
}

void writeQASMheader(const string& outFilename, ofstream& outFile, const string& execName) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    cout << "Writing QASM output to file: " + outFilename + "\n\n";
    outFile << "//QASM netlist compiled by " << execName << "\n";
    outFile << "//Author: Nandakishore Santhi <nsanthi@lanl.gov>\n";
    outFile << "//Date: " << getDateTime() << "\n\n";
    outFile << "OPENQASM 2.0;\n\n";
    outFile << "include \"qelib1.inc\";\n";
}

void writeQASMnetlist(QNet_t& Net, ofstream& outFile) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    cout << "Linking QASM code compiled from given operator" << endl;

    if (!useInline) {
        outFile << "gate operator ";
        for (uint8_t i=0; i<Net.QubitSize; i++) {
            if (i) outFile << ",";
            outFile << signal("q", i);
        }
        outFile << " {\n";
    }

    vector<uint64_t> Counter(Net.Gates.size(), 0); //One counter for each qubit timeline
    while (true) {
        int16_t line=-1; //Dominating qubit-line, ie., least serial# gate among all qubit-lines
        for (uint8_t i=0; i<Net.Gates.size(); i++) { //Get the first gate from among all qubit timelines
            if (Counter[i] >= Net.Gates[i].size()) continue;
            if (line >= 0) if (Net.Gates[i][Counter[i]].serial >= Net.Gates[line][Counter[line]].serial) continue;
            line = i;
        }
        if (line < 0) break;
        QG_t& gate = Net.Gates[line][Counter[line]];
        Counter[line]++;

#ifdef VERY_VERBOSE
        if (Verbose) {
            cout << "Line: " << int(line) << ", ";
            printGate(gate);
        }
#endif
        if (gate.type == Gctl) continue; //This is a control-tap place-holder, so ignore

        if (gate.arg.size() == 2) { //Only controlled gate at this stage is a CNOT gate
            outFile << indent() << "cx " << signal("q", gate.arg[0]) << ", " << signal("q", gate.arg[1]) << ";\n";
            gate.type = Gcx;
        }
        else {
            if (l2norm(gate.W - X) < QC_TOL) { //Use common X gate notation
                outFile << indent() << "x " << signal("q", gate.arg[0]) << ";\n";
                gate.type = Gx;
            }
            else if (l2norm(gate.W - H) < QC_TOL) { //Use common H gate notation
                outFile << indent() << "h " << signal("q", gate.arg[0]) << ";\n";
                gate.type = Gh;
            }
            else { //Determine if we need u1, u2 or u3
                double alpha, beta, gamma, delta;
                //u3(gamma, beta, delta) and phase=alpha (we will ignore the phase); see pages 174,176 of N/C
                //u1=U(0, 0, alpha)
                param4simple(gate.W, alpha, beta, gamma, delta);
                bool isBeta0=(abs(beta) < QC_TOL), isGamma0=(abs(gamma) < QC_TOL), isGammaPi2=(abs(gamma-(M_PI/2.0)) < QC_TOL), isDelta0=(abs(delta) < QC_TOL);

                if (isGamma0 && isBeta0 && !isDelta0) {
                    outFile << indent() << "u1(" << delta << ") " << signal("q", gate.arg[0]) << ";\n";
                    gate.type = Gu1;
                }
                else if (isGammaPi2 && !(isBeta0 && isDelta0)) {
                    outFile << indent() << "u2(" << rnd0(beta) << ", " << rnd0(delta) << ") " << signal("q", gate.arg[0]) << ";\n";
                    gate.type = Gu2;
                }
                else if (!(isGamma0 && isBeta0 && isDelta0)) {
                    outFile << indent() << "u3(" << rnd0(gamma) << ", " << rnd0(beta) << ", " << rnd0(delta) << ") " << signal("q", gate.arg[0]) << ";\n"; //U3
                    gate.type = Gu3;
                }
                else {
                    gate.type = Gnull; //This gate is redundant
                    cerr << __FUNCTION__ << ": alpha ignored case: (" << signal("q", gate.arg[0]) << ": " << alpha << ", " << isBeta0 << ", " << isGamma0 << ", " << isGammaPi2 << ", " << isDelta0 << ")\n";
                    //exit(1);
                }
            }
        }
    }

    if (!useInline) outFile << "}" << endl;
}

void showStatistics(QNet_t& Net) { //Compute cumulative costs
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    vector<uint64_t> Counter(Net.Gates.size(), 0); //One counter for each qubit timeline
    vector<double> Cost(Net.Gates.size(), 0.0); //One dominating cost for each qubit timeline
    vector<uint64_t> gateCount(Gend, 0);
    while (true) {
        int16_t line=-1; //Dominating qubit-line, ie., least serial# gate among all qubit-lines
        for (uint8_t i=0; i<Net.Gates.size(); i++) { //Get the first gate from among all qubit timelines
            if (Counter[i] >= Net.Gates[i].size()) continue;
            if (line >= 0) if (Net.Gates[i][Counter[i]].serial >= Net.Gates[line][Counter[line]].serial) continue;
            line = i;
        }
        if (line < 0) break;
        QG_t& gate = Net.Gates[line][Counter[line]];
        Counter[line]++;

        gateCount[gate.type]++;
        if (gate.type != Gnull) //We consider Gctl also for cost calculations
            Cost[line] = max(Cost[line], Cost[gate.arg[0]]) + gateCost[gate.type]; //Follow the critical path
    }
    double execTime=Cost[0];
    for (uint8_t i=1; i<Cost.size(); i++)
        if (Cost[i] > execTime) execTime = Cost[i];


    uint64_t numGates=0;
    for (uint8_t i=0; i<Net.Gates.size(); i++) numGates += Net.Gates[i].size();
    numGates -= (gateCount[Gnull] + gateCount[Gctl]);

    if (gateCount[Gnull] != 0) { cerr << "WARNING: NULL gates exist\n"; }
    if (gateCount[Gctl] != gateCount[Gcx]) { cerr << "WARNING: #CTAP does not match #CX gates\n"; }

    cout << "\n**********\n";
    cout << "STATISTICS\n";
    cout << "**********\n";
    cout << "\nGenerated unitary operator netlist with " << numGates << " simple gates\n";

    //For timing the following model is used for IBMQXi machines
    //(FC=FrameChange; GD=GaussianDerivative, GF=GaussianFlattop):
    //FC is virtual, GD typically 80ns, GF [100-400]ns depending on
    //which physical qubit located where on a given physical machine
    //U1: 1FC, U2: 2FC+GD, U3: 3FC+2GD, CX: 1FC+2GD+2GF
    uint64_t fieldWidth = 10;
    cout << "\nTiming considering pipelined/parallel gates, ignoring linked external code, and post optimizations using approximate model (ns):\n";

    cout << "\n" << setfill('-') << setw(fieldWidth*4) << "-\n" << setfill(' ');
    cout << setw(fieldWidth) << "FC" << setw(fieldWidth) << "GD" << setw(fieldWidth) << "GF";
    cout << "\n" << setfill('-') << setw(fieldWidth*4) << "-\n" << setfill(' ');
    cout << setw(fieldWidth) << FC << setw(fieldWidth) << GD << setw(fieldWidth) << GF;
    cout << "\n" << setfill('-') << setw(fieldWidth*4) << "-\n" << setfill(' ') << "\n";

    cout << "\tFC: Frame Changes\n"
         << "\tGD: Gaussian Derivatives\n"
         << "\tGF: Gaussian Flattops\n"
         << "\nExecution time estimate based on most critical (longest) computation path:\n\t" << (execTime/1000.0) << " (us)\n\n"
         << "NOTE: Coherence time is expected to be of the order of 90us for IBM 20-qubit machines\n"
         << "      When the time estimate becomes larger than coherence time of machine, dominant states get impossible to identify\n\n";

    cout << "Number of qubit-lines used in this Quantum circuit:\n";
    cout << int(Net.QubitSize) << ": {";
    for (uint8_t i=0; i<Net.Gates.size(); i++) if (Counter[i] > 0) cout << (i ? ", " : " ") << signal("q", i);
    cout << " }\n\n";

    cout << "Number of primitive gate types used in this Quantum circuit:\n";
    cout << "\n" << setfill('-') << setw(fieldWidth*(Gend+1-Gx)) << "-\n" << setfill(' ');
    for (uint8_t i=Gx; i<Gend; i++) cout << setw(fieldWidth) << gateName[i];
    cout << "\n" << setfill('-') << setw(fieldWidth*(Gend+1-Gx)) << "-\n" << setfill(' ');
    for (uint8_t i=Gx; i<Gend; i++) cout << setw(fieldWidth) << gateCount[i];
    cout << "\n" << setfill('-') << setw(fieldWidth*(Gend+1-Gx)) << "-\n" << setfill(' ') << "\n";
}
