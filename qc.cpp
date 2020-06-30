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
//Author: Nandakishore Santhi <nsanthi@lanl.gov>
//Copyright: Intended to be open-source. *Please* be sure to acknowledge author if you use this code
//Purpose:
//  Decompose any given Unitary matrix operator into an OPENQASM netlist of simple physical q-gates
//Date: 20 Nov, 2017
/*
 * POSSIBLE IMPROVEMENTS FOR FUTURE:
 * (1) One can reduce the gate complexity (upto a constant factor, which may matter for low coherence times and/or large netlists) at the expense of using ancilliary (work) qubits as an alternative to the successive sqrt-gate recursions to reduce control-signals per gate to 1
 * (2) Current code can sometimes generate netlists of O(n*4^n) gate complexity when given Unitary Operators; here, n := #qubits. Gray encoding is already being employed in the Givens decomposition stage, which means we can always achieve O(4^n) in the general case with some more optimization effort on the control-bit masks
 * (3) For approximate computation with bounded errors, we may truncate the successive sqrt recursion to smaller depths. In general, one could use slightly more sophisticated methods, including truncating infinitesimally precise gates from the netlist with little loss in fidelity of final measurements. A very rudimentary truncation is already performed by this compiler
 * (4) Logical qubit placement on physical qubit locations can be optimized automatically, based on machine's qubit connection topology
 * (5) Qubit routing generates gate network using 6n + O(1) time-slots. Mirroring and pipelining can reduce this to 2n + O(1) time-slots
 */
#include "qutils.h"
#include "qgate.h"
#include "qunitary.h"
#include "qperm.h"
#include "qasm.h"

// GLOBAL VARIABLES
const double QC_TOL=1e-5; //Use a small enough value to denote near 0 values (NOTE: very low values tend to give unstable results!)
const uint64_t NO_VAL=((1UL<<62)-1); //Used as no-val in FloydWarshall algorithm. Chosen so that there is no overflow when multiplied by 2
const char* gateName[Gend] = {"NULL", "CTAP", "x", "u1", "u2", "u3", "cx", "h"};
const double FC=7, GD=83+7, GF=400+7;
const double gateCost[Gend] = {0, FC+2*GD+2*GF, 3*FC+2*GD, FC, 2*FC+GD, 3*FC+2*GD, FC+2*GD+2*GF, 2*FC+GD};

bool Verbose=false, useAlternateRecursion=false, useInline=false, useNaturalPivot=false, doStatePrep=false, doOpt=false;
uint64_t indirectPenalty=3;
cx_mat X(2, 2), I2(2, 2), H(2, 2); //Store often used 1-qubit gates: 2x2 Pauli, 2x2 Identity, and 2x2 Hadamard

int main(int argc, char** argv) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    string inFilename="", outFilename="", mapFilename="", topoFilename="";
    vector<string> mainFilename;
    bool unitaryOp = true, linkMain = false;
    uint64_t opPos, mainCount=0;
    for (int i=1; i<argc; i++) { //Simple command line options parser
        if (string(argv[i]) == "-iP") { unitaryOp = false; inFilename = argv[++i]; opPos = mainCount; }
        else if (string(argv[i]) == "-iU") { unitaryOp = true; inFilename = argv[++i]; opPos = mainCount; }
        else if (string(argv[i]) == "-iS") { unitaryOp = doStatePrep = true; inFilename = argv[++i]; opPos = mainCount; }
        else if (string(argv[i]) == "-oq") { outFilename = argv[++i]; }
        else if (string(argv[i]) == "-lq") { linkMain = true; mainFilename.push_back(argv[++i]); mainCount++; }
        else if (string(argv[i]) == "-map") { mapFilename = argv[++i]; }
        else if (string(argv[i]) == "-top") { topoFilename = argv[++i]; }
        else if (string(argv[i]) == "-tW") { indirectPenalty = stoul(argv[++i], nullptr, 10); }
        else if (string(argv[i]) == "-alt") useAlternateRecursion = true;
        else if (string(argv[i]) == "-nat") useNaturalPivot = true;
        else if (string(argv[i]) == "-I") useInline = true;
        else if (string(argv[i]) == "-O") doOpt = true;
        else if (string(argv[i]) == "-v") Verbose = true;
        else usageExit(argv[0]);
    }

    if ((inFilename == "") || (outFilename == "")) usageExit(argv[0]); //Mandatory options missing
    if (!fileExists(inFilename)) { cerr << "Could not open input file " + inFilename + "\n"; exit(1); }
    if (linkMain && doStatePrep) { cerr << "Incompatible options combination (link and state-preparation)\n"; exit(1); }

    ofstream outFile(outFilename);
    if (!outFile.is_open()) { cerr << "Could not open output file " + outFilename + "\n"; exit(1); }

    if (Verbose) cout << "Using Armadillo version: " << arma_version::as_string() << endl;

    QNet_t Net;

    cout << "\n************************************************\n"
            "        Quantum Netlist Compiler v1.0a\n"
            " Author: Nandakishore Santhi <nsanthi@lanl.gov>\n"
            "************************************************\n\n";

    if (mapFilename != "") {
        if (!fileExists(mapFilename)) { cerr << "Could not open physical map file " + mapFilename + "\n"; exit(1); }
        readPermutation(mapFilename, Net.physMap); //Read the logical to physical qubits mapping
        Net.useMap = true;
    }
    else Net.useMap = false;

    cout << "\n";
    if (topoFilename != "") {
        if (!fileExists(topoFilename)) { cerr << "Could not open machine topology file " + topoFilename + "\n"; exit(1); }
        Net.topology.load(topoFilename, raw_ascii); //Read topology file

        if (Verbose) Net.topology.print("Quantum Machine Topology:");
        if (Net.topology.n_rows != Net.topology.n_cols) { cerr << "Adjacency matrix for topology graph should be square\n"; exit(1); }
        if (Net.useMap && (Net.physMap.size() < Net.topology.n_rows)) { cerr << "Provided logical to physical qubit mapping is incompatible with the machine topology\n"; exit(1); }

        cout << "Special connection topology for physical qubits on a machine given\n";
        cout << "Indirect edges penalty was set to: " << indirectPenalty << "\n";
        prepPaths(Net); //Pre-compute optimal paths for qubit routing on the given topology
    }
    else Net.useTop = false;
    Net.nxtSerial = 0;

    //Show some of the more important requested options for debuging/information
    if (Net.useMap) cout << "Special placement of logical qubits on the architecture requested\n";
    if (useAlternateRecursion) cout << "Alternate recursive structure for C^n-W gates will be used" << "\n";
    else cout << "Normal recursive structure for C^n-W gates will be used" << "\n";
    if ((unitaryOp && useNaturalPivot) || (!unitaryOp)) cout << "Natural pivoting scheme for qubit lines will be used" << "\n";
    else if (unitaryOp && !useNaturalPivot) cout << "Gray code pivoting scheme for qubit lines will be used" << "\n";
    if (doOpt) cout << "Gate level optimizations will be performed" << "\n";
    else cout << "Gate level optimizations will NOT be performed" << "\n";
    cout << "\n";

    I2 << 1 << 0 << endr << 0 << 1; //Initialize the 2x2 identity
    X << 0 << 1 << endr << 1 << 0; //Initialize the Pauli matrix for NOT gate
    H << 1/sqrt(2.0) << 1/sqrt(2.0) << endr << 1/sqrt(2.0) << -1/sqrt(2.0); //Initialize the 2x2 Hadamard matrix

    writeQASMheader(outFilename, outFile, argv[0]);

    if (linkMain) linkMainFiles(mainFilename, outFile, 0, opPos); //Link optional QASM files

    uint8_t numQubitsInOperator;
    if (unitaryOp) numQubitsInOperator = compileUnitaryOp(inFilename, Net, outFile);
    else numQubitsInOperator = compilePermutationOp(inFilename, Net, outFile);
    cout << "Number of qubits in the given operator: " << int(numQubitsInOperator) << "\n";

    lowerPrimitives(Net);

    writeQASMnetlist(Net, outFile);

    if (linkMain) linkMainFiles(mainFilename, outFile, opPos, mainFilename.size()); //Link optional QASM files
    outFile.close();

    showStatistics(Net);
}
