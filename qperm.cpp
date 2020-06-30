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
#include "qutils.h"
#include "qgate.h"
#include "qunitary.h"
#include "qperm.h"

void readPermutation(const string& filename, vector<uint64_t>& P) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    ifstream inFile(filename);

    uint64_t value;
    if (!inFile.is_open()) { cerr << "Could not open input file " + filename + "\n"; exit(1); }
    while (inFile >> value) P.push_back(value); //Read the elements in the file into a vector
    inFile.close();
    if (P.size() == 0) { cerr << "FORMAT ERROR: permutation could not be read: should be in 0-base Cauchy 1-line form\n"; exit(1); }

    if (Verbose) {
        cout << "Permutation map from file '" << filename << "': [ ";
        for (auto i = P.begin(); i != P.end(); i++) cout << *i << ' ';
        cout << "]\n" << endl;
    }
}

uint8_t compilePermutationOp(const string& inFilename, QNet_t& Net, ofstream& outFile) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    vector<uint64_t> P;
    readPermutation(inFilename, P);

    const uint64_t N = P.size();
    const uint64_t lgN = round(log2(N));
    uint64_t QubitSize = lgN;

    if (Net.useMap && (Net.physMap.size() < QubitSize)) { cerr << "Provided logical to physical qubit mapping is incompatible with the permutation operator\n"; exit(1); }
    if (Net.useMap && (Net.physMap.size() > QubitSize)) QubitSize = Net.physMap.size();
    if (Net.useTop && (Net.topology.n_rows < QubitSize)) { cerr << "Provided topology graph is incompatible with the permutation operator\n"; exit(1); }
    if (Net.useTop && (Net.topology.n_rows > QubitSize)) QubitSize = Net.topology.n_rows;

    vector<vector<uint64_t>> VSwaps;
    swapsDecompose(P, VSwaps); //Decompose P into a sequence of at most N swaps

    outFile << "\n//Quantum netlist implementing permutation operator '"
            << inFilename << "' using simple gates:\n";

    Net.Gates = vector<vector<QG_t>>(QubitSize); //One vector for each qubit timeline
    if (Verbose) cout << "\n________________________________\n";
    //\Pi_i G_i == P; so apply the G_i in same order to affect P (post-fix)
    for (int64_t i=0; i<VSwaps.size(); i++) lowerCnG(lgN, X, VSwaps[i], Net, false);
    if (Verbose) cout << "________________________________\n" << endl;

    return lgN;
}

void swapsDecompose(const vector<uint64_t>& P, vector<vector<uint64_t>>& VSwaps) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    //Decompose P into a sequence of at most N swaps
    vector<bool> done;
    for (uint64_t i=0; i<P.size(); i++) done.push_back(false);

    vector<vector<uint64_t>> cycles;
    for (uint64_t i=0; i<P.size(); i++) {
        if (!done[i]) {
            uint64_t j=i, start=i;
            vector<uint64_t> cycle;
            while (cycles.size() < P.size()) { //Find cycles in the permutation operator
                cycle.push_back(j);
                done[j] = true;

                j = P[j];
                if (j == start) break;
            }
            cycles.push_back(cycle);
        }
    }
    if (Verbose) {
        cout << "CYCLES: [ ";
        for (auto i = cycles.begin(); i != cycles.end(); i++) {
            cout << "(";
            for (auto j = (*i).begin(); j != (*i).end(); j++) cout << " " << (*j);
            cout << " ) ";
        }
        cout << "]\n" << endl;
    }

    for (auto i = cycles.begin(); i != cycles.end(); i++) {
        if ((*i).size() > 1) {
            uint64_t start=(*i)[0];
            for (uint64_t j=1; j<(*i).size(); j++) {
                vector<uint64_t> s;
                s.push_back(start);
                s.push_back((*i)[j]);
                VSwaps.push_back(s);
            }
        }
    }
    if (Verbose) {
        cout << "SWAPS: [ ";
        for (auto i=VSwaps.begin(); i!=VSwaps.end(); i++)
            cout << "(" << (*i)[0] << ", " << (*i)[1] << ") ";
        cout << "]\n" << endl;
    }
}
