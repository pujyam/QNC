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
#ifndef __QGATE_H__
#define __QGATE_H__

inline void printGate(const QG_t& gate) { //Human readable form for C-W or W gates
#ifdef VERY_VERBOSE
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
#endif
    cout << "GATE #" << gate.serial << " "
         << ((gate.type == Gnull) ? ((gate.arg.size()==2) ? "C-W" : "W") : gateName[gate.type])
         << " (";
    for (int i=0; i<gate.W.n_rows; i++)
        for (int j=0; j<gate.W.n_cols; j++)
            cout << " " << gate.W(i, j);
    cout << " ) [";
    for (int i=0; i<gate.arg.size(); i++)
        cout << " " << gate.arg[i];
    cout << " ]\n";
}

void getPath(const Mat<uint64_t> X, uint64_t u, uint64_t v, vector<uint64_t>& path);
void floydWarshall(const Mat<int64_t>& A, Mat<int64_t>& M, vector<vector<vector<uint64_t>>>& Paths);
void prepPaths(QNet_t& Net);

void lowerCnW(uint8_t S, const cx_mat& W, uint8_t target, uint64_t cMask, QNet_t& Net);

void addGate(QNet_t& Net, const cx_mat& W, const vector<uint64_t>& args);
void mapTopology(QNet_t& Net, const QG_t& gate);
void mergeGate(QNet_t& Net, const QG_t& gate);
void appendGate(QNet_t& Net, const QG_t& gate, bool replaceLast=false);

void param4simple(const cx_mat& W, double& alpha, double& beta, double& gamma, double& delta);
void C1Wsimple(const cx_mat& W, cx_mat& A, cx_mat& B, cx_mat& C, cx_mat& D);

void lowerPrimitives(QNet_t& Net);
void mergeDirectIndirectCNOT(QNet_t& Net, uint64_t c, uint64_t t);

#endif //__QGATE_H__
