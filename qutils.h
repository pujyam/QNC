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
#ifndef __QUTILS_H__
#define __QUTILS_H__

#include <iostream>
#include <fstream>
#include <ctime>
#include <complex>

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html
// We dont need/use BLAS or LAPACK, so only Armadillo headers are the prerequisite, which is included
#include <armadillo>

using namespace std;
using namespace arma;

enum GateType:uint8_t {Gnull, Gctl, Gx, Gu1, Gu2, Gu3, Gcx, Gh, Gend}; //Only define gates we use

extern const double QC_TOL; //Use a small enough value to denote near 0 values (NOTE: very low values tend to give unstable results!)
extern const uint64_t NO_VAL; //Used as no-val in FloydWarshall algorithm. Chosen so that there is no overflow when multiplied by 2
extern const char* gateName[Gend];
extern const double FC, GD, GF;
extern const double gateCost[Gend];

typedef struct QG_t {
    cx_mat W;
    vector<uint64_t> arg;
    uint64_t serial;
    GateType type;
} QG_t;

typedef struct QNet_t {
    bool useTop;
    Mat<uint64_t> topology; //Adjacency matrix for the original directed machine-qubit-connectivity-graph
    vector<vector<vector<uint64_t>>> Paths; //To store the optimal paths from phys-qubit source 'u' to destination 'v'

    bool useMap;
    vector<uint64_t> physMap; //To store the mapping between the operator's logical qubits and the machine's physical qubits

    vector<vector<QG_t>> Gates; //One vector for each qubit timeline
    int64_t nxtSerial; //Next gate's serial number in current netlist
    uint8_t QubitSize; //Size of qreg needed for implementation (computed after taking physMap and topology into account)
} QNet_t;

// GLOBAL VARIABLES
extern bool Verbose, useAlternateRecursion, useInline, useNaturalPivot, doStatePrep, doOpt;
extern uint64_t indirectPenalty;
extern cx_mat X, I2, H; //Store often used 1-qubit gates: 2x2 Pauli, 2x2 Identity, and 2x2 Hadamard

inline double l2norm(const cx_mat& M) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    return abs(cdot(M, M));
}

inline complex<double> diagProd(const cx_mat& M) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    complex<double> d=1.0;
    for (uint64_t i=0; i<M.n_rows; i++) d *= M(i, i);
    return d;
}

inline bool isDiagIdentity(const cx_mat& M) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    double d=0.0;
    for (uint64_t i=0; i<M.n_rows; i++) d += abs(M(i, i) - 1.0);
    return (d < QC_TOL);
}

inline complex<double> det2x2(const cx_mat& M) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if ((M.n_rows != 2) || (M.n_cols != 2)) { cerr << __FUNCTION__ << ": Matrix is not 2x2\n"; exit(1); }
    return (M(0,0)*M(1,1) - M(1,0)*M(0,1));
}

inline cx_mat sqrt2x2(const cx_mat& M) { //Compute a (not-unique) square-root of a given complex 2x2 matrix
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if ((M.n_rows != 2) || (M.n_cols != 2)) { cerr << __FUNCTION__ << ": Matrix is not 2x2\n"; exit(1); }

    cx_mat R(2, 2);
    if ((M(0, 1) == 0.0) || (M(1, 0) == 0.0)) { //Some off-diagonal is 0.0
        complex<double> lambda = sqrt(M(0, 0)) + sqrt(M(1, 1));
        if (lambda == 0.0) lambda = sqrt(M(0, 0)) - sqrt(M(1, 1));
        //lambda is now guaranteed to be non-zero, as det(M) != 0.0 for unitary Givens matrices
        R << sqrt(M(0, 0)) << (M(0, 1)/lambda) << endr << (M(1, 0)/lambda) << sqrt(M(1, 1));
    }
    else { //Some off-diagonal element is guaranteed non-zero
        complex<double> s = sqrt(det2x2(M));
        complex<double> t = sqrt(M(0,0) + M(1,1) + 2.0*s);
        //NOTE: By Cayley-Hamilton theorem, A = \tau.R - \delta.I
        //(since a matrix satisfies its own characteristic equation)
        //where, R = sqrt(A) and \tau := sqrt(t), is the trace of R (only if non-zero)
        //So, if A != k.I (multiple of I), then, t != 0, guaranteed!!

        R << M(0, 0) + s << M(0, 1) << endr << M(1, 0) << M(1, 1) + s;
        R /= t; //Since t != 0 as at least one off-diagonal element of A is non-zero
    }
    if (l2norm(R*R - M) > QC_TOL) { cerr << "Matrix square-root algorithm failed for matrix:\n" << M << "\n"; exit(1); }
    if (l2norm(R*R.t() - I2) > QC_TOL) { cerr << "Matrix square-root not unitary for matrix:\n" << M << "\n"; exit(1); }

    return R;
}

inline double rnd0(double x) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    return ((abs(x) < QC_TOL) ? 0.0 : x);
}

inline bool fileExists(const string& name) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    ifstream f(name);
    return f.good();
}

inline string getDateTime(void) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y %I:%M:%S",timeinfo);
    return string(buffer);
}

inline string indent(void) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    return (useInline) ? "" : "  ";
}

inline string signal(const string& name, uint64_t idx) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    return (useInline) ? (name + "[" + to_string(idx) + "]") : (name + to_string(idx));
}

inline void usageExit(const string& execName) {
    cerr << "\n************************************************\n"
            "        Quantum Netlist Compiler v1.0a\n"
            "  Build:   " __DATE__ " " __TIME__ "\n"
            " Author: Nandakishore Santhi <nsanthi@lanl.gov>\n"
            "************************************************\n\n"
            "Purpose: Given a quantum unitary operator, generate a topology constrained quantum gate netlist,\n"
            "         targeting a specific machine architecture. Currently, the only quantum gates that this\n"
            "         compiler makes use of are: u1, u2, u3, and cx (plus derived gates x, and H).\n"
            "         The compiler tries to get the number of quantum gates to be roughly linear in the size\n"
            "         of the operator notation used, without making use of ancilla-qubits.\n\n"
            "Usage: " << execName << " -oq out.qasm [-lq file{0..m}.qasm] -i[U|P|S] operator.op [-lq file{(m+1)..N}.qasm] [-nat] [-alt] [-I] [-v] [-map M] [-top T [-tW w]]\n"
            "\t-iU:  Reads a Unitary Matrix Operator from given file.\n"
            "\t       The input file format is a space seperated list of complex or real numbers, a matrix-row per line.\n"
            "\t       Complex numbers are represented as a tuple within brackets: eg: (1.0, 2) stands for 1.0 + 2.0i\n"
            "\t       Unitary operators have a worst case gate complexity ~O(4^n), where n=#qubits.\n"
            "\t       (one of -iU, -iP, or iS is mandatory)\n"
            "\t-iP:  Reads a Permutation Operator from the given file (in Cauchy's single line notation).\n"
            "\t       The input file format is the single line Cauchy permutation operator form;\n"
            "\t       0 based-index-integers seperated by spaces.\n"
            "\t       Permutation operators have a worst case gate complexity of ~O(n^2.2^n), where n=#qubits.\n"
            "\t       This complexity is often better than for general Unitary operators.\n"
            "\t       (one of -iU, -iP, or iS is mandatory)\n"
            "\t-iS:  Reads a initial state column vector from the given file.\n"
            "\t       In this mode, the compiler produces a quantum circuit which performs a 'state preparation'\n"
            "\t       to get the initial state starting from the default machine qubit-state |0>.\n"
            "\t       (one of -iU, -iP, or iS is mandatory)\n"
            "\t-oq:  Output filename to write an OPENQASM 2.0 netlist (mandatory).\n"
            "\t-lq:  Reads any number of optional OPENQASM-2.0 program files to link with the\n"
            "\t       generated QASM function, in the given order, either before or after the operator,\n"
            "\t       which is inserted where -iP|U is. (optional)\n"
            "\t-map: Uses the given permutation map M to map logical qubits to physical qubits. (optional)\n"
            "\t-top: Constrains the controlled-gates to given physical qubit-topology. Default is an ideal clique topology.\n"
            "\t       (T=adj-matrix, optional)\n"
            "\t-tW:  Penalty factor used in topology based qubit routing for indirect edges. Due to the combinatorial\n"
            "\t       algorithm used for routing, small changes in 'w' can dramatically change gate counts.\n"
            "\t       (w=unsigned-integer, optional, default is 3)\n"
            "\t-alt: Uses an alternate recursive structure to breakdown C^n-U gates. (optional)\n"
            "\t-nat: Use natural pivoting instead of defaul Gray-code based pivoting for Givens decompositions. (optional)\n"
            "\t       This has any effect only when either iU or iS is specified. Permutation operators always use natural encoding.\n"
            "\t-I:   Uses inline operator. Default is to use a function/procedure for the operator gate. (optional)\n"
            "\t-O:   Perform optimizations by merging and removing gates whenever possible. (optional)\n"
            "\t-v:   Sets verbose mode for printing lots of debug information. (optional)\n\n"
            "Copyright: Open-source BSD-3 license (LANL LACC# C18030 - QNC: Quantum Netlist Compiler).\n"
            "           If you found this useful, please cite this work, and also the original author.\n"
            "           This is the only means by which the author gets credit.\n"
            "           Again, please remember to cite this work if you find the code helpful in any way.\n\n";
    exit(0);
}

#endif //__QUTILS_H__
