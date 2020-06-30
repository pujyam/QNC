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

inline cx_mat Ry(double x) { //pg174 N/C
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    cx_mat M(2, 2);
    M << cos(x/2) << -sin(x/2) << endr << sin(x/2) << cos(x/2);

    return M;
}

inline cx_mat Rz(double x) { //pg174 N/C
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    cx_mat M(2, 2);
    M << polar(1.0, -x/2) << 0 << endr << 0 << polar(1.0, x/2);

    return M;
}

inline cx_mat zD(double x) { //Controlled phase-shift-x is same as 1-qubit zD(x) on control-qubit (pg180 N/C)
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    cx_mat M(2, 2);
    M << 1 << 0 << endr << 0 << polar(1.0, x);

    return M;
}

void getPath(const Mat<uint64_t> X, uint64_t u, uint64_t v, vector<uint64_t>& path) { //Get a least cost path u --> v
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if (X(u, v) == NO_VAL) return;

    uint64_t w=u;
    path.push_back(w);
    while (w != v) {
        w = X(w, v);
        path.push_back(w);
        if (Verbose) cerr << " " << w;
    }
}

void floydWarshall(const Mat<int64_t>& A, Mat<int64_t>& M, vector<vector<vector<uint64_t>>>& Paths) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    uint64_t N=A.n_cols;
    M.set_size(N, N);
    Mat<uint64_t> X(N, N);
    for (uint64_t i=0; i<N; i++) {
        for (uint64_t j=0; j<N; j++) {
            X(i, j) = (A(i, j) == 0) ? NO_VAL : j;
            M(i, j) = ((i != j) && (A(i, j) == 0)) ? NO_VAL : A(i, j);
        }
    }

    for (uint64_t k=0; k<N; k++) {
        for (uint64_t i=0; i<N; i++) {
            for (uint64_t j=0; j<N; j++) {
                if (M(i, j) > (M(i, k) + M(k, j))) {
                    M(i, j) = M(i, k) + M(k, j);
                    X(i, j) = X(i, k);
                }
            }
        }
    }

    for (uint64_t i=0; i<N; i++) {
        Paths.push_back(vector<vector<uint64_t>>(N));
        for (uint64_t j=0; j<N; j++) {
            if (Verbose) cerr << "\t(" << i << " --> " << j << "):";
            getPath(X, i, j, Paths[i][j]);
            if (Verbose) cerr << endl;
        }
    }

    if (Verbose) {
        for (uint64_t i=0; i<N; i++) {
            cout << "\nMinimum Cost With Respect to Node:" << i << endl;
            for (uint64_t j=0; j<N; j++) {
                cout << "  (" << i << "->" << j << "): " << M(i, j) << " : [ ";
                for (uint64_t k=0; k<Paths[i][j].size(); k++) cout << Paths[i][j][k] << " ";
                cout << "]\n";
            }
        }
    }
}

void prepPaths(QNet_t& Net) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    Mat<uint64_t>& topology=Net.topology;
    if (topology.n_rows) { //Consider topology constraints
        Net.useTop = true;

        vector<vector<vector<uint64_t>>>& Paths=Net.Paths;
        Mat<int64_t> Costs;
        Mat<int64_t> A(topology.n_rows, topology.n_cols);
        for (uint64_t i=0; i<topology.n_rows; i++) {
            for (uint64_t j=0; j<topology.n_cols; j++) {
                if (topology(i, j) == 0) {
                    if (topology(j, i) != 0) A(i, j) = indirectPenalty*topology(j, i); //H, CNOT, H; so 3x cost is usually good
                    else A(i, j) = 0;
                }
                else A(i, j) = topology(i, j);
            }
        }
        if (Verbose) A.print("A:");

        cout << "Finding an optimal set of qubit routings within the given machine topology ...\n";
        floydWarshall(A, Costs, Paths);
    }
}

void lowerCnW(uint8_t S, const cx_mat& W, uint8_t target, uint64_t cMask, QNet_t& Net) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if (Verbose)
        cout << "Gate (" << W(0, 0) << ", " << W(0, 1) << "; " << W(1, 0) << ", " << W(1, 1)
             << ";) on qubit " << int(target) << " with control mask " << cMask << endl;

    if (!cMask) {
        if (Verbose) cout << "No control mask ...\n";
        addGate(Net, W, {target});
        return;
    }

    //Take the last active control bit from the cMask
    uint8_t b;
    for (b=0; b<S; b++) if (cMask & (1<<b)) break;
    uint64_t nxtCMask = cMask ^ (1<<b);

    if (!nxtCMask) addGate(Net, W, {b, target}); //cMask has only 1 active control qubit, so handle differently
    else { //Recursively reduce the control-set size
        //
        //===/n==*========     ==/n=====*=====*===*== [nxtCMask]
        //       |          =>          |     |   |
        //-------*--------     ------*--+--*--+------ [b]
        //       |          =>       |     |      |
        //-------W--------     ------R-----R'-----R-- [target]
        //
        cx_mat R = sqrt2x2(W);
        if (useAlternateRecursion) { //We could also do:
            lowerCnW(S, X, b, nxtCMask, Net); //C(nxtCMask) NOT
            lowerCnW(S, R.t(), target, (1<<b), Net); //C(b)-R'
            lowerCnW(S, X, b, nxtCMask, Net); //C(nxtCMask) NOT
            lowerCnW(S, R, target, (1<<b), Net); //C(b)-R
            lowerCnW(S, R, target, nxtCMask, Net); //C(nxtCMask)-R
        }
        else { //Do roughly according to N/C pg182
            lowerCnW(S, R, target, (1<<b), Net); //C(b)-R
            lowerCnW(S, X, b, nxtCMask, Net); //C(nxtCMask) NOT
            lowerCnW(S, R.t(), target, (1<<b), Net); //C(b)-R'
            lowerCnW(S, X, b, nxtCMask, Net); //C(nxtCMask) NOT
            lowerCnW(S, R, target, nxtCMask, Net); //C(nxtCMask)-R
        }
    }
}

void addGate(QNet_t& Net, const cx_mat& W, const vector<uint64_t>& args) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    vector<uint64_t> mappedArgs; //Map to physical, then to topology
    for (uint8_t i=0; i<args.size(); i++) mappedArgs.push_back(Net.physMap.size() ? Net.physMap[args[i]] : args[i]);

    //Add gate to structure Gates_t, while optimizing out adjacent gates to simpler ones
    mapTopology(Net, {W, mappedArgs, NO_VAL, Gnull}); //Consider topology constraints for CNOT and CU3 gates
}

void mapTopology(QNet_t& Net, const QG_t& gate) { //Consider topology constraints for CNOT and CU3 gates
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if ((!Net.useTop) || (gate.arg.size() == 1)) { //If not a restrictive topology or if a simple 1-qubit gate
        mergeGate(Net, gate);
        return;
    }

    uint64_t oc=gate.arg[0], t=gate.arg[1], c, l=Net.Paths[oc][t].size(); //control and target qubits
    if (Verbose) cout << "Route and Place qubit " << oc << " adjacent to target qubit " << t << "\n";

    if (l == 0) { cerr << "No path exists in given topology for qubit " << oc << " to control qubit " << t << "\n"; printGate(gate); exit(1); }
    if (l == 1) { cerr << "Qubit " << oc << " cannot control qubit " << t << " (same?)\n"; printGate(gate); exit(1); }

    if (l > 2) { //1-hop path does not exist between control and target, swap qubits to take 'oc' near 't'
        //Paths[c][t][0] == 'oc'; Paths[c][t][l-1] == 't'; Use swaps to take qubit 'oc' to 'Paths[c][t][l-2]'
        c = oc;
        for (uint64_t i=1; i<(l-1); i++) {
            uint64_t x = Net.Paths[oc][t][i];
            if (Verbose) cout << "Forward Swap q" << c << " q" << x << "\n";
            //Swap q[c] and q[x]
            mergeGate(Net, {X, {c, x}, NO_VAL, Gnull});
            mergeGate(Net, {X, {x, c}, NO_VAL, Gnull});
            mergeGate(Net, {X, {c, x}, NO_VAL, Gnull});
            c = x;
        }
    }
    else c = oc; //l == 2 case; where 'oc' is either directly or indirectly 1-hop from 't'

    //1-hop path exists between control and target, no need to swap any qubits
    mergeGate(Net, {gate.W, {c, t}, NO_VAL, Gnull});

    if (l > 2) { //1-hop path does not exist between control and target, swap qubits to take 'c' back to 'oc'
        //Paths[c][t][0] == 'oc'; Paths[c][t][l-1] == 't'; Use swaps to take qubit 'c' back to 'oc'
        for (uint64_t i=l-2; i>0; i--) {
            uint64_t x = Net.Paths[oc][t][i-1];
            if (Verbose) cout << "Reverse Swap q" << c << " q" << x << "\n";
            //Swap q[c] and q[x]
            mergeGate(Net, {X, {c, x}, NO_VAL, Gnull});
            mergeGate(Net, {X, {x, c}, NO_VAL, Gnull});
            mergeGate(Net, {X, {c, x}, NO_VAL, Gnull});
            c = x;
        }
    }
}

void mergeGate(QNet_t& Net, const QG_t& gate) { //Optimize and merge gate to netlist. TODO: Verify logic thoroughly
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;

    if (gate.type == Gctl) {
        cout << "Ignored a control-tap gate\n";
        return;
    }

    uint8_t s = gate.arg.size();
    if (s > 2) { cerr << "Gates should be 1-qubit with either 1 or 0 control qubits\n"; exit(1); }
    if (l2norm(I2 - gate.W) < QC_TOL) {
        if (Verbose) cout << "Ignoring Id gate\n";
        return;
    }

    if (!doOpt) { //Simply add the gate
        appendGate(Net, gate);
        return;
    }

    uint64_t t=gate.arg[s-1], c=gate.arg[0]; //target, control qubits
    if (Net.Gates[t].size() == 0) appendGate(Net, gate); //First qubit on this target qubit's timeline
    else { //See if we can merge with previous gate
        QG_t& gatePrv = Net.Gates[t].back();
        uint8_t sPrv = gatePrv.arg.size();

        bool canMerge = (sPrv == s) && (gatePrv.type != Gctl);
        if (canMerge) { //See if a merging is appropriate
            for(uint8_t i=0; i<s; i++) {
                if (gatePrv.arg[i] != gate.arg[i]) { //Arguments dont match?
                    canMerge = false;
                    break;
                }
            }
        }

        if (canMerge) {
            if (Net.Gates[c].size() == 0) canMerge = false; //Control qubit-line's first gate
            else if (Net.Gates[c].back().serial > gatePrv.serial) canMerge = false; //Control qubit has meanwhile transformed
        }

        if (canMerge) {
            QG_t newGate = gatePrv;
            newGate.W = gate.W*newGate.W; //Paying attention to the order of the operators
#ifdef VERY_VERBOSE
            if (Verbose) {
                cout << "gatePrv: " << gatePrv.W << "\n";
                cout << "gate: " << gate.W << "\n";
                cout << "newGate: " << newGate.W << "\n";
            }
#endif
            if (l2norm(I2 - newGate.W) > QC_TOL) {
                if (Verbose) cout << "Merging to last gate in line " << t << endl;
                appendGate(Net, newGate, true); //Do these gates cancel out?
            }
            else {
                if (Verbose) cout << "Removing last gate in line " << t << endl;
                Net.Gates[t].pop_back(); //Remove old gate
                if (s==2) Net.Gates[c].pop_back(); //Remove old CTAP
            }
        }
        else appendGate(Net, gate); //Merging not possible
    }
}

void appendGate(QNet_t& Net, const QG_t& gate, bool replaceLast) { //Append gate to netlist, and update gateCount, and totalTime
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    uint64_t serial;
    uint8_t s=gate.arg.size();
    uint64_t t=gate.arg[s-1], c=gate.arg[0]; //target, control qubits
    if (replaceLast) {
        serial = Net.Gates[t].back().serial; //Reuse old serial#
        Net.Gates[t].back() = {gate.W, gate.arg, serial, gate.type}; //Replace old gate
    }
    else {
        for(uint8_t i=0; i<gate.arg.size(); i++)
            if (Net.QubitSize < (gate.arg[i]+1)) Net.QubitSize = gate.arg[i]+1; //Update qreg size needed for implementation
        serial = Net.nxtSerial++; //New serial#
#ifdef VERY_VERBOSE
        if (Verbose) {
            cout << "#Qubit lines in Net.Gates: " << Net.Gates.size() << "\n";
            cout << "W:\n" << gate.W << "\n[ " << t << " " << serial << " " << gateName[gate.type] << " ]\n";
        }
#endif
        Net.Gates[t].push_back({gate.W, gate.arg, serial, gate.type});
        if (s==2) {
#ifdef VERY_VERBOSE
            if (Verbose) cout << "W:\n" << gate.W << "\n[ " << c << " " << serial << " " << gateName[gate.type] << " ]\n";
#endif
            Net.Gates[c].push_back({gate.W, gate.arg, serial, Gctl}); //Place-holder to indicate that we have tapped off this control-line
        }
    }
}

void param4simple(const cx_mat& W, double& alpha, double& beta, double& gamma, double& delta) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;

    const complex<double> j(0, 1);
    //This function is to be used to get the rotation angles from a 1-qubit operator for IBM OPENQASM
    //Get alpha, beta, gamma, delta to decompose unitary W (We follow notation on pg176 of Nielson/Chuang)
    //W00 = exp(i(alpha-(beta/2)-(delta/2))cos(gamma/2)
    //W01 = -exp(i(alpha-(beta/2)+(delta/2))sin(gamma/2))
    //W10 = exp(i(alpha+(beta/2)-(delta/2))sin(gamma/2)
    //W11 = exp(i(alpha+(beta/2)+(delta/2))cos(gamma/2)
    //
    //phi := beta; lambda = delta; theta = gamma; exp(j.alpha).U3(theta, phi, lambda) == exp(j.alpha).U3(gamma, beta, delta)
    double R00=abs(W(0, 0)), R10=abs(W(1, 0));
    gamma = 2*atan2(R10, R00); //Dont miss the phase/sign-factor!

    double CosGam2=cos(gamma/2), SinGam2=sin(gamma/2);
#ifdef VERY_VERBOSE
    if (Verbose) cout << "R00: " << R00 << ", R10: " << R10 << ", CosGam2: " << CosGam2 << ", SinGam2: " << SinGam2 << "\n";
#endif

    complex<double> D = sqrt(det2x2(W)); //exp(i*alpha); sqrt() means there is a PI ambiguity in angle alpha
    alpha = arg(D);

    //Factor out the common phase, and make V \in SU(2)
    cx_mat V = W/exp(j*alpha); //(z, -w; w', z';), z := exp(i(-beta-delta)/2)cos(gamma/2); w := exp(i((-beta+delta)/2))sin(gamma/2)
    if (abs(V(0, 0) - conj(V(1, 1))) > QC_TOL) { cerr << "V not SU(2)\n"; exit(1); }
    if (abs(V(0, 1) + conj(V(1, 0))) > QC_TOL) { cerr << "V not SU(2)\n"; exit(1); }
#ifdef VERY_VERBOSE
    if (Verbose) {
        cout << "V:\n" << V << "\n";
        cout << "Phase (exp(i.alpha)): " << exp(j*alpha) << "\n";
    }
#endif

    if (abs(CosGam2) < QC_TOL) {
        beta = arg(V(1, 0)); //We know only (beta-delta)
        delta = -beta;
    }
    else {
        if (abs(SinGam2) < QC_TOL) {
            beta = arg(V(1, 1)); //We know only (beta+delta)
            delta = beta;
        }
        else {
            beta = -arg(-V(0, 0)*V(0, 1)/(CosGam2*SinGam2));
            delta = -arg(V(0, 0)*V(1, 0)/(CosGam2*SinGam2));
            if (Verbose) { //Compare to alternate expressions for beta/delta
                double beta1 = arg(V(1, 0)*V(1, 1)/(CosGam2*SinGam2));
                if (abs(beta - beta1) > QC_TOL) cout << "beta mismatch (" << beta << ", " << beta1 << ")\n";
                double delta1 = arg(-V(0, 1)*V(1, 1)/(CosGam2*SinGam2));
                if (abs(delta - delta1) > QC_TOL) cout << "delta mismatch (" << delta << ", " << delta1 << ")\n";
            }
        }
    }

    //Verify if we got the parameterization correct
    cx_mat WW(2, 2);
    WW << cos(gamma/2)*exp(j*(alpha-(beta/2)-(delta/2)))
       << -sin(gamma/2)*exp(j*(alpha-(beta/2)+(delta/2))) << endr
       << sin(gamma/2)*exp(j*(alpha+(beta/2)-(delta/2)))
       << cos(gamma/2)*exp(j*(alpha+(beta/2)+(delta/2)));
    if (l2norm(WW - W) > QC_TOL) {
        if (Verbose) {
            cout << "Global phase had to be adjusted with a PI rotation\n";
            cout << "\tW:\n" << W << "\n";
        }
        alpha += M_PI;
    }
#if VERY_VERBOSE
    else if (Verbose) cout << "Global phase NOT adjusted\n";
#endif
    WW << cos(gamma/2)*exp(j*(alpha-(beta/2)-(delta/2)))
       << -sin(gamma/2)*exp(j*(alpha-(beta/2)+(delta/2))) << endr
       << sin(gamma/2)*exp(j*(alpha+(beta/2)-(delta/2)))
       << cos(gamma/2)*exp(j*(alpha+(beta/2)+(delta/2)));
    if (l2norm(WW - W) > QC_TOL) {
        cerr << __FUNCTION__ << " failed:\n\tW:\n" << W << "\n\tW.W':" << W*W.t() << "\n\tWW:\n" << WW << "\n\tWW.WW':" << WW*WW.t() << "\n";
        cout << "**** <R00: " << R00 << ", R10: " << R10 << ">\n";
        cout << "**** <CosGam2: " << CosGam2 << ", SinGam2: " << SinGam2 << ">\n";
        cout << "**** <phi00: " << arg(W(0, 0)) << ", phi01: " << arg(-W(0, 1)) << ", phi10: " << arg(W(1, 0)) << ", phi11: " << arg(W(1, 1)) << ">\n";
        cout << "**** <phi00': " << (alpha-(beta/2)-(delta/2)) << ", phi01': " << (alpha-(beta/2)+(delta/2)) << ", phi10': " << (alpha+(beta/2)-(delta/2)) << ", phi11': " << (alpha+(beta/2)+(delta/2)) << ">\n";
        cout << "**** <alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << ", delta: " << delta << ">\n";
        exit(1);
    }

    if (Verbose) cout << "**** <alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << ", delta: " << delta << ">\n";
}

void C1Wsimple(const cx_mat& W, cx_mat& A, cx_mat& B, cx_mat& C, cx_mat& D) { //(UNUSED)
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;

    const complex<double> j(0, 1);
    //This function is of pedantic interest, as IBM OPENQASM does not take unitary matrix as direct inputs
    //Decompose Controlled-W into 1x each 1-qubit A, B, C, D and 2x 2-qubit CNOT gates
    //Step 1: Get alpha, beta, gamma, delta
    double alpha, beta, gamma, delta;
    param4simple(W, alpha, beta, gamma, delta);

    //Return A, B, C, D; st., A.B.C = I, and, D.(A.X.B.X.C) = W
    //
    //-------*--------     -----*-----*--D--
    //       |          =>      |     |
    //-------W--------     --C--+--B--+--A--
    //
    //U3(theta, phi, lambda) == U3(gamma, beta, delta)
    D = zD(alpha);
    A = Rz(beta)*Ry(gamma/2); //u3(gamma/2, beta, 0)
    B = Ry(-gamma/2)*Rz(-(delta+beta)/2); //u3(-gamma/2, 0, -(beta+lambda)/2)
    C = Rz((delta-beta)/2); //u1((delta-beta)/2)

    cx_mat ABC(A*B*C);
    cx_mat dAXBXC(exp(j*alpha)*A*X*B*X*C);
    if ((l2norm(I2 - ABC) > QC_TOL) || (l2norm(W - dAXBXC) > QC_TOL)) {
        cout << __FUNCTION__ << " failed:\n\tW:\n" << W << "\n\tCBA:\n" << ABC << "\n\tdCXBXA:\n" << dAXBXC << "\n";
        exit(1);
    }
}

void lowerPrimitives(QNet_t& Net) { //Append gate to netlist, and update gateCount, and totalTime
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    //Gates in Net are either single-target only or single-qubit with 1-control
    //We will lower them to 1-target and CNOT
    //We will also handle indirect-edge C-W
    vector<vector<QG_t>> Gates(Net.Gates.size()); //One vector for each qubit timeline
    Gates.swap(Net.Gates);
    Net.nxtSerial = 0;

    vector<uint64_t> Counter(Gates.size(), 0); //One counter for each qubit timeline

    while (true) {
        int16_t line=-1; //Dominating qubit-line, ie., least serial# gate among all qubit-lines
        for (uint8_t i=0; i<Gates.size(); i++) { //Get the first gate from among all qubit timelines
            if (Counter[i] >= Gates[i].size()) continue;
            if (line >= 0) if (Gates[i][Counter[i]].serial >= Gates[line][Counter[line]].serial) continue;
            line = i;
        }
        if (line < 0) break;
        QG_t& gate = Gates[line][Counter[line]];
        Counter[line]++;

        if (gate.type == Gctl) continue; //This is a control-tap place-holder, so ignore

        uint8_t s = gate.arg.size();
        if (s > 2) { cerr << "Gates should be 1-qubit with either 1 or 0 control qubits\n"; exit(1); }

        uint64_t t=gate.arg[s-1], c=gate.arg[0]; //target, control qubits
        if (s == 0) { cerr << "Gate has no target!\n"; exit(1); }
        else if (s == 1) mergeGate(Net, gate); //No control-qubits
        else {
            if (l2norm(X - gate.W) < QC_TOL) {
                if (Verbose) cout << "Convert indirect C-NOT to direct C-NOT and simple gates\n";
                mergeDirectIndirectCNOT(Net, c, t); //If the gate is X
            }
            else {
                if (Verbose) cout << "Convert C-W to simple gates\n";
                //Get A, B, C, D; st., A.B.C = I, and, D.(A.X.B.X.C) = W
                //
                //-------*--------     -----*-----*--D-- [c]
                //       |          =>      |     |
                //-------W--------     --C--+--B--+--A-- [t]
                //
                cx_mat A, B, C, D;
                C1Wsimple(gate.W, A, B, C, D);
                mergeGate(Net, {C, {t}, NO_VAL, Gnull});
                mergeDirectIndirectCNOT(Net, c, t);
                mergeGate(Net, {B, {t}, NO_VAL, Gnull});
                mergeDirectIndirectCNOT(Net, c, t);
                mergeGate(Net, {A, {t}, NO_VAL, Gnull});
                mergeGate(Net, {D, {c}, NO_VAL, Gnull});
            }
        }
    }
}

void mergeDirectIndirectCNOT(QNet_t& Net, uint64_t c, uint64_t t) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if (Net.useTop) {
        if (Net.topology(c, t) != 0) //Direct link exists, ie., c -> t
            mergeGate(Net, {X, {c, t}, NO_VAL, Gnull});
        else if (Net.topology(t, c) != 0) { //Indirect link exists, ie., t -> c
            mergeGate(Net, {H, {t}, NO_VAL, Gnull});
            mergeGate(Net, {H, {c}, NO_VAL, Gnull});
            mergeGate(Net, {X, {t, c}, NO_VAL, Gnull}); //Swap roles of c, t
            mergeGate(Net, {H, {t}, NO_VAL, Gnull});
            mergeGate(Net, {H, {c}, NO_VAL, Gnull});
        }
        else { cerr << "No 1-hop connection (" << c << " <=> " << t << ")\n"; exit(1); } //Should not happen at this stage (ie., after topology Mapping)
    }
    else mergeGate(Net, {X, {c, t}, NO_VAL, Gnull});
}
