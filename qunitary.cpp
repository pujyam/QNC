#include "qutils.h"
#include "qgate.h"
#include "qunitary.h"

uint8_t compileUnitaryOp(const string& inFilename, QNet_t& Net, ofstream& outFile) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    cout << "File: " << inFilename << "\n";
    cx_mat U;
    U.load(inFilename, raw_ascii);
    if ((U.n_rows == U.n_cols) && (U.n_rows == 0)) { cerr << "FORMAT ERROR: matrix could not be read\n"; exit(1); }
    if (Verbose) U.print("U:");

    const uint64_t N = U.n_rows;
    const uint64_t lgN = round(log2(N));
    uint64_t QubitSize = lgN;

    if (Net.useMap && (Net.physMap.size() < QubitSize)) { cerr << "Provided logical to physical qubit mapping is incompatible with the unitary operator\n"; exit(1); }
    if (Net.useMap && (Net.physMap.size() > QubitSize)) QubitSize = Net.physMap.size();
    if (Net.useTop && (Net.topology.n_rows < QubitSize)) { cerr << "Provided topology graph is incompatible with the unitary operator\n"; exit(1); }
    if (Net.useTop && (Net.topology.n_rows > QubitSize)) QubitSize = Net.topology.n_rows;

    if (doStatePrep) {
        if (1 != U.n_cols) { cerr << "Initial-state s should be a column-vector\n"; exit(1); }
        double E = l2norm(U);
        if (abs(E - 1.0) > QC_TOL) { cerr << "Initial-state s should have unit energy. Got energy = " << E << "\n"; exit(1); }
    }
    else {
        if (N != U.n_cols) { cerr << "Matrix U should be square\n"; exit(1); }
        if (l2norm(U*U.t() - eye<cx_mat>(arma::size(U))) > QC_TOL) { cerr << "Matrix U should be unitary\n"; exit(1); }
    }

    vector<uint64_t> gamma(N);
    getRowEncoding(gamma);

    vector<cx_mat> VG;
    vector<vector<uint64_t>> VGcontrols;
    GivensDecompose(U, gamma, VG, VGcontrols); //Givens decompose inverse of U

    outFile << "\n//Quantum netlist implementing unitary operator '"
            << inFilename << "' using simple gates:\n";

    Net.Gates = vector<vector<QG_t>>(QubitSize); //One vector for each qubit timeline
    if (Verbose) cout << "\n________________________________\n";
    //\Pi_i G_i' == U; so apply the G_i' in reverse order to affect U (post-fix)
    for (int64_t i=VG.size()-1; i>=0; i--) lowerCnG(lgN, VG[i].t(), VGcontrols[i], Net, !useNaturalPivot);
    if (Verbose) cout << "________________________________\n" << endl;

    return lgN;
}

void getRowEncoding(vector<uint64_t>& gamma) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    //Form a standard Gray code (not unique) or use Natural encoding
    for(uint64_t i=0; i<gamma.size(); i++) gamma[i] = (useNaturalPivot ? i : (i ^ (i/2)));

#ifdef VERY_VERBOSE
    if (Verbose) {
        cout << "Gamma: [ ";
        for (auto i = gamma.begin(); i != gamma.end(); i++) cout << *i << ' ';
        cout << "]\n" << endl;
    }
#endif
}

inline void adjustGivensRows(cx_mat& U, const cx_mat& G, uint64_t gJ, uint64_t gJ1) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    //Adjust the rows gJ, gJ1 of U
    //k,c: u'_k,i * u_k,c + u'_j,i * u_j,c / Norm (== 1 when c==i)
    //j,c: -u_j,i * u_k,c + u_k,i * u_j,c / Norm (== 0 when c==i)
    for (uint64_t c=0; c<U.n_cols; c++) {
        complex<double> tmp = G(0, 0)*U(gJ1, c) + G(0, 1)*U(gJ, c);
        U(gJ, c) = G(1, 0)*U(gJ1, c) + G(1, 1)*U(gJ, c);
        U(gJ1, c) = tmp;
    }
}

void GivensDecompose(const cx_mat& UU, const vector<uint64_t>& gamma, vector<cx_mat>& VG, vector<vector<uint64_t>>& VGcontrols) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    // 2 use cases:
    //
    // (A) Unitary Operator:
    // [(\Pi_{i=1..m} G_i).U = I] => [(\Pi_{i=m..1} G_i') = U]
    //
    // (B) State Preparation:
    // [(\Pi_{i=1..m} G_i).s = e0] => [(\Pi_{i=m..1} G_i').e0 = s]
    //
    //So, in either case, we need the quantum-gate implementation for the product: (\Pi_{i=m..1} G_i')
    cx_mat U = UU; //Work on a copy. NOTE: Matrix U gets clobbered to I in the end
    const uint64_t N = U.n_rows;

    for (uint64_t i=0; i<(doStatePrep ? 1 : (N-1)); i++) { //Column loop
        uint64_t gI = gamma[i];
        for (uint64_t j=(N-1); j>i; j--) { //Row loop
            uint64_t gJ = gamma[j], gJ1 = gamma[j-1];
            cx_mat G(2, 2);
            G << conj(U(gJ1, gI)) << conj(U(gJ, gI)) << endr << -U(gJ, gI) << U(gJ1, gI);
            complex<double> dG = det2x2(G);
            if (abs(dG) < QC_TOL) continue;
            G /= sqrt(dG); //So now G is unitary
            if (l2norm(G - I2) > QC_TOL) { //Null (gJ, gI) using U(gJ1, gI)
                if (Verbose)
                    cout << "\n________________________________\n"
                         << "Objective: Null U(" << gJ << ", " << gI << ") using U(" << gJ1 << ", " << gI << ")\n";

                VG.push_back(G); //Multiply elements in V in reverse to get back U
                vector<uint64_t> control(3);
                control[0] = gJ; control[1] = gJ1; control[2] = gI;
                VGcontrols.push_back(control);
                if (Verbose) {
                    cout << "Givens " << gI << "_G_(" << gJ << "," << gJ1 << ")" << endl;
                    cout << G << endl;
                }

                //Adjust the rows gJ, gJ1 of U
                adjustGivensRows(U, G, gJ, gJ1);

                if (Verbose) {
                    U.print("U:");
                    cout << "________________________________\n";
                }
            }
            if (abs(U(gJ, gI)) > QC_TOL) { cerr << "(" << gJ << "," << gI << ") : " << U(gJ, gI) << ": Matrix was not Givens deomposed correctly\n"; U.print("U:"); exit(1); }
        }
        if (abs(U(gI, gI)-1.0) > QC_TOL) { cerr << "(" << gI << "," << gI << ") : " << U(gI, gI) << ": Matrix was not Givens deomposed correctly\n"; exit(1); }
    }

    complex<double> dU = (doStatePrep ? U(0, 0) : diagProd(U)); //abs(dU) == 1, guaranteed as either ||s||_2 = 1 or |det(U)| = 1
    if (abs(dU-1.0) > QC_TOL) {
        //Only the gamma[(N-1)]^th diagonal element can be non-unity for unitary operators
        //Similarly, only the gamma[0]^th element can be non-unity for state-preparation
        uint64_t gI = gamma[doStatePrep ? 0 : (N-1)], gJ = gamma[doStatePrep ? 1 : 0], gJ1 = gamma[doStatePrep ? 0 : (N-1)];
        //We use the (gJ, gI)^th entry (ie., by extension, the entire gJ^th row) to make (gJ1, gI)^th entry 1
        cx_mat G(2, 2);
        G << 1.0/dU << 0.0 << endr << 0.0 << 1.0;
        VG.push_back(G); //Multiply elements in V in reverse to get back U
        vector<uint64_t> control(3);
        control[0] = gJ; control[1] = gJ1; control[2] = gI;
        VGcontrols.push_back(control);

        adjustGivensRows(U, G, gJ, gJ1);

        bool done = (doStatePrep ? true : isDiagIdentity(U));
        if (!done) { //Elaborate check that the Givens decomposition was correct (Skip for state-preparations)
            cx_mat P=eye<cx_mat>(arma::size(UU));
            for (int64_t i=VG.size()-1; i>=0; i--) {
                cx_mat Qi=eye<cx_mat>(arma::size(UU));

                uint64_t J1=VGcontrols[i][1], J=VGcontrols[i][0];
                Qi(J1, J1) = VG[i](0, 0);
                Qi(J1, J) = VG[i](0, 1);
                Qi(J, J1) = VG[i](1, 0);
                Qi(J, J) = VG[i](1, 1);

                P *= Qi;
            }
            P *= UU;

            if (l2norm(P - eye<cx_mat>(arma::size(P))) > QC_TOL) {
                cerr << "Givens decomposition was not exact!\n\n";
                if (Verbose) cerr << "[(\\Pi_i G_i) * Op]:\n" << P << "\n\n";
                exit(1);
            }
        }
    }
}

void lowerCnG(uint64_t S, const cx_mat& G, const vector<uint64_t>& control, QNet_t& Net, bool useGrayCode) {
    if (Verbose) cout << "In function " << __FUNCTION__ << endl;
    if (Verbose) {
        if (control.size() == 3) //This is for Givens lowering
            cout << "Givens " << control[2] << "_G_(" << control[0] << "," << control[1] << "):\n";
        else //This is for Swap lowering (Natural encoding only)
            cout << "Swap S_(" << control[0] << "," << control[1] << "):\n";
        cout << G << endl;
    }

    uint64_t first=control[0], last=control[1], all=((1<<S)-1);
    if (useGrayCode) {
        if (control.size() != 3) { cerr << "Swaps cannot be lowered with Gray coded qubit sequence\n"; exit(1); }
        if (Verbose) cout << "Using Gray encoded qubit lines\n";
        //Form single qubit NOT gate network to take the mask to all 1 for the two rows on which G should act
        const uint64_t toggle = (first ^ last); //Bits which need to be '0' to act as control
        const uint8_t target = round(log2(toggle));
        const uint64_t mask = (~first) & all; //Bits which need to be '0' to act as control
        const uint64_t cMask = (~toggle) & all;
#ifdef VERY_VERBOSE
        if (Verbose) cout << "Mask: " << mask << ", Toggle: " << toggle << endl;
#endif

        //Prep mask bits to be all '1'
        if (cMask) for (uint8_t b=0; b<S; b++)
            if (mask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act

        //Convert the C^(S-1)G gate to a network of simpler single qubit and CNOT gates in a recursive manner
        lowerCnW(S, G, target, cMask, Net); //1-qubit gate, control-qubits, target-qubit

        //Reverse-prep mask bits to be what they were originally
        if (cMask) for (uint8_t b=0; b<S; b++)
            if (mask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act
    }
    else {
        if (Verbose) cout << "Using Naturally encoded qubit lines\n";
        //NOTE: Does not matter that we always take c1 as first. If c2 were to be 1st, #bit-flips remain same
        //Setup so that last is all-1's
        const uint64_t gmask = (~last) & all; //Bits to be negated to bring last-index to all-one's
        last = all;
        first ^= gmask; //After transforming last <= all; this is new first

#ifdef VERY_VERBOSE
        if (Verbose) cout << "First: " << first << ", Last: " << last << ", Gmask: " << gmask << "\n";
#endif
        for (uint8_t b=0; b<S; b++)
            if (gmask & (1<<b)) {
                addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act
                if (Verbose) cout << "X(" << int(b) << ")\n";
            }

        //Walk from first to last-but-one with a Gray code, flipping 1 bit at a time
        uint64_t current=first, tmask;
        int8_t target;
        for (target=0; target<S; target++) {
            tmask = (1<<target);
            if (!(first & tmask)) {
                uint64_t next = current^tmask;
                if (next == last) break; //When this loop breaks, target := last-bit-to-flip-towards-all-1's
                current = next;

                uint64_t mask = ((~current)&last) ^ tmask;
                for (uint8_t b=0; b<S; b++)
                    if (mask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act

                //Convert the C^(S-1)X gate to a network of simpler single qubit and CNOT gates in a recursive manner
                lowerCnW(S, X, target, last^tmask, Net); //1-qubit gate, control-qubits, target-qubit

                for (uint8_t b=0; b<S; b++)
                    if (mask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act
            }
        }

        //Convert the C^(S-1)X gate to a network of simpler single qubit and CNOT gates in a recursive manner
        lowerCnW(S, G, target, all^tmask, Net); //1-qubit gate, control-qubits, target-qubit

        //Walk back from last-but-one to first with same Gray code, flipping 1 bit at a time
        for (target=target-1; target>=0; target--) {
            tmask = (1<<target);
            if (!(first & tmask)) {
                current ^= tmask;

                uint64_t mask = ((~current)&last) ^ tmask;
                for (uint8_t b=0; b<S; b++)
                    if (mask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act

                //Convert the C^(S-1)X gate to a network of simpler single qubit and CNOT gates in a recursive manner
                lowerCnW(S, X, target, last^tmask, Net); //1-qubit gate, control-qubits, target-qubit

                for (uint8_t b=0; b<S; b++)
                    if (mask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act

                if (current == first) break;
            }
        }

        //Unwind the gmask, by reversing the NOTs
        for (uint8_t b=0; b<S; b++)
            if (gmask & (1<<b)) addGate(Net, X, {b}); //If this bit needs to be '0' for the control-G to act
    }
}
