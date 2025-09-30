#include "qec.h"

#include <iostream>
#include <random>
#include <stdexcept>

Qec::Qec() {
    
}

void Qec::init(u_int32_t num_qubits) {
    n_qubits = num_qubits;
    x = std::vector<std::vector<bool>>(2 * n_qubits + 1, std::vector<bool>(n_qubits, false));
    z = std::vector<std::vector<bool>>(2 * n_qubits + 1, std::vector<bool>(n_qubits, false));
    r = std::vector<bool>(2 * n_qubits + 1, false);
    rng = std::mt19937(std::random_device{}());
    
    // Initialize to identity matrix tableau
    // Destabilizer generators rows 0 to n-1
    // Stabilizer generators rows n to 2n-1
    // Row 2n is scratch space
    
    for (size_t i = 0; i < n_qubits; i++) {
        // Destabilizer rows x_ij = delta_ij, z_ij = 0
        x[i][i] = true;
        
        // Stabilizer rows x_ij = 0, z_ij = delta_{(i-n)j}
        z[n_qubits + i][i] = true;
    }
}

int Qec::g(bool x1, bool z1, bool x2, bool z2) const {
    if (!x1 && !z1) return 0;  // I case
    if (x1 && z1) return z2 - x2;  // Y case
    if (x1 && !z1) return z2 * (2 * x2 - 1);  // Z case
    if (!x1 && z1) return x2 * (1 - 2 * z2);  // X case
    return 0;
}

void Qec::rowsum(size_t h, size_t i) {
    int sum = 0;
    
    // phase sum
    for (size_t j = 0; j < n_qubits; j++) {
        sum += g(x[i][j], z[i][j], x[h][j], z[h][j]);
    }
    
    // updarte phase bit r_h
    int total = 2 * r[h] + 2 * r[i] + sum;
    r[h] = (total % 4 == 2);
    
    // update Paulis
    for (size_t j = 0; j < n_qubits; j++) {
        x[h][j] = x[h][j] ^ x[i][j];
        z[h][j] = z[h][j] ^ z[i][j];
    }
}

void Qec::cnot(size_t control, size_t target) {
    if (control >= n_qubits || target >= n_qubits) {
        throw std::out_of_range("Qubit index out of range");
    }
    
    for (size_t i = 0; i < 2 * n_qubits; i++) {
        // Update phase using r_i = r_i tensor x_ia * z_ib * (x_ib tensor z_ia tensor 1)
        bool x_ia = x[i][control];
        bool z_ia = z[i][control];
        bool x_ib = x[i][target];
        bool z_ib = z[i][target];
        
        r[i] = r[i] ^ (x_ia && z_ib && (x_ib ^ z_ia ^ 1));
        
        // Update Pauli components
        x[i][target] = x_ib ^ x_ia;
        z[i][control] = z_ia ^ z_ib;
    }
}

void Qec::hadamard(size_t qubit) {
    if (qubit >= n_qubits) {
        throw std::out_of_range("Qubit index out of range");
    }
    
    for (size_t i = 0; i < 2 * n_qubits; i++) {
        // Update phase r_i = r_i tensor x_ia * z_ia
        r[i] = r[i] ^ (x[i][qubit] && z[i][qubit]);
        
        // Swap x_ia and z_ia
        std::swap(x[i][qubit], z[i][qubit]);
    }
}

void Qec::phase(size_t qubit) {
    if (qubit >= n_qubits) {
        throw std::out_of_range("Qubit index out of range");
    }
    
    for (size_t i = 0; i < 2 * n_qubits; i++) {
        // Update phase: r_i = r_i tensor x_ia * z_ia
        r[i] = r[i] ^ (x[i][qubit] && z[i][qubit]);
        
        // Update z_ia = z_ia tensor x_ia
        z[i][qubit] = z[i][qubit] ^ x[i][qubit];
    }
}

bool Qec::measure(size_t qubit) {
    if (qubit >= n_qubits) {
        throw std::out_of_range("Qubit index out of range");
    }
    
    // Find first stabilizer with x_pa = 1
    size_t p = 2 * n_qubits; // invalid idx
    for (size_t i = n_qubits; i < 2 * n_qubits; i++) {
        if (x[i][qubit]) {
            p = i;
            break;
        }
    }

    // case 1 rndom measurement
    if (p < 2 * n_qubits) {
        for (size_t i = 0; i < 2 * n_qubits; i++) {
            if (i != p && x[i][qubit]) {
                rowsum(i, p);
            }
        }
        
        // setting (p-n)th destabilizer row equal to p-th stabilizer row
        size_t destab_index = p - n_qubits;
        for (size_t j = 0; j < n_qubits; j++) {
            x[destab_index][j] = x[p][j];
            z[destab_index][j] = z[p][j];
        }
        r[destab_index] = r[p];
        
        // resetting p-th stabilizer row
        for (size_t j = 0; j < n_qubits; j++) {
            x[p][j] = false;
            z[p][j] = false;
        }
        z[p][qubit] = true;
        
        // random outcome
        std::uniform_int_distribution<int> dist(0, 1);
        r[p] = dist(rng);
        
        return r[p];
    } else { // case 2 deterministic
        // resetting scratch row (index 2n)
        for (size_t j = 0; j < n_qubits; j++) {
            x[2 * n_qubits][j] = false;
            z[2 * n_qubits][j] = false;
        }
        r[2 * n_qubits] = false;

        // sum correct stabilizer rows
        for (size_t i = 0; i < n_qubits; i++) {
            if (x[i][qubit]) {
                rowsum(2 * n_qubits, n_qubits + i);
            }
        }
        return r[2 * n_qubits];
    }
}

void Qec::print_tableau() const {
    std::cout << "Destabilizer generators:\n";
    for (size_t i = 0; i < n_qubits; i++) {
        std::cout << (r[i] ? "- " : "+ ");
        for (size_t j = 0; j < n_qubits; j++) {
            if (x[i][j] && z[i][j]) std::cout << "Y";
            else if (x[i][j]) std::cout << "Z";
            else if (z[i][j]) std::cout << "X";
            else std::cout << "I";
        }
        std::cout << "\n";
    }
    
    std::cout << "Stabilizer generators:\n";
    for (size_t i = n_qubits; i < 2 * n_qubits; i++) {
        std::cout << (r[i] ? "- " : "+ ");
        for (size_t j = 0; j < n_qubits; j++) {
            if (x[i][j] && z[i][j]) std::cout << "Y";
            else if (x[i][j]) std::cout << "Z";
            else if (z[i][j]) std::cout << "X";
            else std::cout << "I";
        }
        std::cout << "\n";
    }
}

void Qec::_bind_methods() {
    ClassDB::bind_method(D_METHOD("init", "num_qubits"), &Qec::init);

	ClassDB::bind_method(D_METHOD("cnot", "control", "target"), &Qec::cnot);
	ClassDB::bind_method(D_METHOD("hadamard", "qubit"), &Qec::hadamard);
	ClassDB::bind_method(D_METHOD("phase", "qubit"), &Qec::phase);

	ClassDB::bind_method(D_METHOD("measure", "qubit"), &Qec::measure);
	ClassDB::bind_method(D_METHOD("print_tableau"), &Qec::print_tableau);
}
