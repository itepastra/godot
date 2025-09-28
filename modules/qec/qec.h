#pragma once

#include "core/object/ref_counted.h"
#include <vector>
#include <cstddef>
#include <cstdint>
#include <random>


class Qec : public RefCounted {
	GDCLASS(Qec, RefCounted);
private:
    size_t n_qubits;
    std::vector<std::vector<bool>> x; // 2n x n matrix for X
    std::vector<std::vector<bool>> z; // 2n x n matrix for Z  
    std::vector<bool> r; // 2n phase bits
    std::mt19937 rng;
    
    // Helper function for Pauli mult
    int g(bool x1, bool z1, bool x2, bool z2) const;
	void rowsum(size_t h, size_t i);

protected:
	static void _bind_methods();

public:
	Qec(size_t n_qubits);

    // operations
    void cnot(size_t control, size_t target);
    void hadamard(size_t qubit);
    void phase(size_t qubit);
    
    // measurement to return outcome and update state
    bool measure(size_t qubit);
    
    // helpers
    size_t num_qubits() const { return n_qubits; }
    void print_tableau() const; // for debug
};
