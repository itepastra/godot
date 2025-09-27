#pragma once

#include "core/object/ref_counted.h"
#include <cstdint>
#include <memory>
#include <vector>

#define BLOCK_BITS 6

class Qec : public RefCounted {
	GDCLASS(Qec, RefCounted);

	bool initialized;
	size_t n;
	std::vector<std::vector<uint64_t>> x_stabilizers; // 2n + 1 by n/4 array (the longs inside are packed bits)
	std::vector<std::vector<uint64_t>> z_stabilizers; // same but for z stabilizers
	std::vector<uint_fast8_t> phases;

protected:
	// godot helper functions
	static void _bind_methods();

	// tableau helper functions
	void rowcopy(size_t i, size_t k); // set row i equal to row k
	void rowswap(size_t i, size_t k); // swap row i and row k
	void rowset(size_t i, size_t b); // set row i to the b-th observable

public:
	// qubit gates for godot to use
	void cnot(size_t control, size_t target); // CNOT
	void hadamard(size_t target); // H
	void phase(size_t target); // S
	void phase_dag(size_t target); // S^+ = S S S
	void xgate(size_t target); // X = H S S H
	void ygate(size_t target); // Y = X Y = H S S H S S
	void zgate(size_t target); // Z = S S
	void cphase(size_t control, size_t target); // CZ = H_target CNOT H_target

	// initialisation
	Qec();
	void init(size_t qubit_amount);
};

const uint64_t powers[1 << BLOCK_BITS] = {
	1ull << 0,
	1ull << 1,
	1ull << 2,
	1ull << 3,
	1ull << 4,
	1ull << 5,
	1ull << 6,
	1ull << 7,
	1ull << 8,
	1ull << 9,
	1ull << 10,
	1ull << 11,
	1ull << 12,
	1ull << 13,
	1ull << 14,
	1ull << 15,
	1ull << 16,
	1ull << 17,
	1ull << 18,
	1ull << 19,
	1ull << 20,
	1ull << 21,
	1ull << 22,
	1ull << 23,
	1ull << 24,
	1ull << 25,
	1ull << 26,
	1ull << 27,
	1ull << 28,
	1ull << 29,
	1ull << 30,
	1ull << 31,
	1ull << 32,
	1ull << 33,
	1ull << 34,
	1ull << 35,
	1ull << 36,
	1ull << 37,
	1ull << 38,
	1ull << 39,
	1ull << 40,
	1ull << 41,
	1ull << 42,
	1ull << 43,
	1ull << 44,
	1ull << 45,
	1ull << 46,
	1ull << 47,
	1ull << 48,
	1ull << 49,
	1ull << 50,
	1ull << 51,
	1ull << 52,
	1ull << 53,
	1ull << 54,
	1ull << 55,
	1ull << 56,
	1ull << 57,
	1ull << 58,
	1ull << 59,
	1ull << 60,
	1ull << 61,
	1ull << 62,
	1ull << 63,
};
