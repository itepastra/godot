#include "qec.h"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <random>


Qec::Qec() {
	this->initialized = false;
	this->rng = std::mt19937(12345);
	this->bit_dist = std::uniform_int_distribution<int>(0, 1);
}

void Qec::init(uint32_t qubit_amount) {
	if (this->initialized) {
		return;
	}
	// uint32_t whole_blocks = (qubit_amount >> BLOCK_BITS) + 1; // shrink by amount of qubits that fit in 64 bits
	uint32_t whole_blocks = (qubit_amount + 63) >> BLOCK_BITS;
	uint32_t rows = 2 * qubit_amount + 1;
	this->n = qubit_amount;

	this->phases = std::vector<uint_fast8_t>(2 * qubit_amount + 1, 0);
	this->x_stabilizers.reserve(rows);
	this->z_stabilizers.reserve(rows);

	for (uint32_t i = 0; i < rows; i++) {
		std::vector<uint64_t> line_x(whole_blocks, 0);
		std::vector<uint64_t> line_z(whole_blocks, 0);

		if (i < qubit_amount) {
			line_x[i >> BLOCK_BITS] = powers[i & 63]; // we create a diagonal
		} else if (i < rows - 1) {
			line_z[(i - qubit_amount) >> BLOCK_BITS] = powers[(i - qubit_amount) & 63]; // same with z, but on a different spot
		}
		this->x_stabilizers.push_back(line_x);
		this->z_stabilizers.push_back(line_z);
	}

	this->initialized = true;
}

uint64_t Qec::get_x_stab(uint32_t i, uint32_t j) {
	uint32_t block = j >> BLOCK_BITS;
	uint64_t inner = powers[j & 63];

	return this->x_stabilizers[i][block] & inner;
}

PackedInt64Array Qec::x_stabs() {
	PackedInt64Array array;
	uint32_t whole_blocks = (this->n + 63) >> BLOCK_BITS;
	for (uint32_t i = 0; i < 2 * this->n; i++) {
		for (uint32_t j = 0; j < whole_blocks; j++) {
 		   array.append(this->x_stabilizers[i][j]);
		}
		// for (uint32_t j = 0; j < (2 * this->n >> BLOCK_BITS) + 1; j++) {
		// 	array.append(this->x_stabilizers[i][j]);
		// }
	}
	return array;
}

uint64_t Qec::get_z_stab(uint32_t i, uint32_t j) {
	uint32_t block = j >> BLOCK_BITS;
	uint64_t inner = powers[j & 63];

	return this->z_stabilizers[i][block] & inner;
}

PackedInt64Array Qec::z_stabs() {
	PackedInt64Array array;
	uint32_t whole_blocks = (this->n + 63) >> BLOCK_BITS;
	for (uint32_t i = 0; i < 2 * this->n; i++) {
		for (uint32_t j = 0; j < whole_blocks; j++) {
    		array.append(this->z_stabilizers[i][j]);
		}
		// for (uint32_t j = 0; j < (2 * this->n >> BLOCK_BITS) + 1; j++) {
		// 	array.append(this->z_stabilizers[i][j]);
		// }
	}
	return array;
}

int Qec::measure(uint32_t target) {
    if (!this->initialized) {
        return 0;
    }

    uint32_t block = target >> BLOCK_BITS;
    uint64_t inner = powers[target & 63];

    // Find first stabilizer row p with X acting on the measured qubit
    uint32_t p = 2 * this->n;
    for (uint32_t i = this->n; i < 2 * this->n; i++) {
        if (this->x_stabilizers[i][block] & inner) {
            p = i;
            break;
        }
    }

    if (p < 2 * this->n) {
        // Case 1 Random measurement outcome (some stabilizer anticommutes)
        for (uint32_t i = 0; i < 2 * this->n; i++) {
            if (i != p && (this->x_stabilizers[i][block] & inner)) {
                this->rowmult(i, p);
            }
        }
        this->rowcopy(p - this->n, p);       // copy stabilizer to corresponding destabilizer
        this->rowset(p, target + this->n);   // set stabilizer row p to Z_target

        // Random outcome: 0 or 1
        int outcome = bit_dist(rng);
        this->phases[p] = (outcome == 1 ? 2 : 0);
        return outcome;
    } else {
        // Case 2 Deterministic measurement outcome
        // Row 2n is reserved as a scratch row for measurement
        this->phases[2 * this->n] = 0;
        uint32_t whole_blocks = (this->n + 63) >> BLOCK_BITS;
        for (uint32_t j = 0; j < whole_blocks; j++) {
            this->x_stabilizers[2 * this->n][j] = 0;
            this->z_stabilizers[2 * this->n][j] = 0;
        }

        // Multiply in stabilizers that anticommute with Z_target
        for (uint32_t i = 0; i < this->n; i++) {
            if (this->x_stabilizers[i][block] & inner) {
                this->rowmult(2 * this->n, this->n + i);
            }
        }

        // Outcome is determined by the resulting phase of the scratch row
        return (this->phases[2 * this->n] ? 1 : 0);
    }
}

uint_fast16_t Qec::get_phase(uint32_t i) {
	if (!this->initialized) {
		return 0;
	}
	return this->phases[i];
}

void Qec::cnot(uint32_t control, uint32_t target) {
	if (!this->initialized) {
		return;
	}
	uint32_t control_block = control >> BLOCK_BITS;
	uint32_t target_block = target >> BLOCK_BITS;
	uint64_t control_inner = powers[control & 63];
	uint64_t target_inner = powers[target & 63];

    for (uint32_t i = 0; i < 2 * this->n; i++) {
        bool x_control = (this->x_stabilizers[i][control_block] & control_inner) != 0;
        bool z_target  = (this->z_stabilizers[i][target_block]  & target_inner) != 0;
        bool x_target  = (this->x_stabilizers[i][target_block]  & target_inner) != 0;
        bool z_control = (this->z_stabilizers[i][control_block] & control_inner) != 0;

        // Phase update from Aaronson–Gottesman rule
        // If a row has X on control and Z on target
        // and NOT (X on target and Z on control simultaneously) then flip the phase.
        // if (x_control && z_target && !(x_target && z_control)) {
        //     this->phases[i] ^= 2; // toggle phase between 0 and 2
        // }
		if (x_control && z_target && ((x_target ^ z_control) == 0)) {
    		this->phases[i] ^= 2;
		}

        // Stabilizer update---
        // X(control) to X(control) X(target)
        if (x_control) {
            this->x_stabilizers[i][target_block] ^= target_inner;
        }
        // Z(target) to Z(control) Z(target)
        if (z_target) {
            this->z_stabilizers[i][control_block] ^= control_inner;
        }
    }
}

void Qec::cphase(uint32_t control, uint32_t target) {
	if (!this->initialized) {
		return;
	}
	this->hadamard(target);
	this->cnot(control, target);
	this->hadamard(target);
} // CZ = H_target CNOT H_target
void Qec::hadamard(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	uint32_t block = target >> BLOCK_BITS;
	uint64_t inner = powers[target & 63];
	uint64_t tmp;

	for (uint32_t i = 0; i < 2 * this->n; i++) {
		tmp = this->x_stabilizers[i][block];
		this->x_stabilizers[i][block] ^= (this->x_stabilizers[i][block] ^ this->z_stabilizers[i][block]) & inner;
		this->z_stabilizers[i][block] ^= (this->z_stabilizers[i][block] ^ tmp) & inner;
		if ((this->x_stabilizers[i][block] & inner) && (this->z_stabilizers[i][block] & inner)) {
			this->phases[i] ^= 0b10;
		}
	}
} // H
void Qec::phase(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	uint32_t block = target >> BLOCK_BITS;
	uint64_t inner = powers[target & 63];
	for (uint32_t i = 0; i < 2 * this->n; i++) {
		if ((this->x_stabilizers[i][block] & inner) && (this->z_stabilizers[i][block] & inner)) {
			this->phases[i] ^= 0b10; // flip sign
		}
		this->z_stabilizers[i][block] ^= this->x_stabilizers[i][block] & inner;
	}
} // S
void Qec::phase_dag(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	this->phase(target);
	this->phase(target);
	this->phase(target);
} // S^+ = S S S
void Qec::xgate(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	this->hadamard(target);
	this->phase(target);
	this->phase(target);
	this->hadamard(target);
} // X = H S S H
void Qec::ygate(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	this->hadamard(target);
	this->phase(target);
	this->phase(target);
	this->hadamard(target);
	this->phase(target);
	this->phase(target);
} // Y = Z X = S S H S S H      NOTE: work up to global phase, i think S_dag X S might be better
void Qec::zgate(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	this->phase(target);
	this->phase(target);
} // Z = S S


void Qec::rowcopy(uint32_t i, uint32_t k)
// copy row k to row i
{
	this->x_stabilizers[i].assign(this->x_stabilizers[k].begin(), this->x_stabilizers[k].end());
	this->z_stabilizers[i].assign(this->z_stabilizers[k].begin(), this->z_stabilizers[k].end());
	this->phases[i] = this->phases[k];
}

void Qec::rowmult(uint32_t i, uint32_t k) {
    // Make row i := row k + row i  (multiply Pauli_k * Pauli_i)
    int total = this->phases[i] + this->phases[k];

    for (uint32_t j = 0; j < this->n; ++j) {
        uint32_t block = j >> BLOCK_BITS;
        uint64_t bit   = powers[j & 63];

        // FIRST operand = source row k (paper's i), SECOND = dest row i (paper's h)
        bool x1 = (this->x_stabilizers[k][block] & bit) != 0;
        bool z1 = (this->z_stabilizers[k][block] & bit) != 0;
        bool x2 = (this->x_stabilizers[i][block] & bit) != 0;
        bool z2 = (this->z_stabilizers[i][block] & bit) != 0;

        int g = 0;
        if (!x1 && !z1) {
            g = 0;
        } else if (x1 && z1) {           // Y
            g = (z2 ? 1 : 0) - (x2 ? 1 : 0);     // z2 - x2
        } else if (x1 && !z1) {          // X
            g = z2 ? (x2 ? 1 : -1) : 0;          // z2*(2*x2 - 1)
        } else { // !x1 && z1         Z
            g = x2 ? (z2 ? -1 : 1) : 0;          // x2*(1 - 2*z2)
        }
        total += g;
    }

    // Normalize and set phase with correct g, m from {0,2}
    int m = ((total % 4) + 4) % 4;
    this->phases[i] = (m == 0) ? 0 : 2;

    // XOR the bit patterns: i := i ⊕ k
    uint32_t whole_blocks = (this->n + 63) >> BLOCK_BITS;
    for (uint32_t b = 0; b < whole_blocks; ++b) {
        this->x_stabilizers[i][b] ^= this->x_stabilizers[k][b];
        this->z_stabilizers[i][b] ^= this->z_stabilizers[k][b];
    }
}

void Qec::rowset(uint32_t i, uint32_t k) {
    uint32_t whole_blocks = (this->n + 63) >> BLOCK_BITS;
    uint32_t block = k >> BLOCK_BITS;
    uint64_t inner = powers[k & 63];
	
	for (uint32_t j = 0; j < whole_blocks; j++) {
		this->x_stabilizers[i][j] = 0;
		this->z_stabilizers[i][j] = 0;
	}
    this->phases[i] = 0;
    if (k < this->n) {
        this->x_stabilizers[i][block] = inner; // X_k
    } else {
        this->z_stabilizers[i][block] = inner; // Z_{k-n}
    }
}

void Qec::_bind_methods() {
	ClassDB::bind_method(D_METHOD("init", "qubit_amount"), &Qec::init);

	ClassDB::bind_method(D_METHOD("cnot", "control", "target"), &Qec::cnot);
	ClassDB::bind_method(D_METHOD("cphase", "control", "target"), &Qec::cphase);
	ClassDB::bind_method(D_METHOD("hadamard", "qubit"), &Qec::hadamard);
	ClassDB::bind_method(D_METHOD("phase", "qubit"), &Qec::phase);
	ClassDB::bind_method(D_METHOD("phase_dag", "qubit"), &Qec::phase_dag);
	ClassDB::bind_method(D_METHOD("xgate", "qubit"), &Qec::xgate);
	ClassDB::bind_method(D_METHOD("ygate", "qubit"), &Qec::ygate);
	ClassDB::bind_method(D_METHOD("zgate", "qubit"), &Qec::zgate);

	ClassDB::bind_method(D_METHOD("measure", "qubit"), &Qec::measure);

	ClassDB::bind_method(D_METHOD("get_x_stab", "i", "j"), &Qec::get_x_stab);
	ClassDB::bind_method(D_METHOD("get_z_stab", "i", "j"), &Qec::get_z_stab);
	ClassDB::bind_method(D_METHOD("x_stabs"), &Qec::x_stabs);
	ClassDB::bind_method(D_METHOD("z_stabs"), &Qec::z_stabs);
	ClassDB::bind_method(D_METHOD("get_phase", "i"), &Qec::get_phase);
}