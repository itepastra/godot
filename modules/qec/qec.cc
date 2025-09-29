#include "qec.h"
#include <algorithm>
#include <cstdint>
#include <cstdlib>

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

uint64_t Qec::get_x_stab(uint32_t i, uint32_t j) {
	uint32_t block = j >> BLOCK_BITS;
	uint64_t inner = powers[j & 63];

	return this->x_stabilizers[i][block] & inner;
}

PackedInt64Array Qec::x_stabs() {
	PackedInt64Array array;
	for (uint32_t i = 0; i < 2 * this->n; i++) {
		for (uint32_t j = 0; j < (2 * this->n >> BLOCK_BITS) + 1; j++) {
			array.append(this->x_stabilizers[i][j]);
		}
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
	for (uint32_t i = 0; i < 2 * this->n; i++) {
		for (uint32_t j = 0; j < (2 * this->n >> BLOCK_BITS) + 1; j++) {
			array.append(this->z_stabilizers[i][j]);
		}
	}
	return array;
}

int Qec::measure(uint32_t target) {
	if (!this->initialized) {
		return 0;
	}
	uint32_t block = target >> BLOCK_BITS;
	uint64_t inner = powers[target & 63];

	for (uint32_t other = 0; other < this->n; other++) {
		if (this->x_stabilizers[other + this->n][block] & inner) {
			this->rowcopy(other, other + this->n);
			this->rowset(other + this->n, target + this->n);
			this->phases[other + this->n] = 2 * (rand() & 0b1); // NOTE: I don't trust the p + n index
			for (uint32_t i = 0; i < 2 * this->n; i++) {
				if ((i != other) && (this->x_stabilizers[i][block] & inner)) {
					this->rowmult(i, other);
				}
			}
			if (this->phases[other + this->n]) {
				return 3;
			} else {
				return 2;
			}
		}
	}

	uint32_t m;
	// the outcome is determinate
	for (m = 0; m < this->n; m++) {
		if (this->x_stabilizers[m][block] & inner) {
			break;
		}
	}
	this->rowcopy(2 * this->n, m + this->n);
	for (uint32_t i = m + 1; i < this->n; i++) {
		if (this->x_stabilizers[i][block] & inner) {
			this->rowmult(2 * this->n, i + this->n);
		}
	}
	if (this->phases[2 * this->n]) {
		return 1;
	} else {
		return 0;
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
		// Flip the x stabilizer on the target, if the control contains 1
		if (this->x_stabilizers[i][control_block] & control_inner) {
			this->x_stabilizers[i][target_block] ^= target_inner;
		}
		// Flip the z stabilizer on the control, if the target contains +
		if (this->z_stabilizers[i][target_block] & target_inner) {
			this->z_stabilizers[i][control_block] ^= control_inner;
		}
		// do the phase kickback when all the stabilizers are 1
		if ((this->x_stabilizers[i][control_block] & control_inner) &&
				(this->z_stabilizers[i][target_block] & target_inner) &&
				(this->x_stabilizers[i][target_block] & target_inner) &&
				(this->z_stabilizers[i][control_block] & control_inner)) {
			// TODO: check if doing `phase ^= 0b10` is equivalent and faster
			this->phases[i] ^= 0b10;
		}
		// do the phase kickback in another case TODO: figure out what case this represents
		if ((this->x_stabilizers[i][control_block] & control_inner) &&
				(this->z_stabilizers[i][target_block] & target_inner) &&
				!(this->x_stabilizers[i][target_block] & target_inner) &&
				!(this->z_stabilizers[i][control_block] & control_inner)) {
			this->phases[i] ^= 0b10;
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
			this->phases[i] ^= 0b10;
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
} // Y = Z X = S S H S S H
void Qec::zgate(uint32_t target) {
	if (!this->initialized) {
		return;
	}
	this->phase(target);
	this->phase(target);
} // Z = S S

void Qec::init(uint32_t qubit_amount) {
	if (this->initialized) {
		return;
	}
	uint32_t whole_blocks = (qubit_amount >> BLOCK_BITS) + 1; // shrink by amount of qubits that fit in 64 bits
	uint32_t rows = 2 * qubit_amount + 1;
	this->n = qubit_amount;

	this->phases = std::vector<uint_fast8_t>(2 * qubit_amount + 1, 0);
	this->x_stabilizers.reserve(rows);
	this->z_stabilizers.reserve(rows);

	for (uint32_t i = 0; i < rows; i++) {
		std::vector<uint64_t> line_x(whole_blocks, 0);
		std::vector<uint64_t> line_z(whole_blocks, 0);
		// TODO: what is this for?? gotta either figure it out or ask scott
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

Qec::Qec() {
	this->initialized = false;
}

void Qec::rowcopy(uint32_t i, uint32_t k)
// copy row k to row i
{
	this->x_stabilizers[i].assign(this->x_stabilizers[k].begin(), this->x_stabilizers[k].end());
	this->z_stabilizers[i].assign(this->z_stabilizers[k].begin(), this->z_stabilizers[k].end());
}

void Qec::rowswap(uint32_t i, uint32_t k)
// swap row i and row k
{
	// NOTE: currently I'm just swapping the pointers around, I don't know if this works as expected
	std::vector<uint64_t> tmp_x = this->x_stabilizers[i];
	this->x_stabilizers[i] = this->x_stabilizers[k];
	this->x_stabilizers[k] = tmp_x;
	std::vector<uint64_t> tmp_z = this->z_stabilizers[i];
	this->z_stabilizers[i] = this->z_stabilizers[k];
	this->z_stabilizers[k] = tmp_z;
}

int_fast8_t Qec::clifford(uint32_t i, uint32_t k) {
	int_fast32_t e = 0; // NOTE: check on how it overflows
	for (uint32_t j = 0; j < this->n; j++) {
		uint32_t block = j >> BLOCK_BITS;
		uint64_t inner = powers[j & 63];
		switch (
				(((this->x_stabilizers[k][block] & inner) != 0) << 3) |
				(((this->z_stabilizers[k][block] & inner) != 0) << 2) |
				(((this->x_stabilizers[i][block] & inner) != 0) << 1) |
				(((this->z_stabilizers[i][block] & inner) != 0))) {
			case 0b1011:
			case 0b1101:
			case 0b0110:
				e++;
				break;
			case 0b1001:
			case 0b1110:
			case 0b0111:
				e--;
				break;
		}
	}

	return e & 0b11; // return e mod 4 always positive
}

void Qec::rowmult(uint32_t i, uint32_t k) {
	uint32_t whole_blocks = (this->n >> BLOCK_BITS) + 1;
	this->phases[i] = this->clifford(i, k);
	for (uint32_t j = 0; j < whole_blocks; j++) {
		this->x_stabilizers[i][j] ^= this->x_stabilizers[k][j];
		this->z_stabilizers[i][j] ^= this->z_stabilizers[k][j];
	}
}

void Qec::rowset(uint32_t i, uint32_t k) {
	uint32_t whole_blocks = (this->n >> BLOCK_BITS) + 1;
	uint32_t block = k >> BLOCK_BITS;
	uint64_t inner = powers[k & 63];
	for (uint32_t j = 0; j < whole_blocks; j++) {
		this->x_stabilizers[i][j] = 0;
		this->z_stabilizers[i][j] = 0;
	}
	this->phases[i] = 0;
	if (k < this->n) {
		this->x_stabilizers[i][block] = inner;
	} else {
		this->z_stabilizers[i][block] = inner;
	}
}

uint32_t Qec::gaussian() {
	uint32_t i = this->n; // index to swap the next found row to
	for (uint32_t j = 0; j < this->n; j++)
	// Go over the x_stabilizers and check for any generators with an X below the triangle
	{
		uint32_t block = j >> BLOCK_BITS;
		uint64_t inner = powers[j & 63];
		uint32_t k = 0;
		for (k = this->n; k < 2 * this->n; k++) {
			if (this->x_stabilizers[k][block] & inner) {
				break;
			}
		}
		if (k < 2 * this->n) {
			// there was a generator X in any column
			this->rowswap(i, k); // move the row k to the first row
			this->rowswap(i - this->n, k - this->n); // swap row 0 with row k of the top part
			for (uint32_t a = i + 1; a < 2 * this->n; a++) {
				if (this->x_stabilizers[a][block] & inner) {
					this->rowmult(a, i);
					this->rowmult(i - this->n, a - this->n);
				}
			}
			i++;
		}
	}
	uint32_t g = i - this->n; // log2 of nonzero basis states
	// Do the same for the Z_stabilizers
	for (uint32_t j = 0; j < this->n; j++) {
		uint32_t block = j >> BLOCK_BITS;
		uint64_t inner = powers[j & 63];
		uint32_t k = 0;
		for (k = i; k < 2 * this->n; k++) {
			if (this->z_stabilizers[k][block] & inner) {
				break;
			}
		}
		if (k < 2 * this->n) {
			this->rowswap(i, k);
			this->rowswap(i - this->n, k - this->n);
			for (uint32_t a = i + 1; a < 2 * this->n; a++) {
				if (this->z_stabilizers[a][block] & inner) {
					this->rowmult(a, i);
					this->rowmult(i - this->n, a - this->n);
				}
			}
			i++;
		}
	}

	return g;
}

void Qec::seed(uint32_t log_amount) {
	// Wipe the scratch space
	this->phases[2 * this->n] = 0;
	std::fill(this->x_stabilizers[2 * this->n].begin(), this->x_stabilizers[2 * this->n].end(), 0);
	std::fill(this->z_stabilizers[2 * this->n].begin(), this->z_stabilizers[2 * this->n].end(), 0);

	uint32_t min;
	for (uint32_t i = 2 * this->n - 1; i >= this->n + log_amount; i--) {
		uint_fast8_t tmp_phase = this->phases[i];
		for (uint32_t j = this->n - i; j >= 0; j--) {
			uint32_t block = j >> BLOCK_BITS;
			uint64_t inner = powers[j & 63];
			if (this->z_stabilizers[i][block] & inner) {
				min = j;
				if (this->x_stabilizers[2 * this->n][block] & inner) {
					tmp_phase ^= 0b10;
				}
			}
		}
		if (tmp_phase == 2) {
			uint32_t block = min >> BLOCK_BITS;
			uint64_t inner = powers[min & 63];
			this->x_stabilizers[2 * this->n][block] ^= inner;
		}
	}
}
