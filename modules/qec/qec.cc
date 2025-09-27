#include "qec.h"
#include <vector>

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
}

void Qec::cnot(size_t control, size_t target) {
	if (!this->initialized) {
		return;
	}
	size_t control_block = control >> BLOCK_BITS;
	size_t target_block = target >> BLOCK_BITS;
	uint64_t control_inner = powers[control & 63];
	uint64_t target_inner = powers[target & 63];

	for (size_t i = 0; i < 2 * this->n; i++) {
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
			this->phases[i] = (this->phases[i] + 2) % 4;
		}
		// do the phase kickback in another case TODO: figure out what case this represents
		if ((this->x_stabilizers[i][control_block] & control_inner) &&
				(this->z_stabilizers[i][target_block] & target_inner) &&
				!(this->x_stabilizers[i][target_block] & target_inner) &&
				!(this->z_stabilizers[i][control_block] & control_inner)) {
			this->phases[i] = (this->phases[i] + 2) % 4;
		}
	}
}
void Qec::cphase(size_t control, size_t target) {
	if (!this->initialized) {
		return;
	}
	this->hadamard(target);
	this->cnot(control, target);
	this->hadamard(target);
} // CZ = H_target CNOT H_target
void Qec::hadamard(size_t target) {
	if (!this->initialized) {
		return;
	}
	size_t block = target >> BLOCK_BITS;
	uint64_t inner = powers[target & 63];
	uint64_t tmp;

	for (size_t i = 0; i < 2 * this->n; i++) {
		tmp = this->x_stabilizers[i][block];
		this->x_stabilizers[i][block] ^= (this->x_stabilizers[i][block] ^ this->z_stabilizers[i][block]) & inner;
		this->z_stabilizers[i][block] ^= (this->z_stabilizers[i][block] ^ tmp) & inner;
		if ((this->x_stabilizers[i][block] & inner) && (this->z_stabilizers[i][block] & inner)) {
			this->phases[i] = (this->phases[i] + 2) % 4;
		}
	}
} // H
void Qec::phase(size_t target) {
	if (!this->initialized) {
		return;
	}
	size_t block = target >> BLOCK_BITS;
	uint64_t inner = powers[target & 63];
	for (size_t i = 0; i < 2 * this->n; i++) {
		if ((this->x_stabilizers[i][block] & inner) && (this->z_stabilizers[i][block] & inner)) {
			this->phases[i] = (this->phases[i] + 2) % 4;
		}
		this->z_stabilizers[i][block] ^= this->x_stabilizers[i][block] & inner;
	}
} // S
void Qec::phase_dag(size_t target) {
	if (!this->initialized) {
		return;
	}
	this->phase(target);
	this->phase(target);
	this->phase(target);
} // S^+ = S S S
void Qec::xgate(size_t target) {
	if (!this->initialized) {
		return;
	}
	this->hadamard(target);
	this->phase(target);
	this->phase(target);
	this->hadamard(target);
} // X = H S S H
void Qec::ygate(size_t target) {
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
void Qec::zgate(size_t target) {
	if (!this->initialized) {
		return;
	}
	this->phase(target);
	this->phase(target);
} // Z = S S

void Qec::init(size_t qubit_amount) {
	if (this->initialized) {
		return;
	}
	size_t whole_blocks = (qubit_amount >> BLOCK_BITS) + 1; // shrink by amount of qubits that fit in 64 bits
	size_t rows = 2 * qubit_amount + 1;
	this->n = qubit_amount;

	for (size_t i = 0; i < rows; i++) {
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
		this->phases.push_back(0);
	}

	this->initialized = true;
}

Qec::Qec() {
	this->initialized = false;
}
