#include "qec.h"
#include <cstdint>
#include <cstdlib>

void Qec::_bind_methods() {
	ClassDB::bind_method(D_METHOD("init", "qubit_amount"), &Qec::init);
	//
	// ClassDB::bind_method(D_METHOD("cnot", "control", "target"), &Qec::cnot);
	// ClassDB::bind_method(D_METHOD("cphase", "control", "target"), &Qec::cphase);
	// ClassDB::bind_method(D_METHOD("hadamard", "qubit"), &Qec::hadamard);
	// ClassDB::bind_method(D_METHOD("phase", "qubit"), &Qec::phase);
	// ClassDB::bind_method(D_METHOD("phase_dag", "qubit"), &Qec::phase_dag);
	// ClassDB::bind_method(D_METHOD("xgate", "qubit"), &Qec::xgate);
	// ClassDB::bind_method(D_METHOD("ygate", "qubit"), &Qec::ygate);
	// ClassDB::bind_method(D_METHOD("zgate", "qubit"), &Qec::zgate);
	//
	// ClassDB::bind_method(D_METHOD("measure", "qubit"), &Qec::measure);
	//
	// ClassDB::bind_method(D_METHOD("get_x_stab", "i", "j"), &Qec::get_x_stab);
	// ClassDB::bind_method(D_METHOD("get_z_stab", "i", "j"), &Qec::get_z_stab);
	// ClassDB::bind_method(D_METHOD("x_stabs"), &Qec::x_stabs);
	// ClassDB::bind_method(D_METHOD("z_stabs"), &Qec::z_stabs);
	// ClassDB::bind_method(D_METHOD("get_phase", "i"), &Qec::get_phase);
}

Qec::Qec() {
	this->initialized = false;
}

void Qec::init(node_idx qubit_amount) {
	this->nodes.resize(qubit_amount);
	for (node_idx i = 0; i < qubit_amount; i++) {
		this->hadamard(i); // run hadamard on all qubits to generate the |0...0> state
	}
}

void Qec::hadamard(node_idx target) {
	this->nodes[target].vop = vop_table[yc][this->nodes[target].vop];
}

void Qec::phase(node_idx target) {
	this->nodes[target].vop = vop_table[xb][this->nodes[target].vop];
}

void Qec::phase_dag(node_idx target) {
	this->nodes[target].vop = vop_table[yb][this->nodes[target].vop];
}

void Qec::xgate(node_idx target) {
	this->nodes[target].vop = vop_table[xa][this->nodes[target].vop];
}

void Qec::ygate(node_idx target) {
	this->nodes[target].vop = vop_table[ya][this->nodes[target].vop];
}

void Qec::zgate(node_idx target) {
	this->nodes[target].vop = vop_table[za][this->nodes[target].vop];
}

void Qec::cphase(node_idx control, node_idx target) {
	if ((this->nodes[control].adjacent.size() > 1) ||
			(this->nodes[control].adjacent.size() == 1 && this->nodes[control].adjacent[0] != target)) {
		// remove_VOP(control,target)
	}
	if ((this->nodes[target].adjacent.size() > 1) ||
			(this->nodes[target].adjacent.size() == 1 && this->nodes[target].adjacent[0] != control)) {
		// remove_VOP(target,control)
	}

	if ((this->nodes[control].adjacent.size() > 1) ||
			(this->nodes[control].adjacent.size() == 1 && this->nodes[control].adjacent[0] != target)) {
		// remove_VOP(control,target)
	}

	return;

	// Case 1: both control and target are in {I,Z,S,S^}
	// toggle the adjacency
	{
		if (erase(this->nodes[control].adjacent, target)) {
			// target was in adjacency, also erase the other way
			erase(this->nodes[target].adjacent, control);
		} else {
			this->nodes[control].adjacent.push_back(target);
			this->nodes[target].adjacent.push_back(control);
		}
	}
	// Case 2: at least one of control and target is not in there
}

void Qec::remove_VOP(node_idx a, node_idx b) {
	node_idx c = b;
	for (uint32_t i = 0; i < this->nodes[a].adjacent.size(); ++i) {
		if (this->nodes[a].adjacent[i] != b) {
			c = i;
			break;
		}
	}
	std::vector<uint8_t> decomp = decompositions[this->nodes[a].vop];

	for (uint32_t i = decomp.size(); i > 0; --i) {
		if (decomp[i] == 0) { // decomp[i] == sqrt(-iX)
			this->local_complementation(a);
		} else { // decomp[i] == sqrt(iZ)
			this->local_complementation(b);
		}
	}
}

void Qec::local_complementation(node_idx a) {
	for (uint32_t i = 0; i < this->nodes[a].adjacent.size(); i++) {
		for (uint32_t j = 0; j < this->nodes[a].adjacent.size(); j++) {
			if (i < j) {
				if (erase(this->nodes[i].adjacent, j)) {
					// target was in adjacency, also erase the other way
					erase(this->nodes[j].adjacent, i);
				} else {
					this->nodes[i].adjacent.push_back(j);
					this->nodes[j].adjacent.push_back(i);
				}
			}
		}
		this->nodes[i].vop = vop_table[this->nodes[i].vop][xb];
		this->nodes[a].vop = vop_table[this->nodes[a].vop][zd];
	}
}
