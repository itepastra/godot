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
	//
	// ClassDB::bind_method(D_METHOD("measure", "qubit"), &Qec::measure);
	//
	ClassDB::bind_method(D_METHOD("get_vop", "node"), &Qec::get_vop);
	ClassDB::bind_method(D_METHOD("get_adjacent", "node"), &Qec::get_adjacent);
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
		remove_VOP(control, target);
	}
	if ((this->nodes[target].adjacent.size() > 1) ||
			(this->nodes[target].adjacent.size() == 1 && this->nodes[target].adjacent[0] != control)) {
		remove_VOP(target, control);
	}

	if ((this->nodes[control].adjacent.size() > 1) ||
			(this->nodes[control].adjacent.size() == 1 && this->nodes[control].adjacent[0] != target)) {
		remove_VOP(control, target);
	}

	uint8_t controlvop = this->nodes[control].vop;
	uint8_t targetvop = this->nodes[target].vop;
	bool had_edge = std::find(this->nodes[control].adjacent.begin(),
							this->nodes[control].adjacent.end(), target) !=
			this->nodes[control].adjacent.end();

	if (cphase_table[had_edge][controlvop][targetvop][0]) {
		this->nodes[control].adjacent.push_back(target);
		this->nodes[target].adjacent.push_back(control);
	} else {
		erase(this->nodes[control].adjacent, target);
		erase(this->nodes[target].adjacent, control);
	}

	this->nodes[control].vop = cphase_table[had_edge][controlvop][targetvop][1];
	this->nodes[target].vop = cphase_table[had_edge][controlvop][targetvop][2];
}

void Qec::cnot(node_idx control, node_idx target) {
	this->hadamard(target);
	this->cphase(control, target);
	this->hadamard(target);
}

void Qec::remove_VOP(node_idx a, node_idx b) {
	if (this->nodes[a].adjacent.size() == 0) {
		// This should never be called without any adjacent on the a side
		abort();
	}

	// if necessary we'll use b, but otherwise try using the first other adjacency
	node_idx c = b;
	for (uint32_t i = 0; i < this->nodes[a].adjacent.size(); i++) {
		if (this->nodes[a].adjacent[i] != b) {
			c = i;
			break;
		}
	}
	std::vector<uint8_t> decomp = decompositions[this->nodes[a].vop];

	for (uint32_t i = decomp.size(); i > 0; --i) {
		if (decomp[i] == 0) { // 0 == U
			this->local_complementation(a);
		} else { // 1 == V
			this->local_complementation(c);
		}
	}
}

void Qec::local_complementation(node_idx a) {
	uint32_t size = this->nodes[a].adjacent.size();
	for (uint32_t i = 0; i < size; i++) {
		for (uint32_t j = i + 1; j < size; j++) {
			if (erase(this->nodes[i].adjacent, j)) {
				// target was in adjacency, also erase the other way
				erase(this->nodes[j].adjacent, i);
			} else {
				this->nodes[i].adjacent.push_back(j);
				this->nodes[j].adjacent.push_back(i);
			}
		}
		this->nodes[i].vop = vop_table[this->nodes[i].vop][yb];
	}
	this->nodes[a].vop = vop_table[this->nodes[a].vop][yd];
}

uint8_t Qec::get_vop(node_idx node) {
	return this->nodes[node].vop;
}

PackedInt32Array Qec::get_adjacent(node_idx node) {
	PackedInt32Array array;
	for (uint32_t i = 0; i < this->nodes[node].adjacent.size(); i++) {
		array.append(this->nodes[node].adjacent[i]);
	};
	return array;
}
