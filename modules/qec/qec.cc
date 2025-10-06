#include "qec.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <utility>
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
	//
	ClassDB::bind_method(D_METHOD("mx", "qubit"), &Qec::mx);
	ClassDB::bind_method(D_METHOD("my", "qubit"), &Qec::my);
	ClassDB::bind_method(D_METHOD("mz", "qubit"), &Qec::mz);
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
	printf("performing hadamard on %d, went from vop %d", target, this->nodes[target].vop);
	this->nodes[target].vop = vop_table[yc][this->nodes[target].vop];
	printf(" to vop %d\n", this->nodes[target].vop);
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
	printf("before cphase between %d and %d, the vops are %d and %d\n", control, target, this->nodes[control].vop, this->nodes[target].vop);
	if ((this->nodes[control].adjacent.size() > 1) ||
			(this->nodes[control].adjacent.size() == 1 && this->nodes[control].adjacent[0] != target)) {
		printf("removing Vertex Operator of %d, avoiding %d\n", control, target);
		remove_VOP(control, target);
		printf("the vops for control (%d) and target (%d) are %d and %d\n", control, target, this->nodes[control].vop, this->nodes[target].vop);
	}
	if ((this->nodes[target].adjacent.size() > 1) ||
			(this->nodes[target].adjacent.size() == 1 && this->nodes[target].adjacent[0] != control)) {
		printf("removing Vertex Operator of %d, avoiding %d\n", target, control);
		remove_VOP(target, control);
		printf("the vops for control (%d) and target (%d) are %d and %d\n", control, target, this->nodes[control].vop, this->nodes[target].vop);
	}

	uint8_t vop = this->nodes[control].vop;
	if (((this->nodes[control].adjacent.size() > 1) ||
				(this->nodes[control].adjacent.size() == 1 && this->nodes[control].adjacent[0] != target)) &&
			!(vop == ia || vop == za || vop == yb || vop == xb)) {
		printf("removing Vertex Operator of %d, avoiding %d\n", control, target);
		remove_VOP(control, target);
		printf("the vops for control (%d) and target (%d) are %d and %d\n", control, target, this->nodes[control].vop, this->nodes[target].vop);
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

	printf("looking in cphase_table for %b, %d, %d\n", had_edge, controlvop, targetvop);
	this->nodes[control].vop = cphase_table[had_edge][controlvop][targetvop][1];
	this->nodes[target].vop = cphase_table[had_edge][controlvop][targetvop][2];
	printf("after cphase the vops for control (%d) and target (%d) are %d and %d\n", control, target, this->nodes[control].vop, this->nodes[target].vop);
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
			c = this->nodes[a].adjacent[i];
			break;
		}
	}
	printf("ended up choosing neighbor %d\n", c);
	std::vector<uint8_t> decomp = decompositions[this->nodes[a].vop];

	for (uint32_t i = decomp.size() - 1; i >= 0; i--) {
		if (i > 10) { // There is a very weird bug I feel, but maybe it has a reason, anyways, this fixes it
			break;
		}
		if (decomp[i] == 0) { // 0 == U
			this->local_complementation(a);
		} else { // 1 == V
			this->local_complementation(c);
		}
	}

	printf("the vops of a, c after removal are %d, %d\n", this->nodes[a].vop, this->nodes[c].vop);
}

void Qec::local_complementation(node_idx a) {
	std::vector<node_idx> neighbors = this->nodes[a].adjacent;

	printf("locally complementing %d\n", a);
	for (auto i = neighbors.begin(); i != neighbors.end(); i++) {
		for (auto j = i; j != neighbors.end(); j++) {
			if (*i != *j) {
				printf("toggling edge between %d and %d", *i, *j);
				toggle_edge(*i, *j);
			}
		}
		this->nodes[*i].vop = vop_table[this->nodes[*i].vop][6];
		printf("changed %d vop to %d\n", *i, this->nodes[*i].vop);
	}
	this->nodes[a].vop = vop_table[this->nodes[a].vop][14];
	printf("changed %d vop to %d\n", a, this->nodes[a].vop);
}

uint8_t rand_bool() {
	return rand() & 0b1;
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

void Qec::erase_connection(node_idx a, node_idx b) {
	erase(this->nodes[a].adjacent, b);
	erase(this->nodes[b].adjacent, a);
}

struct edge_hash {
	static uint32_t hash(const std::pair<unsigned int, unsigned int> &p_key) {
		return uint32_t(p_key.first ^ (p_key.second << 16));
	}
};

uint8_t mx(node_idx node) {
	this->hadamard(node);
	uint8_t res = this->measure(node);
	this->hadamard(node);
	return res;
}
uint8_t my(node_idx node) {
	this->phase_dag(node);
	this->hadamard(node);
	uint8_t res = this->measure(node);
	this->hadamard(node);
	this->phase(node);
	return res;
}
uint8_t mz(node_idx node) {
	uint8_t res = this->measure(node);
	return res;
}

// measures a qubit in the Z basis
uint8_t Qec::measure(node_idx node) {
	uint8_t original_basis = za;
	uint8_t zeta;
	uint8_t real_basis;

	{
		uint8_t nop = this->nodes[node].vop;
		if ((nop & 0x03) == 0 || (nop & 0x03) == za) {
			if (nop >= 4 && nop <= 15) {
				zeta = 2;
			} else {
				zeta = 0;
			}
		} else {
			if (nop >= 4 && nop <= 15) {
				zeta = 0;
			} else {
				zeta = 2;
			}
		}

		real_basis = measurement_conj_table[original_basis - xa][nop];
	}

	uint8_t res;
	switch (real_basis) {
		case 1:
			res = this->measure_x(node);
			break;
		case 2:
			res = this->measure_y(node);
			break;
		case 3:
			res = this->measure_z(node);
			break;
		default:
			abort();
	}
	if (zeta == 2) {
		res ^= 0b01;
	}
	return res;
}

// collapses the graph to what would happen if a measurement in X happened
uint8_t Qec::measure_x(node_idx node) {
	if (this->nodes[node].adjacent.size() == 0) {
		return 0b10;
	}

	uint8_t res = rand_bool();

	node_idx other = this->nodes[node].adjacent[0];

	std::vector<node_idx> node_neighbors = this->nodes[node].adjacent;
	std::vector<node_idx> other_neighbors = this->nodes[other].adjacent;

	if (res) {
		// measured a |->
		this->nodes[other].vop = vop_table[this->nodes[other].vop][xc];
		this->nodes[node].vop = vop_table[this->nodes[node].vop][za];

		for (uint32_t i = 0; i < this->nodes[other].adjacent.size(); i++) {
			node_idx waa = this->nodes[other].adjacent[i];
			if (waa != node && !contains(this->nodes[node].adjacent, waa)) {
				this->nodes[waa].vop = vop_table[this->nodes[waa].vop][za];
			}
		}
	} else {
		// measured a |+>
		this->nodes[other].vop = vop_table[this->nodes[other].vop][zc];

		for (uint32_t i = 0; i < this->nodes[node].adjacent.size(); i++) {
			node_idx waa = this->nodes[node].adjacent[i];
			if (waa != other && !contains(this->nodes[other].adjacent, waa)) {
				this->nodes[waa].vop = vop_table[this->nodes[waa].vop][za];
			}
		}
	}

	{
		HashSet<std::pair<node_idx, node_idx>, edge_hash> procd_edges;
		for (auto i = node_neighbors.begin(); i != node_neighbors.end(); i++) {
			for (auto j = other_neighbors.begin(); j != other_neighbors.end(); j++) {
				if ((*i != *j) && (procd_edges.find(std::pair(*i, *j)) == procd_edges.end())) {
					procd_edges.insert(std::pair(*i, *j));
					this->toggle_edge(*i, *j);
				}
			}
		}
	}

	std::vector<node_idx> intersection;
	for (uint32_t i = 0; i < node_neighbors.size(); i++) {
		if (contains(other_neighbors, node_neighbors[i])) {
			intersection.push_back(i);
		}
	}

	for (uint32_t i = 0; i < intersection.size(); i++) {
		for (uint32_t j = i + 1; j < intersection.size(); j++) {
			toggle_edge(intersection[i], intersection[j]);
		}
	}

	for (uint32_t i = 0; i < node_neighbors.size(); i++) {
		if (node_neighbors[i] != other) {
			toggle_edge(other, node_neighbors[i]);
		}
	}

	return res;
}

void Qec::toggle_edge(node_idx i, node_idx j) {
	printf("toggling edge between %d and %d...", i, j);
	if (erase(this->nodes[i].adjacent, j)) {
		printf(" removing\n");
		// target was in adjacency, also erase the other way
		erase(this->nodes[j].adjacent, i);
	} else {
		printf(" adding\n");
		this->nodes[i].adjacent.push_back(j);
		this->nodes[j].adjacent.push_back(i);
	}
}

bool contains(std::vector<node_idx> vec, node_idx val) {
	for (uint32_t i = 0; i < vec.size(); i++) {
		if (vec[i] == val) {
			return true;
		}
	}
	return false;
}

// collapses the graph to what would happen if a measurement in Y happened
uint8_t Qec::measure_y(node_idx node) {
	uint8_t res = rand_bool();
	std::vector<node_idx> neighbors = this->nodes[node].adjacent;
	for (uint32_t i = 0; i < neighbors.size(); i++) {
		if (res) {
			this->nodes[neighbors[i]].vop = vop_table[this->nodes[neighbors[i]].vop][xb]; // S
		} else {
			this->nodes[neighbors[i]].vop = vop_table[this->nodes[neighbors[i]].vop][yb]; // Sdag
		}
	}
	neighbors.push_back(node);
	for (uint32_t i = 0; i < neighbors.size(); i++) {
		for (uint32_t j = i + 1; j < neighbors.size(); j++) {
			if (erase(this->nodes[neighbors[i]].adjacent, neighbors[j])) {
				// target was in adjacency, also erase the other way
				erase(this->nodes[neighbors[j]].adjacent, neighbors[i]);
			} else {
				this->nodes[neighbors[i]].adjacent.push_back(neighbors[j]);
				this->nodes[neighbors[j]].adjacent.push_back(neighbors[i]);
			}
		}
	}
	if (res) {
		this->nodes[node].vop = vop_table[this->nodes[node].vop][yb];
	} else {
		this->nodes[node].vop = vop_table[this->nodes[node].vop][xb];
	}
	return res;
}

// collapses the graph to what would happen if a measurement in Z happened
uint8_t Qec::measure_z(node_idx node) {
	uint8_t res = rand_bool(); // start with a random bit
	std::vector<node_idx> neighbors = this->nodes[node].adjacent;
	// remove entanglement with all neighbors depending on res
	for (uint32_t i = 0; i < neighbors.size(); i++) {
		this->erase_connection(node, neighbors[i]);
		if (res) {
			this->nodes[neighbors[i]].vop = vop_table[this->nodes[neighbors[i]].vop][za];
		}
	}
	if (res) {
		this->nodes[node].vop = vop_table[this->nodes[node].vop][yc];
	} else {
		this->nodes[node].vop = vop_table[vop_table[this->nodes[node].vop][xa]][yc];
	}
	return res;
}
