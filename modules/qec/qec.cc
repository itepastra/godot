#include "qec.h"

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

	ClassDB::bind_method(D_METHOD("relax", "qubit"), &Qec::relax);

	ClassDB::bind_method(D_METHOD("peek_measurement_random", "qubits"), &Qec::peek_measure_random);

	ClassDB::bind_method(D_METHOD("get_entanglement_group", "seed"), &Qec::get_entanglement_group);
	ClassDB::bind_method(D_METHOD("snapshot_entanglement_group", "seed"), &Qec::snapshot_entanglement_group);
	ClassDB::bind_method(D_METHOD("restore_entanglement_group", "snapshot"), &Qec::restore_entanglement_group);
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

	uint8_t vop = this->nodes[control].vop;
	if (((this->nodes[control].adjacent.size() > 1) ||
				(this->nodes[control].adjacent.size() == 1 && this->nodes[control].adjacent[0] != target)) &&
			!(vop == ia || vop == za || vop == yb || vop == xb)) {
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
			c = this->nodes[a].adjacent[i];
			break;
		}
	}
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
}

void Qec::local_complementation(node_idx a) {
	std::vector<node_idx> neighbors = this->nodes[a].adjacent;

	for (auto i = neighbors.begin(); i != neighbors.end(); i++) {
		for (auto j = i; j != neighbors.end(); j++) {
			if (*i != *j) {
				toggle_edge(*i, *j);
			}
		}
		this->nodes[*i].vop = vop_table[this->nodes[*i].vop][6];
	}
	this->nodes[a].vop = vop_table[this->nodes[a].vop][14];
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

bool Qec::has_edge(node_idx a, node_idx b) const {
	for (node_idx v : this->nodes[a].adjacent) {
		if (v == b) {
			return true;
		}
	}
	return false;
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

uint8_t Qec::mx(node_idx node) {
	this->hadamard(node);
	uint8_t res = this->measure(node);
	this->hadamard(node);
	return res;
}
uint8_t Qec::my(node_idx node) {
	this->phase_dag(node);
	this->hadamard(node);
	uint8_t res = this->measure(node);
	this->hadamard(node);
	this->phase(node);
	return res;
}
uint8_t Qec::mz(node_idx node) {
	uint8_t res = this->measure(node);
	return res;
}

void Qec::relax(node_idx node) {
	this->measure(node);
	this->nodes[node].vop = 10;
}

uint8_t Qec::peek_determinism(node_idx node) {
	uint8_t nop = this->nodes[node].vop;
	uint8_t hnop = adj_tbl[nop];
	uint8_t bx = measurement_conj_table[xa - xa][hnop];
	uint8_t by = measurement_conj_table[ya - xa][hnop];
	uint8_t bz = measurement_conj_table[za - xa][hnop];

	if (bx == 1 && this->nodes[node].adjacent.size() == 0) {
		// measuring in the X basis is deterministic
		return 0;
	} else if (by == 1 && this->nodes[node].adjacent.size() == 0) {
		// measuring in the Y basis is deterministic
		return 1;
	} else if (bz == 1 && this->nodes[node].adjacent.size() == 0) {
		// measuring in the Z basis is deterministic
		return 2;
	} else {
		// no measurement is deterministic
		return 3;
	}
}

PackedByteArray Qec::peek_measure_random(PackedInt32Array meas_nodes) {
	std::vector<QecNode> before_state = this->nodes;
	PackedByteArray res;

	uint8_t dir = (rand() % 3) + 1;

	for (auto node = meas_nodes.begin(); node != meas_nodes.end(); ++node) {
		uint8_t det = this->peek_determinism(*node);
		uint8_t resp;

		switch (det) {
			case 0:
				this->measure(*node, xa);
				resp = this->nodes[*node].vop | 0b0100000;
				res.append(resp & 0b11111);
				break;
			case 1:
				this->measure(*node, ya);
				resp = this->nodes[*node].vop | 0b1000000;
				res.append(resp & 0b11111);
				break;
			case 2:
				this->measure(*node, za);
				resp = this->nodes[*node].vop | 0b1100000;
				res.append(resp & 0b11111);
				break;
			case 3:
				this->measure(*node, dir);
				resp = this->nodes[*node].vop | (dir << 5);
				res.append(resp & 0b11111);
				break;
		}
	}

	this->nodes = before_state;
	return res;
}

// measures a qubit in the Z basis
uint8_t Qec::measure(node_idx node, uint8_t basis) {
	uint8_t original_basis = basis;
	uint8_t zeta;
	uint8_t real_basis;

	{
		uint8_t nop = this->nodes[node].vop;
		if ((nop & 0b11) == 0 || (nop & 0b11) == original_basis) {
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

		real_basis = measurement_conj_table[original_basis - xa][adj_tbl[nop]];
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
	if (erase(this->nodes[i].adjacent, j)) {
		// target was in adjacency, also erase the other way
		erase(this->nodes[j].adjacent, i);
	} else {
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
		uint8_t right = vop_table[xa][yc];
		this->nodes[node].vop = vop_table[this->nodes[node].vop][right];
	} else {
		this->nodes[node].vop = vop_table[this->nodes[node].vop][yc];
	}
	return res;
}

const char *Qec::phase_to_str(uint8_t code) {
	switch (code) {
		case PLUS:
			return "+";
		case MINUS:
			return "-";
		case PLUS_I:
			return "+i";
		case MINUS_I:
			return "-i";
		case ONE:
			return "1";
		case ZERO:
			return "0";
	}
	return "?";
}

PackedInt32Array Qec::get_entanglement_group(node_idx seed) const {
	std::vector<node_idx> group;
	compute_entanglement_group_vec(seed, group);
	return pack_vector(group);
}

void Qec::compute_entanglement_group_vec(node_idx seed, std::vector<node_idx> &queue) const {
	const node_idx N = (node_idx)this->nodes.size();

	// reuse seen (thread-unsafe)
	bfs_seen_.assign(N, 0);

	queue.clear();
	queue.reserve(N); // worst case
	queue.push_back(seed);
	bfs_seen_[seed] = 1;

	for (size_t i = 0; i < queue.size(); ++i) {
		node_idx u = queue[i];
		const std::vector<node_idx> &adj = this->nodes[u].adjacent;
		for (node_idx v : adj) {
			if (!bfs_seen_[v]) {
				bfs_seen_[v] = 1;
				queue.push_back(v);
			}
		}
	}

	std::sort(queue.begin(), queue.end());
}

PackedInt32Array Qec::pack_vector(const std::vector<node_idx> &v) {
	PackedInt32Array out;
	const int n = (int)v.size();
	out.resize(n);
	for (int i = 0; i < n; ++i) {
		out.set(i, (int32_t)v[i]);
	}
	return out;
}

Dictionary Qec::snapshot_entanglement_group(node_idx seed) const {
	Dictionary snap;

	// nodes in component
	PackedInt32Array group = get_entanglement_group(seed);
	const int M = group.size();

	// map: node_id -> compact index [0..M)
	// (vector of size total_n works, -1 means not in group)
	std::vector<int> idx_map(this->nodes.size(), -1);
	for (int i = 0; i < M; ++i) {
		idx_map[(node_idx)group[i]] = i;
	}

	// capture vops in same order
	PackedByteArray vops;
	vops.resize(M);
	for (int i = 0; i < M; ++i) {
		node_idx u = (node_idx)group[i];
		vops.set(i, this->nodes[u].vop);
	}

	// capture edges internal to the group (store as list of pairs u,v with u < v)
	PackedInt32Array edges; // flattened [u0,v0,u1,v1,...] using absolute node indices
	std::vector<std::pair<node_idx, node_idx>> pairs;
	for (int i = 0; i < M; ++i) {
		node_idx u = (node_idx)group[i];
		for (node_idx v : this->nodes[u].adjacent) {
			if (idx_map[v] >= 0 && u < v) {
				pairs.emplace_back(u, v);
			}
		}
	}
	edges.resize((int)(pairs.size() * 2));
	for (int i = 0; i < (int)pairs.size(); ++i) {
		edges.set(2 * i + 0, pairs[i].first);
		edges.set(2 * i + 1, pairs[i].second);
	}

	snap["nodes"] = group;
	snap["vops"] = vops;
	snap["edges"] = edges;
	return snap;
}

void Qec::restore_entanglement_group(const Dictionary &snapshot) {
	ERR_FAIL_COND_MSG(!snapshot.has("nodes") || !snapshot.has("vops") || !snapshot.has("edges"),
			"restore_entanglement_group: snapshot missing fields");

	PackedInt32Array group = snapshot["nodes"];
	PackedByteArray vops = snapshot["vops"];
	PackedInt32Array edges = snapshot["edges"];

	const int M = group.size();
	ERR_FAIL_COND_MSG(vops.size() != M, "restore_entanglement_group: vops size mismatch");

	// build membership map for fast checks
	std::vector<uint8_t> in_group(this->nodes.size(), 0);
	for (int i = 0; i < M; ++i) {
		in_group[(node_idx)group[i]] = 1;
	}

	// restore vops
	for (int i = 0; i < M; ++i) {
		node_idx u = (node_idx)group[i];
		this->nodes[u].vop = vops[i];
	}

	// remove all current edges incident to group nodes (both inside and to outside),
	// so we can add back exactly what the snapshot had
	//  must iterate over a copy because we mutate adjacency
	for (int i = 0; i < M; ++i) {
		node_idx u = (node_idx)group[i];
		std::vector<node_idx> current = this->nodes[u].adjacent; // copy
		for (node_idx v : current) {
			// remove edge (u,v) once, erase_connection does symmetric remove
			this->erase_connection(u, v);
		}
	}

	// add back internal edges from snapshot
	ERR_FAIL_COND_MSG(edges.size() % 2 != 0, "restore_entanglement_group: edges length not even");
	for (int i = 0; i < edges.size(); i += 2) {
		node_idx a = (node_idx)edges[i + 0];
		node_idx b = (node_idx)edges[i + 1];
		// only add if both endpoints are in the group
		if (!in_group[a] || !in_group[b]) {
			continue;
		}

		// add edge (a,b) if it's not already present
		if (!this->has_edge(a, b)) {
			this->nodes[a].adjacent.push_back(b);
			this->nodes[b].adjacent.push_back(a);
		}
	}
}
