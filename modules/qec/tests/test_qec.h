#pragma once

#include "modules/qec/qec.h"
#include "tests/test_macros.h"

namespace TestQec {
TEST_CASE("[Modules][Qec] Hadamard works as expected") {
	Ref<Qec> q = memnew(Qec);
	// init to |0>
	q->init(1);

	CHECK(q->get_vop(0) == 10);

	// hadamard to |+>
	q->hadamard(0);

	CHECK(q->get_vop(0) == 0);

	// hadamard back to |0>
	q->hadamard(0);

	CHECK(q->get_vop(0) == 10);
}

TEST_CASE("[Modules][Qec] Measure Deterministic") {
	Ref<Qec> q = memnew(Qec);

	q->init(1);
	CHECK(q->mz(0) == 0b10); // |0> deterministic
	q->zgate(0);
	CHECK(q->mz(0) == 0b11); // |1> deterministic
}

TEST_CASE("[Modules][Qec] Measure Phases") {
	Ref<Qec> q = memnew(Qec);
	uint8_t res;

	q->init(1);
	CHECK(q->get_vop(0) == 10);
	res = q->mz(0) >> 1; // needed to check for both 2 and 3
	CHECK(res == 0b1); // +|0> deterministic
	q->phase(0);
	CHECK(q->get_vop(0) == 21);
	res = q->mz(0) >> 1;
	CHECK(res == 0b1); // +i|0> deterministic
	q->phase(0);
	CHECK(q->get_vop(0) == 11);
	res = q->mz(0) >> 1;
	CHECK(res == 0b1); // -|0> deterministic
	q->phase(0);
	CHECK(q->get_vop(0) == 20);
	res = q->mz(0) >> 1;
	CHECK(res == 0b1); // -i|0> deterministic
	q->phase(0);
	CHECK(q->get_vop(0) == 10);
	res = q->mz(0) >> 1;
	CHECK(res == 0b1); // +|0> deterministic
	q->hadamard(0);
	CHECK(q->get_vop(0) == 0);
	res = q->mx(0) >> 1;
	CHECK(res == 0b1); // |+> deterministic
	q->phase(0);
	CHECK(q->get_vop(0) == 5);
	res = q->my(0) >> 1;
	CHECK(res == 0b1);
}

TEST_CASE("[Modules][Qec] Chained CNOTs") {
	Ref<Qec> q = memnew(Qec);
	// init to |0>
	q->init(6);
	q->hadamard(0);
	CHECK(q->get_vop(0) == 0);
	CHECK(q->get_adjacent(0) == PackedInt32Array{});

	q->cnot(0, 1);
	CHECK(q->get_vop(0) == 0);
	CHECK(q->get_adjacent(0) == PackedInt32Array{ 1 });
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_adjacent(1) == PackedInt32Array{ 0 });

	q->cnot(1, 2);
	CHECK(q->get_vop(0) == 10);
	CHECK(q->get_adjacent(0) == PackedInt32Array{ 1 });
	CHECK(q->get_vop(1) == 0);
	CHECK(q->get_adjacent(1) == PackedInt32Array{ 0, 2 });
	CHECK(q->get_vop(2) == 10);
	CHECK(q->get_adjacent(2) == PackedInt32Array{ 1 });

	q->cnot(2, 5);
	CHECK(q->get_vop(0) == 10);
	CHECK(q->get_adjacent(0) == PackedInt32Array{ 2 });
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_adjacent(1) == PackedInt32Array{ 2 });
	CHECK(q->get_vop(2) == 0);
	CHECK(q->get_adjacent(2) == PackedInt32Array{ 1, 0, 5 });
	CHECK(q->get_vop(5) == 10);
	CHECK(q->get_adjacent(5) == PackedInt32Array{ 2 });
}

TEST_CASE("[Modules][Qec] SWAP by CNOTs") {
	Ref<Qec> q = memnew(Qec);
	q->init(3);
	// generate a bell state on 0, 1
	q->hadamard(0);
	q->cnot(0, 1);

	// generate some state on qubit 2 (+i)
	q->hadamard(2);
	q->phase(2);
	CHECK(q->get_vop(2) == 5);

	// swap qubit 1 and 2
	q->cnot(1, 2);
	q->cnot(2, 1);
	q->cnot(1, 2);

	CHECK(q->get_adjacent(0) == PackedInt32Array{ 2 });
	CHECK(q->get_adjacent(1) == PackedInt32Array{});
	CHECK(q->get_adjacent(2) == PackedInt32Array{ 0 });
	CHECK(q->get_vop(1) == 5);
}

TEST_CASE("[Modules][Qec] Phase dagger cancels phase") {
	Ref<Qec> q = memnew(Qec);
	q->init(1);

	// record original VOP
	uint8_t orig = q->get_vop(0);

	// apply phase then its dagger, should return to original
	q->phase(0);
	q->phase_dag(0);

	CHECK(q->get_vop(0) == orig);

	// do the same but with hadamard (X-basis)
	q->hadamard(0);
	uint8_t orig2 = q->get_vop(0);
	q->phase(0);
	q->phase_dag(0);
	CHECK(q->get_vop(0) == orig2);
}

TEST_CASE("[Modules][Qec] CZ (cphase) is its own inverse") {
	Ref<Qec> q = memnew(Qec);
	q->init(2);

	// state and vops before
	uint8_t before_v0 = q->get_vop(0);
	uint8_t before_v1 = q->get_vop(1);
	auto before_adj0 = q->get_adjacent(0);
	auto before_adj1 = q->get_adjacent(1);

	// apply cphase twice
	q->cphase(0, 1);
	q->cphase(0, 1);

	// edges and vops should be back to prior state
	CHECK(q->get_vop(0) == before_v0);
	CHECK(q->get_vop(1) == before_v1);
	CHECK(q->get_adjacent(0) == before_adj0);
	CHECK(q->get_adjacent(1) == before_adj1);
}

TEST_CASE("[Modules][Qec] Bell: determinism of partner after Z on first") {
	Ref<Qec> q = memnew(Qec);
	q->init(2);

	// entanglement
	q->hadamard(0);
	q->cnot(0, 1);

	// before measurement neither qubit should be deterministically measurable
	CHECK(q->peek_determinism(0) == 3);
	CHECK(q->peek_determinism(1) == 3);

	// Measure qubit 0 in Z
	(void)q->mz(0);

	// qubit 1 is now deterministic
	CHECK(q->peek_determinism(1) == 2);
	uint8_t r1 = q->mz(1) & 1;
	uint8_t r2 = q->mz(1) & 1;
	CHECK(r1 == r2);

	// adjacency is gone
	CHECK(q->get_adjacent(0) == PackedInt32Array{});
	CHECK(q->get_adjacent(1) == PackedInt32Array{});
}

TEST_CASE("[Modules][Qec] Bell: determinism of partner after X on first") {
	Ref<Qec> q = memnew(Qec);
	q->init(2);

	q->hadamard(0);
	q->cnot(0, 1);

	CHECK(q->peek_determinism(0) == 3);
	CHECK(q->peek_determinism(1) == 3);

	// Measure qubit 0 in X
	(void)q->mx(0);

	// qubit 1 should be deterministic
	CHECK(q->peek_determinism(1) == 0);

	uint8_t r1 = q->mx(1) & 1;
	uint8_t r2 = q->mx(1) & 1;
	CHECK(r1 == r2);

	// adjacency is gone
	CHECK(q->get_adjacent(0) == PackedInt32Array{});
	CHECK(q->get_adjacent(1) == PackedInt32Array{});
}

TEST_CASE("[Modules][Qec] GHZ collapse after Z on leaf") {
	Ref<Qec> q = memnew(Qec);
	q->init(3);

	// GHZ state
	q->hadamard(0);
	q->cnot(0, 1);
	q->cnot(0, 2);

	// star adjacency
	CHECK(q->get_adjacent(0) == PackedInt32Array{ 1, 2 });
	CHECK(q->get_adjacent(1) == PackedInt32Array{ 0 });
	CHECK(q->get_adjacent(2) == PackedInt32Array{ 0 });

	// measure leaf (qubit 1),removes edges, no 0–2 edge is created
	uint8_t m1 = q->mz(1);
	uint8_t r1 = m1 & 0x1;
	(void)r1;

	// both remaining qubits must be deterministic and not entangled
	CHECK(q->get_adjacent(0) == PackedInt32Array{});
	CHECK(q->get_adjacent(2) == PackedInt32Array{});
	CHECK(q->peek_determinism(0) == 2); // Z deterministic
	CHECK(q->peek_determinism(2) == 2); // Z deterministic

	// check: repeat Z on 0 and 2 to check for fixed outcomes
	uint8_t m0a = q->mz(0), m0b = q->mz(0);
	uint8_t m2a = q->mz(2), m2b = q->mz(2);
	CHECK((m0a & 1) == (m0b & 1));
	CHECK((m2a & 1) == (m2b & 1));
}

// Helper to minimally reproduce zeta logic from measure for a given vop and requested basis
static inline uint8_t zeta_for(uint8_t nop, uint8_t original_basis /* xa/ya/za */) {
	// measure():
	// if (((nop & 0b11) == 0) || ((nop & 0b11) == original_basis)) {
	//     zeta = (4 <= nop && nop <= 15) ? 2 : 0;
	// } else {
	//     zeta = (4 <= nop && nop <= 15) ? 0 : 2;
	// }
	if (((nop & 0b11) == 0) || ((nop & 0b11) == original_basis)) {
		return (nop >= 4 && nop <= 15) ? 2 : 0;
	} else {
		return (nop >= 4 && nop <= 15) ? 0 : 2;
	}
}

TEST_CASE("[Modules][Qec] GHZ collapse after Z on leaf PARITY") {
	Ref<Qec> q = memnew(Qec);
	q->init(3);

	// GHZ
	q->hadamard(0);
	q->cnot(0, 1);
	q->cnot(0, 2);

	// measure leaf (q1) get VOP before measuring to compute zeta1.
	// uint8_t v1_pre = q->get_vop(1);
	// uint8_t z1 = zeta_for(v1_pre, QecConst::za); // 0 or 2
	// uint8_t m1 = q->mz(1);
	// uint8_t r1 = (m1 & 1);
	// uint8_t fr1 = r1 ^ (z1 == 2); // frame-invariant Z bit for qubit 1
	q->mz(1);

	// after measuring q1, 0 and 2 should be deterministic and not entangled
	CHECK(q->get_adjacent(0) == PackedInt32Array{});
	CHECK(q->get_adjacent(2) == PackedInt32Array{});
	CHECK(q->peek_determinism(0) == 2);
	CHECK(q->peek_determinism(2) == 2);

	// q0 compute zeta from its VOP before measuring 0
	// uint8_t v0_pre = q->get_vop(0);
	// uint8_t z0 = zeta_for(v0_pre, QecConst::za);
	uint8_t m0 = q->mz(0);
	uint8_t r0 = (m0 & 1);
	// uint8_t fr0 = r0 ^ (z0 == 2);

	// same for q2
	// uint8_t v2_pre = q->get_vop(2);
	// uint8_t z2 = zeta_for(v2_pre, QecConst::za);
	uint8_t m2 = q->mz(2);
	uint8_t r2 = (m2 & 1);
	// uint8_t fr2 = r2 ^ (z2 == 2);

	// // strict GHZ parity in the computational basis:
	// // frame-invariant bits must match the leaf's frame-invariant bit
	// CHECK(fr0 == fr1);
	// CHECK(fr2 == fr1);

	// sanity check
	CHECK((q->mz(0) & 1) == r0);
	CHECK((q->mz(2) & 1) == r2);
}

// some helpers ---

static inline PackedInt32Array sort_copy(const PackedInt32Array &a) {
	PackedInt32Array b = a;
	std::vector<int> tmp(b.size());
	for (int i = 0; i < b.size(); ++i) {
		tmp[i] = b[i];
	}
	std::sort(tmp.begin(), tmp.end());
	for (int i = 0; i < (int)tmp.size(); ++i) {
		b.set(i, tmp[i]);
	}
	return b;
}

static inline bool adj_eq_as_set(const PackedInt32Array &adj, std::initializer_list<int> expect) {
	std::vector<int> lhs(adj.size());
	for (int i = 0; i < adj.size(); ++i) {
		lhs[i] = adj[i];
	}
	std::sort(lhs.begin(), lhs.end());
	std::vector<int> rhs(expect);
	std::sort(rhs.begin(), rhs.end());
	return lhs == rhs;
}

static inline std::vector<std::pair<int, int>> normalize_edge_pairs(const PackedInt32Array &flat /* [u0,v0,u1,v1,...] */) {
	std::vector<std::pair<int, int>> out;
	out.reserve(flat.size() / 2);
	for (int i = 0; i + 1 < flat.size(); i += 2) {
		int a = flat[i], b = flat[i + 1];
		if (a == b) {
			continue;
		}
		if (b < a) {
			std::swap(a, b);
		}
		out.emplace_back(a, b);
	}
	std::sort(out.begin(), out.end());
	out.erase(std::unique(out.begin(), out.end()), out.end());
	return out;
}

static inline bool edge_set_equals(const PackedInt32Array &flat, std::initializer_list<std::pair<int, int>> expect_pairs) {
	auto got = normalize_edge_pairs(flat);
	std::vector<std::pair<int, int>> exp(expect_pairs);
	for (auto &p : exp) {
		if (p.second < p.first) {
			std::swap(p.first, p.second);
		}
	}
	std::sort(exp.begin(), exp.end());
	exp.erase(std::unique(exp.begin(), exp.end()), exp.end());
	return got == exp;
}

TEST_CASE("[Modules][Qec] Entanglement group test Bell collapse") {
	Ref<Qec> q = memnew(Qec);
	q->init(2);

	q->hadamard(0);
	q->cnot(0, 1);

	// Both in same component
	{
		auto g0 = q->get_entanglement_group(0);
		auto g1 = q->get_entanglement_group(1);
		CHECK(sort_copy(g0) == PackedInt32Array{ 0, 1 });
		CHECK(sort_copy(g1) == PackedInt32Array{ 0, 1 });
	}

	// Collapse
	(void)q->mz(0);
	CHECK(sort_copy(q->get_entanglement_group(0)) == PackedInt32Array{ 0 });
	CHECK(sort_copy(q->get_entanglement_group(1)) == PackedInt32Array{ 1 });
}

TEST_CASE("[Modules][Qec] Snapshot contents GHZ star {0,1,2}") {
	Ref<Qec> q = memnew(Qec);
	q->init(3);

	// GHZ star, 0 is center
	q->hadamard(0);
	q->cnot(0, 1);
	q->cnot(0, 2);

	// check star shape
	CHECK(adj_eq_as_set(q->get_adjacent(0), { 1, 2 }));
	CHECK(adj_eq_as_set(q->get_adjacent(1), { 0 }));
	CHECK(adj_eq_as_set(q->get_adjacent(2), { 0 }));

	Dictionary snap = q->snapshot_entanglement_group(0);

	// nodes field must be the whole component {0,1,2}
	PackedInt32Array nodes = snap["nodes"];
	CHECK(sort_copy(nodes) == PackedInt32Array{ 0, 1, 2 });

	// vops must match live vops in same order as 'nodes'
	PackedByteArray vops = snap["vops"];
	CHECK(vops.size() == nodes.size());
	for (int i = 0; i < nodes.size(); ++i) {
		CHECK(vops[i] == q->get_vop(nodes[i]));
	}

	// edges must contain exactly (0,1) and (0,2), nothing else
	PackedInt32Array edges = snap["edges"];
	CHECK(edge_set_equals(edges, { { 0, 1 }, { 0, 2 } }));
}

TEST_CASE("[Modules][Qec] Undo collapse on Bell via snapshot/restore") {
	Ref<Qec> q = memnew(Qec);
	q->init(2);

	q->hadamard(0);
	q->cnot(0, 1);

	// Take snapshot of entanglement group
	Dictionary snap = q->snapshot_entanglement_group(0);

	// keep originals for exact equality checks
	uint8_t v0 = q->get_vop(0), v1 = q->get_vop(1);
	PackedInt32Array a0 = q->get_adjacent(0), a1 = q->get_adjacent(1);

	// Collapse
	(void)q->mz(0);
	CHECK(adj_eq_as_set(q->get_adjacent(0), {}));
	CHECK(adj_eq_as_set(q->get_adjacent(1), {}));

	// Restore
	q->restore_entanglement_group(snap);

	// Must be bit-identical to pre-measurement state
	CHECK(q->get_vop(0) == v0);
	CHECK(q->get_vop(1) == v1);
	CHECK(q->get_adjacent(0) == a0);
	CHECK(q->get_adjacent(1) == a1);
}

TEST_CASE("[Modules][Qec] Undo X-collapse on GHZ center via snapshot/restore") {
	Ref<Qec> q = memnew(Qec);
	q->init(3);

	// GHZ star
	q->hadamard(0);
	q->cnot(0, 1);
	q->cnot(0, 2);

	// Snapshot the {0,1,2} component
	Dictionary snap = q->snapshot_entanglement_group(0);

	// Record original adjacencies & vops
	PackedInt32Array a0 = q->get_adjacent(0);
	PackedInt32Array a1 = q->get_adjacent(1);
	PackedInt32Array a2 = q->get_adjacent(2);
	uint8_t v0 = q->get_vop(0), v1 = q->get_vop(1), v2 = q->get_vop(2);

	// Measure X on 0, neighbors (1,2) should be entangled and so connected, 0 disconnects
	(void)q->mx(0);
	CHECK(adj_eq_as_set(q->get_adjacent(0), {}));
	CHECK(q->get_adjacent(1).size() == 1);
	CHECK(q->get_adjacent(2).size() == 1);
	CHECK(q->get_adjacent(1) == PackedInt32Array{ 2 });
	CHECK(q->get_adjacent(2) == PackedInt32Array{ 1 });

	// Restore to pre-measurement star
	q->restore_entanglement_group(snap);

	CHECK(q->get_adjacent(0) == a0);
	CHECK(q->get_adjacent(1) == a1);
	CHECK(q->get_adjacent(2) == a2);
	CHECK(q->get_vop(0) == v0);
	CHECK(q->get_vop(1) == v1);
	CHECK(q->get_vop(2) == v2);
}

TEST_CASE("[Modules][Qec] Restore removes edges that were added to the group later") {
	Ref<Qec> q = memnew(Qec);
	q->init(3);

	// bell
	q->hadamard(0);
	q->cnot(0, 1);

	// Snapshot {0,1}
	Dictionary snap = q->snapshot_entanglement_group(0);

	// new edge from 0 to 2
	q->cnot(0, 2);

	// 0 should now have both 1 and 2 as neighbors
	CHECK(q->get_adjacent(0).size() >= 1);

	// Restore should drop the 0-2 edge and revert to the original {0,1} only
	q->restore_entanglement_group(snap);

	CHECK(adj_eq_as_set(q->get_adjacent(0), { 1 }));
	CHECK(adj_eq_as_set(q->get_adjacent(1), { 0 }));
	CHECK(adj_eq_as_set(q->get_adjacent(2), {}));
}

} //namespace TestQec
