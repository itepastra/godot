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

} //namespace TestQec
