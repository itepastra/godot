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
	CHECK(q->measure(0) == 0b10); // |0> deterministic
	q->zgate(0);
	CHECK(q->measure(0) == 0b11); // |1> deterministic
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
} //namespace TestQec
