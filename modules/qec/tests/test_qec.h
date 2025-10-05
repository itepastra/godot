#pragma once

#include "modules/qec/qec.h"
#include "tests/test_macros.h"

namespace TestQec {
TEST_CASE("[Modules][Qec] Chained CNOT Gates") {
	Ref<Qec> q = memnew(Qec);
	q->init(4);

	// check if all qubits at |0>
	CHECK(q->get_vop(0) == 10);
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_vop(2) == 10);
	CHECK(q->get_vop(3) == 10);

	// rotate qubit to |+>
	q->hadamard(0);
	CHECK(q->get_vop(0) == 0);

	// entangle qubit 0 and 1
	q->cnot(0, 1);
	CHECK(q->get_vop(0) == 0);
	CHECK(q->get_adjacent(0) == PackedInt32Array{ 1 });
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_adjacent(1) == PackedInt32Array{ 0 });

	q->cnot(1, 2);

	q->cnot(2, 3);

	CHECK(q->get_vop(0) == 10);
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_vop(2) == 10);
	CHECK(q->get_vop(3) == 10);
}
} //namespace TestQec
