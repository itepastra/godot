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
} //namespace TestQec
