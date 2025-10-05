#pragma once

#include "modules/qec/qec.h"
#include "tests/test_macros.h"

namespace TestQec {
TEST_CASE("[Modules][Qec] Chained CNOT Gates") {
	Ref<Qec> q = memnew(Qec);
	q->init(4);

	CHECK(q->get_vop(0) == 10);
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_vop(2) == 10);
	CHECK(q->get_vop(3) == 10);

	q->cnot(0, 1);
	q->cnot(1, 2);
	q->cnot(2, 3);
	q->cnot(3, 4);

	CHECK(q->get_vop(0) == 10);
	CHECK(q->get_vop(1) == 10);
	CHECK(q->get_vop(2) == 10);
	CHECK(q->get_vop(3) == 10);
}
} //namespace TestQec
