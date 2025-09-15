#include "qec.h"

void Qec::add(int p_value) {
	count += p_value;
}

void Qec::reset() {
	count = 0;
}

int Qec::get_total() const {
	return count;
}

void Qec::_bind_methods() {
	ClassDB::bind_method(D_METHOD("add", "value"), &Qec::add);
	ClassDB::bind_method(D_METHOD("reset"), &Qec::reset);
	ClassDB::bind_method(D_METHOD("get_total"), &Qec::get_total);
}

Qec::Qec() {
	count = 0;
}
