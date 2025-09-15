#include "qec.h"

void QEC::add(int p_value) {
	count += p_value;
}

void QEC::reset() {
	count = 0;
}

int QEC::get_total() const {
	return count;
}

void QEC::_bind_methods() {
	ClassDB::bind_method(D_METHOD("add", "value"), &QEC::add);
	ClassDB::bind_method(D_METHOD("reset"), &QEC::reset);
	ClassDB::bind_method(D_METHOD("get_total"), &QEC::get_total);
}

QEC::QEC() {
	count = 0;
}
