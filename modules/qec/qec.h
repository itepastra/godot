#pragma once

#include "core/object/ref_counted.h"

class Qec : public RefCounted {
	GDCLASS(Qec, RefCounted);

	int count;

protected:
	static void _bind_methods();

public:
	void add(int p_value);
	void reset();
	int get_total() const;

	Qec();
};
