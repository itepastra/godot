#pragma once

#include "core/object/ref_counted.h"
#include "stim.h"

class QEC : public RefCounted {
	GDCLASS(QEC, RefCounted);

	int count;

protected:
	static void _bind_methods();

public:
	void add(int p_value);
	void reset();
	int get_total() const;

	QEC();
};
