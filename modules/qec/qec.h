#pragma once

#include "core/object/ref_counted.h"
#include "stim.h"
#include "core/variant/array.h"

class Qec : public RefCounted {
	GDCLASS(Qec, RefCounted);

private:
	int count;

	stim::Circuit circuit;
    size_t num_qubits;

protected:
	static void _bind_methods();

public:
	void add(int p_value);
	void reset();
	int get_total() const;

	void init_circuit(int qubit_count);
    void add_x_gate(int qubit);
    void add_cx_gate(int control, int target);
    void add_measurement(int qubit);
    
    Array simulate_measurements(int shot_count);
    Array simulate_detector_samples(int shot_count);

    int get_qubit_count() const;

	Qec();
};
