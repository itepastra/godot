#include "qec.h"
#include "core/error/error_macros.h"
#include "stim/simulators/tableau_simulator.h"

Qec::Qec() {
	count = 0;
	num_qubits = 0;
}

void Qec::add(int p_value) {
	count += p_value;
}

void Qec::reset() {
	count = 0;
}

int Qec::get_total() const {
	return count;
}

void Qec::init_circuit(int qubit_count) {
	circuit = stim::Circuit();
	num_qubits = qubit_count;
}

void Qec::add_x_gate(int qubit) {
	if (qubit < 0 || qubit >= num_qubits) {
        ERR_FAIL_MSG("Qubit index out of range");
    }
    circuit.safe_append_u("X", {(uint32_t)qubit});
}

void Qec::add_cx_gate(int control, int target) {
    if (control < 0 || control >= num_qubits || target < 0 || target >= num_qubits) {
        ERR_FAIL_MSG("Qubit index out of range");
    }
    circuit.safe_append_u("CX", {(uint32_t)control, (uint32_t)target});
}

void Qec::add_measurement(int qubit) {
    if (qubit < 0 || qubit >= num_qubits) {
        ERR_FAIL_MSG("Qubit index out of range");
    }
    circuit.safe_append_u("M", {(uint32_t)qubit});
}

Array Qec::simulate_measurements(int shot_count) {
    if (shot_count <= 0) {
        ERR_FAIL_V_MSG(Array(), "Shot count must be positive");
    }
    
    if (num_qubits == 0) {
        ERR_FAIL_V_MSG(Array(), "Circuit not initialized. Call init_circuit first.");
    }
    
    Array results_array;
    results_array.resize(shot_count);
    
    std::random_device rd;
    std::mt19937_64 rng(rd());
    
    for (int i = 0; i < shot_count; i++) {
        auto measurement_results = stim::TableauSimulator<stim::MAX_BITWORD_WIDTH>::sample_circuit(circuit, rng);
        
        // convert to Godot arr
        Array shot_results;
        size_t num_measurements = circuit.count_measurements();
        shot_results.resize(num_measurements);
        
        for (size_t j = 0; j < num_measurements; j++) {
            shot_results[j] = (bool)measurement_results[j];
        }
        
        results_array[i] = shot_results;
    }
    
    return results_array;
}


int Qec::get_qubit_count() const {
    return (int)num_qubits;
}


void Qec::_bind_methods() {
	ClassDB::bind_method(D_METHOD("add", "value"), &Qec::add);
	ClassDB::bind_method(D_METHOD("reset"), &Qec::reset);
	ClassDB::bind_method(D_METHOD("get_total"), &Qec::get_total);

	ClassDB::bind_method(D_METHOD("init_circuit", "qubit_count"), &Qec::init_circuit);
	ClassDB::bind_method(D_METHOD("add_x_gate", "qubit"), &Qec::add_x_gate);
	ClassDB::bind_method(D_METHOD("add_cx_gate", "control", "target"), &Qec::add_cx_gate);
	ClassDB::bind_method(D_METHOD("add_measurement", "qubit"), &Qec::add_measurement);
	ClassDB::bind_method(D_METHOD("simulate_measurements", "shot_count"), &Qec::simulate_measurements);
	ClassDB::bind_method(D_METHOD("get_qubit_count"), &Qec::get_qubit_count);

}
