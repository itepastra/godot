#include "qec.h"

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
    circuit.safe_append_u("X", {qubit});
}

void Qec::add_cx_gate(int control, int target) {
    if (control < 0 || control >= num_qubits || target < 0 || target >= num_qubits) {
        ERR_FAIL_MSG("Qubit index out of range");
    }
    circuit.safe_append_u("CX", {control, target});
}

void Qec::add_measurement(int qubit) {
    if (qubit < 0 || qubit >= num_qubits) {
        ERR_FAIL_MSG("Qubit index out of range");
    }
    circuit.safe_append_u("M", {qubit});
}

Array Qec::simulate_measurements(int shot_count) {
	Array results;
    
    try {
        // sample measurements from the circuit
        auto samples = circuit.compile_sampler().sample(shot_count);
        
        // convert to godot array
        for (size_t shot = 0; shot < samples.size(); shot++) {
            Array shot_results;
            for (size_t qubit = 0; qubit < samples[shot].size(); qubit++) {
                shot_results.push_back((int)samples[shot][qubit]);
            }
            results.push_back(shot_results);
        }
    } catch (const std::exception& e) {
        ERR_FAIL_V_MSG(Array(), String("Stim simulation error: ") + e.what());
    }
    
    return results;
}


Array Qec::simulate_detector_samples(int shot_count) {
    Array results;
    
    try {
        // detector samples detecting errors
        auto detector_samples = circuit.compile_detector_sampler().sample(shot_count);
        
        // to godot arr
        for (size_t shot = 0; shot < detector_samples.size(); shot++) {
            Array shot_results;
            for (size_t detector = 0; detector < detector_samples[shot].size(); detector++) {
                shot_results.push_back((int)detector_samples[shot][detector]);
            }
            results.push_back(shot_results);
        }
    } catch (const std::exception& e) {
        ERR_FAIL_V_MSG(Array(), String("Stim detector sampling error: ") + e.what());
    }
    
    return results;
}

int get_qubit_count() const {
	    return num_qubits;
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
	ClassDB::bind_method(D_METHOD("simulate_detector_samples", "shot_count"), &Qec::simulate_detector_samples);
	ClassDB::bind_method(D_METHOD("get_qubit_count"), &Qec::get_qubit_count);

}
