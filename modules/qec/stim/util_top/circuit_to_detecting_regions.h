#ifndef _STIM_UTIL_TOP_CIRCUIT_TO_DETECTING_REGIONS
#define _STIM_UTIL_TOP_CIRCUIT_TO_DETECTING_REGIONS

#include "../circuit/circuit.h"
#include "../dem/dem_instruction.h"
#include "../stabilizers/flex_pauli_string.h"

namespace stim {

std::map<DemTarget, std::map<uint64_t, FlexPauliString>> circuit_to_detecting_regions(
    const Circuit &circuit,
    std::set<stim::DemTarget> included_targets,
    std::set<uint64_t> included_ticks,
    bool ignore_anticommutation_errors);

}  // namespace stim

#endif
