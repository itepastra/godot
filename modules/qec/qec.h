#pragma once

#include "core/object/ref_counted.h"
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#define BLOCK_BITS 6
#define node_idx uint32_t
#define ia 0
#define ib 1
#define ic 2
#define id 3
#define ie 4
#define il 5
#define xa 6
#define xb 7
#define xc 8
#define xd 9
#define xe 10
#define xf 11
#define ya 12
#define yb 13
#define yc 14
#define yd 15
#define ye 16
#define yf 17
#define za 18
#define zb 19
#define zc 20
#define zd 21
#define ze 22
#define zf 23

struct Node {
	uint8_t vop = ia;
	std::vector<node_idx> adjacent;
};

class Qec : public RefCounted {
	GDCLASS(Qec, RefCounted);

	bool initialized;
	node_idx qubit_count;
	std::vector<Node> nodes;

protected:
	// godot helper functions
	static void _bind_methods();

	void remove_VOP(node_idx a, node_idx b);
	void local_complementation(node_idx a);

public:
	// qubit gates for godot to use
	void cnot(node_idx control, node_idx target); // CNOT
	void hadamard(node_idx target); // H
	void phase(node_idx target); // S
	void phase_dag(node_idx target); // S^+ = S S S
	void xgate(node_idx target); // X = H S S H
	void ygate(node_idx target); // Y = X Y = H S S H S S
	void zgate(node_idx target); // Z = S S
	void cphase(node_idx control, node_idx target); // CZ = H_target CNOT H_target

	// initialisation
	Qec();
	void init(node_idx qubit_amount);
};

const uint32_t SYMMETRIES = 24;
// I = Ia
// H = Yc
// S^ = Yb
// S = Xb
// Z = S S Ia = Za
// X = H S S H Ia = Xa
// Y = Ya

// This table is copied from page 31 of https://archive.org/download/thesis-anders/Thesis_Anders.pdf
const uint8_t vop_table[SYMMETRIES][SYMMETRIES] = {
	{ ia, ib, ic, id, ie, il, xa, xb, xc, xd, xe, xf, ya, yb, yc, yd, ye, yf, za, zb, zc, zd, ze, zf }, // Ia
	{ ib, ia, il, ie, id, ic, xb, xa, xf, xe, xd, xc, yb, ya, zf, ye, yd, yc, zb, za, zf, ze, zd, zc }, // Ib
	{ ic, ie, ia, il, ib, id, xc, xe, xa, xf, xb, xd, yc, ye, ya, yf, yb, yd, zc, ze, za, zf, zb, zd }, // Ic
	{ id, il, ie, ia, ic, ib, xd, xf, xe, xa, xc, xb, yd, yf, ye, ya, yc, yb, zd, zf, ze, za, zc, zb }, // Id
	{ ie, ic, id, ib, il, ia, xe, xc, xd, xb, xf, xa, ye, yc, yd, yb, yf, ya, ze, zc, zd, zb, zf, za }, // Ie
	{ il, id, ib, ic, ia, ie, xf, xd, xb, xc, xa, xe, yf, yd, yb, yc, ya, ye, zf, zd, zb, zc, za, ze }, // If
	// ------------------------------------------------------------------------------------------------------
	{ xa, yb, zc, xd, ze, yf, ia, zb, yc, id, ye, zf, za, ib, xc, zd, xe, il, ya, xb, ic, yd, ie, xf }, // Xa
	{ xb, ya, zf, xe, zd, yc, ib, za, yf, ie, yd, zc, zb, ia, xf, ze, xd, ic, yb, xa, il, ye, id, xc }, // Xb = S
	{ xc, ye, za, xf, zb, yd, ic, ze, ya, il, yb, zd, zc, ie, xa, zf, xb, id, yc, xe, ia, yf, ib, xd }, // Xc
	{ xd, yf, ze, xa, zc, yb, id, zf, ye, ia, yc, zb, zd, il, xe, za, xc, ib, yd, xf, ie, ya, ic, xb }, // Xd
	{ xe, yc, zd, xb, zf, ya, ie, zc, yd, ib, yf, za, ze, ic, xd, zb, xf, ia, ye, xc, id, yb, il, xa }, // Xe
	{ xf, yd, zb, xc, za, ye, il, zd, yb, ic, ya, ze, zf, id, xb, zc, xa, ie, yf, xd, ib, yc, ia, xe }, // Xf
	// ------------------------------------------------------------------------------------------------------
	{ ya, xb, yc, zd, xe, zf, za, ib, zc, yd, ie, yf, ia, zb, ic, xd, ze, xf, xa, yb, xc, id, ye, il }, // Ya
	{ yb, xa, yf, ze, xd, zc, zb, ia, zf, ye, id, yc, ib, za, il, xe, zd, xc, xb, ya, xf, ie, yd, ic }, // Yb = S^dagger
	{ yc, xe, ya, zf, xb, zd, zc, ie, za, yf, ib, yd, ic, ze, ia, xf, zb, xd, xc, ye, xa, il, yb, id }, // Yc = Hadamard
	{ yd, xf, ye, za, xc, zb, zd, il, ze, ya, ic, yb, id, zf, ie, xa, zc, xb, xd, yf, xe, ia, yc, ib }, // Yd
	{ ye, xc, yd, zb, xf, za, ze, ic, zd, yb, il, ya, ie, zc, id, xb, zf, xa, xe, yc, xd, ib, yf, ia }, // Ye
	{ yf, xd, yb, zc, xa, ze, zf, id, zb, yc, ia, ye, il, zd, ib, xc, za, xe, xf, yd, xb, ic, ya, ie }, // Yf
	// ------------------------------------------------------------------------------------------------------
	{ za, zb, xc, yd, ye, xf, ya, yb, ic, zd, ze, il, xa, xb, zc, id, ie, zf, ia, ib, yc, xd, xe, yf }, // Za
	{ zb, za, xf, ye, yd, xc, yb, ya, il, ze, zd, ic, xb, xa, zf, ie, id, zc, ib, ia, yf, xe, xd, yc }, // Zb
	{ zc, ze, xa, yf, yb, xd, yc, ye, ia, zf, zb, id, xc, xe, za, il, ib, zd, ic, ie, ya, xf, xb, yd }, // Zc
	{ zd, zf, xe, ya, yc, xb, yd, yf, ie, za, zc, ib, xd, xf, ze, ia, ic, zb, id, il, ye, xa, xc, yb }, // Zd
	{ ze, zc, xd, yb, yf, xa, ye, yc, id, zb, zf, ia, xe, xc, zd, ib, il, za, ie, ic, yd, xb, xf, ya }, // Ze
	{ zf, zd, xb, yc, ya, xe, yf, yd, ib, zc, za, ie, xf, xd, zb, ic, ia, ze, il, id, yb, xc, xa, ye } // Zf
};

const std::vector<uint8_t> decompositions[SYMMETRIES] = {
	{},
	{ 0, 0 },
	{},
	{ 1, 1 },
	{},
	{},
	// -----------------------------
	{},
	{},
	{},
	{},
	{},
	{},
	// -----------------------------
	{},
	{},
	{},
	{},
	{},
	{},
	// -----------------------------
	{},
	{},
	{},
	{},
	{},
	{},
};
