#ifndef PHOTOSPLINE_DETAIL_CONTRACT_H
#define PHOTOSPLINE_DETAIL_CONTRACT_H

#include <cstddef>
#include <cstdint>
#include <memory>

namespace photospline{ namespace detail{

// Blocks with more rows than this keep their scratch space on the heap, sparing thread stacks.
constexpr unsigned contract_stack_rows = 1024;

template<typename V, typename F>
inline void splat(V& out, F x){ out = x - V{}; }

///Sums coefficient * (product of one basis value per dimension) over the block of
///(order+1)^ndim coefficients that support the point described by centers, for NV lanes of
///type V at once (V is a scalar for plain evaluation, a SIMD vector for gradients).
///
///The block is processed as rows along the last dimension. Each row is a short dot product
///with the last dimension's basis values, weighted by the outer product of the other
///dimensions' values, which is precomputed together with each row's offset so the loop has no
///branching index arithmetic. Rows go to four running sums in turn, so successive rows do not
///each wait for the previous addition to finish.
///
///Plain and gradient evaluation both use this function, so they add in the same order and
///agree exactly.
///
///\tparam Row the last dimension's order+1 if known at compile time, otherwise 0
///\param basis basis(d, i) returns a pointer to the NV lanes of basis function i in dimension d
///\param result NV lanes to which the sum is added
template<unsigned Row, unsigned NV, typename V, typename Basis>
void contract(const float* coefficients, const uint64_t* strides, const uint32_t* order,
              unsigned ndim, const int* centers, const Basis& basis, V* result)
{
	const unsigned last = ndim - 1;
	const unsigned row = Row ? Row : order[last] + 1;

	int64_t start = 0;
	for (unsigned d = 0; d < ndim; d++)
		start += (centers[d] - int64_t(order[d]))*int64_t(strides[d]);
	const float* block = coefficients + start;
	const V* last_basis[row];
	for (unsigned i = 0; i < row; i++)
		last_basis[i] = basis(last, i);

	V acc[4][NV];
	for (unsigned j = 0; j < 4; j++)
		for (unsigned k = 0; k < NV; k++)
			splat(acc[j][k], 0.f);
	auto add_row = [&](const float* c, const V* weight, V* sum){
		V dot[NV];
		for (unsigned k = 0; k < NV; k++)
			splat(dot[k], 0.f);
		for (unsigned i = 0; i < row; i++) {
			V coefficient;
			splat(coefficient, c[i]);
			for (unsigned k = 0; k < NV; k++)
				dot[k] += last_basis[i][k]*coefficient;
		}
		for (unsigned k = 0; k < NV; k++)
			sum[k] += weight[k]*dot[k];
	};
	auto finish = [&]{
		for (unsigned k = 0; k < NV; k++)
			result[k] += ((acc[0][k] + acc[1][k]) + acc[2][k]) + acc[3][k];
	};

	// With at most two weight factors, form them on the fly; the products and the assignment of
	// rows to running sums are the same as on the general path below, so results are too.
	if (last <= 2) {
		const unsigned n0 = last > 0 ? order[0] + 1 : 1;
		const unsigned n1 = last > 1 ? order[1] + 1 : 1;
		const int64_t s0 = last > 0 ? int64_t(strides[0]) : 0;
		const int64_t s1 = last > 1 ? int64_t(strides[1]) : 0;
		unsigned r = 0;
		for (unsigned i0 = 0; i0 < n0; i0++) {
			for (unsigned i1 = 0; i1 < n1; i1++, r++) {
				V weight[NV];
				for (unsigned k = 0; k < NV; k++) {
					if (last == 0)
						splat(weight[k], 1.f);
					else if (last == 1)
						weight[k] = basis(0, i0)[k];
					else
						weight[k] = basis(0, i0)[k]*basis(1, i1)[k];
				}
				add_row(block + i0*s0 + i1*s1, weight, acc[r % 4]);
			}
		}
		finish();
		return;
	}

	unsigned rows = 1;
	for (unsigned d = 0; d < last; d++)
		rows *= order[d] + 1;
	const bool on_stack = rows <= contract_stack_rows;
	V stack_weight[on_stack ? rows*NV : 1];
	int64_t stack_offset[on_stack ? rows : 1];
	V* weight = stack_weight;
	int64_t* offset = stack_offset;
	std::unique_ptr<char[]> heap;
	if (!on_stack) {
		size_t space = rows*(NV*sizeof(V) + sizeof(int64_t)) + alignof(V);
		heap.reset(new char[space]);
		void* buf = heap.get();
		weight = static_cast<V*>(std::align(alignof(V), rows*NV*sizeof(V), buf, space));
		offset = reinterpret_cast<int64_t*>(weight + rows*NV);
	}

	// Expand the outer product one dimension at a time, in place and from the back so each
	// entry is read before it is overwritten. Dimension 0 seeds it directly (1*x is exact).
	unsigned n = order[0] + 1;
	for (unsigned i = 0; i < n; i++) {
		const V* b = basis(0, i);
		for (unsigned k = 0; k < NV; k++)
			weight[i*NV+k] = b[k];
		offset[i] = int64_t(i)*int64_t(strides[0]);
	}
	for (unsigned d = 1; d < last; d++) {
		const unsigned len = order[d] + 1;
		for (unsigned a = n; a-- > 0;) {
			V w[NV];
			for (unsigned k = 0; k < NV; k++)
				w[k] = weight[a*NV+k];
			const int64_t o = offset[a];
			for (unsigned i = len; i-- > 0;) {
				const V* b = basis(d, i);
				for (unsigned k = 0; k < NV; k++)
					weight[(a*len+i)*NV+k] = w[k]*b[k];
				offset[a*len+i] = o + int64_t(i)*int64_t(strides[d]);
			}
		}
		n *= len;
	}

	unsigned r = 0;
	for (; r + 4 <= n; r += 4) {
		add_row(block + offset[r], weight + r*NV, acc[0]);
		add_row(block + offset[r+1], weight + (r+1)*NV, acc[1]);
		add_row(block + offset[r+2], weight + (r+2)*NV, acc[2]);
		add_row(block + offset[r+3], weight + (r+3)*NV, acc[3]);
	}
	for (unsigned j = 0; r < n; r++, j++)
		add_row(block + offset[r], weight + r*NV, acc[j]);
	finish();
}

}} //namespace photospline::detail

#endif //PHOTOSPLINE_DETAIL_CONTRACT_H
