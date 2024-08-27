import heapq

import numpy as np

from . import symbolic_variable
from .symbolic_variable import Variable


def const_fp_bits(x):
    "Number of fp bits needed to represent x exactly."
    l, h = -55, 55  # noqa: E741
    while l < h:  # noqa: E741
        v = (l + h) // 2
        _v = x * 2.0**v
        if np.round(_v) == _v:
            h = v
        else:
            l = v + 1  # noqa: E741
    return l


def binary_shift_decompose(x: int | float | np.number):
    '''Given a number, return positive and negative bitshifts to represent another number multiplied by it:
    returns pos, neg such that `N*x = sum(x << pos) - sum(x << neg)` for some N
    '''
    x, sign = abs(x), np.sign(x)
    f = const_fp_bits(x)
    x *= 2**f
    x = int(x)
    binary = np.array(list(bin(x)[2:])[::-1], dtype=np.int8)
    binary_m1 = np.array(list(bin(abs(x - 1))[2:])[::-1], dtype=np.int8)

    idx_pos = (np.where(binary)[0] - f).astype(np.int8)
    idx_neg = (np.where(binary_m1 == 0)[0] - f).astype(np.int8)

    if len(idx_neg) + 1 < len(idx_pos):
        r = (np.array([len(binary_m1) - f], dtype=np.int8), idx_neg)
    else:
        r = (idx_pos, np.array([], dtype=np.int8))
    if sign < 0:
        r = r[::-1]
    return r


def get_mat_shift_mask(arr: np.ndarray):
    '''Given an array, return the bitshifts and combined mask to represent the array as a sum of powers of 2:

    ```
    row_scale = 2.**shifts[:, None]
    pow2 = np.sum((2.**np.arange(mask.shape[2]))[None, None, :]
    arr == row_scale * pow2 * mask, axis=2)
    ```
    '''
    if arr.ndim == 1:
        arr = arr[:, None]
    assert arr.ndim == 2, "Only 2D arrays are supported"  # ch_in, ch_out
    shape = arr.shape
    arr_flat = arr.ravel()
    shifts = np.full(shape[0], 127, dtype=np.int8)
    poss, negs = [], []

    for i, x in enumerate(arr_flat):
        i0 = np.unravel_index(i, shape)[0]
        p, n = binary_shift_decompose(x)
        shifts[i0] = min(*p[:1], *n[:1], shifts[i0], 127)
        poss.append(p)
        negs.append(n)

    mask_width = 0
    for i, (p, n) in enumerate(zip(poss, negs)):
        s = shifts[i % shape[0]]
        mask_width = max(*(p[-1:] - s), *(n[-1:] - s), mask_width, -128)
    mask_width += 1
    assert mask_width > 0
    mask = np.zeros((len(arr_flat), mask_width), dtype=np.int8)
    for i, (p, n) in enumerate(zip(poss, negs)):
        s = shifts[i % shape[0]]
        mask[i, p - s] = 1
        mask[i, n - s] = -1
    shifts[shifts == 127] = 0  # 127 means the input value is exactly 0
    mask = mask.reshape(shape + (mask_width,))
    return shifts, mask


def get_bit_reduction_loc(combination_mask: np.ndarray, bit_mask: list[int] | None = None):
    '''
    Get a bit location to reduce in the next iteration.
    Args:
        `combination_mask`: (..., n_bits), represents the binary mask of a vector in binary form in a 2-d array. Sign of each bit is represented by 1 and -1. np.sum([[1,2,4,...]] * combination_mask) gives the original vector.
        `bit_mask`: list of bits to consider. If None, all bits are considered.

    Returns:
        the bit location to reduce. Currently, the bit with the least number of non-zero entries is returned.
    '''  # noqa: E501
    bit_mask = bit_mask or list(range(combination_mask.shape[-1]))
    assert len(bit_mask) > 0, "bit_mask must be non-empty"
    idx = int(np.argmin(np.sum(combination_mask[..., bit_mask] != 0, axis=0)))
    return bit_mask[idx]


def bit_reduction(combination_mask: np.ndarray, bit_mask: list[int], bit_loc: int):
    '''
    Gather indices of entries that has non-zero component at `bit_loc`, and reduce the bit_loc from the combination_mask (remove it from bit_mask). Returns the necessary operations to perform this operation.
    Args:
        `combination_mask`: (..., n_bits), represents the binary mask of a vector in binary form in a 2-d array. Sign of each bit is represented by 1 and -1. np.sum([[1,2,4,...]] * combination_mask) gives the original vector.
        `bit_mask`: list of bits to consider. If None, all bits are considered.
        `bit_loc`: the bit location to reduce.

    Returns:
        `extract_from`: the summed value for `bit_loc` should be computed from the values at these indices.
        `gather_to`: where the values should go for the next iteration. If -1, the value is not used for the next iteration and should be discarded.
        `idx`: the unique indices of the combination_mask.
    '''  # noqa: E501

    pos_bits = np.where(combination_mask[..., bit_loc] == 1)[0]
    neg_bits = np.where(combination_mask[..., bit_loc] == -1)[0]
    extract_from = np.concatenate([pos_bits, -neg_bits - 1])  # extra -1 to break ties in the case of 0

    bit_mask = [i for i in bit_mask if i != bit_loc]
    combination_mask = combination_mask[..., bit_mask]
    _, idx, inverse_idx = np.unique(combination_mask, axis=0, return_index=True, return_inverse=True)
    gather_to = inverse_idx

    zeros = np.all(combination_mask == 0, axis=-1)
    gather_to[zeros] = -1
    return extract_from, gather_to, idx


def to_operations(arr: np.ndarray):
    """For a 2d array as linear operator, decompose it as a series of operations.
    y = arr @ v is equivalent to:
    Returns:
        `shift`, `gather_tos`, `extract_froms`, `bit_extract_order`
    ```
    v << shifts elementwise
    buffer0 = v
    init r = empty
    for gather_to, extract_from in zip(gather_tos, extract_froms):
        init buffer1 = zeros
        for i, _to in enumerate(gather_to):
            if _to != -1:
                buffer1[_to] += buffer0[i]
        for i, _from in enumerate(extract_from):
            if _from >= 0:
                r[*bit_extract_order[i]] += buffer1[_from]
            else:
                r[*bit_extract_order[i]] -= buffer1[-_from+1]
        buffer0 = buffer1
    y = np.sum(r * 2**np.arange(r.shape[-1]), axis=-1)
    return y
    ```

    """
    shift, combination_mask = get_mat_shift_mask(arr)
    ch_in, ch_out, bw = combination_mask.shape

    # hackey way to support 2d arrays. Turned out to be NOT optimal at all. should break into dotp for better performance
    combination_mask = combination_mask.reshape(ch_in, ch_out * bw)

    bit_mask = list(range(combination_mask.shape[-1]))

    gather_tos, extract_froms = [], []
    bit_extract_order = []

    _, gather_to, idx = bit_reduction(combination_mask, bit_mask, -1)
    gather_tos.append(gather_to)
    combination_mask = combination_mask[idx]

    for _ in range(ch_out * bw):
        loc = get_bit_reduction_loc(combination_mask, bit_mask)
        bit_mask.remove(loc)
        extract_from, gather_to, idx = bit_reduction(combination_mask, bit_mask, loc)
        ch, bit_loc = np.unravel_index(loc, (ch_out, bw))
        bit_extract_order.append((ch, bit_loc))
        extract_froms.append(extract_from)
        gather_tos.append(gather_to)
        combination_mask = combination_mask[idx]

    gather_tos = gather_tos[:-1]

    return shift, gather_tos, extract_froms, bit_extract_order


def balanced_reduction(vec: list):
    if len(vec) == 0:
        return 0
    vec = [x for x in vec if isinstance(x, Variable)]
    bias = sum(x for x in vec if not isinstance(x, Variable))
    with symbolic_variable.order_metrics(('depth', 'i')), symbolic_variable.fuse_associative_ops(False):
        heapq.heapify(vec)
        while len(vec) > 1:
            v1, v2 = heapq.heappop(vec), heapq.heappop(vec)
            heapq.heappush(vec, v1 + v2)
    return vec[0] + bias if vec else bias  # type: ignore


def _compile_dense(kernel: np.ndarray, inp: np.ndarray):
    shifts, gather_tos, extract_froms, bit_extract_order = to_operations(kernel)
    ch_in, ch_out = kernel.shape
    buf0 = inp * 2.0**shifts
    r: list[float | Variable | list[Variable]] = np.empty((ch_out, 0), dtype=object).tolist()
    with symbolic_variable.fuse_associative_ops(False):
        for i, (gather_to, extract_from) in enumerate(zip(gather_tos, extract_froms)):
            _buf1 = [[] for _ in range(np.max(gather_to) + 1)]
            for v, to in zip(buf0, gather_to):
                if to >= 0:
                    _buf1[to].append(v)
            buf1 = [balanced_reduction(v) for v in _buf1]
            bit_loc = bit_extract_order[i][1]
            _r = [(buf1[_from] if _from >= 0 else -buf1[-_from - 1]) * 2.0**bit_loc for _from in extract_from]
            r[bit_extract_order[i][0]].extend(_r)
            buf0 = buf1
        for i, _r in enumerate(r):
            x = balanced_reduction(_r)  # type: ignore
            if isinstance(x, Variable):
                x.fix_precision(recursive=True)
            r[i] = x
    return r


def compile_dense(kernel: np.ndarray, inp: np.ndarray):
    out = []
    for _kernel in kernel.T:  # ch_in, 1
        out.append(_compile_dense(_kernel[:, None], inp)[0])
    return np.array(out).T


def compile_conv(kernel: np.ndarray, inp: list | np.ndarray):
    """
    Apply a single kernel to the input x
    """
    *_ch_in, ch_out = kernel.shape
    ch_in = int(np.prod(_ch_in))
    inp = np.reshape(inp, ch_in)
    kernel = np.reshape(kernel, (ch_in, ch_out))
    return compile_dense(kernel, inp)
