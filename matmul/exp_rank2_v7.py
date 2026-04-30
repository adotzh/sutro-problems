#!/usr/bin/env python3
"""
Exploring address layout permutations to reduce cost below 73,602.

Key insight: Current layout (sA_cache=1, tmp=2, sB@3-6, sC@7-38) has:
  addr 1 (sA_cache): 4096 reads at cost 1 = 4,096
  addr 2 (tmp): 3840 reads at cost 2 = 7,680
  addr 3 (sB[0]): 1024 reads at cost 2 = 2,048
  addr 4 (sB[1]): 1024 reads at cost 2 = 2,048
  addr 5 (sB[2]): 1024 reads at cost 3 = 3,072
  addr 6 (sB[3]): 1024 reads at cost 3 = 3,072
  sC@7-38: 128 reads each, total 20,736

The read counts are FIXED (determined by Ti, Tj, N):
  sA_cache: 4096 reads
  tmp: 3840 reads
  sB: 1024 reads each (Tj=4 cells)
  sC: 128 reads each (Ti*Tj=32 cells)

We want to minimize sum(reads[i] * cost(addr[i])).

With reads fixed, this is a WEIGHTED ASSIGNMENT problem:
  Assign 39 items to 39 distinct addresses, minimizing weighted cost.
  Items: sA(4096), tmp(3840), sB(1024)*4, sC(128)*32
  Optimal: assign highest-read item to lowest-cost address!

Sorted by reads (descending): sA(4096) > tmp(3840) > sB*4(1024 each) > sC*32(128 each)
Sorted by address cost (ascending): addr 1(cost 1), addr 2(cost 2), addr 3(cost 2), addr 4(cost 2),...

OPTIMAL ASSIGNMENT:
  addr 1 (cost 1): sA (4096 reads) → 4096
  addr 2 (cost 2): tmp (3840 reads) → 7680
  addr 3 (cost 2): sB[0] (1024 reads) → 2048
  addr 4 (cost 2): sB[1] (1024 reads) → 2048
  addr 5 (cost 3): sB[2] (1024 reads) → 3072
  addr 6 (cost 3): sB[3] (1024 reads) → 3072
  addr 7 (cost 3): sC[0] (128 reads) → 384
  ... sC cells at addr 7-38 ...

This IS the current layout! The current layout is already OPTIMAL for the given
read count structure, up to the constraints (sC must be contiguous? No -- sC
only needs each cell to have a UNIQUE address; they don't need to be contiguous).

The question is: can we BREAK the sC contiguity assumption and assign sC cells
to addresses optimally (interleaved with sB etc.)?

The IR generator uses sC_base + ii*Tj + jj, which IS contiguous. But there's no
fundamental reason sC must be contiguous -- we could put each sC(ii,jj) at
any unused address.

KEY QUESTION: If we allow non-contiguous sC, what's the optimal layout?

With 39 items (sA:1, tmp:1, sB:4, sC:32) and 38 distinct addresses (addr 1..39):
Wait -- actually we need any 39 distinct addresses in 1..N for some N.

The optimal assignment (by reads descending → address cost ascending) gives:
  addr 1 (cost 1): sA (4096)
  addr 2,3,4 (cost 2): tmp (3840), sB[0] (1024), sB[1] (1024)
  addr 5,6,...,9 (cost 3): sB[2] (1024), sB[3] (1024), sC[...] x3 (128 each)
  etc.

This DOES match the current layout since:
  sA=1, tmp=2, sB@3-6, sC@7 has 7 items using addr 1-7 optimally.

Can we fit MORE into cost=2 tier (addr 2-4)?
  Addr 2: cost 2 (4 addresses: 2,3,4)
  We have: tmp(3840), sB[0](1024), sB[1](1024) → 3 items at cost 2.
  Actually cost(2)=cost(3)=cost(4)=2. So addr 2,3,4 all have cost 2.

Items at cost 2: tmp(3840), sB[0](1024), sB[1](1024) = 3 items at 3 addresses. Perfect.
Items at cost 3 (addr 5..9): sB[2](1024), sB[3](1024), sC[0..2](128 each) = 5 items, 5 addresses. Perfect.

So the current layout IS the optimal weighted assignment!

CONCLUSION: No further improvement possible from address permutation.

NEW IDEA: Can we REDUCE the number of items (reduce reads)?

Currently:
  sA: 4096 reads = Ti * nbk * Tj * nbi * nbj = 8*16*4*2*4 = 4096 (FIXED by N,Ti,Tj)
  tmp: 3840 reads = Ti * (nbk-1) * Tj * nbi * nbj = 8*15*4*2*4 = 3840 (FIXED by N,Ti,Tj)
  sB: 1024 reads per cell = Ti * nbk * nbi * nbj = 8*16*2*4 = 1024 (FIXED)
  sC: 128 reads per cell = nbk * nbi * nbj = 16*2*4 = 128 (FIXED)

The only way to reduce reads is to change the computation structure.

Can we use a different multiply-accumulate pattern?

IDEA: UNROLLED ACCUMULATION without tmp
For bk=0: mul sC, sA, sB  (initialize sC)
For bk=1: we need sC += sA*sB. Currently: mul tmp, sA, sB; add sC, tmp.
What if instead: add sC, sA; mul sC, sB; sub sC, sB (WRONG ALGEBRA)

Actually: fused multiply-add would be perfect: sC += sA*sB in one instruction.
But the IR only has: add dest,src (in-place) and mul dest,s1,s2 (3-op).
There's no fused MAC instruction.

Is there a way to implement MAC in 1 instruction? Not with current IR.

IDEA: What if we eliminate tmp by using a DIFFERENT ACC PATTERN?
For example:
  Copy sC to tmp: copy tmp, sC  (then tmp = sC, but writes are free, so no cost)
  Wait: "copy tmp, sC" reads sC (costs sC read) and writes tmp (free).
  Then: mul sC, sA, sB  (overwrites sC with product)
  add sC, tmp  (sC = product + old_sC = correct!)
  But this reads sC ONCE (in the copy) and tmp ONCE (in the add).
  Total reads: sC(1) + sA(1) + sB(1) + sC(0, already have product) + tmp(1) = 4 reads
  vs current: sA(1) + sB(1) [in mul] + sC(1) + tmp(1) [in add] = 4 reads
  SAME! But now we're reading sC instead of tmp. sC costs more than tmp (3+ vs 2).
  WORSE.

IDEA: Can we avoid the tmp read in 'add sC, tmp' by pre-computing something?
  add sC, [expression that evaluates to sA*sB without reading it from memory]
  No -- we already computed sA*sB and stored in tmp; must read it.

COMPLETELY DIFFERENT IDEA: What if we increase Ti dramatically?
For Ti=16, Tj=4: nbi=1, nbj=4, nbk=16.
  sA_cache reads = 16*16*4*1*4 = 4096 (same)
  tmp reads = 16*15*4*1*4 = 3840 (same)
  sB reads = 16*16*1*4 = 1024 per cell (same)
  sC: 64 cells at addr 7-70.
  sC cost = 128 * sum_{k=0}^{63} cost(7+k) > current 20,736

For Ti=16, nbi=1: A is loaded only nbi*nbj=4 times (vs nbi*nbj=8 for Ti=8).
But nbi*nbj = 1*4 = 4 (vs 2*4 = 8), so A reads = Ti*nbk*nbi*nbj... wait:

sC reads per cell = nbk * nbi * nbj = 16 * 1 * 4 = 64 (vs 16*2*4=128 for Ti=8)?
No: for each (bi,bj,bk,ii,jj): sC is read once (for bk>0) and once for copy-out.
Total per sC cell = (nbk-1)*nbi*nbj + nbi*nbj = nbk*nbi*nbj.
For Ti=16, Tj=4: 16*1*4 = 64. For Ti=8, Tj=4: 16*2*4 = 128.
DIFFERENT! Ti=16 has FEWER sC reads per cell but MORE cells (64 vs 32).
Total sC reads: Ti*Tj*nbk*nbi*nbj = 16*4*16*1*4 = 4096 (always same, as noted).

For Ti=16, Tj=4, sC cells 7..70 (64 cells):
  cost sum = sum_{k=0}^{63} cost(7+k) = cost(7)+...+cost(70)
  cost(7)=3,...,cost(9)=3, cost(10)=4,...,cost(16)=4, cost(17)=5,...,cost(25)=5,
  cost(26)=6,...,cost(36)=6, cost(37)=7,...,cost(49)=7, cost(50)=8,...,cost(64)=8, cost(65)=9,...

  Let me compute: cost(n) = isqrt(n-1)+1
  cost(7)=3(6→2),cost(8)=3,cost(9)=3,cost(10)=4(9→3),cost(16)=4(15→3),cost(17)=5(16→4),
  cost(25)=5(24→4),cost(26)=6(25→5),cost(36)=6(35→5),cost(37)=7(36→6),cost(49)=7(48→6),
  cost(50)=8(49→7),cost(64)=8(63→7),cost(65)=9,...

  addr 7-9 (3 cells): cost 3
  addr 10-16 (7 cells): cost 4
  addr 17-25 (9 cells): cost 5
  addr 26-36 (11 cells): cost 6
  addr 37-49 (13 cells): cost 7
  addr 50-64 (15 cells): cost 8
  addr 65-70 (6 cells): cost 9

  Sum = 3*3 + 7*4 + 9*5 + 11*6 + 13*7 + 15*8 + 6*9
      = 9 + 28 + 45 + 66 + 91 + 120 + 54 = 413

  sC cost for Ti=16, Tj=4: 64 * avg_cost = but unequal reads!
  Each cell has 64 reads (nbk*nbi*nbj = 16*1*4 = 64).
  sC cost = 64 * sum = 64 * 413 = 26,432
  vs Ti=8, Tj=4: 128 * 162 = 20,736.
  MUCH WORSE for Ti=16.

So Ti=8 is better because it has fewer sC cells AND more reads per cell (lower avg cost).

DEEP INSIGHT: Total sC cost = (N/Ti) * (N/Tj) * nbk * nbi * nbj * [sum of costs of Ti*Tj cells]
= N^2/Ti/Tj * nbk * nbi * nbj * sum_{k=0}^{Ti*Tj-1} cost(sC_base+k)
= 256/Ti/Tj * N * (N/Ti) * (N/Tj) * sum_of_costs
= Hmm let me redo:
Total sC reads = Ti * Tj * nbk * nbi * nbj = Ti*Tj*N*(N/Ti)*(N/Tj) = Ti*Tj*N*N/Ti/Tj * N = N^3? No.
= Ti*Tj * N * (N/Ti) * (N/Tj) = N * N/Ti/Tj * Ti*Tj * N? Let me be careful:
nbk=N, nbi=N/Ti, nbj=N/Tj.
Total sC reads = Ti*Tj * nbk * nbi * nbj = Ti*Tj * N * (N/Ti) * (N/Tj) = Ti*Tj * N^3/(Ti*Tj) = N^3 = 4096.
So total sC reads = 4096 always, but the distribution across cells determines cost.

We want sum_{k=0}^{Ti*Tj-1} (4096/(Ti*Tj)) * cost(sC_base + k) to be minimal.
= (4096/Ti/Tj) * sum_{k=0}^{Ti*Tj-1} cost(sC_base + k)

With sC_base = 3+Tj (after sA_cache, tmp, sB@3..(Tj+2)), the sum over Ti*Tj cells starting at 3+Tj.

To minimize: use FEWER cells (larger Ti*Tj doesn't help if each cell has fewer reads).
With 4096 total reads split across Ti*Tj cells, the key is:
  If Ti*Tj = 1: 4096 reads at addr sC_base (say addr 8): 4096 * cost(8) = 12,288
  If Ti*Tj = 2: 2048 reads each at addr sC_base, sC_base+1: 2048*(3+3) = 12,288 (addr 7-8)
  Wait: sC_base = 3+Tj = 3+1 = 4 for Tj=1. So addr 4 (cost 2).
  Ti*Tj=1, Tj=1, Ti=16: sC at addr 4 (cost 2): 4096*2 = 8192.

  But wait: nbk*nbi*nbj = 16*1*16 = 256 reads per sC cell? Let me recheck.
  For Ti=16, Tj=1: nbi=1, nbj=16, nbk=16.
  Reads per sC cell = nbk * nbi * nbj = 16*1*16 = 256? No!
  Oh wait: sC(ii,jj) is read once per bk>0 (15 times) + once for copy-out = 16 times per (bi,bj).
  Over all (bi,bj): nbi*nbj = 1*16 = 16 pairs.
  So reads per sC cell = 16 * 16 = 256. But Ti*Tj * 256 = 16*1*256 = 4096. ✓

  For Ti*Tj=1: sC_base = 4 (Tj=1, addr = 3+1+ii*1+jj = 4).
  sC cost = 16 * 256 * cost(4) = 4096 * 2 = 8,192!

Hmm wait: for Ti=16, Tj=1: sC is Ti*Tj=16 cells! Addr 4..19 (Ti=16 cells).
Each sC cell has 256 reads.
sC cost = sum_{ii=0}^{15} 256 * cost(4+ii) = 256 * sum_{k=4}^{19} cost(k)
= 256 * (2+2+3+3+3+3+3+4+4+4+4+4+4+5+5+5) = 256 * 58 = 14,848

For Ti=8, Tj=4: sC is 32 cells at addr 7-38, 128 reads each.
sC cost = 128 * 162 = 20,736.

For Ti=16, Tj=1: sC cost = 14,848 (much better!).
But Ti=16, Tj=1 overall score was 78,291 — because of much higher bulk costs.

So sC is better for Ti=16,Tj=1 but bulk is worse. The trade-off favors Ti=8,Tj=4 overall.

THE KEY LEVERAGE: Is there any way to reduce the sC cost for Ti=8,Tj=4 specifically?

For Ti=8, Tj=4 with sA_cache=1, tmp=2, sB@3-6, sC@7-38:
  sC cells: 32 cells at addr 7-38, 128 reads each.
  sC cost = 128 * sum_{k=7}^{38} cost(k) = 128 * 162 = 20,736.

Could we put sC at LOWER addresses? The lowest available is addr 7 (after sA,tmp,sB).
Not possible to go lower without moving sB or removing other cells.

What if we use FEWER sC cells with MORE reads per cell?
E.g., Tj=2 with Ti=8: sC is 16 cells at addr 5-20, 256 reads each.
sC cost = 256 * sum_{k=5}^{20} cost(k) = 256 * (3+3+3+4+4+4+4+4+4+5+5+5+5+5+5+5) = 256 * 68 = 17,408.
vs Ti=8,Tj=4: 20,736. BETTER sC, but overall Ti=8,Tj=2 score was 78,294 (worse than 73,602).

So reducing sC helps but increasing nbj (more bj blocks → more bulk B reads) hurts more.

The ONLY remaining lever I haven't exploited: INSTRUCTION REUSE.

For a specific (bi,bj) pair, can we reuse computed products from a neighboring (bi,bj) pair?
No — each tile is independent and must be computed from scratch.

WAIT: The mul instruction for bk=0 writes directly to sC. For bk>0, writes to tmp then adds to sC.
What if for the VERY FIRST (bk=0, ii=0, jj=0): we write to sC(0,0) = addr 7.
What if sC(0,0) IS addr 2 (tmp)? Then for bk>0, writing to tmp overwrites sC(0,0)!
We'd need to save sC(0,0) before using tmp. Adding a copy reads sC(0,0) at cost 7. Worse.

OK I think 73,602 might be near-optimal for the parametric family we're exploring.

Let me try a few more creative structural approaches:
1. Process two bj columns simultaneously, sharing sA reads
2. Different loop orders (bj > bi > bk > ii > jj)
3. What if we precompute A*B partial sums in a different order?
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 73_602
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


def generate_two_bj_parallel(Ti: int) -> str | None:
    """
    Process TWO bj columns simultaneously, sharing sA reads.
    For Tj=4: nbj=4. Process bj=0 and bj=2 together (or 0&1, etc.)

    Layout: sA_cache=1, tmp=2, sB_left@3-6, sB_right@7-10, sC_left@11-42, sC_right@43-74
    For each (bi, bk): load sA once for ALL ii.
      But sB requires TWO separate loads (sB_left and sB_right).
      For each ii: read sA_cache, mul with sB_left, mul with sB_right.

    This doubles sB loads but halves sA loads (load A[i,bk] once, not twice).

    With Ti=8, Tj=4, doing 2 bj blocks at once:
      nbi=2, nbj=4 (now process 2 at a time → nbj_effective=2 outer iterations)
      nbk=16

      sA_cache reads: Ti * nbk * Tj * nbi * nbj_effective * 2 = 8*16*4*2*2*2 = 4096
      No wait: for each (bi, bj_pair, bk, ii): read sA_cache (Tj+Tj) = 8 times (for 2 Tj groups).
      So sA reads = Ti * nbk * (2*Tj) * nbi * nbj_effective = 8*16*8*2*2 = 4096. SAME!

    sA reads are the same regardless. No benefit.

    Actually the POINT of two-column is: can we eliminate one sA load by sharing?
    Current: for each (bi, bj, bk, ii): copy sA_cache, A[...]; then Tj reads from sA_cache.
    If two-column: for each (bi, bj_pair, bk, ii): copy sA_cache once, then 2*Tj reads from sA_cache.
    Copies: nbi * (nbj/2) * nbk * Ti = 2 * 2 * 16 * 8 = 512 (vs 2*4*16*8 = 1024 current).
    READS from sA_cache: 2*Ti * nbk * Tj * nbi * (nbj/2) = 2*8*16*4*2*2 = 4096 (same).

    The COPIES of sA (writes to addr 1) are free! So fewer copies = no cost change.
    But the READS from sA_cache remain the same: 4096.

    HOWEVER: the tmp reads might change!
    For each (bi, bj_pair, bk, ii, jj, col): tmp reads for bk>0.
    Total tmp reads still = Ti * (nbk-1) * (2*Tj) * nbi * (nbj/2) = 8*15*8*2*2 = 3840. SAME.

    So absolutely no benefit. The cost depends only on reads, not copies (writes).
    """
    return None  # No benefit expected


def generate_sA1_tmp2_bj_bi_swap(Ti: int, Tj: int) -> str | None:
    """
    Try loop order: bi > bk > bj > ii > jj (moving bk before bj).

    Impact: A[bi*Ti+ii, bk] is read once per (bi,bk,ii) for ALL bj.
    But then sB must be reloaded for each bj under the same bk.
    Effectively same reads as before? Yes. But is sC handled differently?

    For (bi,bk,bj,ii,jj): accumulate into sC(ii,jj).
    But with bk outer: sC must persist across all bj iterations (can't reuse sC for different bj).
    This means we need FULL sC for the tile... same as before!

    Unless: we use a SEPARATE sC for each bj and accumulate directly to bulk C.
    That's the "accumulate to bulk" approach (very expensive, as analyzed).

    NOT USEFUL.
    """
    return None


def generate_with_different_sC_assignment(Ti: int, Tj: int) -> str | None:
    """
    Use non-contiguous sC addresses: assign sC cells to OPTIMAL positions.

    For Ti=8, Tj=4: 32 sC cells, each with 128 reads.
    With sA_cache=1, tmp=2, sB@3-6: lowest 6 addrs used.
    Next available: addr 7, 8, 9, ...
    All 32 sC cells have the same read count (128), so we just want them at
    the 32 lowest available addresses: 7-38. This IS the current layout.

    But what if we allow sC cells to share addresses with sB or tmp?
    Not possible (writes would corrupt data).

    CONCLUSION: Contiguous sC@7-38 is already optimal.
    """
    return None


def generate_fused_bk(Ti: int, Tj: int, fuse_k: int = 2) -> str | None:
    """
    Fuse fuse_k consecutive bk iterations to avoid tmp for the first.

    For fuse_k=2: process bk=0,1 together without tmp for bk=0.
    Sequence for ii=0, jj=0:
      mul sC(0,0), sA_0, sB_0(0)   [bk=0: direct init]
      copy sA_1, A[bi*Ti, 1]        [load A for bk=1]
      mul tmp, sA_1, sB_1(0)
      add sC(0,0), tmp               [bk=1: accumulate]
    vs current:
      mul sC(0,0), sA_0, sB_0(0)   [bk=0]
      mul tmp, sA_1, sB_1(0)        [bk=1]
      add sC(0,0), tmp

    Actually: the reads are IDENTICAL. The difference is only in code structure, not reads.

    What if we combine mul+add for bk=0 and bk=1 into a single mul_add (impossible in IR)?

    What if we process two consecutive ii at once?
    For ii=0: copy sA_0, A[...]; for jj: mul sC(0,jj), sA_0, sB(jj)
    For ii=1: copy sA_1, A[...]; for jj: mul sC(1,jj), sA_1, sB(jj)

    Could we save reads by sharing sA between ii=0 and ii=1? No -- they're different A cells.

    NOT USEFUL.
    """
    return None


def generate_transposed_B(Ti: int, Tj: int) -> str | None:
    """
    What if B is stored transposed in bulk? B_transposed[j,k] = B[k,j].

    For loading sB: copy sB(jj), B_T[bj*Tj+jj, bk] instead of B[bk, bj*Tj+jj].
    The ADDRESS of B_T[j,k] = B_T_base + j*N + k.

    B_transposed row: j from 0..15 (= Tj*nbj addresses first dimension).
    If B_transposed is stored at higher addresses:
    B_T_base = A_base + 256 = 39 + 256 = 295 (same as B currently).
    B_T[j,k] = 295 + j*16 + k.
    B[k,j] = 295 + k*16 + j.

    The ADDRESSES are different! B_T storage means the j-major order.
    For sB loading: B_T[bj*Tj+jj, bk] = 295 + (bj*Tj+jj)*16 + bk.
    For fixed bj,bk: addresses are 295 + bj*Tj*16 + jj*16 + bk.
    These are NOT sequential (stride 16 vs stride 1).
    The ACTUAL ADDRESSES accessed for each (bj, bk): spread out across B region.

    Does transposing change the cost? The cost depends only on the address value.
    With transposed B: for bj=0, bk=0: B_T addresses = 295+0, 295+16, 295+32, 295+48.
    With normal B: for bj=0, bk=0: B addresses = 295, 295+1, 295+2, 295+3.
    Normal B addresses are LOWER (295-298 vs 295, 311, 327, 343).

    So transposed B is WORSE (higher addresses for the same data).
    Current B layout is optimal: row-major so bk varies slowly.
    """
    return None


def generate_alternative_sB_layout(Ti: int, Tj: int, sB_offset: int) -> str | None:
    """
    Try sB at different starting addresses, allowing non-standard layout.
    E.g., sB@1, 2, 3, 4 (before sA_cache). But then sA_cache needs addr 5.
    sA_cache at addr 5 (cost 3) vs addr 1 (cost 1): penalty = 4096*(3-1) = 8192. Much worse.

    What about: sA_cache=1, sB@2-5 (shift sB down by 1), tmp=6?
    sB at addr 2-5: cost 2,2,3,3 (same as addr 3-6).
    tmp at addr 6 (cost 3) vs addr 2 (cost 2): penalty = 3840*(3-2) = 3840. WORSE.
    """
    return None


def generate_split_sA(Ti: int, Tj: int) -> str | None:
    """
    Split A into multiple cached registers.

    For Ti=8, Tj=4: what if we process ii=0..3 in one pass and ii=4..7 in another?
    Each half-pass: Ti/2=4 sA cells at addr 1-4, tmp at addr 5, sB@6-9, sC@10-25 (16 cells).

    Two half-passes: process bi=0,ii=0..3 then bi=0,ii=4..7.
    Each has nbi_half=4 bi blocks (since we're splitting Ti).

    Hmm, this is equivalent to Ti=4 tiling! Which we already know gives cost 79,520.
    No benefit.
    """
    return None


def generate_reorder_jj_first(Ti: int, Tj: int) -> str | None:
    """
    Alternative: load ALL Tj sB elements, then for each ii do all jj.
    This is the current approach. What if jj is outermost?

    Loop: bi > bj > bk > jj > ii
    For each (bk, jj): sB_cache = B[bk, bj*Tj+jj] at addr 3 (cost 2).
    For each ii: sA_cell(ii) at addr 4..11 (Ti=8 cells, cost 2-4).
                mul into sC(ii, jj).

    sB_cache at addr 3: cost 2, 4096 reads.
    sA at addr 4-11: 3840 reads each? No.
    sA[ii] read for ALL jj (but here jj is outer, so sA read once per (bk,jj,ii)):
    Actually: sA_cell(ii) reads = Tj * nbk * nbi * nbj = 4*16*2*4 = 512 per ii-cell.
    8 ii cells × 512 reads = 4096 total sA reads (same).

    But now sA is spread across 8 addresses (3..10 or 4..11) at cost 2-4.
    vs current sA at addr 1 (4096 reads, cost 1).
    sA cost = 512*(2+2+2+3+3+3+3+4) [addr 4-11] = 512*22 = 11,264
    vs current sA_cache cost = 4096*1 = 4096.
    WORSE by 7,168.

    NO BENEFIT.
    """
    return None


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    print("=== Analytical verification: current layout is optimal ===")

    # Compute optimal assignment cost for Ti=8, Tj=4
    Ti, Tj = 8, 4
    N_val = 16
    nbi, nbj, nbk = N_val // Ti, N_val // Tj, N_val

    # Read counts (fixed)
    read_counts = {
        "sA_cache": Ti * nbk * Tj * nbi * nbj,   # = 4096
        "tmp": Ti * (nbk - 1) * Tj * nbi * nbj,   # = 3840
    }
    for jj in range(Tj):
        read_counts[f"sB[{jj}]"] = Ti * nbk * nbi * nbj   # = 1024 each

    for ii in range(Ti):
        for jj in range(Tj):
            read_counts[f"sC[{ii},{jj}]"] = nbk * nbi * nbj  # = 128 each

    # Sort by read count descending
    items = sorted(read_counts.items(), key=lambda x: -x[1])

    print(f"Items sorted by read count:")
    for name, reads in items[:10]:
        print(f"  {name}: {reads} reads")
    print(f"  ... ({len(items)} total items)")

    # Assign to lowest-cost addresses in order
    total_cost = 0
    addr = 1
    print(f"\nOptimal assignment:")
    for name, reads in items:
        c = addr_cost(addr)
        total_cost += reads * c
        addr += 1
    print(f"  Optimal scratchpad cost = {total_cost:,}")

    # Add bulk costs (A,B,C) - same regardless of scratch layout
    scratch_end = len(items) + 1
    A_base = scratch_end
    B_base = A_base + N_val * N_val
    C_base = B_base + N_val * N_val

    # A read nbi*nbj times each (= 8 times)
    bulk_A_cost = sum(nbi * nbj * addr_cost(A_base + i) for i in range(N_val * N_val))
    # B read nbi times each (= 2 times)
    bulk_B_cost = sum(nbi * addr_cost(B_base + i) for i in range(N_val * N_val))
    # C read 1 time each
    bulk_C_cost = sum(addr_cost(C_base + i) for i in range(N_val * N_val))

    print(f"  Bulk A cost = {bulk_A_cost:,}")
    print(f"  Bulk B cost = {bulk_B_cost:,}")
    print(f"  Bulk C cost = {bulk_C_cost:,}")
    print(f"  Total = {total_cost + bulk_A_cost + bulk_B_cost + bulk_C_cost:,}")
    print(f"  (vs actual 73,602)")

    # Now compute current layout cost analytically
    print(f"\n=== Current layout analysis ===")
    current_scratch_cost = (
        4096 * addr_cost(1) +   # sA_cache=1
        3840 * addr_cost(2) +   # tmp=2
        sum(1024 * addr_cost(3 + jj) for jj in range(Tj)) +   # sB@3-6
        sum(128 * addr_cost(7 + ii * Tj + jj) for ii in range(Ti) for jj in range(Tj))   # sC@7-38
    )
    print(f"  Current scratchpad cost = {current_scratch_cost:,}")
    print(f"  Optimal scratchpad cost = {total_cost:,}")
    print(f"  Difference = {current_scratch_cost - total_cost:,}")

    # What IS the optimal assignment?
    print(f"\n=== Optimal vs current assignment (first few) ===")
    addr = 1
    for name, reads in items[:8]:
        c = addr_cost(addr)
        print(f"  addr={addr} (cost={c}): {name} ({reads} reads)")
        addr += 1
