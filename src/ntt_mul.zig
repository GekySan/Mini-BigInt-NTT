const std = @import("std");
const Allocator = std.mem.Allocator;

const kBaseExp: u8 = 19;
const kBase: u64 = 10_000_000_000_000_000_000;

const Limb = u64;
const dLimb = u128;

const kFmtLimb1 = "{d}";
const kFmtLimb = "{d:0>19}";

const BigInt = std.ArrayList(Limb);

/// Supprime les blocs (limbs) de poids fort qui sont à zéro d'un BigInt.
fn BigTrim(self: *BigInt) void {
    while (self.items.len > 0 and self.items[self.items.len - 1] == 0) {
        _ = self.pop();
    }
}

/// Affiche la représentation décimale d'un BigInt dans le writer fourni.
fn BigPrint(self: BigInt, writer: anytype) !void {
    if (self.items.len == 0) {
        try writer.print("0", .{});
        return;
    }
    try writer.print(kFmtLimb1, .{self.items[self.items.len - 1]});
    var i: usize = self.items.len - 1;
    while (i > 0) {
        i -= 1;
        try writer.print(kFmtLimb, .{self.items[i]});
    }
}

/// Analyse une chaîne de caractères décimale et crée un BigInt.
fn BigFromDecimalString(allocator: Allocator, s: []const u8) !BigInt {
    if (s.len == 0) return BigInt.init(allocator);

    const chunks = (s.len + kBaseExp - 1) / kBaseExp;
    var res = try BigInt.initCapacity(allocator, chunks);
    try res.resize(chunks);

    var i: usize = 0;
    while (i < chunks) : (i += 1) {
        const end = s.len - i * kBaseExp;
        const start = if (end >= kBaseExp) end - kBaseExp else 0;
        const chunkStr = s[start..end];
        res.items[i] = try std.fmt.parseUnsigned(u64, chunkStr, 10);
    }
    BigTrim(&res);
    return res;
}

/// Calcule le plafond du logarithme en base 2 d'un entier u64.
inline fn CeilLog2U64(a: u64) u32 {
    return if (a <= 1) 0 else 64 - @clz(a - 1);
}

/// Effectue une multiplication modulaire (a * b) % m.
inline fn MulMod(a: u64, b: u64, m: u64) u64 {
    return @truncate((@as(dLimb, a) * @as(dLimb, b)) % @as(dLimb, m));
}

/// Calcule (a ^ e) % m en utilisant l'exponentiation modulaire.
inline fn PowMod(a: u64, e: u32, m: u64) u64 {
    var res: u64 = 1;
    var base = a;
    var exp = e;
    while (exp > 0) {
        if (exp & 1 == 1) res = MulMod(res, base, m);
        base = MulMod(base, base, m);
        exp >>= 1;
    }
    return res;
}

const kNbMods = 3;
const kNttProot2Exp = 53;

const kNttMods = [_]Limb{ 0x280000000000001, 0x660000000000001, 0xbe0000000000001 };
const kNttProot = [_][kNbMods]Limb{
    .{ 0xcfd41b9100000, 0x11ce6ee1f1b68ff, 0x71fba0101cb9d56 },
    .{ 0x1034f97971ea68a, 0x4442bba3b520d83, 0x835e9c6fae891c2 },
};
const kNttModsCr = [_]Limb{ 0x170842108421086, 0xcaaaaaaaaaaaac, 0x3822e8ba2e8ba31 };

var s_nttProotPow: [2][kNbMods][kNttProot2Exp]Limb = undefined;

/// Effectue une addition modulaire (a + b) % m.
inline fn AddMod(a: Limb, b: Limb, m: Limb) Limb {
    const r = a +% b -% m;
    return if (r > a) r +% m else r;
}

/// Effectue une soustraction modulaire (a - b) % m.
inline fn SubMod(a: Limb, b: Limb, m: Limb) Limb {
    const r = a -% b;
    return if (r > a) r +% m else r;
}

const NttLimb = Limb;

/// Effectue une Transformée Numérique Théorique (NTT) en place.
fn NttFft(outBuf: []NttLimb, inBuf: []NttLimb, k: u32, proot: Limb, m: Limb) void {
    const n = @as(Limb, 1) << @as(u6, @intCast(k));
    var nbBlocks: Limb = n;
    var fftPerBlock: Limb = 1;
    const strideIn: Limb = n / 2;
    var tabIn = inBuf;
    var tabOut = outBuf;
    var cMul = proot;

    while (nbBlocks != 1) {
        nbBlocks >>= 1;
        if (nbBlocks == 1) tabOut = outBuf;
        var p: Limb = 0;
        var idx: Limb = 0;
        var c: Limb = 1;
        var b_idx: Limb = 0;
        while (b_idx < nbBlocks) : (b_idx += 1) {
            var j: Limb = 0;
            while (j < fftPerBlock) : (j += 1) {
                const a0 = tabIn[idx + j];
                const a1 = tabIn[idx + j + strideIn];
                const b0 = AddMod(a0, a1, m);
                var b1 = SubMod(a0, a1, m);
                b1 = MulMod(b1, c, m);
                tabOut[p + j] = b0;
                tabOut[p + j + fftPerBlock] = b1;
            }
            c = MulMod(c, cMul, m);
            idx += fftPerBlock;
            p += 2 * fftPerBlock;
        }
        cMul = MulMod(cMul, cMul, m);
        fftPerBlock <<= 1;
        std.mem.swap([]NttLimb, &tabIn, &tabOut);
    }
}

/// Effectue une multiplication point par point et met à l'échelle pour l'iFFT.
fn NttVecMul(tab1: []NttLimb, tab2: []const NttLimb, k: u32, m: Limb) void {
    const cInv2 = (m + 1) / 2;
    const cInv = PowMod(cInv2, k, m);
    const n = @as(Limb, 1) << @as(u6, @intCast(k));
    for (0..n) |i| {
        const a = MulMod(tab1[i], tab2[i], m);
        tab1[i] = MulMod(a, cInv, m);
    }
}

/// Prépare un BigInt pour la NTT.
fn LimbToNtt(dst: []NttLimb, fft_len: u64, src: []const Limb, midx: usize) void {
    const m = kNttMods[midx];
    for (src, 0..) |s, i| dst[i] = s % m;
    @memset(dst[src.len..@intCast(fft_len)], 0);
}

/// Reconstruit un BigInt à partir des résultats de la NTT via CRT.
fn NttToLimb(out: []Limb, rLen: u64, buf: []const NttLimb, k: u32) void {
    const fftLen = @as(u64, 1) << @as(u6, @intCast(k));
    var carry: [kNbMods]Limb = .{0} ** kNbMods;

    for (0..rLen) |i| {
        var y: [kNbMods]Limb = undefined;
        for (0..kNbMods) |j| y[j] = buf[@intCast(i + fftLen * j)];

        var l: usize = 0;
        var j: usize = 0;
        while (j < kNbMods - 1) : (j += 1) {
            var t: usize = j + 1;
            while (t < kNbMods) : (t += 1) {
                const m = kNttMods[t];
                const dj = SubMod(y[t], y[j], m);
                y[t] = MulMod(dj, kNttModsCr[l], m);
                l += 1;
            }
        }

        var u: [kNbMods]Limb = undefined;
        u[0] = y[kNbMods - 1];
        l = 1;
        var j_signed: isize = kNbMods - 2;
        while (j_signed >= 0) : (j_signed -= 1) {
            const currentJ = @as(usize, @intCast(j_signed));
            var R = y[currentJ];
            for (0..l) |t| {
                const p: dLimb = @as(dLimb, u[t]) * @as(dLimb, kNttMods[currentJ]) + @as(dLimb, R);
                u[t] = @truncate(p % kBase);
                R = @truncate(p / kBase);
            }
            u[l] = R;
            l += 1;
        }

        var c: u64 = 0;
        for (0..kNbMods) |modIdx| {
            const sum: dLimb = @as(dLimb, u[modIdx]) + @as(dLimb, carry[modIdx]) + @as(dLimb, c);
            u[modIdx] = @truncate(sum % kBase);
            c = @truncate(sum / kBase);
        }
        out[@intCast(i)] = u[0];
        for (0..kNbMods - 1) |modIdx| carry[modIdx] = u[modIdx + 1];
    }
}

/// Précalcule les tables pour la NTT.
fn NttInit() void {
    for ([_]bool{ false, true }, 0..) |_, inv| {
        for (0..kNbMods) |j| {
            var c = kNttProot[inv][j];
            for (0..kNttProot2Exp) |i| {
                s_nttProotPow[inv][j][i] = c;
                c = MulMod(c, c, kNttMods[j]);
            }
        }
    }
}

/// Multiplication «école», complexité (O(n^2)).
fn BigMulSchoolbook(r: *BigInt, a: BigInt, b: BigInt, allocator: Allocator) !void {
    _ = allocator;
    if (a.items.len == 0 or b.items.len == 0) {
        r.shrinkRetainingCapacity(0);
        return;
    }
    const n = a.items.len;
    const m = b.items.len;
    const rl = n + m + 1;
    try r.resize(rl);
    @memset(r.items, 0);

    for (a.items, 0..) |aLimb, i| {
        var carry: u64 = 0;
        for (b.items, 0..) |bLimb, j| {
            var p: dLimb = @as(dLimb, aLimb) * @as(dLimb, bLimb);
            p += @as(dLimb, r.items[i + j]);
            p += @as(dLimb, carry);
            r.items[i + j] = @truncate(p % kBase);
            carry = @truncate(p / kBase);
        }
        r.items[i + m] += carry;
    }
    var c: u64 = 0;
    for (0..rl) |k| {
        const sum: dLimb = @as(dLimb, r.items[k]) + @as(dLimb, c);
        r.items[k] = @truncate(sum % kBase);
        c = @truncate(sum / kBase);
    }
    BigTrim(r);
}

/// Multiplication via NTT.
fn BigMulNtt(r: *BigInt, a: BigInt, b: BigInt, allocator: Allocator) !void {
    if (a.items.len == 0 or b.items.len == 0) {
        r.shrinkRetainingCapacity(0);
        return;
    }
    const outLen: u64 = @intCast(a.items.len + b.items.len);
    const k = CeilLog2U64((outLen * kBaseExp + kBaseExp - 1) / kBaseExp);
    std.debug.assert(k <= kNttProot2Exp);
    const n: u64 = @as(u64, 1) << @as(u6, @intCast(k));

    var buf1 = try allocator.alloc(NttLimb, @intCast(n * kNbMods));
    defer allocator.free(buf1);
    const buf2 = try allocator.alloc(NttLimb, @intCast(n));
    defer allocator.free(buf2);
    const buf3 = try allocator.alloc(NttLimb, @intCast(n));
    defer allocator.free(buf3);

    const rootIdx = kNttProot2Exp - k;
    for (0..kNbMods) |j| {
        const buf1Slice = buf1[n * @as(u64, @intCast(j)) .. n * @as(u64, @intCast(j + 1))];
        LimbToNtt(buf1Slice, n, a.items, j);
        LimbToNtt(buf2, n, b.items, j);
        const m = kNttMods[j];
        NttFft(buf3, buf1Slice, k, s_nttProotPow[0][j][rootIdx], m);
        NttFft(buf1Slice, buf2, k, s_nttProotPow[0][j][rootIdx], m);
        NttVecMul(buf3, buf1Slice, k, m);
        NttFft(buf1Slice, buf3, k, s_nttProotPow[1][j][rootIdx], m);
    }

    try r.resize(@intCast(outLen + 1));
    @memset(r.items, 0);
    NttToLimb(r.items, outLen, buf1, k);
    BigTrim(r);
}

// Exemple, même si que cet exemple est mauvais.
pub fn main() !void {
    NttInit();
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const stdout = std.io.getStdOut().writer();

    const stringA = "1234567890123456789345678901234567891234567890123456789123456789012345678912345678901234567891234567890123456789123456789017893";
    const stringB = "9876543210987654321321098765432198765432109876543219876543210987654321987654321098765432198765432109876543219876543210987654321";

    const A = try BigFromDecimalString(allocator, stringA);
    const B = try BigFromDecimalString(allocator, stringB);
    var R = BigInt.init(allocator);

    try stdout.print("NTT Multiplication Result :\n", .{});
    try BigMulNtt(&R, A, B, allocator);
    try BigPrint(R, stdout);
    try stdout.print("\n", .{});
}
