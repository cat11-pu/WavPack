#!/usr/bin/env python3
import os, random, struct, subprocess, sys, shutil, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED_DIR = os.path.join(ROOT, "fuzzing", "regression")
CORPUS = os.path.join(ROOT, "fuzz-local", "corpus")
CRASH = os.path.join(ROOT, "fuzz-local", "crashes")
HARNESS = os.path.join(ROOT, "fuzz-local", "harness-asan")
os.makedirs(CORPUS, exist_ok=True); os.makedirs(CRASH, exist_ok=True)

def parse_blocks(data):
    blocks = []
    i = 0
    while i + 32 <= len(data):
        j = data.find(b"wvpk", i)
        if j < 0 or j + 32 > len(data): break
        cksize = struct.unpack_from("<I", data, j + 4)[0]
        end = j + 8 + cksize
        if end > len(data):
            blocks.append((j, len(data), cksize, True))
            break
        blocks.append((j, end, cksize, False))
        i = end
    return blocks

def parse_metadata(block, hstart):
    # returns list of (id_byte_pos, len_field_pos, data_start, data_end)
    cksize = struct.unpack_from("<I", block, hstart + 4)[0]
    end = hstart + 8 + min(cksize, len(block) - hstart - 8)
    p = hstart + 32
    out = []
    while p + 2 <= end:
        idpos = p
        mid = block[p]; l1 = block[p+1]; p += 2
        length = l1 << 1
        if mid & 0x80:  # ID_LARGE
            if p + 2 > end: break
            length += block[p] << 9; length += block[p+1] << 17; p += 2
        dstart = p
        dend = dstart + length + (length & 1)
        if dend > end: break
        out.append((idpos, dstart, dend, mid & 0x3f, length))
        p = dend
    return out, end

def mutate(data, rng):
    d = bytearray(data)
    blocks = parse_blocks(d)
    if not blocks:
        pos = rng.randrange(len(d)); d[pos] ^= 1 << rng.randrange(8); return bytes(d)
    mode = rng.randrange(10)
    b = rng.choice(blocks)
    h, end, cksize, truncated = b
    if mode == 0:  # flip header field byte
        field = rng.choice([4,5,6,7, 8,9,10,11, 12,13,14,15, 16,17, 20,21,22,23, 24,25,26,27])
        if h + field < len(d): d[h+field] = rng.randrange(256)
    elif mode == 1:  # cksize extremes (oversize allocation)
        val = rng.choice([0,1,2,23,24,25,31,32,255, 0xffff, 0xffffff, 0xfffffffe, 0xffffffff])
        struct.pack_into("<I", d, h+4, val)
    elif mode == 2:  # block_samples extremes
        val = rng.choice([0,1,2,255,256, 0x7fff,0x8000,0xffff, 0x7fffffff,0x80000000,0xffffffff])
        struct.pack_into("<I", d, h+12, val)
    elif mode == 3:  # truncate at random point in/around block
        cut = rng.randrange(max(1,h), min(len(d), end+1) if end>h else len(d))
        d = d[:cut]
    elif mode == 4:  # flags/total samples/index dwords
        off = rng.choice([8, 16, 20, 24])
        val = rng.choice([0,1,2,3, 0xff,0xffff,0xffffff,0xffffffff,0x80000000,0x7fffffff])
        struct.pack_into("<I", d, h+off, val)
    elif mode == 5:  # metadata-aware mutation
        blk = bytes(d)
        metas, bend = parse_metadata(blk, h)
        if metas:
            (idpos, ds, de, mid, ln) = rng.choice(metas)
            sm = rng.randrange(4)
            if sm == 0: d[idpos] = rng.randrange(256)
            elif sm == 1 and ln:
                d[rng.randrange(ds, de)] = rng.randrange(256)
            elif sm == 2:
                d[idpos+1] = rng.choice([0,1,2,255])
            elif ln >= 1 and ds < de:
                d[ds] = rng.choice([0,1,2,255])
        else:
            hi = min(len(d), bend)
            if h + 32 < hi:
                pos = rng.randrange(h+32, hi); d[pos] ^= 0xff
            else:
                d[rng.randrange(len(d))] ^= 0xff
    elif mode == 6:  # insert garbage bytes inside block
        pos = rng.randrange(h, min(len(d), end))
        d[pos:pos] = os.urandom(rng.randrange(1, 17))
    elif mode == 7:  # duplicate / delete a block
        if len(blocks) > 1:
            bb = rng.choice(blocks)
            chunk = bytes(d[bb[0]:bb[1]])
            if rng.random() < 0.5:
                d[bb[0]:bb[1]] = b""
            else:
                d[bb[0]:bb[0]] = chunk
    elif mode == 8:  # random byte flip
        for _ in range(rng.randrange(1,5)):
            pos = rng.randrange(len(d)); d[pos] ^= 1 << rng.randrange(8)
    elif mode == 9:  # interesting values at random 4-byte aligned spots
        pos = rng.randrange(0, max(1,len(d)-4))
        struct.pack_into("<I", d, pos, rng.choice([0,1,2,0x7f,0x80,0xff,0x100,0xffff,0x10000,0xffffffff]))
    return bytes(d)

def run_batch(files, env_extra=None):
    env = dict(os.environ)
    env["ASAN_OPTIONS"] = "detect_leaks=0:abort_on_error=1:symbolize=0"
    env["UBSAN_OPTIONS"] = "halt_on_error=1:print_stacktrace=0:symbolize=0"
    if env_extra: env.update(env_extra)
    try:
        p = subprocess.run([HARNESS] + files, capture_output=True, timeout=120, env=env)
        return p.returncode, p.stderr
    except subprocess.TimeoutExpired:
        return -99, b"TIMEOUT"

def main():
    iters = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    batch = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    seeds = [os.path.join(SEED_DIR, f) for f in os.listdir(SEED_DIR)]
    rng = random.Random(0x3619)
    pending = []
    crashes = 0; timeouts = 0; ran = 0
    t0 = time.time()
    for i in range(iters):
        seed = rng.choice(seeds)
        with open(seed, "rb") as f: data = f.read()
        # sometimes stack several mutations
        rounds = rng.choice([1,1,1,2,3])
        for _ in range(rounds): data = mutate(data, rng)
        path = os.path.join(CORPUS, f"mut_{i:06d}.wv")
        with open(path, "wb") as f: f.write(data)
        pending.append(path)
        if len(pending) >= batch or i == iters - 1:
            rc, err = run_batch(pending)
            ran += len(pending)
            if rc != 0:
                if b"TIMEOUT" in err:
                    timeouts += len(pending)
                    for fp in pending: shutil.copy(fp, os.path.join(CRASH, "timeout-" + os.path.basename(fp)))
                else:
                    # bisect which file crashed by running individually
                    for fp in pending:
                        rc2, err2 = run_batch([fp])
                        if rc2 != 0:
                            crashes += 1
                            name = os.path.basename(fp)
                            shutil.copy(fp, os.path.join(CRASH, name))
                            with open(os.path.join(CRASH, name + ".err"), "wb") as ef: ef.write(err2[:4000])
                            print(f"[CRASH] {name} rc={rc2}", flush=True)
            for fp in pending:
                try: os.unlink(fp)
                except OSError: pass
            pending = []
            if i % 500 < batch:
                print(f"progress {ran}/{iters} crashes={crashes} timeouts={timeouts} elapsed={time.time()-t0:.0f}s", flush=True)
    print(f"DONE ran={ran} crashes={crashes} timeouts={timeouts}")

if __name__ == "__main__":
    main()
