#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "wavpack.h"
#include <string.h>

typedef struct {
    unsigned char ungetc_char, ungetc_flag;
    unsigned char *sptr, *dptr, *eptr;
} RC;

static int32_t rb (void *id, void *data, int32_t bcount) {
    RC *r = id; unsigned char *o = data;
    while (bcount) {
        if (r->ungetc_flag) { *o++ = r->ungetc_char; r->ungetc_flag = 0; bcount--; }
        else {
            size_t n = r->eptr - r->dptr;
            if (!n) break;
            if (n > (size_t)bcount) n = bcount;
            memcpy (o, r->dptr, n); r->dptr += n; o += n; bcount -= (int32_t)n;
        }
    }
    return (int32_t)(o - (unsigned char *)data);
}
static int32_t wb (void *id, void *data, int32_t b) { (void)id;(void)data; return b; }
static int64_t gp (void *id) { RC *r = id; return r->dptr - r->sptr; }
static int spa (void *id, int64_t p) { RC *r = id; if (p < 0 || r->sptr + p > r->eptr) return 1; r->dptr = r->sptr + p; return 0; }
static int spr (void *id, int64_t d, int m) {
    RC *r = id; unsigned char *ref = m==SEEK_SET?r->sptr:m==SEEK_CUR?r->dptr:r->eptr;
    if (ref + d < r->sptr || ref + d > r->eptr) return 1; r->dptr = ref + d; return 0;
}
static int pb (void *id, int c) { RC *r = id; r->ungetc_char = (unsigned char)c; r->ungetc_flag = 1; return c; }
static int64_t gl (void *id) { RC *r = id; return r->eptr - r->sptr; }
static int cs (void *id) { (void)id; return 1; }
static int cc (void *id) { (void)id; return 0; }

static WavpackStreamReader64 reader = { rb, wb, gp, spa, spr, pb, gl, cs, NULL, cc };

int LLVMFuzzerTestOneInput (const uint8_t *data, size_t size) {
    int flags = OPEN_TAGS | OPEN_EDIT_TAGS | OPEN_WRAPPER | OPEN_DSD_AS_PCM | OPEN_NO_CHECKSUM | OPEN_NORMALIZE;
    RC raw; WavpackContext *wpc; char error[128];
    memset (&raw, 0, sizeof raw);
    raw.dptr = raw.sptr = (unsigned char *)data; raw.eptr = raw.dptr + size;
    wpc = WavpackOpenFileInputEx64 (&reader, &raw, NULL, error, flags, 15);
    if (!wpc) return 0;
    int nc = WavpackGetNumChannels (wpc);
    if (nc && nc <= 256) {
        int32_t *buf = malloc (1024 * nc * sizeof (int32_t));
        unsigned char md5[16]; uint32_t r;
        do { r = WavpackUnpackSamples (wpc, buf, 1024); } while (r);
        WavpackGetMD5Sum (wpc, md5);
        free (buf);
    }
    WavpackCloseFile (wpc);
    return 0;
}

int main (int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        FILE *f = fopen (argv[i], "rb");
        if (!f) continue;
        fseek (f, 0, SEEK_END); long n = ftell (f); fseek (f, 0, SEEK_SET);
        if (n <= 0 || n > 4*1024*1024) { fclose(f); continue; }
        unsigned char *b = malloc ((size_t)n);
        if (fread (b, 1, (size_t)n, f) == (size_t)n) LLVMFuzzerTestOneInput (b, (size_t)n);
        free (b); fclose (f);
    }
    return 0;
}
