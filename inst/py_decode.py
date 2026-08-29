#!/usr/bin/env python3
"""py_decode.py -- helper for mzML binary data array decoding.
Called from R via system2() to handle base64 + zlib decompression.

Usage:
  python3 py_decode.py <dtype> <byteorder>  (base64 string read from stdin)

  dtype:    '32' (float32) or '64' (float64)
  byteorder: 'little' or 'big'

The base64 payload is read from stdin rather than argv: a single
high-resolution profile scan's m/z array can encode to hundreds of KB,
well past typical OS command-line length limits (e.g. ~32K on Windows).

Outputs decoded numbers as newline-separated text to stdout.
"""
import sys
import zlib
import base64
import struct

def decode(b64_str, dtype='64', byteorder='little'):
    raw = base64.b64decode(b64_str.strip())
    # Try zlib decompress; if it fails, assume uncompressed
    try:
        raw = zlib.decompress(raw)
    except zlib.error:
        pass  # not compressed
    fmt_char = '<' if byteorder == 'little' else '>'
    if dtype == '32':
        fmt = fmt_char + ('f' * (len(raw) // 4))
        vals = struct.unpack(fmt, raw[:len(raw) // 4 * 4])
    else:
        fmt = fmt_char + ('d' * (len(raw) // 8))
        vals = struct.unpack(fmt, raw[:len(raw) // 8 * 8])
    return vals

if __name__ == '__main__':
    b64 = sys.stdin.read()
    dtype = sys.argv[1] if len(sys.argv) > 1 else '64'
    byteorder = sys.argv[2] if len(sys.argv) > 2 else 'little'
    vals = decode(b64, dtype, byteorder)
    for v in vals:
        sys.stdout.write(f"{v}\n")
