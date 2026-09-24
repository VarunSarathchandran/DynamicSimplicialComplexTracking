#!/usr/bin/env python3
"""Report embedding and subsetting of every font in a PDF.

IEEE PDF eXpress requires all fonts embedded (Type 3 discouraged); subsetting
is expected. Usage:  python3 checkfonts.py paper.pdf [more.pdf ...]
"""
import re, sys, zlib

def all_text(data):
    parts = [data]
    for m in re.finditer(rb'stream\r?\n', data):
        s = m.end(); e = data.find(b'endstream', s)
        try: parts.append(zlib.decompress(data[s:e]))
        except Exception: pass
    return b'\n'.join(parts)

def dicts(text, key):
    """Yield the text of every << ... >> that contains `key`, depth-aware."""
    for m in re.finditer(key, text):
        i = m.start()
        # walk back to the opening << of this dictionary
        depth = 0; j = i
        while j > 0:
            if text[j:j+2] == b'>>': depth += 1; j -= 2; continue
            if text[j:j+2] == b'<<':
                if depth == 0: break
                depth -= 1; j -= 2; continue
            j -= 1
        start = j
        depth = 0; k = start
        while k < len(text):
            if text[k:k+2] == b'<<': depth += 1; k += 2; continue
            if text[k:k+2] == b'>>':
                depth -= 1; k += 2
                if depth == 0: break
                continue
            k += 1
        yield text[start:k]

def check(path):
    data = open(path, 'rb').read()
    text = all_text(data)
    # Descriptors tell us embedding; fonts tell us names/types.
    desc = {}
    for d in dicts(text, rb'/Type\s*/FontDescriptor'):
        n = re.search(rb'/FontName\s*/([^\s/>\[\]]+)', d)
        if n:
            desc[n.group(1)] = any(t in d for t in (b'/FontFile', b'/FontFile2', b'/FontFile3'))
    rows = {}
    for d in dicts(text, rb'/Type\s*/Font\b'):
        sub = re.search(rb'/Subtype\s*/(\w+)', d)
        base = re.search(rb'/BaseFont\s*/([^\s/>\[\]]+)', d)
        if not sub or not base: continue
        st, bn = sub.group(1).decode(), base.group(1)
        if st == 'Type0': continue          # composite: descendant carries descriptor
        name = bn.decode(errors='replace')
        subset = bool(re.match(rb'^[A-Z]{6}\+', bn))
        if st == 'Type3': emb = 'type3'
        else: emb = desc.get(bn, None)
        rows[name] = (st, emb, subset)
    print(f"\n{path}")
    print(f"  {'font':45} {'type':10} {'embedded':20} {'subset'}")
    bad = 0
    for name, (st, emb, subset) in sorted(rows.items()):
        e = {True:'yes', False:'NO', None:'NO (no descriptor)', 'type3':'type3'}[emb]
        flag = ''
        if emb is not True: flag = '   <-- PROBLEM' ; bad += 1
        elif not subset: flag = '   <-- not subset (usually OK, but check)'
        print(f"  {name:45} {st:10} {e:20} {'yes' if subset else 'no'}{flag}")
    if not rows: print("  (no simple fonts found — may be all Type0/CID; inspect with pdffonts)")
    print(f"  => {'PASS' if bad == 0 and rows else 'FAIL'}: {bad} font(s) not embedded")
    return bad

if __name__ == '__main__':
    results = [check(p) for p in sys.argv[1:]]   # list, not any(): check every file
    sys.exit(1 if any(results) else 0)
