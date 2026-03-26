#!/usr/bin/env python3
"""Parse Cosmed Fitmate v2.4 databases and dump RMR test results.

Expects these files in the same directory:
  ANAGRAFE.DBF / .fpt  (patients)
  SESSIONI.DBF / .fpt  (sessions)
  TEST.DBF     / .fpt  (test results with binary blobs)

Usage:
  python fitmate_dump.py                # list all patients
  python fitmate_dump.py 1              # dump RMR results for patient 1
  python fitmate_dump.py "DAYAH"        # dump RMR results by last name
"""

import struct, sys, os, statistics
from pathlib import Path

try:
    from dbfread import DBF, FieldParser
except ImportError:
    sys.exit("pip install dbfread")

TIPO_RMR = 26
FIO2 = 20.93  # atmospheric O2 %

class CB4FieldParser(FieldParser):
    """Handle CodeBase4 'W' (wide/UTF-16LE) fields."""
    def parseW(self, field, data):
        if not data or data == b'\x00' * len(data):
            return None
        if isinstance(data, bytes):
            # Inline W fields: raw UTF-16LE bytes
            try:
                return data.decode('utf-16-le').rstrip('\x00') or None
            except Exception:
                pass
            # Fallback: hex-encoded UTF-16LE (seen in some memo contexts)
            try:
                h = data.decode('ascii', errors='ignore').strip().rstrip('\x00')
                if h:
                    return bytes.fromhex(h).decode('utf-16-le').rstrip('\x00') or None
            except Exception:
                pass
        return None


def open_table(directory, name):
    return DBF(os.path.join(directory, name),
               ignore_missing_memofile=True,
               parserclass=CB4FieldParser,
               encoding='latin-1')


def parse_rmr_blob(blob):
    """Parse a Tipo=26 (RMR) binary blob from T_TEST memo field.

    Header layout (504 bytes):
      0x34  i32   test subtype flag (1 = standard RMR)
      0x68  i64   total breath count
      0x72  u16   barometric pressure (Pb) in mmHg
      0x80  f64   assumed RQ (always 0.85)
      0x88  f64   measured RQ (STPD-corrected average)
      0x90  i32   breaths to skip before averaging
      0xB0  i32   RMR in kcal/day

    Per-breath records start at offset 0x200, each 104 bytes (13 × f64):
      d[0]  Rf    respiratory frequency (breaths/min)
      d[1]  HR    heart rate (0 if not connected)
      d[2]  VE    minute ventilation (L/min)
      d[3]  FeO2  expired O2 fraction (%)
      d[4]  (reserved)
      d[5]  Vt    tidal volume (L)
      d[6..12]    (reserved, always 0 on Fitmate; may be used by Fitmate Pro)
    """
    if isinstance(blob, str):
        raw = blob.encode('latin-1')
    elif isinstance(blob, bytes):
        raw = blob
    else:
        return None

    if len(raw) < 0x200:
        return None

    n_breaths = struct.unpack_from('<q', raw, 0x68)[0]
    pb        = struct.unpack_from('<H', raw, 0x72)[0]
    rq_assumed= struct.unpack_from('<d', raw, 0x80)[0]
    rq        = struct.unpack_from('<d', raw, 0x88)[0]
    skip      = struct.unpack_from('<i', raw, 0x90)[0]
    rmr       = struct.unpack_from('<i', raw, 0xB0)[0]

    breaths = []
    for b in range(n_breaths):
        off = 0x200 + b * 13 * 8
        if off + 13 * 8 > len(raw):
            break
        rf   = struct.unpack_from('<d', raw, off + 0*8)[0]
        hr   = struct.unpack_from('<d', raw, off + 1*8)[0]
        ve   = struct.unpack_from('<d', raw, off + 2*8)[0]
        feo2 = struct.unpack_from('<d', raw, off + 3*8)[0]
        vt   = struct.unpack_from('<d', raw, off + 5*8)[0]
        breaths.append({'rf': rf, 'hr': hr, 've': ve, 'feo2': feo2, 'vt': vt})

    # Compute averages over the averaging window (skip initial breaths)
    avg_window = breaths[skip:] if skip < len(breaths) else breaths
    avg = {}
    if avg_window:
        avg['rf']   = statistics.mean(b['rf']   for b in avg_window)
        avg['ve']   = statistics.mean(b['ve']   for b in avg_window)
        avg['feo2'] = statistics.mean(b['feo2'] for b in avg_window)
        avg['vt']   = statistics.mean(b['vt']   for b in avg_window)
        avg['hr']   = statistics.mean(b['hr']   for b in avg_window)

    total_minutes = n_breaths / avg['rf'] if avg.get('rf', 0) > 0 else 0

    return {
        'n_breaths': n_breaths,
        'pb_mmhg': pb,
        'rq': rq,
        'rmr_kcal': rmr,
        'skip': skip,
        'avg': avg,
        'breaths': breaths,
        'total_min': total_minutes,
    }


def main():
    dbf_dir = os.path.dirname(os.path.abspath(__file__))

    # Load patients
    patients = {}
    for rec in open_table(dbf_dir, 'ANAGRAFE.DBF'):
        d = dict(rec)
        pid = d['PROGRESS']
        patients[pid] = {
            'id': pid,
            'last': d.get('A_LASTNAME') or '',
            'first': d.get('A_FRSTNAME') or '',
            'member': d.get('A_MEMBERID') or '',
            'dob': d.get('A_BIRTHDT'),
            'sex': 'F' if d.get('A_SEX') else 'M',
        }

    # Load sessions
    sessions = {}
    for rec in open_table(dbf_dir, 'SESSIONI.DBF'):
        d = dict(rec)
        sid = d['PROGRESS']
        sessions[sid] = {
            'id': sid,
            'anag': d['S_ANAG'],
            'opened': d.get('S_OPENEDON'),
            'closed': d.get('S_CLOSEDON'),
            'height': d.get('S_HEIGHT') or 0,
            'weight': d.get('S_WEIGHT') or 0,
        }

    # No argument: list patients
    if len(sys.argv) < 2:
        print("Patients:")
        for p in sorted(patients.values(), key=lambda x: x['id']):
            n_sess = sum(1 for s in sessions.values() if s['anag'] == p['id'])
            print(f"  [{p['id']}] {p['last']}, {p['first']}  "
                  f"DOB={p['dob']}  Sex={p['sex']}  Sessions={n_sess}")
        print(f"\nUsage: {sys.argv[0]} <patient# or lastname>")
        return

    # Find patient
    query = sys.argv[1]
    target = None
    if query.isdigit():
        target = int(query)
    else:
        for p in patients.values():
            if p['last'].upper() == query.upper():
                target = p['id']
                break
    if target is None or target not in patients:
        sys.exit(f"Patient not found: {query}")

    pat = patients[target]
    print(f"Patient: {pat['last']}, {pat['first']}  DOB={pat['dob']}  Sex={pat['sex']}")
    print()

    # Load RMR tests for this patient
    tests = []
    for rec in open_table(dbf_dir, 'TEST.DBF'):
        d = dict(rec)
        if d['T_ANAG'] != target or d['T_TIPO'] != TIPO_RMR:
            continue
        blob = d.get('T_TEST')
        sess = sessions.get(d['T_SESSION'], {})
        tests.append({
            'date': d['T_DATE'],
            'session': d['T_SESSION'],
            'height': sess.get('height', 0),
            'weight': sess.get('weight', 0),
            'blob': blob,
        })

    if not tests:
        print("No RMR tests found.")
        return

    print(f"{'Date':>12} {'RMR':>6} {'RQ':>6} {'VE':>5} {'FeO2':>6} "
          f"{'Rf':>5} {'Vt':>5} {'Pb':>5} {'Wt':>6} {'Dur':>5} {'#br':>4}")
    print('-' * 80)

    for t in sorted(tests, key=lambda x: x['date']):
        parsed = parse_rmr_blob(t['blob']) if t['blob'] else None
        if parsed is None:
            print(f"{t['date']}  (no data)")
            continue

        a = parsed['avg']
        dur = f"{parsed['total_min']:.1f}m"
        print(f"{t['date']}  "
              f"{parsed['rmr_kcal']:5d} "
              f"{parsed['rq']:5.3f}  "
              f"{a.get('ve',0):4.1f}  "
              f"{a.get('feo2',0):5.2f} "
              f"{a.get('rf',0):5.1f} "
              f"{a.get('vt',0):4.2f} "
              f"{parsed['pb_mmhg']:5d} "
              f"{t['weight']:5.1f} "
              f"{dur:>5} "
              f"{parsed['n_breaths']:4d}")

    # Detailed per-breath dump if -v
    if len(sys.argv) > 2 and sys.argv[2] == '-v':
        print("\n\nPer-breath detail (last test):")
        t = sorted(tests, key=lambda x: x['date'])[-1]
        parsed = parse_rmr_blob(t['blob'])
        if parsed:
            print(f"{'#':>4} {'Rf':>6} {'VE':>6} {'FeO2':>6} {'Vt':>6} {'HR':>5}")
            for i, b in enumerate(parsed['breaths']):
                marker = ' *' if i == parsed['skip'] else ''
                print(f"{i:4d} {b['rf']:6.2f} {b['ve']:6.2f} "
                      f"{b['feo2']:6.2f} {b['vt']:6.3f} {b['hr']:5.0f}{marker}")
            print(f"\n  * = averaging begins at breath {parsed['skip']}")


if __name__ == '__main__':
    main()
