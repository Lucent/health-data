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

import csv, struct, sys, os, statistics
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

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


def _stpd_factor(pb_mmhg):
    """BTPS-to-STPD correction factor (body temp 37°C, PH2O=47 mmHg)."""
    return (pb_mmhg - 47) / 760.0 * 273.15 / 310.15


def _vo2_stpd(ve, feo2, rq, pb_mmhg):
    """Per-breath VO2 in mL/min STPD using Haldane transformation.

    ve:   minute ventilation (L/min, BTPS)
    feo2: expired O2 fraction (%)
    rq:   respiratory quotient (assumed)
    pb_mmhg: barometric pressure
    """
    feco2 = (FIO2 - feo2) * rq
    fen2 = 100.0 - feo2 - feco2
    fin2 = 100.0 - FIO2 - 0.04  # 0.04% atmospheric CO2
    vi = ve * fen2 / fin2  # Haldane: inspired volume
    vo2_btps = vi * FIO2 / 100.0 - ve * feo2 / 100.0  # L/min BTPS
    return vo2_btps * _stpd_factor(pb_mmhg) * 1000.0  # mL/min STPD


def _minute_bin_cv(times_min, values):
    """CV% computed from 1-minute bin averages (Fitmate method)."""
    if len(values) < 2:
        return 0.0
    t0 = times_min[0]
    total = int(times_min[-1] - t0) + 1
    bins = []
    for m in range(total):
        mask = (times_min >= t0 + m) & (times_min < t0 + m + 1)
        if np.any(mask):
            bins.append(np.mean(values[mask]))
    if len(bins) < 2:
        return 0.0
    bins = np.array(bins)
    return float(np.std(bins, ddof=1) / np.mean(bins) * 100.0)


def compute_rmr_stats(parsed):
    """Compute summary stats from a parsed RMR blob for rmr.csv columns."""
    breaths = parsed.get('breaths', [])
    n = len(breaths)
    skip = parsed.get('skip', 0)
    rq = parsed.get('rq', 0.85)
    pb = parsed.get('pb_mmhg', 760)
    if n == 0 or skip >= n:
        return {}

    window = breaths[skip:]
    ve   = np.array([b['ve']   for b in window])
    rf   = np.array([b['rf']   for b in window])
    feo2 = np.array([b['feo2'] for b in window])
    vo2  = np.array([_vo2_stpd(b['ve'], b['feo2'], rq, pb) for b in window])

    # Time axis for the averaging window
    times = []
    t = 0.0
    for b in window:
        times.append(t)
        if b['rf'] > 0:
            t += 60.0 / b['rf']
    times = np.array(times) / 60.0

    return {
        'vo2_mL_min': round(float(np.mean(vo2)), 1),
        've_L_min': round(float(np.mean(ve)), 2),
        'rf_br_min': round(float(np.mean(rf)), 1),
        'feo2_pct': round(float(np.mean(feo2)), 2),
        'cv_ve_pct': round(_minute_bin_cv(times, ve), 1),
        'cv_vo2_pct': round(_minute_bin_cv(times, vo2), 1),
        'n_breaths': n,
        'skip': skip,
        'duration_min': round(parsed.get('total_min', 0), 1),
    }


def plot_rmr_test(parsed, date, out_path):
    """Plot breath-by-breath RMR test data to PNG."""
    breaths = parsed['breaths']
    n = len(breaths)
    if n == 0:
        return

    # Cumulative time axis from per-breath Rf
    times = []
    t = 0.0
    for b in breaths:
        times.append(t)
        if b['rf'] > 0:
            t += 60.0 / b['rf']  # seconds per breath
    times = np.array(times) / 60.0  # convert to minutes

    rf   = np.array([b['rf']   for b in breaths])
    ve   = np.array([b['ve']   for b in breaths])
    feo2 = np.array([b['feo2'] for b in breaths])
    vt   = np.array([b['vt']   for b in breaths])
    hr   = np.array([b['hr']   for b in breaths])

    # Per-breath VO2 (mL/min STPD) and RMR (kcal/day) via Haldane + STPD
    rq = parsed['rq']
    pb = parsed['pb_mmhg']
    vo2 = np.array([_vo2_stpd(v, f, rq, pb) for v, f in zip(ve, feo2)])
    rmr_breath = vo2 / 1000.0 * (3.941 + 1.106 * rq) * 1440.0

    skip = parsed['skip']
    rmr = parsed['rmr_kcal']

    has_hr = np.any(hr > 0)
    n_panels = 5 if has_hr else 4
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 2.2 * n_panels + 0.8),
                             sharex=True)

    panels = [
        (rmr_breath, 'RMR', 'kcal/day', 'tab:red'),
        (ve,         'VE',  'L/min',    'tab:blue'),
        (feo2,       'FeO₂', '%',       'tab:green'),
        (rf,         'Rf',  'br/min',   'tab:purple'),
    ]
    if has_hr:
        panels.append((hr, 'HR', 'bpm', 'tab:orange'))

    for ax, (y, label, unit, color) in zip(axes, panels):
        ax.plot(times, y, color=color, linewidth=0.8, alpha=0.85)
        if skip > 0 and skip < n:
            ax.axvline(times[skip], color='gray', linestyle='--', linewidth=0.7,
                       alpha=0.6)
        # Averaging window mean
        avg_vals = y[skip:] if skip < n else y
        avg_mean = np.mean(avg_vals)
        t_start = times[skip] if skip < n else times[0]
        ax.hlines(avg_mean, t_start, times[-1], color=color,
                  linestyle=':', linewidth=1.0, alpha=0.7)
        ax.set_ylabel(f'{label} ({unit})', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.2)

    axes[-1].set_xlabel('Time (min)', fontsize=9)

    # Summary stats annotation
    stats = compute_rmr_stats(parsed)
    if stats:
        text = (
            f"VO₂  {stats['vo2_mL_min']:.0f} mL/min  (CV {stats['cv_vo2_pct']:.1f}%)\n"
            f"VE   {stats['ve_L_min']:.1f} L/min  (CV {stats['cv_ve_pct']:.1f}%)\n"
            f"FeO₂ {stats['feo2_pct']:.2f}%\n"
            f"Rf   {stats['rf_br_min']:.1f} br/min\n"
            f"{stats['n_breaths']} breaths, skip {stats['skip']}, {stats['duration_min']:.1f} min"
        )
        fig.text(0.99, 0.99, text, transform=fig.transFigure,
                 fontsize=8, fontfamily='monospace',
                 verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                           edgecolor='gray', alpha=0.9))

    fig.suptitle(f'RMR Test — {date}    RMR={rmr} kcal/day    RQ={rq:.3f}',
                 fontsize=11, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 0.78, 0.96])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  → {out_path}")


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

    # Load sessions (optional — file may not be present)
    sessions = {}
    sess_path = os.path.join(dbf_dir, 'SESSIONI.DBF')
    if os.path.exists(sess_path):
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

        png_path = os.path.join(dbf_dir, f"{t['date']}.png")
        plot_rmr_test(parsed, t['date'], png_path)

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
