"""Refit the sky path of the galstreams li2021 AAU tracks to the Li+2021 spectroscopic members.

Keeps each track's phi1 sampling and its distance / pm / RV columns; replaces only (ra, dec) with
phi2 = poly(phi1) fit to the members of that piece. Frame: pole from the li2021 summary files,
phi1 = 0 at RA = 31 deg (ATLAS: phi1 < 0, Aliqa Uma: phi1 > 0).

Member table: Li et al. 2021, ApJ 911, 149, Table 1 (downloaded from MRT_URL unless a local path is given).

Usage: python scripts/refit_aau_li2021_tracks.py [--members PATH] [--outdir DIR]
       (default outdir: galstreams/tracks, i.e. overwrites the li2021 AAU track files in place)
"""
import os, argparse, urllib.request, tempfile, numpy as np
import astropy.units as u, astropy.coordinates as ac
from astropy.table import Table
import gala.coordinates as gc

MRT_URL = "https://content.cld.iop.org/journals/0004-637X/911/2/149/revision2/apjabeb18t1_mrt.txt"
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "galstreams", "tracks") + os.sep

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--members", help="local copy of the Li+2021 Table 1 MRT file (default: download)")
ap.add_argument("--outdir", default=T, help="where to write the track files (default: galstreams/tracks)")
args = ap.parse_args()
OUT = args.outdir
os.makedirs(OUT, exist_ok=True)
if args.members is None:
    args.members = os.path.join(tempfile.mkdtemp(), "li2021_table1_mrt.txt")
    urllib.request.urlretrieve(MRT_URL, args.members)
DEG = {"AAU-ATLAS": 3, "AAU-AliqaUma": 2}

pole = ac.SkyCoord(74.25906739590538*u.deg, 47.93312434806639*u.deg)
fr = gc.GreatCircleICRSFrame.from_pole_ra0(pole, 31.*u.deg)
c = ac.SkyCoord([20, 30]*u.deg, [-26.8, -33]*u.deg).transform_to(fr)
SGN = 1 if np.diff(c.phi1.wrap_at(180*u.deg).deg)[0] > 0 else -1
def to_phi(ra, dec):
    c = ac.SkyCoord(ra=np.asarray(ra)*u.deg, dec=np.asarray(dec)*u.deg).transform_to(fr)
    return SGN*c.phi1.wrap_at(180*u.deg).deg, c.phi2.deg
def to_radec(p1, p2):
    c = ac.SkyCoord(phi1=SGN*np.asarray(p1)*u.deg, phi2=np.asarray(p2)*u.deg, frame=fr).icrs
    return c.ra.deg, c.dec.deg

def rewrite_radec(src, dst, cols):
    """Copy an ECSV verbatim, replacing only the given columns' values in the data rows."""
    lines = open(src).read().splitlines(keepends=True)
    hdr = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]
    names = body[0].strip().split(",")
    idx = {k: names.index(k) for k in cols}
    out = hdr + [body[0]]
    for i, l in enumerate(body[1:]):
        f = l.rstrip("\n").split(",")
        for k, j in idx.items():
            f[j] = repr(float(np.atleast_1d(cols[k])[i]))
        out.append(",".join(f) + "\n")
    open(dst, "w").write("".join(out))

li = np.genfromtxt(args.members, skip_header=23, usecols=(1, 2))
m1, m2 = to_phi(li[:, 0], li[:, 1])

for name, piece in [("AAU-ATLAS", m1 < 0), ("AAU-AliqaUma", m1 >= 0)]:
    fn = f"track.st.{name}.li2021.ecsv"
    t = Table.read(T + fn)          # read originals before any in-place rewrite
    t1, t2 = to_phi(t["ra"], t["dec"])
    x, y = m1[piece], m2[piece]
    coef = np.polyfit(x, y, DEG[name])
    for it in range(3):                       # 3-sigma clip
        r = y - np.polyval(coef, x); k = np.abs(r) < 3*np.std(r)
        coef = np.polyfit(x[k], y[k], DEG[name])
    r = y - np.polyval(coef, x)
    old = y - np.interp(x, *zip(*sorted(zip(t1, t2))))
    print(f"{name}: N={len(x)} deg={DEG[name]} coef={np.array2string(coef, precision=5)}")
    print(f"   member dphi2 vs OLD track: median {np.median(old):+.3f}  rms {np.sqrt(np.mean(old**2)):.3f}")
    print(f"   member dphi2 vs NEW track: median {np.median(r):+.3f}  rms {np.sqrt(np.mean(r**2)):.3f}")
    print(f"   track phi1 range {t1.min():.2f} .. {t1.max():.2f};  members {x.min():.2f} .. {x.max():.2f}")
    ra, dec = to_radec(t1, np.polyval(coef, t1))
    rewrite_radec(T + fn, os.path.join(OUT, fn), {"ra": ra, "dec": dec})

    # summary: move end/mid points onto the new path (only if they lie within this piece)
    sfn = fn.replace(".ecsv", ".summary.ecsv"); s = Table.read(T + sfn); new = {}
    for key in ["end_o", "end_f", "mid"]:
        p1k, _ = to_phi(s[f"{key}.ra"], s[f"{key}.dec"])
        if t1.min() - 0.01 <= p1k[0] <= t1.max() + 0.01:
            rk, dk = to_radec(p1k, np.polyval(coef, p1k))
            new[f"{key}.ra"], new[f"{key}.dec"] = rk, dk
        else:
            print(f"   summary {key} at phi1={p1k[0]:.2f} is outside this piece; left unchanged")
    rewrite_radec(T + sfn, os.path.join(OUT, sfn), new)
print("wrote", OUT)
