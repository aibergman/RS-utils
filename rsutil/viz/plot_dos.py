#!/usr/bin/env python3
"""Plot total and species/site projected DOS from RS-LMTO-ASA output.

Spin-polarized (default): magneticdos.out + <species>_dos.out, up positive, down negative.
Non-spin-polarized (--no-spin): totaldos.out + <species>_dos.out (total column).

Examples:
    plot_dos.py                          # run in a calculation folder
    plot_dos.py path/to/calc --no-spin
    plot_dos.py -s Ni-S Ni-S-1 --fermi 0.7 --xlim -0.5 0.2 -o dos.png
    plot_dos.py --layout overlay --smooth 5
"""
import argparse
import glob
import math
import os
import sys

import numpy as np


def load(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def smooth(y, n):
    if n <= 1:
        return y
    k = np.ones(n) / n
    return np.convolve(y, k, mode="same")


def read_fermi_energy(folder):
    """Read the Fermi energy in Ry from ``folder/input.nml`` if available."""
    input_path = os.path.join(folder, "input.nml")
    if not os.path.isfile(input_path):
        return None

    try:
        import f90nml

        document = f90nml.read(input_path)
    except (OSError, ValueError) as exc:
        raise ValueError(f"could not read Fermi energy from {input_path}: {exc}") from exc

    energy = document.get("energy", {})
    fermi = energy.get("fermi") if hasattr(energy, "get") else None
    if fermi is None:
        return None
    try:
        fermi = float(fermi)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{input_path} has a non-numeric &energy fermi value") from exc
    if not math.isfinite(fermi):
        raise ValueError(f"{input_path} has a non-finite &energy fermi value")
    return fermi


def find_species_files(folder, names=None):
    files = sorted(
        f for f in glob.glob(os.path.join(folder, "*_dos.out"))
        if not f.endswith("orbital_dos.out")
    )
    out = []
    for f in files:
        name = os.path.basename(f)[: -len("_dos.out")]
        if names is None or name in names:
            out.append((name, f))
    return out


def enable_overlay_legend_toggle(fig, legend, groups):
    """Make each overlay legend entry toggle all artists for that DOS file."""
    legend_artists = {}
    legend_lines = legend.get_lines()
    legend_texts = legend.get_texts()

    for legend_line, legend_text in zip(legend_lines, legend_texts):
        group = groups.get(legend_text.get_text())
        if group is None:
            continue
        legend_line.set_picker(5)
        legend_text.set_picker(True)
        legend_artists[legend_line] = (group, legend_line, legend_text)
        legend_artists[legend_text] = (group, legend_line, legend_text)

    def on_pick(event):
        picked = legend_artists.get(event.artist)
        if picked is None:
            return
        group, legend_line, legend_text = picked
        visible = not any(artist.get_visible() for artist in group)
        for artist in group:
            artist.set_visible(visible)
        alpha = 1.0 if visible else 0.2
        legend_line.set_alpha(alpha)
        legend_text.set_alpha(alpha)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("pick_event", on_pick)


def get_curves(name, path, spin):
    """Return (energy, {label: y}) for one file."""
    d = load(path)
    e = d[:, 0]
    if name == "System":
        if spin:  # magneticdos.out: energy dos_up dos_dw ...
            return e, {"up": d[:, 1], "dw": -d[:, 2]}
        return e, {"total": d[:, 1]}
    # species: energy site_dos_total site_dos_up site_dos_dw
    if spin:
        return e, {"up": d[:, 2], "dw": -d[:, 3]}
    return e, {"total": d[:, 1]}


def add_plot_dos_arguments(parser):
    """Add DOS file plotting options to an argparse parser."""
    parser.add_argument("folder", nargs="?", default=".", help="calculation folder (default: .)")
    parser.add_argument("--no-spin", action="store_true", help="plot non-spin-polarized DOS (totaldos.out)")
    parser.add_argument("-s", "--species", nargs="+", help="only these species (default: all)")
    parser.add_argument("--no-system", action="store_true", help="skip full-system DOS")
    parser.add_argument("--no-species", action="store_true", help="skip species DOS")
    parser.add_argument("--layout", choices=["stack", "overlay"], default="overlay",
                        help="stack: one panel per curve set; overlay: all in one axes (default: overlay)")
    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy in Ry (default: read &energy fermi from input.nml)")
    parser.add_argument("--xlim", nargs=2, type=float, metavar=("MIN", "MAX"))
    parser.add_argument("--ylim", nargs=2, type=float, metavar=("MIN", "MAX"))
    parser.add_argument("--smooth", type=int, default=1, help="moving-average window in points")
    fill = parser.add_mutually_exclusive_group()
    fill.add_argument("--fill", dest="fill", action="store_true", help="fill under curves (default)")
    fill.add_argument("--no-fill", dest="fill", action="store_false", help="do not fill under curves")
    parser.set_defaults(fill=True)
    parser.add_argument("--xlabel", default=None)
    parser.add_argument("-o", "--output", help="save figure instead of showing")
    parser.add_argument("--dpi", type=int, default=200)
    return parser


def plot_dos(args):
    """Plot DOS files using parsed arguments and return a process-style status."""
    if args.smooth < 1:
        raise ValueError("--smooth must be at least 1")

    if args.fermi is None:
        fermi = read_fermi_energy(args.folder)
        has_fermi_reference = fermi is not None
        if fermi is None:
            print(
                f"warning: no &energy fermi found in {os.path.join(args.folder, 'input.nml')}; "
                "using zero energy shift",
                file=sys.stderr,
            )
            fermi = 0.0
    else:
        fermi = args.fermi
        has_fermi_reference = True

    spin = not args.no_spin
    entries = []
    if not args.no_system:
        sysfile = os.path.join(args.folder, "magneticdos.out" if spin else "totaldos.out")
        if os.path.exists(sysfile):
            entries.append(("System", sysfile))
        else:
            print(f"warning: {sysfile} not found", file=sys.stderr)
    if not args.no_species:
        entries += find_species_files(args.folder, args.species)
    if not entries:
        raise ValueError("No DOS files found.")

    import matplotlib.pyplot as plt

    n = len(entries)
    if args.layout == "stack":
        fig, axes = plt.subplots(n, 1, sharex=True, figsize=(6, 2.4 * n + 0.6), squeeze=False)
        axes = axes[:, 0]
    else:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        axes = [ax] * n

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    artist_groups = {}
    for i, (name, path) in enumerate(entries):
        e, curves = get_curves(name, path, spin)
        e = e - fermi
        ax = axes[i]
        c = colors[i % len(colors)]
        group = artist_groups.setdefault(name, [])
        for label, y in curves.items():
            y = smooth(y, args.smooth)
            ls = "--" if label == "dw" and args.layout == "overlay" else "-"
            (line,) = ax.plot(e, y, color=c, ls=ls, lw=1.2, label=name if label in ("up", "total") else None)
            group.append(line)
            if args.fill:
                fill = ax.fill_between(e, 0, y, color=c, alpha=0.25)
                group.append(fill)
        if args.layout == "stack":
            ax.set_title(name, fontsize=9, loc="right")
            ax.set_ylabel("DOS (states/Ry)")
        ax.axvline(0, color="k", lw=0.6, ls=":")
        if spin:
            ax.axhline(0, color="k", lw=0.6)
        if args.xlim:
            ax.set_xlim(*args.xlim)
        if args.ylim:
            ax.set_ylim(*args.ylim)

    if args.layout == "overlay":
        axes[0].set_ylabel("DOS (states/Ry)")
        legend = axes[0].legend(frameon=False)
        if args.output is None:
            legend.set_title("Click to toggle")
            enable_overlay_legend_toggle(fig, legend, artist_groups)
    axes[-1].set_xlabel(args.xlabel or ("E - E$_F$ (Ry)" if has_fermi_reference else "Energy (Ry)"))
    fig.tight_layout()
    if args.output:
        fig.savefig(args.output, dpi=args.dpi)
        print(f"saved {args.output}")
    else:
        plt.show()
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_plot_dos_arguments(p)
    args = p.parse_args(argv)
    try:
        return plot_dos(args)
    except (ImportError, OSError, ValueError) as exc:
        p.error(str(exc))


if __name__ == "__main__":
    main()
