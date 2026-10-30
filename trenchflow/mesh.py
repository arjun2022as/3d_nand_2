#!/usr/bin/env python3
"""Write a 2D Fluent ASCII mesh (.msh) of a rinse channel over one rectangular trench.

Copied from 3d_nand/ansys_busnaina/make_trench_mesh.py (commit 141e41b); the only change is
three options, all off by default so the baseline mesh is byte-identical:
  --top-type      symmetry (baseline) | velocity-inlet (transverse inflow, Cases B/C/R)
  --left-type     velocity-inlet named "inlet" (baseline) | pressure-outlet named "left" (open side, Case R)
  --split-cavity  split the trench into "cavity" (upper half) and "cavity-lower" (lower half)
  --grade-ratio   geometric cell growth away from the trench (taller/longer domains, Case S)

Geometry follows Lin, Busnaina & Suni (2002), Fig. 1. y = 0 is the wafer surface.

     y=H  +------------------ top (symmetry = free surface) -----------------+
          |                                                                  |
   inlet  |                    fluid zone "channel"                          |  outlet
          |                                                                  |
     y=0  +------ wafer (wall) ------+-- mouth --+------ wafer (wall) -------+
                                     |           |
                                     |  "cavity" |  trench-walls (wall)
                                     |           |
    y=-D                             +-----------+
         x=-Lu                      x=0         x=W                       x=W+Ld

Uniform square cells of size h = W / n_per_w. Pure Python, no dependencies.

Fluent 2D face convention (ASCII format, section 13): each face is "n0 n1 c0 c1"
with c0 on the LEFT of the direction n0 -> n1, and c1 = 0 on boundaries. This is
the opposite of the 3D right-hand rule; the mesh is checked against it before writing.
"""
import argparse
import json
import sys

# Fluent bc-type codes (section 13) and zone type names (section 45)
BC = {"interior": 2, "wall": 3, "pressure-outlet": 5, "symmetry": 7, "velocity-inlet": 10}

FACE_ZONES = [  # (id, name, type) - written in this order, ids contiguous per zone
    (4, "interior-channel", "interior"),
    (5, "interior-cavity", "interior"),
    (6, "mouth", "interior"),
    (7, "inlet", "velocity-inlet"),
    (8, "outlet", "pressure-outlet"),
    (9, "top", "symmetry"),
    (10, "wafer", "wall"),
    (11, "trench-walls", "wall"),
    (13, "cavity-split", "interior"),      # only with --split-cavity (empty zones are not written)
]
CELL_ZONES = [(2, "channel"), (3, "cavity"), (12, "cavity-lower")]
NODE_ZONE_ID = 1


def graded(L, h, core, ratio, hmax):
    """Cell sizes from the trench outward over length L: uniform h up to `core`, then growing by
    `ratio` per cell up to hmax. The outermost cells are rescaled so the sizes sum to L exactly."""
    sizes, total = [], 0.0
    while total < L - 1e-12 * L:
        s = h if not sizes or total < core - 1e-12 * L else min(hmax, sizes[-1] * ratio)
        sizes.append(s)
        total += s
    over = total - L
    if over > 0:                        # shrink the outer graded cells, or the last one if all uniform
        k = next((i for i, s in enumerate(sizes) if s > h * (1 + 1e-9)), len(sizes) - 1)
        tail = sum(sizes[k:])
        sizes[k:] = [s * (tail - over) / tail for s in sizes[k:]]
    return sizes


def build(W, D, H, Lu, Ld, n_per_w, split_cavity=False, grade=None):
    """grade = None: uniform square cells (baseline). grade = dict(ratio, hmax, core_x, core_y) in units
    of W: uniform cells within core_x of the trench and below core_y, geometric growth outside."""
    h = W / n_per_w
    nw, ndp = n_per_w, max(1, round(D / h))
    if grade:
        g = grade
        up = graded(Lu, h, g["core_x"] * W, g["ratio"], g["hmax"] * W)
        dn = graded(Ld, h, g["core_x"] * W, g["ratio"], g["hmax"] * W)
        hy = graded(H, h, g["core_y"] * W, g["ratio"], g["hmax"] * W)
        xs_up = [-sum(up[:len(up) - k]) for k in range(len(up))]      # from x = -Lu toward the trench
        xs_dn = [W + sum(dn[:k]) for k in range(len(dn) + 1)]
        ys_ch = [sum(hy[:k]) for k in range(len(hy) + 1)]
        nu, nd, nh = len(up), len(dn), len(hy)
    else:
        nu, nd, nh = max(1, round(Lu / h)), max(1, round(Ld / h)), max(1, round(H / h))
        xs_up = [-Lu + Lu * i / nu for i in range(nu)]
        xs_dn = [W + Ld * i / nd for i in range(nd + 1)]
        ys_ch = [H * j / nh for j in range(nh + 1)]
    NX = nu + nw + nd

    xs = xs_up + [W * i / nw for i in range(nw)] + xs_dn
    ys = [-D + D * j / ndp for j in range(ndp)] + ys_ch

    # cells as lattice (i, j) of their lower-left corner; channel first, then cavity
    cells = [(i, j, 2) for j in range(ndp, ndp + nh) for i in range(NX)]
    j_split = ndp // 2 if split_cavity else 0          # rows below j_split form "cavity-lower"
    cells += [(i, j, 3) for j in range(j_split, ndp) for i in range(nu, nu + nw)]
    cells += [(i, j, 12) for j in range(j_split) for i in range(nu, nu + nw)]

    used = set()
    for i, j, _ in cells:
        used.update({(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)})
    node_id, nodes = {}, []
    for ij in sorted(used, key=lambda t: (t[1], t[0])):
        node_id[ij] = len(nodes) + 1
        nodes.append((xs[ij[0]], ys[ij[1]]))

    # faces: walk each cell counter-clockwise, so the cell is on the left of p -> q (= c0)
    faces = {}
    for c, (i, j, _) in enumerate(cells, start=1):
        corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        for k in range(4):
            p, q = node_id[corners[k]], node_id[corners[(k + 1) % 4]]
            key = (min(p, q), max(p, q))
            if key in faces:
                if faces[key][3]:
                    raise RuntimeError(f"face {key} shared by more than two cells")
                faces[key][3] = c
            else:
                faces[key] = [p, q, c, 0]

    zone_of_cell = [None] + [z for _, _, z in cells]
    tol = h * 1e-6
    by_zone = {zid: [] for zid, _, _ in FACE_ZONES}
    for f in faces.values():
        p, q, c0, c1 = f
        if c1:
            z0, z1 = zone_of_cell[c0], zone_of_cell[c1]
            pair = {z0, z1}
            by_zone[4 if pair == {2} else 5 if pair in ({3}, {12}) else 13 if pair == {3, 12} else 6].append(f)
            continue
        xm = 0.5 * (nodes[p - 1][0] + nodes[q - 1][0])
        ym = 0.5 * (nodes[p - 1][1] + nodes[q - 1][1])
        if abs(xm + Lu) < tol:
            by_zone[7].append(f)
        elif abs(xm - (W + Ld)) < tol:
            by_zone[8].append(f)
        elif abs(ym - H) < tol:
            by_zone[9].append(f)
        elif abs(ym) < tol:
            by_zone[10].append(f)
        else:
            by_zone[11].append(f)

    info = dict(h=h, n_cells=len(cells), n_nodes=len(nodes), n_faces=len(faces),
                n_channel=sum(1 for c in cells if c[2] == 2),
                n_cavity=sum(1 for c in cells if c[2] in (3, 12)),
                n_cavity_lower=sum(1 for c in cells if c[2] == 12),
                cells_x=[nu, nw, nd], cells_y=[ndp, nh])
    return nodes, cells, xs, ys, by_zone, info


def check(nodes, cells, xs, ys, by_zone, W, D, H, Lu, Ld):
    """Verify orientation, connectivity, boundary lengths and area before writing."""
    def centroid(c):
        i, j, _ = cells[c - 1]
        return 0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1])

    refs = [0] * (len(cells) + 1)
    for flist in by_zone.values():
        for p, q, c0, c1 in flist:
            (px, py), (qx, qy) = nodes[p - 1], nodes[q - 1]
            for c, sign in ((c0, 1), (c1, -1)):
                if not c:
                    continue
                cx, cy = centroid(c)
                cross = (qx - px) * (cy - py) - (qy - py) * (cx - px)
                if sign * cross <= 0:
                    raise RuntimeError(f"bad orientation: face {p}-{q}, cell {c}")
                refs[c] += 1
    if any(r != 4 for r in refs[1:]):
        raise RuntimeError("a cell does not have exactly 4 faces")

    def length(zid):
        return sum(abs(nodes[q - 1][0] - nodes[p - 1][0]) + abs(nodes[q - 1][1] - nodes[p - 1][1])
                   for p, q, _, _ in by_zone[zid])
    expect = {6: W, 7: H, 8: H, 9: Lu + W + Ld, 10: Lu + Ld, 11: 2 * D + W}
    if by_zone[13]:
        expect[13] = W
    for zid, L in expect.items():
        if abs(length(zid) - L) > 1e-9 * L:
            raise RuntimeError(f"zone {zid} length {length(zid)} != {L}")

    area = sum((xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j]) for i, j, _ in cells)
    expect_area = (Lu + W + Ld) * H + W * D
    if abs(area - expect_area) > 1e-9 * expect_area:
        raise RuntimeError(f"area {area} != {expect_area}")


def write_msh(path, nodes, cells, by_zone, title, top_type="symmetry", left_type="velocity-inlet"):
    nn, nc = len(nodes), len(cells)
    nf = sum(len(v) for v in by_zone.values())
    rename = {"top": ("top", top_type),
              "inlet": ("inlet", "velocity-inlet") if left_type == "velocity-inlet" else ("left", left_type)}
    faces = [(zid, *rename.get(name, (name, ztype))) for zid, name, ztype in FACE_ZONES if by_zone[zid]]
    cell_zones = [(zid, name, sum(1 for c in cells if c[2] == zid)) for zid, name in CELL_ZONES]
    cell_zones = [z for z in cell_zones if z[2]]
    with open(path, "w", newline="\n") as f:
        f.write(f'(0 "{title}")\n(2 2)\n')
        f.write(f"(10 (0 1 {nn:x} 0 2))\n(12 (0 1 {nc:x} 0))\n(13 (0 1 {nf:x} 0))\n")

        f.write(f"(10 ({NODE_ZONE_ID:x} 1 {nn:x} 1 2)(\n")
        f.writelines(f"{x:.15e} {y:.15e}\n" for x, y in nodes)
        f.write("))\n")

        first = 1
        for zid, _, count in cell_zones:          # cells are stored zone by zone (see build)
            f.write(f"(12 ({zid:x} {first:x} {first + count - 1:x} 1 3))\n")
            first += count

        first = 1
        for zid, name, ztype in faces:
            flist = by_zone[zid]
            last = first + len(flist) - 1
            f.write(f"(13 ({zid:x} {first:x} {last:x} {BC[ztype]:x} 2)(\n")
            f.writelines(f"{p:x} {q:x} {c0:x} {c1:x}\n" for p, q, c0, c1 in flist)
            f.write("))\n")
            first = last + 1

        for zid, name, _ in cell_zones:
            f.write(f"(45 ({zid} fluid {name})())\n")
        for zid, name, ztype in faces:
            f.write(f"(45 ({zid} {ztype} {name})())\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--W", type=float, required=True, help="trench width [m]")
    ap.add_argument("--AR", type=float, default=1.0, help="aspect ratio D/W")
    ap.add_argument("--H", type=float, default=2.5, help="channel height / W")
    ap.add_argument("--Lu", type=float, default=2.0, help="upstream length / W")
    ap.add_argument("--Ld", type=float, default=3.0, help="downstream length / W")
    ap.add_argument("--n-per-w", type=int, default=40, help="cells across the trench width")
    ap.add_argument("--top-type", choices=["symmetry", "velocity-inlet"], default="symmetry",
                    help="top boundary: symmetry (baseline) or velocity-inlet (transverse inflow)")
    ap.add_argument("--left-type", choices=["velocity-inlet", "pressure-outlet"], default="velocity-inlet",
                    help="left boundary: velocity inlet 'inlet' (baseline) or open pressure outlet 'left'")
    ap.add_argument("--split-cavity", action="store_true", help="separate upper/lower trench cell zones")
    ap.add_argument("--grade-ratio", type=float, default=0.0,
                    help="> 1: geometric cell growth away from the trench (0 = uniform, baseline)")
    ap.add_argument("--grade-hmax", type=float, default=0.5, help="largest graded cell / W")
    ap.add_argument("--core-x", type=float, default=1.0, help="uniform region either side of the trench / W")
    ap.add_argument("--core-y", type=float, default=2.5, help="uniform region above the wafer / W")
    ap.add_argument("-o", "--out", required=True, help="output .msh path")
    a = ap.parse_args(argv)

    W = a.W
    D, H, Lu, Ld = a.AR * W, a.H * W, a.Lu * W, a.Ld * W
    grade = (dict(ratio=a.grade_ratio, hmax=a.grade_hmax, core_x=a.core_x, core_y=a.core_y)
             if a.grade_ratio > 1 else None)
    nodes, cells, xs, ys, by_zone, info = build(W, D, H, Lu, Ld, a.n_per_w, a.split_cavity, grade)
    check(nodes, cells, xs, ys, by_zone, W, D, H, Lu, Ld)
    title = f"Trench W={W:g} m D={D:g} m H={H:g} m, {a.n_per_w} cells/W"
    write_msh(a.out, nodes, cells, by_zone, title, a.top_type, a.left_type)
    info.update(W=W, D=D, H=H, Lu=Lu, Ld=Ld, msh=a.out, top_type=a.top_type,
                left_type=a.left_type, split_cavity=a.split_cavity, grade=grade)
    print(json.dumps(info))
    return info


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
