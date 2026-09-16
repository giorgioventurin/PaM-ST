"""Deterministic fixtures for the validation scripts."""
import csv, math, os, random

DEV = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(DEV, "data")
REAL = "../../data/PDAC_202602/annotated_cells/LSP31891_P181_rep1_annotated.csv"
W, N, K = 3000.0, 4500, 6
BASE = [0.35, 0.25, 0.15, 0.12, 0.08, 0.05]


def write(name, rows):
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, name), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["X_centroid", "Y_centroid", "Cell_Type"])
        writer.writerows(rows)
    return len(rows)


def pick(rng, weights):
    u, cumulative = rng.random(), 0.0
    for i, weight in enumerate(weights):
        cumulative += weight
        if u < cumulative:
            return i
    return len(weights) - 1


def niches(seed=11, count=40, radius=150.0):
    """Background composition plus two-type niches: fine-scale structure."""
    rng = random.Random(seed)
    centres = [(rng.uniform(0, W), rng.uniform(0, W)) for _ in range(count)]
    rows = []
    for _ in range(N):
        x, y = rng.uniform(0, W), rng.uniform(0, W)
        inside = any((x - cx) ** 2 + (y - cy) ** 2 < radius ** 2 for cx, cy in centres)
        label = rng.choice([3, 5]) if inside else pick(rng, BASE)
        rows.append([x, y, f"T{label}"])
    return rows


def compartments(seed=12):
    """Two sharp compartments: coarse structure only."""
    rng = random.Random(seed)
    tumour = [0.55, 0.05, 0.05, 0.15, 0.15, 0.05]
    stroma = [0.05, 0.50, 0.25, 0.05, 0.05, 0.10]
    rows = []
    for _ in range(N):
        x, y = rng.uniform(0, W), rng.uniform(0, W)
        side = math.sin(x / 700.0) + math.cos(y / 900.0) > 0
        rows.append([x, y, f"T{pick(rng, tumour if side else stroma)}"])
    return rows


def real_crop(half_width=4000.0):
    """An 8000 x 8000 window of LSP31891, centred on the median cell."""
    source = os.path.join(DEV, REAL)
    cells = [c for c in csv.DictReader(open(source)) if c["Cell_Type"].lower() != "unclassified"]
    xs = sorted(float(c["X_centroid"]) for c in cells)
    ys = sorted(float(c["Y_centroid"]) for c in cells)
    cx, cy = xs[len(xs) // 2], ys[len(ys) // 2]
    return [[c["X_centroid"], c["Y_centroid"], c["Cell_Type"]] for c in cells
            if abs(float(c["X_centroid"]) - cx) < half_width
            and abs(float(c["Y_centroid"]) - cy) < half_width]


if __name__ == "__main__":
    print("niches.csv", write("niches.csv", niches()))
    print("compartments.csv", write("compartments.csv", compartments()))
    try:
        print("real_crop.csv", write("real_crop.csv", real_crop()))
    except FileNotFoundError:
        print("real_crop.csv skipped: PDAC data not found at", REAL)
