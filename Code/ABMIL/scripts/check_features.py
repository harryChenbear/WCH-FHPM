#!/usr/bin/env python
import argparse
import sys
from pathlib import Path
ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
from abmil.io import load_features


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pt_dir", required=True)
    parser.add_argument("--expected_dim", type=int, default=1024)
    parser.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()
    pt_dir = Path(args.pt_dir)
    files = sorted(pt_dir.glob("*.pt"))
    print("pt_dir:", pt_dir)
    print("n_pt:", len(files))
    for f in files[:args.limit]:
        x = load_features(f, expected_dim=args.expected_dim)
        print(f.name, tuple(x.shape), x.dtype)


if __name__ == "__main__":
    main()
