<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import hashlib
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import openslide
import pandas as pd
import timm
import torch
from PIL import Image, ImageDraw, ImageFilter, ImageFont
from torchvision import transforms


SLIDE_EXTS = {".svs", ".ndpi", ".tif", ".tiff", ".mrxs", ".scn"}

try:
    RESAMPLING_NEAREST = Image.Resampling.NEAREST
except AttributeError:
    RESAMPLING_NEAREST = Image.NEAREST


def safe_name(s, max_len=120):
    s = str(s)
    b = re.sub(r"[^A-Za-z0-9\u4e00-\u9fff._-]+", "_", s).strip("._-")
    h = hashlib.md5(s.encode("utf-8", errors="ignore")).hexdigest()[:8]
    return f"{b[:max_len]}__{h}"


def find_slides(root):
    return sorted(
        str(p)
        for p in Path(root).rglob("*")
        if p.is_file() and p.suffix.lower() in SLIDE_EXTS
    )


def positive_index(label_mapping, task):
    if not isinstance(label_mapping, dict) or not label_mapping:
        return 1

    def norm(x):
        return str(x).lower().replace("_", "").replace("-", "").replace(" ", "")

    lm = {norm(k): int(v) for k, v in label_mapping.items()}

    if task == "fh":
        if "fh" in lm:
            return lm["fh"]  # exact match; avoids matching nonFH

        for k, v in lm.items():
            if not k.startswith("non") and k in {
                "fhrcc",
                "fhdrcc",
                "fhdrrcc",
                "fhdrc",
            }:
                return v

        return 1

    if task == "tumor":
        for k in ["tumor", "tumour", "cancer", "neoplastic"]:
            if k in lm:
                return lm[k]

        for k, v in lm.items():
            if (
                "normal" not in k
                and not k.startswith("non")
                and "negative" not in k
                and ("tumor" in k or "tumour" in k)
            ):
                return v

        return 1

    return 1


def load_model(path, device):
    ckpt = torch.load(path, map_location="cpu", weights_only=False)

    model_name = ckpt.get("model_name", "vit_base_patch16_224")
    img_size = int(ckpt.get("img_size", 512))
    num_classes = int(ckpt.get("num_classes", 2))

    state = (
        ckpt.get("state_dict")
        or ckpt.get("model_state_dict")
        or ckpt.get("model")
        or ckpt
    )

    state = {
        k.replace("module.", "", 1) if k.startswith("module.") else k: v
        for k, v in state.items()
    }

    model = timm.create_model(
        model_name,
        pretrained=False,
        num_classes=num_classes,
        img_size=img_size,
    )
    model.load_state_dict(state, strict=True)
    model.to(device)
    model.eval()

    prep = ckpt.get("preprocess", {})
    mean = prep.get("mean", prep.get("normalize_mean", [0.485, 0.456, 0.406]))
    std = prep.get("std", prep.get("normalize_std", [0.229, 0.224, 0.225]))

    info = {
        "img_size": img_size,
        "mean": mean,
        "std": std,
        "label_mapping": ckpt.get("label_mapping", {}),
        "model_path": str(path),
    }

    return model, info


def build_tf(info):
    return transforms.Compose(
        [
            transforms.Resize((info["img_size"], info["img_size"])),
            transforms.ToTensor(),
            transforms.Normalize(info["mean"], info["std"]),
        ]
    )


def tissue_fraction(img, sat_thr=20, gray_high=245, gray_low=10):
    a = np.asarray(img.convert("RGB"))

    mx = a.max(2).astype(np.int16)
    mn = a.min(2).astype(np.int16)
    sat = mx - mn
    gray = a.mean(2)

    return float(((sat > sat_thr) & (gray < gray_high) & (gray > gray_low)).mean())


def batch_predict(model, imgs, tf, device, pos, batch_size):
    out = []

    if not imgs:
        return out

    with torch.no_grad():
        for i in range(0, len(imgs), batch_size):
            x = torch.stack([tf(im) for im in imgs[i : i + batch_size]], 0).to(
                device,
                non_blocking=True,
            )

            with torch.cuda.amp.autocast(enabled=(device.type == "cuda")):
                p = torch.softmax(model(x), 1)[:, pos].detach().cpu().numpy()

            out += [float(v) for v in p]

    return out


def get_font(size, bold=False):
    font_paths = [
        (
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
            if bold
            else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
        ),
        (
            "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"
            if bold
            else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"
        ),
    ]

    for p in font_paths:
        if os.path.exists(p):
            return ImageFont.truetype(p, size)

    return ImageFont.load_default()


def interp(c1, c2, t):
    return tuple(
        np.clip(np.array(c1) * (1 - t) + np.array(c2) * t, 0, 255)
        .astype(np.uint8)
        .tolist()
    )


def color_prob(p):
    p = float(np.clip(0.5 + 1.18 * (float(np.clip(p, 0, 1)) - 0.5), 0, 1))

    anchors = [
        (0, (32, 72, 192)),
        (0.25, (92, 140, 230)),
        (0.5, (247, 245, 242)),
        (0.75, (239, 126, 92)),
        (1, (197, 26, 44)),
    ]

    for (p1, c1), (p2, c2) in zip(anchors[:-1], anchors[1:]):
        if p1 <= p <= p2:
            return interp(c1, c2, (p - p1) / (p2 - p1))

    return anchors[-1][1]


def rr_mask(size, r):
    m = Image.new("L", size, 0)
    d = ImageDraw.Draw(m)
    d.rounded_rectangle((0, 0, size[0] - 1, size[1] - 1), radius=r, fill=255)
    return m


def center_text(draw, xy, text, font, fill):
    b = draw.textbbox((0, 0), text, font=font)
    draw.text(
        (xy[0] - (b[2] - b[0]) / 2, xy[1] - (b[3] - b[1]) / 2),
        text,
        font=font,
        fill=fill,
    )


def make_heatmap(
    patch_df,
    final_prob,
    slide_name,
    out_png,
    canvas_w=1600,
    canvas_h=1400,
):
    df = patch_df[(patch_df["is_tumor_associated"] == 1)].copy()
    df["FH_probability"] = pd.to_numeric(df["FH_probability"], errors="coerce")
    df = df.dropna(subset=["FH_probability"])

    if len(df) == 0:
        return None

    xs = np.sort(df.x.unique())
    ys = np.sort(df.y.unique())
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: i for i, y in enumerate(ys)}

    grid = np.full((len(ys), len(xs), 3), 255, dtype=np.uint8)

    for _, r in df.iterrows():
        grid[yi[r.y], xi[r.x]] = color_prob(r.FH_probability)

    raw = Image.fromarray(grid, "RGB")
    bg = Image.new("RGBA", (canvas_w, canvas_h), (248, 248, 250, 255))

    # Card shadow
    shadow = Image.new("RGBA", bg.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle(
        (36, 46, canvas_w - 36, canvas_h - 26),
        radius=34,
        fill=(0, 0, 0, 42),
    )
    shadow = shadow.filter(ImageFilter.GaussianBlur(16))
    bg.alpha_composite(shadow)

    # Card
    card = Image.new("RGBA", (canvas_w - 72, canvas_h - 72), (255, 255, 255, 255))
    bg.paste(card, (36, 36), rr_mask(card.size, 34))

    draw = ImageDraw.Draw(bg)
    title = get_font(34, True)
    probfont = get_font(30, True)
    tickfont = get_font(20, False)

    center_text(
        draw,
        (canvas_w // 2, 92),
        "Heatmap of the Prediction",
        title,
        (25, 30, 40),
    )

    hx, hy, hw, hh = 105, 165, 1130, 980

    panel = Image.new("RGBA", (hw, hh), (252, 252, 253, 255))
    bg.paste(panel, (hx, hy), rr_mask((hw, hh), 22))

    draw.rounded_rectangle(
        (hx, hy, hx + hw, hy + hh),
        radius=22,
        outline=(232, 234, 238),
        width=2,
    )

    aw, ah = hw - 40, hh - 40
    rw, rh = raw.size
    sc = min(aw / rw, ah / rh)

    nw = max(1, int(rw * sc))
    nh = max(1, int(rh * sc))

    fit = raw.resize((nw, nh), RESAMPLING_NEAREST)

    mx = hx + 20 + (aw - nw) // 2
    my = hy + 20 + (ah - nh) // 2

    bg.paste(fit, (mx, my))

    # Colorbar
    cbx, cby, cbw, cbh = 1290, 240, 62, 740

    cb = Image.new("RGB", (cbw, cbh), "white")
    cd = ImageDraw.Draw(cb)

    for j in range(cbh):
        cd.line((0, j, cbw, j), fill=color_prob(1 - j / max(cbh - 1, 1)))

    bg.paste(cb, (cbx, cby), rr_mask((cbw, cbh), 12))
    draw = ImageDraw.Draw(bg)

    for t in [0.2, 0.4, 0.6, 0.8]:
        ty = cby + int((1 - t) * cbh)
        draw.line(
            (cbx + cbw + 8, ty, cbx + cbw + 20, ty),
            fill=(90, 94, 102),
            width=2,
        )
        draw.text(
            (cbx + cbw + 28, ty - 11),
            f"{t:.1f}",
            font=tickfont,
            fill=(88, 92, 98),
        )

    prefix = "Predicted probability of FHdRCC: "
    risk = "high risk" if final_prob >= 0.5 else "low risk"
    val = f"{final_prob:.3f} ({risk})"

    bw = draw.textbbox((0, 0), prefix, font=probfont)[2]
    vw = draw.textbbox((0, 0), val, font=probfont)[2]

    x = (canvas_w - bw - vw) // 2
    y = 1248

    draw.text((x, y), prefix, font=probfont, fill=(18, 28, 44))
    draw.text(
        (x + bw, y),
        val,
        font=probfont,
        fill=(153, 0, 52) if final_prob >= 0.5 else (22, 84, 153),
    )

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    bg.convert("RGB").save(out_png, quality=95)

    return out_png


def process_slide(
    slide_path,
    args,
    tumor_model,
    tumor_info,
    fh_model,
    fh_info,
    device,
):
    name = os.path.basename(slide_path)
    uid = safe_name(Path(name).stem)
    outdir = os.path.join(args.out_dir, uid)
    os.makedirs(outdir, exist_ok=True)

    patch_csv = os.path.join(outdir, "patch_predictions.csv")
    result_csv = os.path.join(outdir, "slide_result.csv")
    heatmap = os.path.join(outdir, "FH_probability_heatmap.png")

    t0 = time.time()

    slide = openslide.OpenSlide(slide_path)
    W, H = slide.dimensions

    tumor_tf = build_tf(tumor_info)
    fh_tf = build_tf(fh_info)

    tumor_pos = positive_index(tumor_info.get("label_mapping", {}), "tumor")
    fh_pos = positive_index(fh_info.get("label_mapping", {}), "fh")

    print(
        f"\n===== Slide: {name} =====\n"
        f"Size: {W} x {H}\n"
        f"tumor_pos_index={tumor_pos}; fh_pos_index={fh_pos}",
        flush=True,
    )

    rows = []
    imgs = []
    metas = []

    n_grid = 0
    n_tissue = 0
    n_failed = 0
    n_tumor = 0

    last = time.time()

    def flush():
        nonlocal imgs, metas, rows, n_tumor

        if not imgs:
            return

        tps = batch_predict(
            tumor_model,
            imgs,
            tumor_tf,
            device,
            tumor_pos,
            args.batch_size,
        )

        tims = []
        trows = []

        for im, meta, tp in zip(imgs, metas, tps):
            is_t = int(tp >= args.tumor_threshold)

            row = dict(meta)
            row.update(
                {
                    "tumor_probability": float(tp),
                    "is_tumor_associated": is_t,
                    "FH_probability": np.nan,
                }
            )

            if is_t:
                tims.append(im)
                trows.append(row)
            else:
                rows.append(row)

        if tims:
            fps = batch_predict(
                fh_model,
                tims,
                fh_tf,
                device,
                fh_pos,
                args.batch_size,
            )

            for row, fp in zip(trows, fps):
                row["FH_probability"] = float(fp)
                rows.append(row)
                n_tumor += 1

        imgs = []
        metas = []

    for y in range(0, max(1, H - args.tile_size + 1), args.step_size):
        for x in range(0, max(1, W - args.tile_size + 1), args.step_size):
            n_grid += 1

            try:
                img = slide.read_region(
                    (int(x), int(y)),
                    0,
                    (args.tile_size, args.tile_size),
                ).convert("RGB")
            except Exception:
                n_failed += 1
                continue

            tfra = tissue_fraction(img)
            now = time.time()

            if n_grid % args.progress_every == 0 or now - last > 60:
                print(
                    f"{name} | "
                    f"grid={n_grid} "
                    f"tissue={n_tissue} "
                    f"tumor={n_tumor} "
                    f"failed={n_failed} "
                    f"elapsed={(now - t0) / 60:.1f}min",
                    flush=True,
                )
                last = now

            if tfra < args.min_tissue_fraction:
                continue

            n_tissue += 1

            imgs.append(img)
            metas.append(
                {
                    "slide_name": name,
                    "x": int(x),
                    "y": int(y),
                    "tissue_fraction": float(tfra),
                }
            )

            if len(imgs) >= args.read_batch_size:
                flush()

    flush()
    slide.close()

    df = pd.DataFrame(rows)
    df.to_csv(patch_csv, index=False, encoding="utf-8-sig")

    if n_tumor > 0:
        probs = (
            df.loc[df["is_tumor_associated"] == 1, "FH_probability"]
            .dropna()
            .astype(float)
            .values
        )
        final = float(np.mean(probs))
        pred = int(final >= args.fh_cutoff)
        status = "OK"
    else:
        final = np.nan
        pred = ""
        status = "NO_TUMOR_PATCHES"

    res = {
        "slide_name": name,
        "slide_path": slide_path,
        "status": status,
        "width": W,
        "height": H,
        "tile_size": args.tile_size,
        "step_size": args.step_size,
        "min_tissue_fraction": args.min_tissue_fraction,
        "tumor_threshold": args.tumor_threshold,
        "fh_cutoff": args.fh_cutoff,
        "n_grid_patches": n_grid,
        "n_tissue_patches": n_tissue,
        "n_tumor_patches": n_tumor,
        "tumor_fraction": n_tumor / max(n_tissue, 1),
        "final_fh_probability": final,
        "final_prediction": "FH" if pred == 1 else "nonFH" if pred == 0 else "",
        "final_pred_label_0.5": pred,
        "patch_predictions_csv": patch_csv,
        "heatmap_png": heatmap,
        "elapsed_minutes": (time.time() - t0) / 60,
    }

    pd.DataFrame([res]).to_csv(result_csv, index=False, encoding="utf-8-sig")

    if status == "OK":
        make_heatmap(df, final, name, heatmap)

    print(
        f"Done: {name} | "
        f"status={status} | "
        f"FH_prob={final} | "
        f"pred={res['final_prediction']} | "
        f"elapsed={res['elapsed_minutes']:.1f}min",
        flush=True,
    )

    return res


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Complete WCH-FHPM workflow: tissue filtering, tumor detector, "
            "FH predictor, heatmap."
        )
    )

    ap.add_argument("--slide", default="")
    ap.add_argument("--slide_dir", default="")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--models_dir", default="")
    ap.add_argument("--tumor_model", default="")
    ap.add_argument("--fh_model", default="")

    ap.add_argument("--tile_size", type=int, default=512)
    ap.add_argument("--step_size", type=int, default=512)
    ap.add_argument("--min_tissue_fraction", type=float, default=0.20)

    ap.add_argument("--tumor_threshold", type=float, default=0.20)
    ap.add_argument("--fh_cutoff", type=float, default=0.50)

    ap.add_argument("--read_batch_size", type=int, default=128)
    ap.add_argument("--batch_size", type=int, default=128)
    ap.add_argument("--progress_every", type=int, default=1000)

    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"])

    args = ap.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    args.models_dir = args.models_dir or os.path.join(script_dir, "models")
    args.tumor_model = args.tumor_model or os.path.join(
        args.models_dir,
        "WCH_TumorDetector_ViTBase512.pth",
    )
    args.fh_model = args.fh_model or os.path.join(
        args.models_dir,
        "WCH_FHPM_ViTBase512.pth",
    )

    os.makedirs(args.out_dir, exist_ok=True)

    device = torch.device(
        "cuda" if args.device == "cuda" and torch.cuda.is_available() else "cpu"
    )

    tumor_model, tumor_info = load_model(args.tumor_model, device)
    fh_model, fh_info = load_model(args.fh_model, device)

    if torch.cuda.device_count() > 1 and device.type == "cuda":
        print(f"Using DataParallel on {torch.cuda.device_count()} GPUs", flush=True)
        tumor_model = torch.nn.DataParallel(tumor_model)
        fh_model = torch.nn.DataParallel(fh_model)

    run_config = {
        "tumor_model": args.tumor_model,
        "fh_model": args.fh_model,
        "tile_size": args.tile_size,
        "step_size": args.step_size,
        "min_tissue_fraction": args.min_tissue_fraction,
        "tumor_threshold": args.tumor_threshold,
        "fh_cutoff": args.fh_cutoff,
        "device": str(device),
    }

    with open(
        os.path.join(args.out_dir, "run_config.json"),
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(run_config, f, indent=2, ensure_ascii=False)

    slides = [args.slide] if args.slide else find_slides(args.slide_dir) if args.slide_dir else []

    if not slides:
        raise ValueError("Please provide --slide or --slide_dir")

    results = []

    for s in slides:
        try:
            results.append(
                process_slide(
                    s,
                    args,
                    tumor_model,
                    tumor_info,
                    fh_model,
                    fh_info,
                    device,
                )
            )
        except Exception as e:
            print(f"[FAILED] {s}: {e}", flush=True)
            results.append(
                {
                    "slide_name": os.path.basename(s),
                    "slide_path": s,
                    "status": "FAILED",
                    "error": str(e),
                }
            )

        pd.DataFrame(results).to_csv(
            os.path.join(args.out_dir, "all_slide_results.csv"),
            index=False,
            encoding="utf-8-sig",
        )

    print("All done.")
    print("Summary:", os.path.join(args.out_dir, "all_slide_results.csv"))


if __name__ == "__main__":
    main()
=======
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import hashlib
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import openslide
import pandas as pd
import timm
import torch
from PIL import Image, ImageDraw, ImageFilter, ImageFont
from torchvision import transforms


SLIDE_EXTS = {".svs", ".ndpi", ".tif", ".tiff", ".mrxs", ".scn"}

try:
    RESAMPLING_NEAREST = Image.Resampling.NEAREST
except AttributeError:
    RESAMPLING_NEAREST = Image.NEAREST


def safe_name(s, max_len=120):
    s = str(s)
    b = re.sub(r"[^A-Za-z0-9\u4e00-\u9fff._-]+", "_", s).strip("._-")
    h = hashlib.md5(s.encode("utf-8", errors="ignore")).hexdigest()[:8]
    return f"{b[:max_len]}__{h}"


def find_slides(root):
    return sorted(
        str(p)
        for p in Path(root).rglob("*")
        if p.is_file() and p.suffix.lower() in SLIDE_EXTS
    )


def positive_index(label_mapping, task):
    if not isinstance(label_mapping, dict) or not label_mapping:
        return 1

    def norm(x):
        return str(x).lower().replace("_", "").replace("-", "").replace(" ", "")

    lm = {norm(k): int(v) for k, v in label_mapping.items()}

    if task == "fh":
        if "fh" in lm:
            return lm["fh"]  # exact match; avoids matching nonFH

        for k, v in lm.items():
            if not k.startswith("non") and k in {
                "fhrcc",
                "fhdrcc",
                "fhdrrcc",
                "fhdrc",
            }:
                return v

        return 1

    if task == "tumor":
        for k in ["tumor", "tumour", "cancer", "neoplastic"]:
            if k in lm:
                return lm[k]

        for k, v in lm.items():
            if (
                "normal" not in k
                and not k.startswith("non")
                and "negative" not in k
                and ("tumor" in k or "tumour" in k)
            ):
                return v

        return 1

    return 1


def load_model(path, device):
    ckpt = torch.load(path, map_location="cpu", weights_only=False)

    model_name = ckpt.get("model_name", "vit_base_patch16_224")
    img_size = int(ckpt.get("img_size", 512))
    num_classes = int(ckpt.get("num_classes", 2))

    state = (
        ckpt.get("state_dict")
        or ckpt.get("model_state_dict")
        or ckpt.get("model")
        or ckpt
    )

    state = {
        k.replace("module.", "", 1) if k.startswith("module.") else k: v
        for k, v in state.items()
    }

    model = timm.create_model(
        model_name,
        pretrained=False,
        num_classes=num_classes,
        img_size=img_size,
    )
    model.load_state_dict(state, strict=True)
    model.to(device)
    model.eval()

    prep = ckpt.get("preprocess", {})
    mean = prep.get("mean", prep.get("normalize_mean", [0.485, 0.456, 0.406]))
    std = prep.get("std", prep.get("normalize_std", [0.229, 0.224, 0.225]))

    info = {
        "img_size": img_size,
        "mean": mean,
        "std": std,
        "label_mapping": ckpt.get("label_mapping", {}),
        "model_path": str(path),
    }

    return model, info


def build_tf(info):
    return transforms.Compose(
        [
            transforms.Resize((info["img_size"], info["img_size"])),
            transforms.ToTensor(),
            transforms.Normalize(info["mean"], info["std"]),
        ]
    )


def tissue_fraction(img, sat_thr=20, gray_high=245, gray_low=10):
    a = np.asarray(img.convert("RGB"))

    mx = a.max(2).astype(np.int16)
    mn = a.min(2).astype(np.int16)
    sat = mx - mn
    gray = a.mean(2)

    return float(((sat > sat_thr) & (gray < gray_high) & (gray > gray_low)).mean())


def batch_predict(model, imgs, tf, device, pos, batch_size):
    out = []

    if not imgs:
        return out

    with torch.no_grad():
        for i in range(0, len(imgs), batch_size):
            x = torch.stack([tf(im) for im in imgs[i : i + batch_size]], 0).to(
                device,
                non_blocking=True,
            )

            with torch.cuda.amp.autocast(enabled=(device.type == "cuda")):
                p = torch.softmax(model(x), 1)[:, pos].detach().cpu().numpy()

            out += [float(v) for v in p]

    return out


def get_font(size, bold=False):
    font_paths = [
        (
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
            if bold
            else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
        ),
        (
            "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"
            if bold
            else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"
        ),
    ]

    for p in font_paths:
        if os.path.exists(p):
            return ImageFont.truetype(p, size)

    return ImageFont.load_default()


def interp(c1, c2, t):
    return tuple(
        np.clip(np.array(c1) * (1 - t) + np.array(c2) * t, 0, 255)
        .astype(np.uint8)
        .tolist()
    )


def color_prob(p):
    p = float(np.clip(0.5 + 1.18 * (float(np.clip(p, 0, 1)) - 0.5), 0, 1))

    anchors = [
        (0, (32, 72, 192)),
        (0.25, (92, 140, 230)),
        (0.5, (247, 245, 242)),
        (0.75, (239, 126, 92)),
        (1, (197, 26, 44)),
    ]

    for (p1, c1), (p2, c2) in zip(anchors[:-1], anchors[1:]):
        if p1 <= p <= p2:
            return interp(c1, c2, (p - p1) / (p2 - p1))

    return anchors[-1][1]


def rr_mask(size, r):
    m = Image.new("L", size, 0)
    d = ImageDraw.Draw(m)
    d.rounded_rectangle((0, 0, size[0] - 1, size[1] - 1), radius=r, fill=255)
    return m


def center_text(draw, xy, text, font, fill):
    b = draw.textbbox((0, 0), text, font=font)
    draw.text(
        (xy[0] - (b[2] - b[0]) / 2, xy[1] - (b[3] - b[1]) / 2),
        text,
        font=font,
        fill=fill,
    )


def make_heatmap(
    patch_df,
    final_prob,
    slide_name,
    out_png,
    canvas_w=1600,
    canvas_h=1400,
):
    df = patch_df[(patch_df["is_tumor_associated"] == 1)].copy()
    df["FH_probability"] = pd.to_numeric(df["FH_probability"], errors="coerce")
    df = df.dropna(subset=["FH_probability"])

    if len(df) == 0:
        return None

    xs = np.sort(df.x.unique())
    ys = np.sort(df.y.unique())
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: i for i, y in enumerate(ys)}

    grid = np.full((len(ys), len(xs), 3), 255, dtype=np.uint8)

    for _, r in df.iterrows():
        grid[yi[r.y], xi[r.x]] = color_prob(r.FH_probability)

    raw = Image.fromarray(grid, "RGB")
    bg = Image.new("RGBA", (canvas_w, canvas_h), (248, 248, 250, 255))

    # Card shadow
    shadow = Image.new("RGBA", bg.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle(
        (36, 46, canvas_w - 36, canvas_h - 26),
        radius=34,
        fill=(0, 0, 0, 42),
    )
    shadow = shadow.filter(ImageFilter.GaussianBlur(16))
    bg.alpha_composite(shadow)

    # Card
    card = Image.new("RGBA", (canvas_w - 72, canvas_h - 72), (255, 255, 255, 255))
    bg.paste(card, (36, 36), rr_mask(card.size, 34))

    draw = ImageDraw.Draw(bg)
    title = get_font(34, True)
    probfont = get_font(30, True)
    tickfont = get_font(20, False)

    center_text(
        draw,
        (canvas_w // 2, 92),
        "Heatmap of the Prediction",
        title,
        (25, 30, 40),
    )

    hx, hy, hw, hh = 105, 165, 1130, 980

    panel = Image.new("RGBA", (hw, hh), (252, 252, 253, 255))
    bg.paste(panel, (hx, hy), rr_mask((hw, hh), 22))

    draw.rounded_rectangle(
        (hx, hy, hx + hw, hy + hh),
        radius=22,
        outline=(232, 234, 238),
        width=2,
    )

    aw, ah = hw - 40, hh - 40
    rw, rh = raw.size
    sc = min(aw / rw, ah / rh)

    nw = max(1, int(rw * sc))
    nh = max(1, int(rh * sc))

    fit = raw.resize((nw, nh), RESAMPLING_NEAREST)

    mx = hx + 20 + (aw - nw) // 2
    my = hy + 20 + (ah - nh) // 2

    bg.paste(fit, (mx, my))

    # Colorbar
    cbx, cby, cbw, cbh = 1290, 240, 62, 740

    cb = Image.new("RGB", (cbw, cbh), "white")
    cd = ImageDraw.Draw(cb)

    for j in range(cbh):
        cd.line((0, j, cbw, j), fill=color_prob(1 - j / max(cbh - 1, 1)))

    bg.paste(cb, (cbx, cby), rr_mask((cbw, cbh), 12))
    draw = ImageDraw.Draw(bg)

    for t in [0.2, 0.4, 0.6, 0.8]:
        ty = cby + int((1 - t) * cbh)
        draw.line(
            (cbx + cbw + 8, ty, cbx + cbw + 20, ty),
            fill=(90, 94, 102),
            width=2,
        )
        draw.text(
            (cbx + cbw + 28, ty - 11),
            f"{t:.1f}",
            font=tickfont,
            fill=(88, 92, 98),
        )

    prefix = "Predicted probability of FHdRCC: "
    risk = "high risk" if final_prob >= 0.5 else "low risk"
    val = f"{final_prob:.3f} ({risk})"

    bw = draw.textbbox((0, 0), prefix, font=probfont)[2]
    vw = draw.textbbox((0, 0), val, font=probfont)[2]

    x = (canvas_w - bw - vw) // 2
    y = 1248

    draw.text((x, y), prefix, font=probfont, fill=(18, 28, 44))
    draw.text(
        (x + bw, y),
        val,
        font=probfont,
        fill=(153, 0, 52) if final_prob >= 0.5 else (22, 84, 153),
    )

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    bg.convert("RGB").save(out_png, quality=95)

    return out_png


def process_slide(
    slide_path,
    args,
    tumor_model,
    tumor_info,
    fh_model,
    fh_info,
    device,
):
    name = os.path.basename(slide_path)
    uid = safe_name(Path(name).stem)
    outdir = os.path.join(args.out_dir, uid)
    os.makedirs(outdir, exist_ok=True)

    patch_csv = os.path.join(outdir, "patch_predictions.csv")
    result_csv = os.path.join(outdir, "slide_result.csv")
    heatmap = os.path.join(outdir, "FH_probability_heatmap.png")

    t0 = time.time()

    slide = openslide.OpenSlide(slide_path)
    W, H = slide.dimensions

    tumor_tf = build_tf(tumor_info)
    fh_tf = build_tf(fh_info)

    tumor_pos = positive_index(tumor_info.get("label_mapping", {}), "tumor")
    fh_pos = positive_index(fh_info.get("label_mapping", {}), "fh")

    print(
        f"\n===== Slide: {name} =====\n"
        f"Size: {W} x {H}\n"
        f"tumor_pos_index={tumor_pos}; fh_pos_index={fh_pos}",
        flush=True,
    )

    rows = []
    imgs = []
    metas = []

    n_grid = 0
    n_tissue = 0
    n_failed = 0
    n_tumor = 0

    last = time.time()

    def flush():
        nonlocal imgs, metas, rows, n_tumor

        if not imgs:
            return

        tps = batch_predict(
            tumor_model,
            imgs,
            tumor_tf,
            device,
            tumor_pos,
            args.batch_size,
        )

        tims = []
        trows = []

        for im, meta, tp in zip(imgs, metas, tps):
            is_t = int(tp >= args.tumor_threshold)

            row = dict(meta)
            row.update(
                {
                    "tumor_probability": float(tp),
                    "is_tumor_associated": is_t,
                    "FH_probability": np.nan,
                }
            )

            if is_t:
                tims.append(im)
                trows.append(row)
            else:
                rows.append(row)

        if tims:
            fps = batch_predict(
                fh_model,
                tims,
                fh_tf,
                device,
                fh_pos,
                args.batch_size,
            )

            for row, fp in zip(trows, fps):
                row["FH_probability"] = float(fp)
                rows.append(row)
                n_tumor += 1

        imgs = []
        metas = []

    for y in range(0, max(1, H - args.tile_size + 1), args.step_size):
        for x in range(0, max(1, W - args.tile_size + 1), args.step_size):
            n_grid += 1

            try:
                img = slide.read_region(
                    (int(x), int(y)),
                    0,
                    (args.tile_size, args.tile_size),
                ).convert("RGB")
            except Exception:
                n_failed += 1
                continue

            tfra = tissue_fraction(img)
            now = time.time()

            if n_grid % args.progress_every == 0 or now - last > 60:
                print(
                    f"{name} | "
                    f"grid={n_grid} "
                    f"tissue={n_tissue} "
                    f"tumor={n_tumor} "
                    f"failed={n_failed} "
                    f"elapsed={(now - t0) / 60:.1f}min",
                    flush=True,
                )
                last = now

            if tfra < args.min_tissue_fraction:
                continue

            n_tissue += 1

            imgs.append(img)
            metas.append(
                {
                    "slide_name": name,
                    "x": int(x),
                    "y": int(y),
                    "tissue_fraction": float(tfra),
                }
            )

            if len(imgs) >= args.read_batch_size:
                flush()

    flush()
    slide.close()

    df = pd.DataFrame(rows)
    df.to_csv(patch_csv, index=False, encoding="utf-8-sig")

    if n_tumor > 0:
        probs = (
            df.loc[df["is_tumor_associated"] == 1, "FH_probability"]
            .dropna()
            .astype(float)
            .values
        )
        final = float(np.mean(probs))
        pred = int(final >= args.fh_cutoff)
        status = "OK"
    else:
        final = np.nan
        pred = ""
        status = "NO_TUMOR_PATCHES"

    res = {
        "slide_name": name,
        "slide_path": slide_path,
        "status": status,
        "width": W,
        "height": H,
        "tile_size": args.tile_size,
        "step_size": args.step_size,
        "min_tissue_fraction": args.min_tissue_fraction,
        "tumor_threshold": args.tumor_threshold,
        "fh_cutoff": args.fh_cutoff,
        "n_grid_patches": n_grid,
        "n_tissue_patches": n_tissue,
        "n_tumor_patches": n_tumor,
        "tumor_fraction": n_tumor / max(n_tissue, 1),
        "final_fh_probability": final,
        "final_prediction": "FH" if pred == 1 else "nonFH" if pred == 0 else "",
        "final_pred_label_0.5": pred,
        "patch_predictions_csv": patch_csv,
        "heatmap_png": heatmap,
        "elapsed_minutes": (time.time() - t0) / 60,
    }

    pd.DataFrame([res]).to_csv(result_csv, index=False, encoding="utf-8-sig")

    if status == "OK":
        make_heatmap(df, final, name, heatmap)

    print(
        f"Done: {name} | "
        f"status={status} | "
        f"FH_prob={final} | "
        f"pred={res['final_prediction']} | "
        f"elapsed={res['elapsed_minutes']:.1f}min",
        flush=True,
    )

    return res


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Complete WCH-FHPM workflow: tissue filtering, tumor detector, "
            "FH predictor, heatmap."
        )
    )

    ap.add_argument("--slide", default="")
    ap.add_argument("--slide_dir", default="")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--models_dir", default="")
    ap.add_argument("--tumor_model", default="")
    ap.add_argument("--fh_model", default="")

    ap.add_argument("--tile_size", type=int, default=512)
    ap.add_argument("--step_size", type=int, default=512)
    ap.add_argument("--min_tissue_fraction", type=float, default=0.20)

    ap.add_argument("--tumor_threshold", type=float, default=0.20)
    ap.add_argument("--fh_cutoff", type=float, default=0.50)

    ap.add_argument("--read_batch_size", type=int, default=128)
    ap.add_argument("--batch_size", type=int, default=128)
    ap.add_argument("--progress_every", type=int, default=1000)

    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"])

    args = ap.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    args.models_dir = args.models_dir or os.path.join(script_dir, "models")
    args.tumor_model = args.tumor_model or os.path.join(
        args.models_dir,
        "WCH_TumorDetector_ViTBase512.pth",
    )
    args.fh_model = args.fh_model or os.path.join(
        args.models_dir,
        "WCH_FHPM_ViTBase512.pth",
    )

    os.makedirs(args.out_dir, exist_ok=True)

    device = torch.device(
        "cuda" if args.device == "cuda" and torch.cuda.is_available() else "cpu"
    )

    tumor_model, tumor_info = load_model(args.tumor_model, device)
    fh_model, fh_info = load_model(args.fh_model, device)

    if torch.cuda.device_count() > 1 and device.type == "cuda":
        print(f"Using DataParallel on {torch.cuda.device_count()} GPUs", flush=True)
        tumor_model = torch.nn.DataParallel(tumor_model)
        fh_model = torch.nn.DataParallel(fh_model)

    run_config = {
        "tumor_model": args.tumor_model,
        "fh_model": args.fh_model,
        "tile_size": args.tile_size,
        "step_size": args.step_size,
        "min_tissue_fraction": args.min_tissue_fraction,
        "tumor_threshold": args.tumor_threshold,
        "fh_cutoff": args.fh_cutoff,
        "device": str(device),
    }

    with open(
        os.path.join(args.out_dir, "run_config.json"),
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(run_config, f, indent=2, ensure_ascii=False)

    slides = [args.slide] if args.slide else find_slides(args.slide_dir) if args.slide_dir else []

    if not slides:
        raise ValueError("Please provide --slide or --slide_dir")

    results = []

    for s in slides:
        try:
            results.append(
                process_slide(
                    s,
                    args,
                    tumor_model,
                    tumor_info,
                    fh_model,
                    fh_info,
                    device,
                )
            )
        except Exception as e:
            print(f"[FAILED] {s}: {e}", flush=True)
            results.append(
                {
                    "slide_name": os.path.basename(s),
                    "slide_path": s,
                    "status": "FAILED",
                    "error": str(e),
                }
            )

        pd.DataFrame(results).to_csv(
            os.path.join(args.out_dir, "all_slide_results.csv"),
            index=False,
            encoding="utf-8-sig",
        )

    print("All done.")
    print("Summary:", os.path.join(args.out_dir, "all_slide_results.csv"))


if __name__ == "__main__":
    main()
>>>>>>> 16e8c11e175de65394360c69a18aa6090b08807d
