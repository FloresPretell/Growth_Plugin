import os
import glob
import math
import re
from pathlib import Path
from PIL import Image, ImageDraw

def parse_case_name(case_path):
    """
    Extract D and V from folder names like:
    D12_V0p5     -> D=12,   V=0.5
    D1p5_V2      -> D=1.5,  V=2
    D1.5_V2.5    -> D=1.5,  V=2.5
    D2.5_V2p5    -> D=2.5,  V=2.5
    """
    name = os.path.basename(case_path.rstrip("/"))

    m = re.match(r"^D(\d+(?:[p.]\d+)?)_V(\d+(?:[p.]\d+)?)$", name)
    if not m:
        return None

    d_str = m.group(1).replace("p", ".")
    v_str = m.group(2).replace("p", ".")

    d_val = float(d_str)
    v_val = float(v_str)

    return d_val, v_val

def build_case_matrix(base_dir):
    """
    Finds all folders matching D*_V* and arranges them into a matrix:
    rows = sorted D values
    cols = sorted V values
    """
    case_dirs = sorted(glob.glob(os.path.join(base_dir, "D*_V*")))

    parsed = []
    d_values = set()
    v_values = set()

    for cdir in case_dirs:
        parsed_case = parse_case_name(cdir)
        if parsed_case is None:
            continue
        d_val, v_val = parsed_case
        parsed.append((d_val, v_val, cdir))
        d_values.add(d_val)
        v_values.add(v_val)

    d_values = sorted(d_values)
    v_values = sorted(v_values)

    # dictionary for quick lookup
    case_dict = {(d, v): path for d, v, path in parsed}

    return d_values, v_values, case_dict


def create_mosaic_gif_matrix(
    base_dir,
    output_file="mosaic_matrix.gif",
    tile_size=(140, 100),
    duration=250,
    keep_last_frame=True,
    png_pattern="*.png",
    png_subdir=os.path.join("png", "u_t")
):
    base_dir = str(Path(base_dir).expanduser().resolve())
    output_file = Path(output_file).expanduser().resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    d_values, v_values, case_dict = build_case_matrix(base_dir)

    if not d_values or not v_values:
        raise ValueError("No valid D*_V* case folders found.")

    nrows = len(d_values)
    ncols = len(v_values)

    # Load PNG lists per case
    all_cases = {}
    max_frames = 0

    for d in d_values:
            for v in v_values:
                case_path = case_dict.get((d, v), None)
                if case_path is None:
                    all_cases[(d, v)] = None
                    continue
    
                png_dir = os.path.join(case_path, png_subdir)
    
                if not os.path.isdir(png_dir):
                    all_cases[(d, v)] = []
                    continue
    
                png_files = sorted(glob.glob(os.path.join(png_dir, png_pattern)))
                if png_files:
                    all_cases[(d, v)] = png_files
                    max_frames = max(max_frames, len(png_files))
                else:
                    all_cases[(d, v)] = []
    
    if max_frames == 0:
            raise ValueError(f"No PNG files found under subdirectory '{png_subdir}'.")


    tile_w, tile_h = tile_size

    # Extra space:
    left_label_w = 70   # space for row labels (D values)
    top_label_h = 30    # space for column labels (V values)
    title_h = 18        # per tile title
    time_h = 15         # per tile time label

    full_tile_h = title_h + tile_h + time_h
    canvas_w = left_label_w + ncols * tile_w
    canvas_h = top_label_h + nrows * full_tile_h

    frames = []

    for frame_idx in range(max_frames):
        canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
        draw = ImageDraw.Draw(canvas)

        # Column headers: V values
        for col, v in enumerate(v_values):
            x = left_label_w + col * tile_w
            draw.text((x + 10, 5), f"V={v:g}", fill="black")

        # Row headers: D values
        for row, d in enumerate(d_values):
            y = top_label_h + row * full_tile_h
            draw.text((10, y + title_h + tile_h // 2 - 5), f"D={d}", fill="black")

        # Fill tiles
        for row, d in enumerate(d_values):
            for col, v in enumerate(v_values):
                x = left_label_w + col * tile_w
                y = top_label_h + row * full_tile_h

                case_name = f"D{d}_V{str(v).replace('.', 'p')}"
                files = all_cases[(d, v)]

                # Small case title
                draw.text((x + 3, y + 1), case_name, fill="black")

                if files is None:
                    # Folder missing entirely
                    blank = Image.new("RGB", (tile_w, tile_h), "lightgray")
                    canvas.paste(blank, (x, y + title_h))
                    draw.text((x + 10, y + title_h + tile_h // 2 - 5), "missing", fill="black")
                    continue

                if len(files) == 0:
                    # Folder exists but no png
                    blank = Image.new("RGB", (tile_w, tile_h), "gainsboro")
                    canvas.paste(blank, (x, y + title_h))
                    draw.text((x + 10, y + title_h + tile_h // 2 - 5), "no png", fill="black")
                    continue

                if frame_idx < len(files):
                    frame_to_use = frame_idx
                    finished = False
                else:
                    if keep_last_frame:
                        frame_to_use = len(files) - 1
                        finished = True
                    else:
                        frame_to_use = None
                        finished = True

                if frame_to_use is not None:
                    img = Image.open(files[frame_to_use]).convert("RGB")
                    img = img.resize((tile_w, tile_h))
                    canvas.paste(img, (x, y + title_h))

                    label = f"t={frame_to_use}"
                    if finished:
                        label += " (last)"
                    draw.text((x + 3, y + title_h + tile_h + 1), label, fill="gray")
                else:
                    blank = Image.new("RGB", (tile_w, tile_h), "lightgray")
                    canvas.paste(blank, (x, y + title_h))
                    draw.text((x + 10, y + title_h + tile_h // 2 - 5), "finished", fill="black")

        frames.append(canvas)

    frames[0].save(
        str(output_file),
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=0
    )

    print(f"Saved GIF to {output_file}")
    print(f"Matrix size: {nrows} rows (D) x {ncols} cols (V)")
    print(f"D values: {d_values}")
    print(f"V values: {v_values}")


# create_mosaic_gif_matrix(
#     base_dir="/scratch/flore0a/Dataset_test1/Cases",
#     output_file="/scratch/flore0a/Dataset_test1/Cases/mosaic_matrix.gif",
#     tile_size=(120, 90),
#     duration=250,
#     keep_last_frame=True,
#     png_pattern="*.png",
#     png_subdir=os.path.join("Simulation2D", "merged", "png", "u_t")
# )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create a mosaic GIF from D*_V* case folders.")
    parser.add_argument("base_dir", help="Base directory containing D*_V* folders")
    parser.add_argument("--output", default="mosaic_matrix.gif", help="Output GIF filename")
    parser.add_argument("--tile_size", type=int, nargs=2, default=(140, 100), help="Tile size (width height)")
    parser.add_argument("--duration", type=int, default=250, help="Frame duration in ms")
    parser.add_argument("--keep_last_frame", action="store_true", help="Keep last frame for cases that finish early")
    parser.add_argument("--png_pattern", default="*.png", help="Pattern to match PNG files in subdirectories")
    parser.add_argument("--png_subdir", default=os.path.join("Simulation2D", "merged", "png", "u_t"), help="Subdirectory under each case folder where PNGs are located")

    args = parser.parse_args()
    create_mosaic_gif_matrix(
        base_dir=args.base_dir,
        output_file=args.output,
        tile_size=tuple(args.tile_size),
        duration=args.duration,
        keep_last_frame=args.keep_last_frame,
        png_pattern=args.png_pattern,
        png_subdir=args.png_subdir,
    )
