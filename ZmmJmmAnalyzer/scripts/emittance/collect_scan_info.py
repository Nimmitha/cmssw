import os
import json
import pandas as pd

# Configuration
INPUT_DIR = "scan_JSON"
OUTPUT_DIR = "scan_summary"
POINTS_REQUIRED = 9  # Change this as needed

os.makedirs(OUTPUT_DIR, exist_ok=True)


def is_valid_scan(content, points_req):
    """Check if all scans have at least points_req"""
    scans = content.get("Scans", [])
    num_points = [scan.get("NumPoints", 0) for scan in scans]
    if len(num_points) != 2:
        print(f"Invalid number of scans: {len(num_points)}")
        return False

    if (num_points[0] == num_points[1] == points_req):
        return True
    else:
        return False


def load_and_validate_json(filepath):
    """Load JSON content and apply early filtering"""
    with open(filepath, "r") as f:
        content = json.load(f)

    if not is_valid_scan(content, POINTS_REQUIRED):
        return None, f"Skipped: Invalid number of points (not {POINTS_REQUIRED})"
    return content, None


def validate_orthogonal_displacement(scan_type, x_disp, y_disp, step_index, filename):
    """Warn if orthogonal displacement is unexpectedly non-zero"""
    if scan_type == "X1" and abs(y_disp) > 1e-6:
        print(f"YDisplacement not zero in X1 scan ({filename}, step {step_index + 1}): {y_disp}")
    elif scan_type == "Y1" and abs(x_disp) > 1e-6:
        print(f"XDisplacement not zero in Y1 scan ({filename}, step {step_index + 1}): {x_disp}")


def extract_scan_rows(content):
    """Extract per-step scan information"""
    rows = []
    fill = content.get("Fill")
    date = content.get("Date")
    run_value = content.get("Run", [None])[0]

    for scan in content.get("Scans", []):
        scan_name = scan.get("ScanName", "Unknown")
        scan_number = scan.get("ScanNumber", 0)
        scan_type = scan_name
        start_times = scan.get("StartTimes", [])
        stop_times = scan.get("StopTimes", [])
        rel_displacements = scan.get("RelativeDisplacements", [])
        x_displacements = scan.get("XDisplacements", [])
        y_displacements = scan.get("YDisplacements", [])
        n_steps = len(start_times)

        if not (len(stop_times) == len(rel_displacements) == n_steps):
            print(f"Mismatch in scan {scan_name} step lengths.")
            continue

        for i in range(n_steps):
            validate_orthogonal_displacement(scan_type, x_displacements[i], y_displacements[i], i, scan_name)

            rows.append({
                "fill": fill,
                "run": run_value,
                "scan": scan_number,
                "type": scan_type,
                "step": i + 1,
                "tStart": start_times[i],
                "tStop": stop_times[i],
                "sep": rel_displacements[i],
            })
    return rows


def extract_summary_row(content):
    """Extract scan summary information"""
    run_list = content.get("Run", [])
    run_value = run_list[0] if run_list else None

    for key in ["Scan_1", "Scan_2", "InputDIPFile", "ParticleTypeB2", "Offset", "ScanNames"]:
        content.pop(key, None)

    row = {
        "Fill": content.get("Fill"),
        "Date": content.get("Date"),
        "Run": run_value,
        "BetaStar": content.get("BetaStar"),
        "Angle": content.get("Angle"),
        "ParticleTypeB1": content.get("ParticleTypeB1"),
        "EnergyB1": content.get("EnergyB1"),
        "EnergyB2": content.get("EnergyB2"),
    }

    for i, scan_type in enumerate(content.get("ScanTypes", [])):
        row[f"ScanType_{i+1}"] = scan_type

    for i, time_window in enumerate(content.get("ScanTimeWindows", [])):
        if len(time_window) == 2:
            row[f"ScanTimeStart_{i+1}"] = time_window[0]
            row[f"ScanTimeStop_{i+1}"] = time_window[1]

    for i, scan in enumerate(content.get("Scans", [])):
        row[f"NumPoints_{i+1}"] = scan.get("NumPoints", 0)

    return row


def save_scan_csv(rows, fill, run):
    """Save per-step scan CSV"""
    df = pd.DataFrame(rows)
    outname = f"scan_{fill}_Run{run}.csv"
    outpath = os.path.join(OUTPUT_DIR, outname)
    df.to_csv(outpath, index=False)
    print(f"Saved scan CSV: {outname}")


def main():
    summary_rows = []

    for filename in os.listdir(INPUT_DIR):
        if not filename.startswith("Scan_") or not filename.endswith(".json"):
            continue

        filepath = os.path.join(INPUT_DIR, filename)
        content, skip_reason = load_and_validate_json(filepath)

        if content is None:
            print(f"{filename}: {skip_reason}")
            continue

        scan_rows = extract_scan_rows(content)
        if scan_rows:
            save_scan_csv(scan_rows, content.get("Fill"), content.get("Run", [None])[0])

        summary_rows.append(extract_summary_row(content))

    # Build and save summary DataFrame
    df_summary = pd.DataFrame(summary_rows)
    for col in df_summary.columns:
        if "ScanTimeStart" in col or "ScanTimeStop" in col:
            df_summary[col] = pd.to_datetime(df_summary[col], unit="s", errors="coerce")

    df_summary = df_summary.sort_values(by=["Fill", "Date"])
    summary_path = os.path.join(OUTPUT_DIR, "scan_summary.csv")
    df_summary.to_csv(summary_path, index=False)
    print(f"\nSaved summary CSV: scan_summary.csv")


if __name__ == "__main__":
    main()
