import ROOT
import pandas as pd


def get_lines_data():
    # Open file
    f = ROOT.TFile.Open("plots/tracker/HitEfficiency_vs_IntLumiRunIII_Layers_run3.root")

    # Get the canvas
    c = f.Get("HitEfficiency_vs_IntLumiRunIII_Layers_run3")

    # List all primitives drawn on the canvas
    primitives = c.GetListOfPrimitives()

    lines = [p for p in primitives if isinstance(p, ROOT.TLine)]
    print("Total TLine objects:", len(lines))

    def is_vertical(L, tol=1e-9):
        return abs(L.GetX1() - L.GetX2()) < tol

    def is_horizontal(L, tol=1e-9):
        return abs(L.GetY1() - L.GetY2()) < tol

    line_rows = []
    for i, L in enumerate(lines):
        row = {
            "i": i,
            "x1": L.GetX1(), "x2": L.GetX2(),
            "y1": L.GetY1(), "y2": L.GetY2(),
            "color": L.GetLineColor(),      # ROOT color index
            "style": L.GetLineStyle(),      # 1=solid, 2=dashed, 3=dotted, ...
            "width": L.GetLineWidth(),
            "orientation": "vertical" if is_vertical(L) else ("horizontal" if is_horizontal(L) else "slanted")
        }
        line_rows.append(row)

    lines_df = pd.DataFrame(line_rows)

    # Example mapping — update once you see actual numbers in lines_df
    STYLE_MAP = {
        ("vertical", 860, 9): "Layer1 HV change",     # dashed, magenta-ish
        ("vertical", 1, 6): "Gain calibration",     # dotted, black
        ("vertical", 810, 9): "Technical stop",     # dashed, orange
        ("vertical", 1, 1): "Year",     # dashed, orange
    }
    def label_line(row):
        return STYLE_MAP.get((row["orientation"], row["color"], row["style"]), "Other")
    lines_df["label"] = lines_df.apply(label_line, axis=1)

    # print(lines_df)
    # lines_df
    return lines_df


def get_charge_plot_FPIX_data():
    # Open file
    f = ROOT.TFile.Open("plots/tracker/AvgOnCluChargeNorm_vs_IntLumiRunIII_Disks_run3.root")

    # Get the canvas
    c = f.Get("AvgOnCluChargeNorm_vs_IntLumiRunIII_Disks_run3")

    # List all primitives drawn on the canvas
    primitives = c.GetListOfPrimitives()

    dfs = {}

    for lay in [1, 2, 3]:
        h = primitives.FindObject(f"AvgOnCluChargeNorm_vs_IntLumiRunIII_Disk{lay}_2021Data")
        assert h, f"Layer {lay} histogram not found"

        xs, ys, yerrs = [], [], []
        for b in range(1, h.GetNbinsX() + 1):
            x = h.GetBinCenter(b)
            y = h.GetBinContent(b)
            e = h.GetBinError(b)
            if y != 0 or e != 0:
                if y > 50:
                    continue
                xs.append(x)
                ys.append(y)
                yerrs.append(e)

        # build a temporary df for this layer
        df = pd.DataFrame({
            "x": xs,
            f"y{lay}_chr_FPIX": ys,
            # f"yerr{lay}": yerrs
        })
        dfs[lay] = df

    # merge all dfs on "x"
    from functools import reduce
    merged = reduce(lambda left, right: pd.merge(left, right, on="x", how="outer"), dfs.values())

    # average the y values and propagate the errors
    ycols = [f"y{lay}_chr_FPIX" for lay in dfs.keys()]
    # ycols = [f"y{lay}" for lay in [1, 2, 3]]
    # yerrcols = [f"yerr{lay}" for lay in dfs.keys()]
    merged["y_chr_FPIX"] = merged[ycols].mean(axis=1)
    # merged["yerr_chr_FPIX"] = (merged[yerrcols]**2).sum(axis=1)**0.5 / len(yerrcols)
    # merged = merged[["x", "y", "yerr"]].sort_values("x")
    merged
    # return merged[["x", "y_chr_FPIX", "yerr_chr_FPIX"]].sort_values("x")
    return merged.sort_values("x")

# get_charge_plot_FPIX_data()


def get_efficiency_plot_data():
    f = ROOT.TFile.Open("plots/tracker/HitEfficiency_vs_IntLumiRunIII_Layers_run3.root")
    # Get the canvas
    c = f.Get("HitEfficiency_vs_IntLumiRunIII_Layers_run3")
    # List all primitives drawn on the canvas
    primitives = c.GetListOfPrimitives()

    dfs = {}

    for lay in [1, 2]:
        h = primitives.FindObject(f"HitEfficiency_vs_IntLumiRunIII_2022Data_Lay{lay}")
        assert h, f"Layer {lay} histogram not found"

        xs, ys, yerrs = [], [], []
        for b in range(1, h.GetNbinsX() + 1):
            x = h.GetBinCenter(b)
            y = h.GetBinContent(b)
            e = h.GetBinError(b)
            if y != 0 or e != 0:
                if y > 0.99:
                    continue
                xs.append(x)
                ys.append(y)
                yerrs.append(e)

        # build a temporary df for this layer
        df = pd.DataFrame({
            "x": xs,
            f"y{lay}_eff": ys,
            # f"yerr{lay}": yerrs
        })
        dfs[lay] = df

    # merge all dfs on "x"
    from functools import reduce
    merged = reduce(lambda left, right: pd.merge(left, right, on="x", how="outer"), dfs.values())

    # average the y values and propagate the errors
    # ycols = [f"y{lay}" for lay in dfs.keys()]
    ycols = [f"y{lay}_eff" for lay in dfs.keys()]
    # yerrcols = [f"yerr{lay}" for lay in dfs.keys()]
    merged["y_eff"] = merged[ycols].mean(axis=1)
    # merged["yerr"] = (merged[yerrcols]**2).sum(axis=1)**0.5 / len(yerrcols)

    # merged["y_eff"] = merged["y"]**2
    # merged["yerr_eff"] = 2 * merged["y"]**0.5 * merged["yerr"]
    # merged = merged[["x", "y", "yerr"]].sort_values("x")
    merged

    # return merged[["x", "y_eff", "yerr_eff"]].sort_values("x")
    return merged.sort_values("x")

# get_efficiency_plot_data()

def get_charge_plot_BPIX_data():
    # Open file
    f = ROOT.TFile.Open("plots/tracker/AvgOnCluChargeNorm_vs_IntLumiRunIII_Layers_run3.root")

    # Get the canvas
    c = f.Get("AvgOnCluChargeNorm_vs_IntLumiRunIII_Layers_run3")

    # List all primitives drawn on the canvas
    primitives = c.GetListOfPrimitives()

    dfs = {}

    for lay in [1, 2, 3, 4]:
        h = primitives.FindObject(f"AvgOnCluChargeNorm_vs_IntLumiRunIII_Lay{lay}_2021Data")
        assert h, f"Layer {lay} histogram not found"

        xs, ys, yerrs = [], [], []
        for b in range(1, h.GetNbinsX() + 1):
            x = h.GetBinCenter(b)
            y = h.GetBinContent(b)
            e = h.GetBinError(b)
            if y != 0 or e != 0:
                if y > 50:
                    continue
                xs.append(x)
                ys.append(y)
                yerrs.append(e)

        # build a temporary df for this layer
        df = pd.DataFrame({
            "x": xs,
            f"y{lay}_chr_BPIX": ys,
            # f"yerr{lay}": yerrs
        })
        dfs[lay] = df

    # merge all dfs on "x"
    from functools import reduce
    merged = reduce(lambda left, right: pd.merge(left, right, on="x", how="outer"), dfs.values())

    # average the y values and propagate the errors
    ycols = [f"y{lay}_chr_BPIX" for lay in dfs.keys()]
    # ycols = [f"y{lay}_chr_BPIX" for lay in [3, 4]]
    # yerrcols = [f"yerr{lay}" for lay in dfs.keys()]
    merged["y_chr_BPIX"] = merged[ycols].mean(axis=1)
    # merged["yerr_chr_BPIX"] = (merged[yerrcols]**2).sum(axis=1)**0.5 / len(yerrcols)
    # merged = merged[["x", "y", "yerr"]].sort_values("x")
    merged
    # return merged[["x", "y_chr_BPIX", "yerr_chr_BPIX"]].sort_values("x")
    return merged.sort_values("x")

# get_charge_plot_BPIX_data()