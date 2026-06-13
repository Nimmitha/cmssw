# ZeeJmm CMSSW Scripts

This directory is for CMSSW-side helpers only. Use it for rootupler or
preselection production tasks that must run inside the ZeeJmm CMSSW setup.

The current script entry point here is:

```bash
./scripts/run_Rootupler.sh
```

Reduced candidate production, candidate diagnostics, and cut scans now live in
the centralized `ml_reboot` workflow:

```bash
cd /home/nimmitha/LPCfiles/run2/reboot/ml_reboot

python3 workflow/selection/make_zeejmm_candidates.py
python3 workflow/selection/diagnose_candidates.py --channel zeejmm
python3 workflow/selection/scan_cuts.py --channel zeejmm --candidate-variant onecand
```

The old Zee-only `scan_cuts.py`, `howto.txt`, and duplicate-candidate inspector
were retired because the reduced candidate files are now produced and inspected
from `ml_reboot`.
