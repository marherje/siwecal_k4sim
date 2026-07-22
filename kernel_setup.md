# Fast Clone of the `siwecal_k4sim` Jupyter Kernel

Use this when the other account is on the same server and can copy files from:

```bash
/home/llr/ilc/shi
```

Assume the target account has the project at:

```bash
~/code/siwecal_k4sim
```

## 1. Copy the Existing Kernel

```bash
mkdir -p ~/.local/share/jupyter/kernels
cp -r /home/llr/ilc/shi/.local/share/jupyter/kernels/siwecal_k4sim_2026_02_01 \
      ~/.local/share/jupyter/kernels/
```

## 2. Replace the Old User Path

```bash
KDIR=~/.local/share/jupyter/kernels/siwecal_k4sim_2026_02_01

sed -i "s#/home/llr/ilc/shi#${HOME}#g" \
    "$KDIR/run.sh" \
    "$KDIR/kernel.json"

chmod +x "$KDIR/run.sh"
```

## 3. Check

```bash
jupyter kernelspec list
```

You should see:

```text
siwecal_k4sim_2026_02_01
```

In Jupyter, select:

```text
Python (siwecal_k4sim key4hep 2026-02-01)
```

## 4. Quick Test in Notebook

```python
import ROOT
print(ROOT.gROOT.GetVersion())
```

Expected:

```text
6.38.00
```

If it fails, check the startup log:

```bash
cat /tmp/siwecal_k4sim_2026_02_01_env.log
```
