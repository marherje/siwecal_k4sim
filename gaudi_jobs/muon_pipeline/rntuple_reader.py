import ROOT
import matplotlib.pyplot as plt
import pandas as pd


def read_file(filename = "ShipHits.root", container = "MTCSciFi"):

    n = ROOT.RNTupleReader.Open(container, filename).GetNEntries()
    if n == 0:
        print(f"[read_file] {container} is empty, returning empty DataFrame")
        return pd.DataFrame(columns=['E', 'x', 'y', 'z', 'bf_x', 'bf_y'])

    df = ROOT.RDF.FromRNTuple(container, filename)
    print(list(df.GetColumnNames()))

    df = pd.DataFrame(df.AsNumpy(['E', 'x', 'y', 'z', 'bf_x', 'bf_y']))
    print(df)
    return df

def plot_columns(df, det = "scint", label = ""):
    for col in df.columns:
        fig, ax = plt.subplots()
        df[col].hist(bins = 100, ax = ax)
        ax.set_xlabel(col)
        fig.savefig(f"plots/{det}_{col}{label}.pdf")


df_scint = read_file(filename = "ShipHits.root", container = "MTCScint")
df_scifi = read_file(filename = "ShipHits.root", container = "MTCSciFi")


out_label = "_added_light_propagation_3_corrected_cartesian_segmentation"
plot_columns(df_scifi, det = "scifi", label = out_label)
plot_columns(df_scint, det = "scint", label = out_label)

