import pandas as pd

yue_df = pd.read_csv("../../merged_ppi_6_angstrom.tsv", sep='\t')

contacts_a_df = yue_df[["protein_I_II", "protein_A", "contact_A", "source"]]
contacts_b_df = yue_df[["protein_I_II", "protein_B", "contact_B", "source"]]
contacts_a_df.columns = contacts_b_df.columns = ["domain_id", "protein_id", "positions", "source"]

contacts_df = pd.concat([contacts_a_df, contacts_b_df], ignore_index=True).drop_duplicates()
contacts_df["domain_id"] = contacts_df["domain_id"] + "_" + contacts_df["protein_id"]


def encode_positions(raw_positions):
    """Re-encode a "[23, 24, 25, ...]" bracketed residue list as a
    compact, lossless run-length "23-25,27,30-31" segment list."""
    positions = sorted(int(p) for p in raw_positions[1:-1].split(", "))
    segments = []
    start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
            continue
        segments.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = pos
    segments.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(segments)


contacts_df["positions"] = contacts_df["positions"].apply(encode_positions)

contacts_df.to_csv("../../reference_data/ppi_binding_sites.tsv", sep='\t', index=False)