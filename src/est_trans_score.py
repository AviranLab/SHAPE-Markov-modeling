"""
Schemes
    --scheme in : I-Paired / NI-Paired.
    --scheme pu : Paired / Unpaired.

Models
    --model 1order   : first-order Markov chain over the whole class
    --model 1dropend : first-order Markov chain estimated per motif with the last nucleotide of each motif dropped
    --model 2order   : second-order Markov chain

Input formats (all three files hold multiple sequences, keyed by ">name")
    --fa     >name / sequence
    --dot    >name / sequence / dot-bracket structure
    --shape  >name / reactivities

Output
    A headerless CSV: 2x2 for 1order/1dropend (rows = current state low/high,
    columns = next state low/high), 4x2 for 2order.

Usage
    python est_trans_score.py \\
        --dot data/all.dot --fa data/all.fa --shape data/all.shape \\
        --model 1order --scheme in --threshold 0.20 \\
        --out score_1order_thd20.csv
"""

import argparse
import itertools
import os
from collections import defaultdict
import numpy as np
import pandas as pd


def identify_function(shape_data_list, dot_bracket, sequence, name):
    stack = []
    pair_map = {}

    for i, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(i)
        elif char == ")":
            if stack:
                j = stack.pop()
                pair_map[j] = i

    sorted_pairs = sorted(pair_map.items())
    group = sorted(sorted_pairs, key=lambda x: x[0], reverse=True)
    list_x0 = [x[0] for x in sorted_pairs]
    list_x1 = [x[1] for x in sorted_pairs]

    df = pd.DataFrame({
        "index": range(len(dot_bracket)),
        "type": None,
        "paired_index": None,
    })

    old_branch_index = None

    for i in range(len(group)):
        if i == 0:
            if group[i][1] - group[i][0] == 1:
                df.loc[group[i][0]] = [group[i][0], "Stacked", None]
                df.loc[group[i][1]] = [group[i][1], "Stacked", None]
            else:
                for j in range(group[i][1] - group[i][0] + 1):
                    if j == 0:
                        df.loc[group[i][0] + j] = [group[i][0] + j, "Helix_end", None]
                    elif j == group[i][1] - group[i][0]:
                        df.loc[group[i][0] + j] = [group[i][0] + j, "Helix_end", None]
                    else:
                        df.loc[group[i][0] + j] = [group[i][0] + j, "Hairpin", None]

        elif i == len(group) - 1:
            df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
            df.loc[group[i][1]] = [group[i][1], "Helix_end", None]

        else:
            if group[i - 1][0] - group[i][0] == 1 and group[i][1] - group[i - 1][1] == 1:
                df.loc[group[i][0]] = [group[i][0], "Stacked", None]
                df.loc[group[i][1]] = [group[i][1], "Stacked", None]

            elif group[i - 1][0] - group[i][0] == 1 and group[i][1] - group[i - 1][1] > 1:
                if all(not (group[i - 1][1] < x < group[i][1]) for x in list_x0):
                    df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                    df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                    df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]
                    for m in range(group[i][1] - group[i - 1][1]):
                        if m == 0:
                            df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                        else:
                            df.loc[group[i][1] - m] = [group[i][1] - m, "Bulge", None]
                else:
                    df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                    df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                    df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                    df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]

            elif group[i - 1][0] - group[i][0] > 1 and group[i][1] - group[i - 1][1] == 1:
                df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]
                for m in range(group[i - 1][0] - group[i][0]):
                    if m == 0:
                        df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                    else:
                        df.loc[group[i][0] + m] = [group[i][0] + m, "Bulge", None]

            elif group[i - 1][0] - group[i][0] > 1 and group[i][1] - group[i - 1][1] > 1:
                if old_branch_index is None:
                    if all(not (group[i - 1][1] < x < group[i][1]) for x in list_x0):
                        df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                        df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]
                        for m in range(group[i - 1][0] - group[i][0]):
                            if m == 0:
                                df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                            else:
                                df.loc[group[i][0] + m] = [group[i][0] + m, "Internal_Loop", f"{name}_Int_Loop_{i}"]
                        for m in range(group[i][1] - group[i - 1][1]):
                            if m == 0:
                                df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                            else:
                                df.loc[group[i][1] - m] = [group[i][1] - m, "Internal_Loop", f"{name}_Int_Loop_{i}"]
                    else:
                        df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                        df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                        df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                        df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]
                else:
                    if group[old_branch_index][1] - group[i][1] > 0:
                        df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                        df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]
                        for m in range(group[i - 1][0] - group[i][0]):
                            if m == 0:
                                df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                            else:
                                df.loc[group[i][0] + m] = [group[i][0] + m, "Internal_Loop", f"{name}_Int_Loop_{i}"]
                        for m in range(group[i][1] - group[i - 1][1]):
                            if m == 0:
                                df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                            else:
                                df.loc[group[i][1] - m] = [group[i][1] - m, "Internal_Loop", f"{name}_Int_Loop_{i}"]
                    else:
                        df.loc[group[i][0]] = [group[i][0], "Helix_end", None]
                        df.loc[group[i][1]] = [group[i][1], "Helix_end", None]
                        old_branch_index = None
                        df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                        df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]

            elif group[i][1] - group[i - 1][0] < 0:
                old_branch_index = i - 1
                if group[i][1] - group[i][0] == 1:
                    df.loc[group[i][0]] = [group[i][0], "Stacked", None]
                    df.loc[group[i][1]] = [group[i][1], "Stacked", None]
                else:
                    for j in range(group[i][1] - group[i][0] + 1):
                        if j == 0:
                            df.loc[group[i][0] + j] = [group[i][0] + j, "Helix_end", None]
                        elif j == group[i][1] - group[i][0]:
                            df.loc[group[i][0] + j] = [group[i][0] + j, "Helix_end", None]
                        else:
                            df.loc[group[i][0] + j] = [group[i][0] + j, "Hairpin", None]
                    df.loc[group[i - 1][0]] = [group[i - 1][0], "Helix_end", None]
                    df.loc[group[i - 1][1]] = [group[i - 1][1], "Helix_end", None]

    for i in range(len(df)):
        if df.loc[i, "type"] is None:
            df.loc[i] = [i, "Free_Single_Stranded", None]
        else:
            break

    for i in range(len(df)):
        idx = len(df) - 1 - i
        if df.loc[idx, "type"] is None:
            df.loc[idx] = [idx, "Free_Single_Stranded", None]
        else:
            break

    for i in range(len(df)):
        if df.loc[i, "type"] is None:
            df.loc[i] = [i, "Multi_Loop", None]

    df["sequence"] = list(sequence)
    df = df[["index", "sequence", "type", "paired_index"]]
    df["index"] = df["index"] + 1
    df["shape_value"] = shape_data_list

    for k in range(len(group)):
        index1 = list_x0[k]
        index2 = list_x1[k]
        df.loc[index1, "paired_index"] = index2 + 1
        df.loc[index2, "paired_index"] = index1 + 1

    df["paired_index"] = df["paired_index"].fillna("Single")
    return df


def _clean_shape_value(token):
    # Missing values (-999, nan) are marked NaN here; the rest is clipped to [0, 3]
    value = float(token)
    if not np.isfinite(value) or value <= -500:
        return float("nan")
    return min(max(value, 0.0), 3.0)


def read_shape_data(path):
    shape_data = {}

    with open(path, "r") as f:
        lines = [x.strip() for x in f if x.strip()]

    if not any(line.startswith(">") for line in lines):
        raise ValueError(
            "For multi-sequence input, the .shape file must use '>name' headers "
            "followed by the reactivities of that sequence."
        )

    records = {}
    name = None
    for line in lines:
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = []
        elif name is not None:
            records[name].append(line.split())

    for name, rows in records.items():
        # Two accepted layouts: one line of values, or one "position value" pair per line.
        if len(rows) > 1 and all(len(r) == 2 for r in rows):
            shape_data[name] = [_clean_shape_value(r[1]) for r in rows]
        else:
            shape_data[name] = [_clean_shape_value(tok) for row in rows for tok in row]

    return shape_data


def read_dot_data(path):
    dot_data = {}

    with open(path, "r") as f:
        lines = [x.strip() for x in f if x.strip()]

    if len(lines) % 3 != 0 or not lines[0].startswith(">"):
        raise ValueError(
            f"{path}: expected repeating triplets of '>name', sequence and "
            f"dot-bracket structure lines (got {len(lines)} non-empty lines)."
        )

    for i in range(0, len(lines), 3):
        name_line, seq_line, dot_line = lines[i], lines[i + 1], lines[i + 2]

        if not name_line.startswith(">"):
            raise ValueError(f"{path}: expected a '>name' line, got: {name_line[:40]}")

        name = name_line[1:].strip()
        dot_data[name] = {
            "sequence": seq_line,
            "structure": dot_line,
        }

    return dot_data


def read_fa_data(path):
    fa_data = {}
    current_name = None
    seq_lines = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    fa_data[current_name] = "".join(seq_lines)
                current_name = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

    if current_name is not None:
        fa_data[current_name] = "".join(seq_lines)

    return fa_data


def build_dataframes(dot_path, fa_path, shape_path):
    shape_data = read_shape_data(shape_path)
    dot_data = read_dot_data(dot_path)
    fa_data = read_fa_data(fa_path)

    df_dict = {}
    name_list = []

    for name, info in dot_data.items():
        if name not in shape_data:
            print(f"[WARN] Missing shape data for {name}; skipped.")
            continue

        if name in fa_data:
            sequence = list(fa_data[name])
        else:
            sequence = list(info["sequence"])

        structure = list(info["structure"])
        shape_values = shape_data[name]

        if len(sequence) != len(structure) or len(sequence) != len(shape_values):
            print(
                f"[WARN] Length mismatch for {name}: "
                f"seq={len(sequence)}, dot={len(structure)}, shape={len(shape_values)}; skipped."
            )
            continue

        df = identify_function(shape_values, structure, sequence, name)
        df["dot_bracket"] = structure
        # Drop missing (-999 -> NaN) 
        df = df.dropna(subset=["shape_value"]).reset_index(drop=True)

        df_dict[name] = df
        name_list.append(name)

    stacked_parts, not_stacked_parts, paired_parts, unpaired_parts = [], [], [], []

    for name in name_list:
        df_temp = df_dict[name].copy()
        df_temp["name"] = name

        is_paired = df_temp["type"].isin(["Stacked", "Helix_end"])
        is_stacked = df_temp["type"] == "Stacked"

        stacked_parts.append(df_temp[is_stacked].copy())
        not_stacked_parts.append(df_temp[~is_stacked].copy())
        paired_parts.append(df_temp[is_paired].copy())
        unpaired_parts.append(df_temp[~is_paired].copy())

    if not name_list:
        raise ValueError("No usable sequence found; check the --dot / --fa / --shape inputs.")

    df_stacked = pd.concat(stacked_parts, ignore_index=True)
    df_not_stacked = pd.concat(not_stacked_parts, ignore_index=True)
    df_paired = pd.concat(paired_parts, ignore_index=True)
    df_unpaired = pd.concat(unpaired_parts, ignore_index=True)

    return df_stacked, df_not_stacked, df_paired, df_unpaired


def _split_by_H_windows(current_motif_rows, base_motif_id, h_value="Helix_end",
                        compress_H=False, drop_empty=True):
    rows = list(current_motif_rows)
    types = [r["type"] for r in rows]

    if compress_H:
        new_rows, new_types = [], []
        for i, t in enumerate(types):
            if t == h_value and new_types and new_types[-1] == h_value:
                continue
            new_rows.append(rows[i])
            new_types.append(types[i])
        rows, types = new_rows, new_types

    h_pos = [i for i, t in enumerate(types) if t == h_value]

    if len(h_pos) == 0:
        return [f"{base_motif_id}_0"], [pd.DataFrame(rows)]

    if len(h_pos) == 1 and (h_pos[0] == 0 or h_pos[0] == len(types) - 1):
        return [f"{base_motif_id}_0"], [pd.DataFrame(rows)]

    if len(h_pos) < 2:
        return [f"{base_motif_id}_0"], [pd.DataFrame(rows)]

    dfs, ids = [], []
    counter = 0

    for a, b in zip(h_pos, h_pos[1:]):
        if drop_empty and b == a + 1:
            continue
        dfs.append(pd.DataFrame(rows[a:b + 1]))
        ids.append(f"{base_motif_id}_{counter}")
        counter += 1

    return ids, dfs


def trans_matrix_two_status(data, threshold):
    data = data.copy()

    data["status"] = np.where(
        pd.isna(data["shape_value"]),
        "missing",
        np.where(data["shape_value"] < threshold, "low", "high")
    )

    status_list = data["status"].tolist()
    states = ["low", "high", "missing"]
    state_mapping = {s: i for i, s in enumerate(states)}

    transition_counts = np.zeros((2, 2))

    data = data.reset_index(drop=True)

    for i in range(len(data) - 1):
        if (
            data.loc[i + 1, "index"] - data.loc[i, "index"] == 1
            and data.loc[i + 1, "name"] == data.loc[i, "name"]
        ):
            current_state = state_mapping[status_list[i]]
            next_state = state_mapping[status_list[i + 1]]

            if current_state in [0, 1] and next_state in [0, 1]:
                transition_counts[current_state, next_state] += 1

    return transition_counts


def cal_score_two_status(data1, data2, threshold):
    count1 = trans_matrix_two_status(data1, threshold)
    count2 = trans_matrix_two_status(data2, threshold)

    transition_d1 = count1 / count1.sum()
    transition_d2 = count2 / count2.sum()

    score = np.zeros((2, 2))

    for i in range(2):
        for j in range(2):
            score[i, j] = np.log(transition_d1[i, j] + 1e-10) - np.log(transition_d2[i, j] + 1e-10)

    return score


def split_motifs_from_df(df):
    df = df.copy().reset_index(drop=True)
    motifs = {}
    current_motif = []
    motif_id = 0

    prev_index = None
    prev_name = None
    prev_type = None

    for _, row in df.iterrows():
        current_index = row["index"]
        current_name = row["name"]
        current_type = row["type"]

        if (
            prev_index is not None
            and current_index == prev_index + 1
            and current_name == prev_name
            and not (current_type == "Helix_end" and prev_type == "Helix_end")
        ):
            current_motif.append(row)
        else:
            if current_motif:
                ids, dfs = _split_by_H_windows(
                    current_motif,
                    motif_id,
                    h_value="Helix_end",
                    compress_H=False,
                    drop_empty=True,
                )
                if len(ids) > 1:
                    for _id, _df in zip(ids, dfs):
                        motifs[_id] = _df
                else:
                    motifs[motif_id] = pd.DataFrame(current_motif)
                motif_id += 1

            current_motif = [row]

        prev_index = current_index
        prev_name = current_name
        prev_type = current_type

    if current_motif:
        ids, dfs = _split_by_H_windows(
            current_motif,
            motif_id,
            h_value="Helix_end",
            compress_H=False,
            drop_empty=True,
        )
        if len(ids) > 1:
            for _id, _df in zip(ids, dfs):
                motifs[_id] = _df
        else:
            motifs[motif_id] = pd.DataFrame(current_motif)

    return motifs


def trans_matrix_two_status_topseq(data, threshold, drop=1):
    df = data.copy()

    df["status"] = np.where(
        pd.isna(df["shape_value"]),
        "missing",
        np.where(df["shape_value"] < threshold, "low", "high")
    )

    motifs = split_motifs_from_df(df)

    count_matrix = np.zeros((2, 2))
    state_mapping = {"low": 0, "high": 1, "missing": 2}

    for _, df_motif in motifs.items():
        df_motif_top = df_motif.iloc[:len(df_motif) - drop].copy()

        for i in range(len(df_motif_top) - 1):
            current_state = state_mapping[df_motif_top.iloc[i]["status"]]
            next_state = state_mapping[df_motif_top.iloc[i + 1]["status"]]

            if current_state in [0, 1] and next_state in [0, 1]:
                count_matrix[current_state, next_state] += 1

    return count_matrix


def cal_score_two_status_topseq(data1, data2, threshold):
    count1 = trans_matrix_two_status_topseq(data1, threshold)
    count2 = trans_matrix_two_status_topseq(data2, threshold)

    transition_d1 = count1 / count1.sum()
    transition_d2 = count2 / count2.sum()

    score = np.zeros((2, 2))

    for i in range(2):
        for j in range(2):
            score[i, j] = np.log(transition_d1[i, j] + 1e-10) - np.log(transition_d2[i, j] + 1e-10)

    return score


def trans_matrix_k_order(data, threshold, k):
    data = data.copy().reset_index(drop=True)

    data["status"] = np.where(
        pd.isna(data["shape_value"]),
        "missing",
        np.where(data["shape_value"] < threshold, "low", "high")
    )

    status_list = data["status"].tolist()
    name_list = data["name"].tolist()
    index_list = data["index"].tolist()

    transition_counts = defaultdict(int)

    for i in range(k, len(data)):
        if (
            name_list[i - k] != name_list[i]
            or any(index_list[j + 1] - index_list[j] != 1 for j in range(i - k, i))
        ):
            continue

        prev_states = tuple(status_list[i - k:i])
        curr_state = status_list[i]

        if curr_state in ["low", "high"] and all(s in ["low", "high"] for s in prev_states):
            transition_counts[(prev_states, curr_state)] += 1

    return transition_counts


def cal_score_k_order(data1, data2, threshold, k):
    states = ("low", "high")

    count1 = trans_matrix_k_order(data1, threshold, k)
    count2 = trans_matrix_k_order(data2, threshold, k)

    total1 = sum(count1.values())
    total2 = sum(count2.values())

    trans1 = {key: val / total1 for key, val in count1.items()}
    trans2 = {key: val / total2 for key, val in count2.items()}

    prefixes = list(itertools.product(states, repeat=k))
    state_idx = {s: j for j, s in enumerate(states)}
    prefix_idx = {p: i for i, p in enumerate(prefixes)}

    score_mat = np.zeros((len(prefixes), len(states)))

    for prefix in prefixes:
        for s in states:
            p1 = max(trans1.get((prefix, s), 0.0), 1e-10)
            p2 = max(trans2.get((prefix, s), 0.0), 1e-10)
            score_mat[prefix_idx[prefix], state_idx[s]] = np.log(p1 / p2)

    return score_mat


def main():
    parser = argparse.ArgumentParser(
        description="Estimate Markov transition score matrices from SHAPE and structure files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("--dot", required=True, help="Input dot-bracket file. Multi-sequence format: >name, sequence, structure.")
    parser.add_argument("--fa", required=True, help="Input FASTA file. Multi-sequence format: >name, sequence.")
    parser.add_argument("--shape", required=True, help="Input SHAPE file. Multi-sequence format: >name, then one line of values or one 'position value' pair per line.")
    parser.add_argument("--model", required=True, choices=["1order", "1dropend", "2order"],
                        help="Markov model. '1dropend' is the Dropend model (published as score_1orddropend_thd{NN}.csv).")
    parser.add_argument("--scheme", choices=["in", "pu"],
                        help="Markov scheme, as in Fold: 'in' = I-Paired/NI-Paired (primary matrices, "
                             "--mc-matrix-a/b); 'pu' = Paired/Unpaired (--mc-pu-matrix-a/b).")
    parser.add_argument("--category", choices=["stacked", "paired"],
                        help="Deprecated alias for --scheme ('stacked' = in, 'paired' = pu).")
    parser.add_argument("--threshold", required=True, type=float,
                        help="SHAPE reactivity threshold separating the 'low' and 'high' states, e.g. 0.20.")
    parser.add_argument("--out", required=True, help="Output CSV path (headerless).")

    args = parser.parse_args()

    scheme = args.scheme or {"stacked": "in", "paired": "pu"}.get(args.category)
    if scheme is None:
        parser.error("--scheme is required (in or pu).")

    df_stacked, df_not_stacked, df_paired, df_unpaired = build_dataframes(
        dot_path=args.dot,
        fa_path=args.fa,
        shape_path=args.shape,
    )

    if scheme == "in":
        data1 = df_not_stacked.reset_index(drop=True)   # NI-Paired
        data2 = df_stacked.reset_index(drop=True)       # I-Paired
    else:
        data1 = df_unpaired.reset_index(drop=True)      # Unpaired
        data2 = df_paired.reset_index(drop=True)        # Paired

    if args.model == "1order":
        score = cal_score_two_status(data1, data2, args.threshold)
    elif args.model == "1dropend":
        score = cal_score_two_status_topseq(data1, data2, args.threshold)
    elif args.model == "2order":
        score = cal_score_k_order(data1, data2, args.threshold, k=2)
    else:
        raise ValueError(f"Unknown model: {args.model}")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    pd.DataFrame(score).to_csv(args.out, index=False, header=False)

    print(f"[OK] Saved score matrix to: {args.out}")
    print(f"[INFO] Model = {args.model}")
    print(f"[INFO] Scheme = {scheme}")
    print(f"[INFO] Threshold = {args.threshold}")


if __name__ == "__main__":
    main()