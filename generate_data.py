import gzip
import pandas as pd

dfam_file = "Dfam-curated_only-1.embl.gz"
output_csv = "te_sequences.csv"
target_classes = ["SINE", "LINE", "LTR", "DNA"]

def parse_dfam_embl_gz(filepath, target_classes):
    records = []
    current_id = None
    current_type = None
    current_seq = []
    in_sq = False
    in_cc_annot = False

    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith("ID "):
                # Save previous record
                if current_id and current_type and current_seq:
                    if current_type in target_classes:
                        records.append({
                            "accession": current_id,
                            "class": current_type,
                            "sequence": "".join(current_seq).upper()
                        })
                # Reset
                current_id = line.split()[1].rstrip(';')
                current_type = None
                current_seq = []
                in_sq = False
                in_cc_annot = False

            # Detect start of RepeatMasker annotations
            elif line.startswith("CC   RepeatMasker Annotations:"):
                in_cc_annot = True

            # Extract Type from CC lines
            elif in_cc_annot and line.startswith("CC        Type:"):
                # Format: "CC        Type: SINE"
                current_type = line.split(":")[1].strip()
                in_cc_annot = False  # Type found, no need to stay in this block

            # Sequence start
            elif line.startswith("SQ "):
                in_sq = True

            # Sequence lines
            elif in_sq:
                if line.startswith("//"):
                    # End of record
                    if current_id and current_type and current_seq:
                        if current_type in target_classes:
                            records.append({
                                "accession": current_id,
                                "class": current_type,
                                "sequence": "".join(current_seq).upper()
                            })
                    in_sq = False
                else:
                    tokens = line.split()
                    for token in tokens:
                        if token.isalpha():
                            current_seq.append(token)
    return records

data = parse_dfam_embl_gz(dfam_file, target_classes)
df = pd.DataFrame(data)
df.to_csv(output_csv, index=False)
print(f"Saved {len(df)} unique records to {output_csv}")
print(df['class'].value_counts())

from sklearn.model_selection import train_test_split

input = "te_sequences.csv"
train_ratio = 0.8
dev_ratio = test_ratio = 0.1
random_state = 42

df = pd.read_csv(input)

print(f"Total records: {len(df)}")
print("Class distribution before split:")
print(df["class"].value_counts())

classes = sorted(df["class"].unique())
class_to_label = {cls: i for i, cls in enumerate(classes)}
df["label"] = df["class"].map(class_to_label)
print("\nLabel mapping:", class_to_label)

train, temp = train_test_split(
    df,
    test_size=(1 - train_ratio),
    stratify=df["label"],
    random_state=random_state
)

dev, test = train_test_split(
    temp,
    test_size=test_ratio / (dev_ratio + test_ratio),
    random_state=random_state
)

train_out = train[["sequence", "label"]]
dev_out = dev[["sequence", "label"]]
test_out = test[["sequence", "label"]]

train_out.to_csv("train.csv", index=False)
dev_out.to_csv("dev.csv", index=False)
test_out.to_csv("test.csv", index=False)

print(f"\nSaved:")
print(f"  train.csv: {len(train_out)} rows")
print(f"  dev.csv:   {len(dev_out)} rows")
print(f"  test.csv:  {len(test_out)} rows")

print("\nLabel distribution in train.csv:")
print(train["class"].value_counts())
