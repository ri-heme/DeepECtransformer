import csv
import logging
import subprocess

BLAST_HEADERS = "qseqid sseqid evalue bitscore score qlen slen length pident"


def run_blastp(target_fasta, blastp_result, db_dir, threads=1):
    logging.info("BLASTp prediction starts on the dataset")
    subprocess.call(
        f"diamond blastp -d {db_dir} -q {target_fasta} -o {blastp_result} "
        f"--threads {threads} --id 50 --outfmt 6 {BLAST_HEADERS}",
        shell=True,
        stderr=subprocess.STDOUT,
    )
    logging.info("BLASTp prediction ended on the dataset")


def read_best_blast_result(blastp_result):
    query_db_set_info = {}
    with open(blastp_result, "r") as fp:
        reader = csv.DictReader(fp, fieldnames=BLAST_HEADERS.split(), delimiter="\t")
        for line in reader:
            query_id = line["qseqid"]
            db_id = line["sseqid"]

            ec_numbers = db_id.split("|")[1].strip()
            score = float(line["score"])
            qlen = float(line["qlen"])
            length = float(line["length"])
            pident = float(line["pident"])

            if pident < 50:
                continue
            coverage = length / qlen
            if coverage >= 0.75:
                if query_id not in query_db_set_info:
                    query_db_set_info[query_id] = [ec_numbers, score]
                else:
                    prev_score = query_db_set_info[query_id][1]
                    if score > prev_score:
                        query_db_set_info[query_id] = [ec_numbers, score]
                    # Merge results if multiple best hits
                    elif score == prev_score:
                        query_db_set_info[query_id][0] += ec_numbers

    # Format EC numbers as string separated with semicolon
    output_dict = {}
    for query_id, (ec_numbers, _) in query_db_set_info.items():
        ec_fmt = set()
        for ec_num in ec_numbers.split(";"):
            if not ec_num.startswith("EC:"):
                ec_num = f"EC:{ec_num}"
            ec_fmt.add(ec_num)
        output_dict[query_id] = ";".join(sorted(ec_fmt))

    return query_db_set_info


def merge_predictions(dl_pred_result, blastp_pred_result, output_dir):
    dl_pred = {}
    with open(dl_pred_result, "r") as f1:
        f1.readline()
        while True:
            line = f1.readline()
            if not line:
                break
            seq_id, ec, _ = line.strip().split("\t")
            if seq_id not in dl_pred:
                dl_pred[seq_id] = ec
            else:
                dl_pred[seq_id] += f";{ec}"

    blastp_pred = {}
    with open(blastp_pred_result, "r") as f2:
        f2.readline()
        while True:
            line = f2.readline()
            if not line:
                break
            seq_id, ec = line.strip().split("\t")
            blastp_pred[seq_id] = ec

    merged = {}
    for seq_id, pred in dl_pred.items():
        if pred == "None":
            if seq_id in blastp_pred:
                merged[seq_id] = blastp_pred[seq_id]
            else:
                merged[seq_id] = "None"
        else:
            merged[seq_id] = dl_pred[seq_id]

    with open(f"{output_dir}/DeepECv2_result.txt", "w") as fp:
        fp.write("sequence_ID\tprediction\n")
        for seq_id, pred in merged.items():
            fp.write(f"{seq_id}\t{pred}\n")
