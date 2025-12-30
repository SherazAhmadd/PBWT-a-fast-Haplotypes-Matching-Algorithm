import numpy as np
import sys
import pandas as pd

class PBWT:
    def __init__(self, matrix):
        self.X = np.array(matrix, dtype=np.int8)
        self.N, self.M = self.X.shape
        self.a = []
        self.d = []

    def build(self):
        a_curr = np.arange(self.N)
        d_curr = np.zeros(self.N, dtype=int)

        self.a.append(a_curr.copy())
        self.d.append(d_curr.copy())

        for k in range(self.M):
            zeros_idx, ones_idx = [], []
            zeros_div, ones_div = [], []

            p = q = k + 1

            for i in range(self.N):
                idx = a_curr[i]
                div = d_curr[i]

                p = max(p, div)
                q = max(q, div)

                if self.X[idx, k] == 0:
                    zeros_idx.append(idx)
                    zeros_div.append(p)
                    p = 0
                else:
                    ones_idx.append(idx)
                    ones_div.append(q)
                    q = 0

            a_curr = np.array(zeros_idx + ones_idx)
            d_curr = np.array(zeros_div + ones_div)

            self.a.append(a_curr)
            self.d.append(d_curr)

    def report_long_matches(self, min_length):
        matches = []
        for k in range(1, self.M + 1):
            for i in range(1, self.N):
                if self.d[k][i] <= k - min_length:
                    h1 = self.a[k][i]
                    h2 = self.a[k][i - 1]
                    start = self.d[k][i]
                    matches.append((h1, h2, start, k))
        return matches


# =========================
# MAIN
# =========================
if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage: python pbwt_csv.py <haplotypes.csv|xlsx> <min_match_length>")
        sys.exit(1)

    file_path = sys.argv[1]
    min_L = int(sys.argv[2])

    print("Loading haplotype file...")

    # Load CSV or Excel
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
    elif file_path.endswith(".xlsx"):
        df = pd.read_excel(file_path)
    else:
        raise ValueError("File must be .csv or .xlsx")

    # Drop haplotype ID column
    matrix = df.iloc[:, 1:].values

    print(f"Haplotypes: {matrix.shape[0]}")
    print(f"Sites: {matrix.shape[1]}")
    print(f"Minimum match length: {min_L}")

    print("\nBuilding PBWT...")
    pbwt = PBWT(matrix)
    pbwt.build()

    print("Searching for matches...")
    results = pbwt.report_long_matches(min_L)

    # =========================
# SAVE RESULTS TO CSV
# =========================
    if results:
        out_df = pd.DataFrame(
            results,
            columns=["hap1", "hap2", "start_site", "end_site"]
        )
        out_df["length"] = out_df["end_site"] - out_df["start_site"]

        output_file = f"pbwt_matches_minL{min_L}.csv"
        out_df.to_csv(output_file, index=False)

        print(f"\nResults saved to: {output_file}")


    print("\n===== RESULTS =====")
    if not results:
        print("No matches found.")
    else:
        print(f"Found {len(results)} matches:\n")
        for h1, h2, start, end in results:
            print(
                f"Haplotype {h1} ↔ {h2} | "
                f"Start: {start}, End: {end}, Length: {end - start}"
            )
