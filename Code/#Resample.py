from pathlib import Path

import pandas as pd


INPUT_CSV = Path(r"C:\Users\Carlos Coronel\Downloads\MaQ.csv")
OUTPUT_CSV = Path(r"C:\Users\Carlos Coronel\Downloads\MaQ_quarterly.csv")
DATE_COL = "observation_date"


def main() -> None:
    df = pd.read_csv(INPUT_CSV)

    if DATE_COL not in df.columns:
        raise KeyError(f"Missing required date column: {DATE_COL}")

    # Parse dd/mm/yyyy dates and drop invalid date rows.
    df[DATE_COL] = pd.to_datetime(df[DATE_COL], format="%d/%m/%Y", errors="coerce")
    df = df.dropna(subset=[DATE_COL]).sort_values(DATE_COL)

    # Force numeric parsing for all non-date columns.
    non_date_cols = [col for col in df.columns if col != DATE_COL]
    for col in non_date_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    numeric_cols = [col for col in non_date_cols if pd.api.types.is_numeric_dtype(df[col])]
    if not numeric_cols:
        raise ValueError("No numeric columns found to resample.")

    df = df[[DATE_COL] + numeric_cols].set_index(DATE_COL)

    # Quarterly average with quarter-start labels.
    quarterly_mean = df[numeric_cols].resample("QS").mean()
    quarterly_count = df[numeric_cols].resample("QS").count()

    # Keep only quarters with 3 valid monthly observations in every numeric column.
    complete_quarters = (quarterly_count == 3).all(axis=1)
    quarterly = quarterly_mean.loc[complete_quarters].reset_index()

    quarterly[DATE_COL] = quarterly[DATE_COL].dt.strftime("%d/%m/%Y")
    quarterly.to_csv(OUTPUT_CSV, index=False)

    print(f"Saved quarterly CSV to: {OUTPUT_CSV}")
    print(f"Rows exported: {len(quarterly)}")


if __name__ == "__main__":
    main()
